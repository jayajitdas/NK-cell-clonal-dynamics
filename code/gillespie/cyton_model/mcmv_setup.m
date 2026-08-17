function mcmv_setup()

% Only using one cell cycle stage for now
n_cycle = 1;
end_t = 8; % 8 days
bI=0.662;
bM=1.441;
dI=0.129;
dM=0.964;
r=0.135;
kI = bI-dI;
kM = bM-dM;

% stoich = [-1 1 0 0 0 0 0
%     0 0 0 -1 1 0 0
%     0 0 0 0 -1 1 0
%     0 0 0 0 0 -1 1
%     0 0 1 1 0 0 -1
%     0 0 0 -1 1 0 0
%     0 0 0 0 -1 1 0
%     0 0 0 0 0 -1 1
%     0 0 1 1 0 0 -1];

% [CD27+, CD27-, n_divisions, "cell cycle stage" (just a 1 for constant division potential)]
stoich = [-1 1 0 0 % differentiation
    0 0 1 0 % immature division
    0 0 1 0]; % mature division
%add extra column for time
stoich(:,end+1) = zeros(length(stoich(:,1)),1);

% rxn_reactants = [1 0 1 0 0 0 0
%     1 0 0 1 0 0 0
%     1 0 0 0 1 0 0
%     1 0 0 0 0 1 0
%     1 0 0 0 0 0 1
%     0 1 0 1 0 0 0
%     0 1 0 0 1 0 0
%     0 1 0 0 0 1 0
%     0 1 0 0 0 0 1];

rxn_reactants = [1 0 1 0 % differentiation
    1 0 0 1 % immature division
    0 1 0 1]; % mature division

% k = [0.126/24/2
%     0.1
%     0.15
%     0.4
%     0.6
%     [0.1
%     0.15
%     0.4
%     0.6]/1.5]';

% Rates for single cell cycle stage that generate negative correlation

% Comment either this or the k = repmat(rate_options...) at the begining of
% the z loop along with the z loop.
k = [r/2 % Differentiation
    kI % immature proliferation
    kM]; % mature proliferation
k = repmat(k,n_cycle,1);

% Data to determine residuals
moments = [1091.1412404642856,4170.573045250178,15726003.519255938,226877736.17836666,10125044.067471344];
errors = [520.413683,1955.93348,14821133.0,173384159,5416863.80];

% Randomly sample some possible rates, we'll iterate them and see if they
% produce small residuals
%rate_options = rand(3,10000);

for z=1%:10000 % iterate over parameter scan
tic;

%k = repmat(rate_options(:,z),n_cycle,1);

for n=1:10000 % iterate over number of clones
% Stochastic simulation for each iteration/clone

% op_m = [1 0 0 1 0 0 0 0];
op_m = [1 0 0 1 0];
finished_cells = [];
n_cell_bool=true;
% while any cell still exists in the simulation
while any(op_m,'all') && n_cell_bool
    
    %determine each cell cycle stage
    state = zeros(length(op_m(:,1)),1);
    for i=1:length(state)
        state(i) = find(op_m(i,end-n_cycle:end-1));
    end
    %calculate all the possible reactions and their propensities
    options = propensities_from_reactants(rxn_reactants,op_m(:,1:end-1),k,state);
    %cumulative sum of propensities so we can effectively choose one
    options = cumsum(options,2);
    %calculate delta t
    dt = log(1./rand(length(op_m(:,1)),1)).*(1./options(:,end));
    %if any cells would progress past end_t, remove them from simulation
    finished_cells = [finished_cells; op_m(op_m(:,end)+dt>end_t,:)];
    op_m(:,end) = op_m(:,end) + dt;
    options(op_m(:,end)>end_t,:) = [];
    op_m(op_m(:,end)>end_t,:) = [];
    %generate random numbers; these will dictate which reaction each cell
    %chooses from its options
    w = rand(length(op_m(:,1)),1);
    %translate random number into a value corresponding to sum of
    %propensities
    choices = w.*(options(:,end));
    %find which reaction corresponds to each choice for each cell
    rxn = zeros(length(op_m(:,1)),1);
    for i=1:length(op_m(:,1))
        rxn(i) = find(choices(i)<options(i,:),1);
    end
    %update each cell according to the stoichiometry of the chosen reaction
    op_m = op_m + stoich(rxn,:);
    
    % Division doesn't impact protein abundances for this simulation
    %for cells that divide, we must randomly assign protein values to each
    %daughter cell
    % initialize
        % first_cell_n_protein = op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:)));
        % % for each cell
        % for i=1:length(first_cell_n_protein(:,1))
        %     % for each protein
        %     for j=1:length(first_cell_n_protein(1,:))
        %         % assign a number of proteins to the first cell
        %         % binomial probability; 50% chance of each cell getting
        %         % each protein
        %         first_cell_n_protein(i,j) = sum(rand(first_cell_n_protein(i,j),1)>0.5);
        %     end
        % end
        % assign remaining to second cell
        % second_cell_n_protein = op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:)))-first_cell_n_protein;
    % for each protein that doesn't obey inheritence, revert all the
    % phospho-protein to unphosphorylated
    % for i=1:length(uninherit_rxn)
    %     first_cell_n_protein(uninherit_rxn(i,:)==1,:) = first_cell_n_protein(uninherit_rxn(i,:)==1,:)+first_cell_n_protein(uninherit_rxn(i,:)==-1,:);
    %     second_cell_n_protein(uninherit_rxn(i,:)==1,:) = second_cell_n_protein(uninherit_rxn(i,:)==1,:)+second_cell_n_protein(uninherit_rxn(i,:)==-1,:);
    %     first_cell_n_protein(uninherit_rxn(i,:)==-1,:) = 0;
    %     second_cell_n_protein(uninherit_rxn(i,:)==-1,:) = 0;
    % end
    % add the daughter cells to op_m
    op_m = [op_m; op_m(rxn==2|rxn==length(stoich(:,1)),:)];

    % If clone size is too big, stop the simulation so we don't spend too
    % much time on it
    % if length(op_m(:,1))+length(finished_cells(:,1)) > 200000
    %     n_cell_bool = false;
    % end
    %op_m(rxn==length(stoich(:,1)),1:length(rxn_stoich(1,:))) = first_cell_n_protein;
    %op_m = [op_m; second_cell_n_protein op_m(rxn==length(stoich(:,1)),length(rxn_stoich(1,:))+1:end)];

end

% Determine n_pos and n_neg for the clone

% For clones that were too big and infeasible (see end of loop above)
% if ~n_cell_bool
%     pos(n) = 200001;
%     neg(n) = 0;
%     clone_size(n) = 0;
%     percent_cd27(n) = 0;
% else
pos(n) = sum(finished_cells(:,1));
neg(n) = sum(finished_cells(:,2));

clone_size(n) = length(finished_cells(:,1));
percent_cd27(n) = mean(finished_cells(:,1))*100;
% end

end

% Remove clones we interrupted
wrong = pos==200001;
pos(wrong) = [];
neg(wrong) = [];
percent_cd27(wrong) = [];
clone_size(wrong) = [];

% Calculate correlation and moments for these parameters
disp(corrcoef(percent_cd27,clone_size))
temp = corrcoef(percent_cd27,clone_size);
corr_cd27(z) = temp(1,2);
mean_pos(z) = mean(pos);
mean_neg(z) = mean(neg);
std_pos(z) = mean(pos.^2);
std_neg(z) = mean(neg.^2);
cov_cd27(z) = mean(pos.*neg);

%disp([mean_pos,mean_neg,std_pos,std_neg,cov_cd27])
%scatter(percent_cd27,clone_size)

toc;

end

% Calculate least squares of all parameters
resid = repmat(moments,z,1)-[mean_pos;mean_neg;std_pos;std_neg;cov_cd27]';
lsq = sum((resid./errors).^2,2);
%disp(mean([mean_pos;mean_neg;std_pos;std_neg;cov_cd27]',1))
%disp(std([mean_pos;mean_neg;std_pos;std_neg;cov_cd27]',1)./mean([mean_pos;mean_neg;std_pos;std_neg;cov_cd27]',1))

end