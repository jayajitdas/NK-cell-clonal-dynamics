function revision_figures()

% 3-stage model results

n = 100000;
n_points = 56;
bI	= 1.064;
bINT	=1.613;
bM	=0.700;
dI	=0.054;
dINT	=0.285;
dM	=1.370 ;
r1	=0.126;
r2	=0.164;

growth_M = [bI 0 0
    r1 bINT 0
    0 r2 bM];

for i=1:n
cells = gillespie_with_birth_and_death([1,0,0],growth_M,[dI,dINT,dM],8);
clone_size(i) = sum(cells);
nI(i) = cells(1);
nM(i) = cells(2)+cells(3);
cd27_pos(i) = nI(i)/clone_size(i)*100;
end

cd27_pos(clone_size==0)=[];
clone_size(clone_size==0)=[];
[temp,p]=corrcoef(cd27_pos(1:n_points),clone_size(1:n_points));
temp=corrcoef(cd27_pos,clone_size);
figure()
scatter(cd27_pos(1:n_points),clone_size(1:n_points),'filled')
set(gca,'yscale','log','fontsize',20)
xlabel('%CD27+', 'fontweight','bold')
ylabel('Clone Size','fontweight','bold')
disp(p(1,2))
disp(temp(1,2))

clear;

% 2-stage model results

bI=0.662;
bM=1.441;
dI=0.129;
dM=0.964;
r=0.135;
kI = bI-dI;
kM = bM-dM;
n = 100000;
n_points = 56;

% Simulation for 2-stage
growth_M = [kI 0
    r kM];

for i=1:n
cells = gillespie_with_birth_and_death([1,0],growth_M,[0,0],8);
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
nI(i) = cells(1);
nM(i) = cells(2);
end

[temp,p]=corrcoef(cd27_pos(1:n_points),clone_size(1:n_points));
temp=corrcoef(cd27_pos,clone_size);
out(1,:) = [temp(1,2) p(1,2) mean(nI) mean(nM) mean(nI.^2) mean(nM.^2) mean(nI.*nM)];
figure()
scatter(cd27_pos(1:n_points),clone_size(1:n_points),'filled')
set(gca,'yscale','log','fontsize',20)
xlabel('%CD27+', 'fontweight','bold')
ylabel('Clone Size','fontweight','bold')
title('Base')

% Simulation for 2-stage with time-dependence
% Parameters of the underlying normal distribution
mu = 1.0;        % mean of the normal distribution
sigma = 0.5;     % standard deviation

% Generate standard normal random numbers and transform
X = mu + sigma * randn(n, 1);  % randn is standard normal
t = exp(X);                    % Y follows lognormal(mu, sigma^2)
t(t>8) = 8;
growth_M = [kI 0
    r kM];

for i=1:n
cells = gillespie_with_birth_and_death([1,0],growth_M,[0,0],8-t(i));
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
nI(i) = cells(1);
nM(i) = cells(2);
end

[temp,p]=corrcoef(cd27_pos(1:n_points),clone_size(1:n_points));
temp=corrcoef(cd27_pos,clone_size);
out(2,:) = [temp(1,2) p(1,2) mean(nI) mean(nM) mean(nI.^2) mean(nM.^2) mean(nI.*nM)];
figure()
%scatter(cd27_pos,clone_size,[],t,'filled')
scatter(cd27_pos(1:n_points),clone_size(1:n_points),'filled')
set(gca,'yscale','log','fontsize',20)
xlabel('%CD27+', 'fontweight','bold')
ylabel('Clone Size','fontweight','bold')
title('Activation Time Distribution')

% Simulation for 2-stage with birth and death
growth_M = [bI 0
    r bM];

for i=1:n
cells = gillespie_with_birth_and_death([1,0],growth_M,[dI,dM],8);
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
nI(i) = cells(1);
nM(i) = cells(2);
end

cd27_pos(clone_size==0)=[];
clone_size(clone_size==0)=[];
[temp,p]=corrcoef(cd27_pos(1:n_points),clone_size(1:n_points));
temp=corrcoef(cd27_pos,clone_size);
out(3,:) = [temp(1,2) p(1,2) mean(nI) mean(nM) mean(nI.^2) mean(nM.^2) mean(nI.*nM)];
figure()
scatter(cd27_pos(1:n_points),clone_size(1:n_points),'filled')
set(gca,'yscale','log','fontsize',20)
xlabel('%CD27+', 'fontweight','bold')
ylabel('Clone Size','fontweight','bold')
title('Birth and Death')

% Simulation for both conditions
growth_M = [bI 0
    r bM];

for i=1:n
cells = gillespie_with_birth_and_death([1,0],growth_M,[dI,dM],8-t(i));
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
nI(i) = cells(1);
nM(i) = cells(2);
end

cd27_pos(clone_size==0)=[];
clone_size(clone_size==0)=[];
[temp,p]=corrcoef(cd27_pos(1:n_points),clone_size(1:n_points));
temp=corrcoef(cd27_pos,clone_size);
out(4,:) = [temp(1,2) p(1,2) mean(nI) mean(nM) mean(nI.^2) mean(nM.^2) mean(nI.*nM)];
figure()
scatter(cd27_pos(1:n_points),clone_size(1:n_points),'filled')
set(gca,'yscale','log','fontsize',20)
xlabel('%CD27+', 'fontweight','bold')
ylabel('Clone Size','fontweight','bold')
title('Activation Time Distribution and Death')



end