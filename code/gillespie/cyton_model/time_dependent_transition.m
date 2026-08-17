function time_dependent_transition()

n=1000;
n_loop=10000;
min_corr = 1;

bI=rand(n_loop,1)*2;
%bM=bI*rand();
bM=rand(n_loop,1)*2;
r1=rand(n_loop,1)*2;

% Parameters of the underlying normal distribution
mu = 1.0;        % mean of the normal distribution
sigma = 0.5;     % standard deviation

% Generate standard normal random numbers and transform
X = mu + sigma * randn(n, 1);  % randn is standard normal
t = exp(X);                    % Y follows lognormal(mu, sigma^2)
t(t>8) = 8;

t= rand(n,1)*8;

for z=1:n_loop
%bI=1.064;
%bM=0.700;
%r1=0.126;



growth_M = [bI(z) 0
    r1(z) bM(z)];

for i=1:n
cells = gillespie_with_birth_and_death([1,0],growth_M,[0,0],t(i));
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
nI = cells(1);
nM = cells(2);
end

temp=corrcoef(cd27_pos,clone_size);
out(z,:) = [temp(1,2) mean(nI) mean(nM) mean(nI.^2) mean(nM.^2) mean(nI.*nM)];

if temp(1,2)<min_corr
    min_corr = temp(1,2);
    min_M = growth_M;
    disp(growth_M)
    disp(temp(1,2))
end

end

end