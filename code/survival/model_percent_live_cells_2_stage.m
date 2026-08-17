function out = model_percent_live_cells_2_stage(vec)

b1 = vec(1);
b2 = vec(2);
d1 = vec(3);
d2 = vec(4);
diff1 = vec(5);
asym = vec(6);

clearance = vec(7);

x0 = [100 100 0 0];
t = linspace(0,7,100);
%t = 4;

M = [b1+asym 0 0 0;
    diff1+asym b2 0 0;
    d1 0 -clearance 0;
    0 d2 0 -clearance];
for i=1:length(M)
    M(i,i) = 2*M(i,i)-sum(M(:,i));
end

for i=1:length(t)
    n_t(i,:) = expm(M*t(i))*x0';
    percent(i,:) = (n_t(i,1:2)./(n_t(i,1:2)+n_t(i,3:4)))*100;
end

out = percent(end,1)>percent(end,2);

end