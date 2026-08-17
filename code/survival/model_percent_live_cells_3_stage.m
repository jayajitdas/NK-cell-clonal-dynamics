function [more_live,percent] = model_percent_live_cells_3_stage(vec)

b1 = vec(1);
b2 = vec(2);
b3 = vec(3);
d1 = vec(4);
d2 = vec(5);
d3 = vec(6);
diff1 = vec(7);
diff2 = vec(8);

clearance = vec(9);

n_mature = vec(10);

x0 = [100 60 70 0 40 30];
t = 4;

M = [b1 0 0 0 0 0;
    diff1 b2 0 0 0 0;
    0 diff2 b3 0 0 0;
    d1 0 0 -clearance 0 0;
    0 d2 0 0 -clearance 0;
    0 0 d3 0 0 -clearance];
for i=1:length(M)
    M(i,i) = 2*M(i,i)-sum(M(:,i));
end

if n_mature==1
    for i=1:length(t)
        n_t(i,:) = expm(M*t(i))*x0';
        percent(i,:) = (n_t(i,1:3)./(n_t(i,1:3)+n_t(i,4:6)))*100;
        out = [(n_t(i,1)+n_t(i,2))/(n_t(i,1)+n_t(i,4)+n_t(i,2)+n_t(i,5)) (n_t(i,3))/(n_t(i,3)+n_t(i,6))];
    end
elseif n_mature==2
    for i=1:length(t)
        n_t(i,:) = expm(M*t(i))*x0';
        percent(i,:) = (n_t(i,1:3)./(n_t(i,1:3)+n_t(i,4:6)))*100;
        %percent(i,:) = 100*[n_t(i,1)/(n_t(i,1)+n_t(i,4)) (n_t(i,3)+n_t(i,2))/(n_t(i,3)+n_t(i,6)+n_t(i,2)+n_t(i,5))];
        out = 100*[n_t(i,1)/(n_t(i,1)+n_t(i,4)) (n_t(i,3)+n_t(i,2))/(n_t(i,3)+n_t(i,6)+n_t(i,2)+n_t(i,5))];
    end
end

more_live = out(1)>out(2);

end