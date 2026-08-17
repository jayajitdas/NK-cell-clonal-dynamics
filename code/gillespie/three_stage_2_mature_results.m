function three_stage_2_mature_results()

bI=1.064;
bINT=1.613;
bM=0.700;
dI=0.054;
dINT=0.285;
dM=1.370;
r1=0.126;
r2=0.164;

growth_M = [bI 0 0
    r1 bINT 0
    0 r2 bM];

for i=1:10000
cells = gillespie_with_birth_and_death([1,0,0],growth_M,[dI,dINT,dM],8);
clone_size(i) = sum(cells);
cd27_pos(i) = cells(1)/clone_size(i)*100;
ly6c_pos(i) = cells(3)/clone_size(i)*100;
end

end