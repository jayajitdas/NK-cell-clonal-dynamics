# This is a script to submit bootstrap sample repetitions
# Input arguments to lmfit_routines_bootstrap.py are for
# the correlation constraint and the repetition number (for
# naming output files)

# Set a correlation constraint
corr=-0.2

# directory work, adjust as necessary
cd asym_trf
mkdir bootstrap
cd ../

cd three_stage_2_mature_death
mkdir bootstrap
cd ../

cd three_stage_death
mkdir bootstrap
cd ../

# for 500 repetitions
for i in {1..500}
do
cd asym_trf/bootstrap
# 'nice' argument sets priority.  I was less interested in asym_division,
# so I set it at a lower priority
sbatch --time=10-0 --nice=2000+$i ../lmfit_routines_bootstrap.py $corr $i
sleep 1
cd ../../

cd three_stage_2_mature_death/bootstrap
# I wanted this and the next set to complete in parallel rather than
# the next to wait for these to complete, thus nic=1000+i
sbatch --time=10-0 --nice=1000+$i ../lmfit_routines_bootstrap.py $corr $i
sleep 1
cd ../../

cd three_stage_death/bootstrap
sbatch --time=10-0 --nice=1000+$i ../lmfit_routines_bootstrap.py $corr $i
sleep 1
cd ../../

done

