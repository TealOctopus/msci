# 1. Generate the pairs + CP input files (quick, fine on the login node)
#    Doing this first means the MOO job reuses the same init.file etc.
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate msci
python ~/MSci/CounterpoiseCode/make_cp_jobs.py Y6 20000 5     # use "5", not "5.0" – it goes into the folder name

# 2. Submit both. Counterpoise (run the qsub line it printed, from HPCFiles/):
cd ~/MSci/CounterpoiseCode/HPCFiles
qsub -J 1-<3*Npairs> -v MOL=Y6,T=20000,CUTOFF=5 batch_cp.pbs
#    MOO:
cd ~/MSci/MOO
qsub -v MOL=Y6,T=20000 RunAllPairs_OneSnapshot.sh

# 3. When the CP jobs are done (qstat empty):
cd ~/MSci/CounterpoiseCode/HPCFiles
bash CheckForFailedJobs.sh Y6_20000ps_5A      # look at $EPHEMERAL/outputs_cp/Y6_20000ps_5A/FailedJob.txt
(cat $EPHEMERAL/outputs_cp/Y6_20000ps_5A/FailedJob.txt checks)
qsub -v MOL=Y6,T=20000,CUTOFF=5 extract_Ts_batch.pbs    # must be submitted from HPCFiles/

# 4. Compare
python ~/MSci/common/compare_Ts.py Y6_20000ps_5A