The python script you need to run is called MakeCPFiles.py
Input structures can be found in the Structures file - I have given five snapshots for two different molecules (Y6 and L8BO)
	Note, you will need to edit line 178 of MakeCPFiles.py to point to the location where you save the structures 
Options:
makeCPfiles: setting this to 1 will mean the code saves the counterpoise jobfiles in savefile_CP
PBC: whether or not to use periodic boundary conditions. keep this set to 1 (legacy)
save_csv: set this to 1 to save the pairs and atoms to csv files you can look at/use later (instead of having to re-calculate them every time)
big_batch: probably more useful later, but it is a way of randomly selecting a subset of the total number of couplings for which to generate job files
Nsamples: the number of samples which are randomly selected per snapshot if you set big_batch = 1
max_dist: how close (in A) two carbon atoms have to be for a pair of molecules to count as neighbours (only generate counterpoise job files for neighbours)
which: either 'Y6' or 'L8BO'

To start, I would suggest looking at one snapshot and setting dist to 5/6 angstroms

The folder HPCFiles should be copied to the HPC. It contains all the code you need to run the counterpoise calculations
The job to submit first is batch_cp.pbs
	You will need to change the line #PBS -J 1-24 so that the final number is the total number of gaussian jobs you need to run (so, if you have 6 pairs, this would be 18 as each pair needs 3 gaussian calculations)
	You will need to change the submission_folder to wherever you saved the counterpoise job files 
The outputs will be saved to ephemeral - note that files left here will be deleted after 30 days
One the Gaussian calculations have finished, you can use CheckForFailedJobs.sh to see if any need to be run again
	You will need to change the output_dir variable to wherever your gaussian output files are
When all Gaussian calculations are finished, you can use extract_Ts_batch to get the couplings
	This file calls CalcJs.py - again, you will need to change the data_file variable to point to wherever your gaussian output files are
