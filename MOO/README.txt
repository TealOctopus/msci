For MOO, you should just be able to run any of the shell scripts on the HPC and do the cacluations. You will need to adjust the following:
EXAMPLE_NAME=name of the file where all the .py fiels etc are located
OUTPUT_DIR=Where you want to save the files
WORKDIR=where you are running the job from

I would start by running the 'RunAllPairs_OneSnapshot.sh' file first to check enerything is working ok before you try and run a batch job. 

Note that you will need to run the PBC and non-PBC job scripts to get all the possible couplings

Once the calculations are finsihed, you can extract the couplings using ExtractTs.py