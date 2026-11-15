# output location
output_dir=$EPHEMERAL/outputs_cp/$1    # run as: bash CheckForFailedJobs.sh Y6_20000ps_5A

# go into the folder from which the job was submitted, and get an array of input files
cd $output_dir


for file in $(find -type f -name "*.log"); do
if ! tail $file -n 1 | grep -q "Normal termination"; then
 echo $file >> FailedJob.txt
fi 
done 