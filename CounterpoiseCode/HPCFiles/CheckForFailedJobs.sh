# output location
output_dir=$EPHEMERAL/"outputs_cp"/_job_files_Ts

# go into the folder from which the job was submitted, and get an array of input files
cd $output_dir


for file in $(find -type f -name "*.log"); do
if ! tail $file -n 1 | grep -q "Normal termination"; then
 echo $file >> FailedJob.log
fi 
done 