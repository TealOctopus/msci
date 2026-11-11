#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=12GB

T=20000
EXAMPLE_NAME="Y6-MD-PBC"
OUT_CSV="$PBS_O_WORKDIR/${EXAMPLE_NAME}-Ts-MOO-${T}ps.csv"

eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate msci

cd $PBS_O_WORKDIR #this is the directory you are working in when you type qsub
python GenxyzFiles_PBC.py $T $HOME/MSci/CounterpoiseCode/Structures/Y6 #find the pairs
#Read in initial and final molecules in the pair
mapfile -t init_array < "init.file"
mapfile -t fin_array < "fin.file"
mapfile -t fin_lattice_ar < "finlattice.file"
mapfile -t dist_array < "dist.txt"
#Get length of the array to know how many pairs there are
NPAIRS=${#init_array[@]}
NPAIRS=$(($NPAIRS - 1))
# select one of the input files based on current array job index
for i in $(seq 0 1 $NPAIRS); do 
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}
	python MakePair_PBC.py $init_mol $fin_mol $fin_lattice 
	#write the pairs, ie combine .xyz files into Mol{i}_Mol{j}.xyz
done

module purge
module load Gaussian

# Start with an empty MOsPBC folder so no orbitals from an earlier run can be used
# mv the other one to old so it isnt just written over
rm -rf MOsPBC_old
[ -d MOsPBC ] && mv MOsPBC MOsPBC_old
mkdir -p MOsPBC

run_zindo () { #define a function
	local name=$1 #$1 means the first argument of run_zindo, cant have any spaces
	#local sets it as a local, rather than global, variable

	# Skip this molecule if we already have an MO for it
	if [ -f "MOsPBC/zindo_MOs_${name}" ]; then
		return
	fi

	cat zindo_MOs_header.gjf MoleculesPBC/${name}.xyz > zindo_job.gjf
	#combine 'header' file with the molecule file of init_mol to make a zindo input file
	#this is because gaussian needs a specific structure with details in this header
	#header contains:
	#Resources required
	# zindo punch=mo nosymm
	#This means write out (ie 'punch') the MOs to the fort.7 file (where fortran writes to by default) (punch=mo) and do not reorientate the molecule (nosymm)
	#Then charge (0) and spin multiplicity (1) are specified
	#charge 0 because we just want the MOs of the neutral molecule
	#1 just means all electrons are paired, so we have a singlet

	echo >> zindo_job.gjf #adds a blank line to the end of the file, which is required for gaussian to run properly

	cp zindo_job.gjf $TMPDIR/ #copy input to tmpdir, a fast dir to work from
	cd $TMPDIR
	rm -f fort.7 #remove any fort.7 files that may be there from a previous run

	# Only keep the result if Gaussian finished without an error AND fort.7 exists and is not empty
	if g16 zindo_job.gjf && [ -s fort.7 ]; then #this line runs the gaussian (g16), and the if also then tests it the command finished successfully
		mv fort.7 zindo_MOs_${name}
		cp zindo_MOs_${name} $PBS_O_WORKDIR/MOsPBC/
	else
		echo "GAUSSIAN FAILED for ${name}" >&2 #>&2 means it goes to the error log, not the output log
		cp zindo_job.log $PBS_O_WORKDIR/FAILED_${name}.log 2>/dev/null #saves the log file of the failed gaussianfor debugging,
		#2>/dev/null means if the file doesn't exist, or there is an error, it will not print anything 
	fi

	cd $PBS_O_WORKDIR
}


for i in $(seq 0 1 $NPAIRS); do 
	run_zindo Mol${init_array[$i]}_Lattice0
	run_zindo Mol${fin_array[$i]}_Lattice${fin_lattice_ar[$i]}

done 

for file in ./MOsPBC/*; do #these are the zindo output files, but also contains other info we don't need
    if [ -f "$file" ]; then
		num=$(grep -n -m 1 "(5D15.8)" $file | cut -d: -f1) 
		#5D15.8 is a tag that marks the start of the coefficient block 
		num2=$(($num - 1))
        sed -i "1,${num2}d" "$file"
		#remove everything above the 5D15.8 line, we dont need it
    fi
done

module purge #clears the modules
module load "tools/prod"
module load "GCC/11.3.0"


cd $HOME/MSci/MOO/src
mkdir -p ../bin
make clean 
make all

cd $PBS_O_WORKDIR

echo "init,fin,fin_lattice,dist,Th_MOO,Te_MOO" > "$OUT_CSV"
#writes the CSV header line to the output file

for i in $(seq 0 1 $NPAIRS); do
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}

	A=Mol${init_mol}_Lattice0
	B=Mol${fin_mol}_Lattice${fin_lattice}
	#names the molecules

	result=$(../bin/abject MoleculesPBC/$A.xyz MOsPBC/zindo_MOs_$A \
	         MoleculesPBC/$B.xyz MOsPBC/zindo_MOs_$B \
	         MoleculesPBC/Mol${init_mol}_Mol${fin_mol}.xyz)
	#this line runs abject with 5 inputs: A's geometry and orbitals
	#B's geometry and orbitals, and the dimer geometry

	#abject is the actual 'MOO' part

	read -r Jh Je extra <<< "$result"
	# splits the output data from abject into Jh, Je and 'extra'
	# extra is generally discarded, unless the calc failed

	if [[ "$Jh" =~ ^-?[0-9] && -n "$Je" ]]; then
		echo "$init_mol,$fin_mol,$fin_lattice,${dist_array[$i]},$Jh,$Je" >> "$OUT_CSV"
		# write to the output CSV file if the calculation was successful
	else
		echo "FAILED Mol${init_mol}_Mol${fin_mol}: $result" >&2
	fi
done

