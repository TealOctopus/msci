#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=1:mem=12GB

T=20000
EXAMPLE_NAME="Y6-MD-PBC"
OUT_CSV="$PBS_O_WORKDIR/${EXAMPLE_NAME}-Ts-MOO-${T}ps.csv"

eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate msci

cd $PBS_O_WORKDIR
python GenxyzFiles_PBC.py $T $HOME/MSci/CounterpoiseCode/Structures/Y6
#Read in initial and final molecules in the pair
mapfile -t init_array < "init.file"
mapfile -t fin_array < "fin.file"
mapfile -t fin_lattice_ar < "finlattice.file"
mapfile -t dist_array < "dist.file"
#Get length of the array to know how many pairs there are
NPAIRS=${#init_array[@]}
NPAIRS=$(($NPAIRS - 1))
# select one of the input files based on current array job index
for i in $(seq 0 1 $NPAIRS); do 
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}
	python MakePair_PBC.py $init_mol $fin_mol $fin_lattice
done

module purge
module load Gaussian

for i in $(seq 0 1 $NPAIRS); do 
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}
	cat zindo_MOs_header.gjf MoleculesPBC/Mol${init_mol}_Lattice0.xyz > zindo_job.gjf	
	echo >> zindo_job.gjf
	
	cp zindo_job.gjf $TMPDIR/
	cd $TMPDIR
	g16 zindo_job.gjf
	mv fort.7 zindo_MOs_Mol${init_mol}_Lattice0
	cp zindo_MOs_Mol${init_mol}_Lattice0 $PBS_O_WORKDIR/MOsPBC/
	cd $PBS_O_WORKDIR
	
	cat zindo_MOs_header.gjf MoleculesPBC/Mol${fin_mol}_Lattice${fin_lattice}.xyz > zindo_job.gjf	
	echo >> zindo_job.gjf	
	
	cp zindo_job.gjf $TMPDIR/
	cd $TMPDIR
	g16 zindo_job.gjf
	mv fort.7 zindo_MOs_Mol${fin_mol}_Lattice${fin_lattice}
	cp zindo_MOs_Mol${fin_mol}_Lattice${fin_lattice} $PBS_O_WORKDIR/MOsPBC/
	cd $PBS_O_WORKDIR
done 

for file in ./MOsPBC/*; do
    if [ -f "$file" ]; then
		num=$(grep -n -m 1 "(5D15.8)" $file | cut -d: -f1)
		num2=$(($num - 1))
        sed -i "1,${num2}d" "$file"
    fi
done

module purge
module load "tools/prod"
module load "GCC/11.3.0"


cd $HOME/MSci/MOO/src
mkdir -p ../bin
make clean
make all

cd $PBS_O_WORKDIR

echo "init,fin,fin_lattice,dist,Th_MOO,Te_MOO" > "$OUT_CSV"

for i in $(seq 0 1 $NPAIRS); do
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}
	A=Mol${init_mol}_Lattice0
	B=Mol${fin_mol}_Lattice${fin_lattice}

	result=$(../bin/abject MoleculesPBC/$A.xyz MOsPBC/zindo_MOs_$A \
	         MoleculesPBC/$B.xyz MOsPBC/zindo_MOs_$B \
	         MoleculesPBC/Mol${init_mol}_Mol${fin_mol}.xyz)
	read -r Jh Je extra <<< "$result"

	if [[ "$Jh" =~ ^-?[0-9] && -n "$Je" ]]; then
		echo "$init_mol,$fin_mol,$fin_lattice,${dist_array[$i]},$Jh,$Je" >> "$OUT_CSV"
	else
		echo "FAILED Mol${init_mol}_Mol${fin_mol}: $result" >&2
	fi
done

