#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -J 2-51
#PBS -R eo

T=$((20000 + 600*($PBS_ARRAY_INDEX - 1)))
EXAMPLE_NAME="Y6-MD-PBC"
OUTPUT_DIR="$PBS_O_WORKDIR/pyMOO_results_${EXAMPLE_NAME}_${T}ps"
mkdir -p "$OUTPUT_DIR"

cp -R $HOME/pyMOO-main $TMPDIR

eval "$(/rds/general/user/ljh3218/home/miniforge3/bin/conda shell.bash hook)"
conda activate general

WORKDIR="$TMPDIR/pyMOO-main/examples/Y6-MD"
cd $WORKDIR
python GenxyzFiles_PBC.py $T $HOME/MD/Y6Crystal_300K_grofiles
cp init.file $OUTPUT_DIR
cp fin.file $OUTPUT_DIR
cp dist.txt $OUTPUT_DIR
cp finlattice.file $OUTPUT_DIR
#Read in initial and final molecules in the pair
mapfile -t init_array < "init.file"
mapfile -t fin_array < "fin.file"
mapfile -t fin_lattice_ar < "finlattice.file"
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
	cp zindo_MOs_Mol${init_mol}_Lattice0 $WORKDIR/MOsPBC/
	cd $WORKDIR
	
	cat zindo_MOs_header.gjf MoleculesPBC/Mol${fin_mol}_Lattice${fin_lattice}.xyz > zindo_job.gjf	
	echo >> zindo_job.gjf	
	
	cp zindo_job.gjf $TMPDIR/
	cd $TMPDIR
	g16 zindo_job.gjf
	mv fort.7 zindo_MOs_Mol${fin_mol}_Lattice${fin_lattice}
	cp zindo_MOs_Mol${fin_mol}_Lattice${fin_lattice} $WORKDIR/MOsPBC/
	cd $WORKDIR
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


cd $TMPDIR/pyMOO-main/src
mkdir -p ../bin
make clean
make all

cd $WORKDIR

for i in $(seq 0 1 $NPAIRS); do 
	init_mol=${init_array[$i]}
	fin_mol=${fin_array[$i]}
	fin_lattice=${fin_lattice_ar[$i]}
	OUTPUT_NAME="Mol${init_mol}_Mol${fin_mol}"
	echo $OUTPUT_NAME
	../../bin/abject MoleculesPBC/Mol${init_mol}_Lattice0.xyz MOsPBC/zindo_MOs_Mol${init_mol}_Lattice0 MoleculesPBC/Mol${fin_mol}_Lattice${fin_lattice}.xyz MOsPBC/zindo_MOs_Mol${fin_mol}_Lattice${fin_lattice} MoleculesPBC/$OUTPUT_NAME.xyz > "$OUTPUT_DIR/$OUTPUT_NAME.txt"
done




