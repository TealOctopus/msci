#!/bin/bash
#PBS -N MOO_one
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=12GB
# Re-run MOO for ONE dimer. New ZINDO orbitals go in MOsPBC_rerun/, so MOsPBC/ and Ts_MOO.csv are untouched.
# It also re-runs abject on the ORIGINAL orbitals, so you can see whether the orbitals or the
# coupling step is what changed.
#
# Submit e.g.:
#   qsub -v MOL=Y6,T=20000,CUTOFF=5,INIT=12,FIN=40,LAT=3 RunOnePair_MOO.sh

RUN_DIR="$HOME/MSci/runs/${MOL}_${T}ps_${CUTOFF}A"
cd "$RUN_DIR"
mkdir -p MOsPBC_rerun

A=Mol${INIT}_Lattice0
B=Mol${FIN}_Lattice${LAT}
DIMER_XYZ=MoleculesPBC/Mol${INIT}_Mol${FIN}_Lattice${LAT}.xyz

module purge
module load Gaussian
for name in $A $B; do
    cat $HOME/MSci/MOO/zindo_MOs_header.gjf MoleculesPBC/${name}.xyz > $TMPDIR/zindo_job.gjf
    echo >> $TMPDIR/zindo_job.gjf
    cd $TMPDIR
    rm -f fort.7
    if g16 zindo_job.gjf && [ -s fort.7 ]; then
        num=$(grep -n -m 1 "(5D15.8)" fort.7 | cut -d: -f1)
        sed -i "1,$((num - 1))d" fort.7          # strip header, same as the main script
        cp fort.7 "$RUN_DIR/MOsPBC_rerun/zindo_MOs_${name}"
    else
        echo "GAUSSIAN FAILED for ${name}"
        cp zindo_job.log "$RUN_DIR/FAILED_rerun_${name}.log"
    fi
    cd "$RUN_DIR"
done

module purge
module load "tools/prod"
module load "GCC/11.3.0"

echo "Pair $A + $B"
echo "Original orbitals (Jh Je): $($HOME/MSci/MOO/bin/abject MoleculesPBC/$A.xyz MOsPBC/zindo_MOs_$A MoleculesPBC/$B.xyz MOsPBC/zindo_MOs_$B $DIMER_XYZ)"
echo "New orbitals      (Jh Je): $($HOME/MSci/MOO/bin/abject MoleculesPBC/$A.xyz MOsPBC_rerun/zindo_MOs_$A MoleculesPBC/$B.xyz MOsPBC_rerun/zindo_MOs_$B $DIMER_XYZ)"
echo "Ts_MOO.csv line:"
grep "^${INIT},${FIN},${LAT}," Ts_MOO.csv
