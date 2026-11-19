"""Calculate Th, Te for one dimer from its three CP logs.

Usage:  python CalcOneJ.py <root folder> <dimer folder name>
  e.g.  python CalcOneJ.py $EPHEMERAL/outputs_cp/Y6_20000ps_5A Mol12_Lattice0_Mol40_Lattice3_20000ps
Expects <root>/<dimer>/<dimer>_part1.log, _part2.log and _pair.log
The result is printed AND saved to <root>/<dimer>/J_result.txt
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from CounterpoiseJarvVersion import CalcJ

root, dimer = sys.argv[1], sys.argv[2]
result = CalcJ(root, dimer)
if result is False:
    line = f"{dimer}: CalcJ failed (basis function counts of part1 and part2 don't match)"
else:
    Th, Te = result
    line = f"{dimer}: Th_CP = {Th:.8f}  Te_CP = {Te:.8f}"
print(line, flush=True)

out_file = os.path.join(root, dimer, 'J_result.txt')
with open(out_file, 'w') as f:
    f.write(line + '\n')
print(f"Saved to {out_file}")
