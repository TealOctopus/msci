# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 14:46:36 2024

@author: ljh3218
"""

from CounterpoiseJarvVersion import CalcJ
import os 
import numpy as np
import re
import sys 

args = sys.argv   
cwd = os.getcwd()
mol = args[1]
PBC = int(args[2])

data_file = os.path.join(os.environ['EPHEMERAL'], 'outputs_cp', 'JobFiles')
subfiles = sorted(d for d in os.listdir(data_file) if os.path.isdir(os.path.join(data_file, d)))

if PBC == 0:
    results = np.zeros((len(subfiles), 4))
elif PBC == 1:
    results = np.zeros((len(subfiles), 5))

for i,file in enumerate(subfiles):
    mols = [int(s) for s in re.findall(r'\d+', file)]
    if PBC == 0:
        results[i,0:2] = mols 
        results[i,2:] = CalcJ(data_file, file)
    elif PBC == 1:
        results[i,0:2] = [mols[0], mols[2]]
        results[i,2] = mols[3]
        results[i,3:] = CalcJ(data_file, file)        
    

if PBC:
    np.savetxt('Ts_{}.csv'.format(mol), results, delimiter = ',', header = 'initial, final, final_lattice, Th, Te')
else:
    np.savetxt('Ts_{}.csv'.format(mol), results, delimiter = ',', header = 'initial, final, Th, Te')
