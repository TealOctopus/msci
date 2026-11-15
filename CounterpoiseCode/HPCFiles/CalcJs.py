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
run_name = args[1] # e.g. Y6_20000ps_5A


data_file = os.path.join(os.environ['EPHEMERAL'], 'outputs_cp', run_name)
subfiles = sorted(d for d in os.listdir(data_file) if os.path.isdir(os.path.join(data_file, d)))

rows = []
  for file in subfiles:
      mols = [int(s) for s in re.findall(r'\d+', file)]   # [init, 0, fin, fin_lattice, t]
      try:
          Th, Te = CalcJ(data_file, file) #prints the error, but keeps running
      except Exception as e:
          print(f"FAILED {file}: {e}")
          continue
      rows.append([mols[0], mols[2], mols[3], mols[4], Th, Te])

  out = os.path.expanduser(f'~/MSci/runs/{run_name}/Ts_CP.csv') #Ts means transfer integrals
  np.savetxt(out, np.array(rows), delimiter=',', comments='',
             header='init,fin,fin_lattice,t_ps,Th_CP,Te_CP')
