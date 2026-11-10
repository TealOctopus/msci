# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:43:37 2026

@author: Lucy
"""

import pandas as pd
import numpy as np

init = np.loadtxt('init.file')
fin = np.loadtxt('fin.file')
seps = np.loadtxt('dist.txt')
#Ts = pd.read_csv('ITIC-Ts.csv', header = 0)
Te_MOO = []
Th_MOO = []
error = []

for i in range(len(init)):
    try:
        values = np.loadtxt('pyMOO_results_Y6-MD-PBC_20000ps/Mol{:.0f}_Mol{:.0f}.txt'.format(init[i],
                                                                         fin[i]))
        Te_MOO.append(values[1])
        Th_MOO.append(values[0])
    except ValueError:
        print('Error for pair Mol {:.0f} and Mol {:.0f}'.format(init[i],fin[i]))
        error.append(i)

init_pd = np.delete(init, error)
fin_pd = np.delete(fin, error)
seps_pd = np.delete(seps, error)

Ts = pd.DataFrame({'init':init_pd, 'fin':fin_pd, 'dist':seps_pd,
                   'Th_MOO': Th_MOO, 'Te_MOO':Te_MOO})
Ts.to_csv('Y6-PBC-Ts-MOO-20000ps.csv')
    
