# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:43:37 2026

@author: Lucy
"""

import pandas as pd
import numpy as np

for t in range(20600,50001,600):
    init = np.loadtxt('./pyMOO_results_Y6-MD-PBC_{}ps/init.file'.format(t))
    fin = np.loadtxt('./pyMOO_results_Y6-MD-PBC_{}ps/fin.file'.format(t))
    seps = np.loadtxt('./pyMOO_results_Y6-MD-PBC_{}ps/dist.txt'.format(t))

    Te_MOO = []
    Th_MOO = []
    error = []

    for i in range(len(init)):
        try:
            values = np.loadtxt('pyMOO_results_Y6-MD-PBC_{}ps/Mol{:.0f}_Mol{:.0f}.txt'.format(t,init[i],
                                                                         fin[i]))
            Te_MOO.append(values[1])
            Th_MOO.append(values[0])
        except ValueError:
            print('Error for pair Mol {:.0f} and Mol {:.0f} at time {} ps'.format(init[i],fin[i]),t)
            error.append(i)

    init_pd = np.delete(init, error)
    fin_pd = np.delete(fin, error)
    seps_pd = np.delete(seps, error)

    Ts = pd.DataFrame({'init':init_pd, 'fin':fin_pd, #'dist':seps_pd,
                   'Th_MOO': Th_MOO, 'Te_MOO':Te_MOO})
    Ts.to_csv('CouplingsPBC/Y6-PBC-Ts-MOO-{}ps.csv'.format(t))
    