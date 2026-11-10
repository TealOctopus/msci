# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:28:56 2026

@author: Lucy
"""

import numpy as np
import pandas as pd
import sys

def make_pair(Mol1, Mol2, FinLattice):
    first = pd.read_csv('MoleculesPBC/Mol{}_Lattice0.xyz'.format(Mol1), 
                       sep = '\t', names = ['Symbol', 'X', 'Y', 'Z'],
                       header = None)
    second = pd.read_csv('MoleculesPBC/Mol{}_Lattice{}.xyz'.format(Mol2, FinLattice), 
                        sep = '\t', names = ['Symbol', 'X', 'Y', 'Z'],
                        header = None)
    pair = pd.concat((first, second), axis = 0,
                         ignore_index = True)
    pair.to_csv('MoleculesPBC/Mol{}_Mol{}.xyz'.format(Mol1, Mol2), 
               sep = '\t',
               float_format = '%.7f',
               index = False,
               lineterminator = '\n',
               header = False)
    
if __name__ == "__main__":
    #sys.argv[0] is the script name
    args = sys.argv
    make_pair(args[1], args[2], args[3])