# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 14:11:01 2024

@author: ljh3218
"""
import numpy as np
from cclib.parser import ccopen

def CalcJ(root, namejob, verbose = 0):
    molA_parser=ccopen("{0}/{1}/{1}_part1.log".format(root,namejob))
    molB_parser=ccopen("{0}/{1}/{1}_part2.log".format(root, namejob))
    molAB_parser=ccopen("{0}/{1}/{1}_pair.log".format(root,namejob))

    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()
    print ("Parsed {}".format(namejob))

    if molA.nbasis!=molB.nbasis:
        print("Count of basis functions doesn't match. Failing.")
        return False
    
    nhomoA = molA.homos
    nhomoB = molB.homos
    nhomoAB = molAB.homos
    
    nlumoA = nhomoA + 1
    nlumoB = nhomoB + 1
    
    #Get molecular orbitals
    MOsA = molA.mocoeffs[0]
    MOsB = molB.mocoeffs[0]
    MOsAB = molAB.mocoeffs[0]
    
    #Get eigenvalues of pai
    EvalsAB = molAB.moenergies[0]
    
    #Calculate the molecular orbitals of A and B in the AB basis set
    SAB = molAB.aooverlaps
    MolAB_Pro = (np.dot(MOsAB, SAB)).T
    PsiA_AB_BS = np.dot(MOsA, MolAB_Pro)
    PsiB_AB_BS = np.dot(MOsB, MolAB_Pro)
    
    #Calculate the matrix of transfer integrals
    JAB=np.dot(np.dot(PsiA_AB_BS, np.diagflat(EvalsAB)), PsiB_AB_BS.T)
    JAA=np.dot(np.dot(PsiA_AB_BS, np.diagflat(EvalsAB)), PsiA_AB_BS.T)
    JBB=np.dot(np.dot(PsiB_AB_BS, np.diagflat(EvalsAB)), PsiB_AB_BS.T)
    S = np.dot(PsiA_AB_BS, PsiB_AB_BS.T)
    
    #Symmetric Lowdin transformation
    #c.f. Baumeier, Kirkpatrick and Adrienko, Phys. Chem. Chem. Phys, 2010, 12, 11103-11113 eqn 10
    J_eff = (JAB - 0.5*(JAA + JBB)*S)/(1-S**2)

    j0h = J_eff[nhomoA, nhomoB][0]
    j0e = J_eff[nlumoA, nlumoB][0]
    
    if verbose == 1:
    
         print('Minimum gap of pair in isolation (eV): {}'.format(min([molB.moenergies[0][nhomoB + 1] - molA.moenergies[0][nhomoA],
                                                           molA.moenergies[0][nhomoA + 1] - molB.moenergies[0][nhomoB],
                                                           molA.moenergies[0][nhomoA + 1] - molA.moenergies[0][nhomoA],
                                                           molB.moenergies[0][nhomoB + 1] - molB.moenergies[0][nhomoB]])),
               '\n',
               'Gap in pair (eV): {}'.format(EvalsAB[nhomoAB + 1] - EvalsAB[nhomoAB]))

    return j0h, j0e