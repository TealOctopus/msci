"""Look in more detail at the CP coupling for one dimer.

Prints:
  1. frontier orbital energies (HOMO-2 ... LUMO+2) of molecule A, molecule B and the pair
  2. the coupling matrix J_eff between A's and B's frontier orbitals (meV)
  3. the overlap S between them
The CalcJ value is the HOMO-HOMO / LUMO-LUMO element of (2).

Usage:  python DiagnoseOneJ.py <root folder> <dimer folder name>
Output is printed AND saved to <root>/<dimer>/J_diagnosis.txt
"""
import os
import sys
import numpy as np
from cclib.io import ccopen


# Names for printing, in the same order as idxA / idxB: H-2, H-1, HOMO, LUMO, L+1, L+2
def make_labels(K):
    labels = []
    for i in range(2 * K):
        if i < K - 1:
            labels.append(f"H-{K-1-i}")      # orbitals below the HOMO
        elif i == K - 1:
            labels.append("HOMO")
        elif i == K:
            labels.append("LUMO")
        else:
            labels.append(f"L+{i-K}")        # orbitals above the LUMO
    return labels

def diagnose(root, dimer, K=3):
    """Return a dictionary of arrays describing the frontier orbitals and couplings of one dimer."""

    # K=  how many orbitals either side of the gap to show
    #ie take K 'below gap' (ie Homo, homo-1... homo- (k-1))
    #and k above (lumo, lumo +1, ..lumo +(k-1))
    d = os.path.join(root, dimer)
    A  = ccopen(f"{d}/{dimer}_part1.log").parse()
    B  = ccopen(f"{d}/{dimer}_part2.log").parse()
    AB = ccopen(f"{d}/{dimer}_pair.log").parse()
    #ie here we read in log files for the first run

    """
    Gaussian lists the MOs in order of increasing energy
    homos is the the index of the highest occupied orbital
    - homos is an array, with one entry per spin. We have spin 0, corrresponding to
    homos[0]
    so homos[0]+1 is entry for LUMO
    """
    hA, hB, hAB = A.homos[0], B.homos[0], AB.homos[0]
    #the index should be the same for hA and hB, since they are same molecule

    idxA = np.arange(hA - K + 1, hA + K + 1)      # indexes for HOMO-2 ... LUMO+2
    idxB = np.arange(hB - K + 1, hB + K + 1)
    idxAB = np.arange(hAB - K + 1, hAB + K + 1)

    # Same maths as CalcJ, but keep the whole frontier block
    MOsA, MOsB, MOsAB = A.mocoeffs[0], B.mocoeffs[0], AB.mocoeffs[0]
    #e.g. MOsA is the matric C. A column 
    # for each MO psi_i: phi_i = C_ia x_a, where x_a are the basis functions
    #basis functions not necessarily orthogonal
    #so each column holds how many of each basis function make up one MO

    E = AB.moenergies[0]           # dimers orbital energy, eV
    proj = (MOsAB @ AB.aooverlaps).T 
    #aooverlaps is the overlap matrix S of the basis functions
    #S_ab = <x_a|x_b>
    PA = MOsA @ proj
    PB = MOsB @ proj
    #essentially, rewrites e.g As orbitals in the pairs basis
    JAB = PA @ np.diag(E) @ PB.T
    JAA = PA @ np.diag(E) @ PA.T
    JBB = PB @ np.diag(E) @ PB.T
    #see MOO notes for explation of maths
    #Jij = <phi_i^A|H|phi_j^B>

    S_MO = PA @ PB.T #overlap of the AOs 
    eA = np.diag(JAA)[:, None]     # site energies of A's orbitals (column)
    eB = np.diag(JBB)[None, :]     # site energies of B's orbitals (row))
    Jeff = (JAB - 0.5 * (eA + eB) * S_MO) / (1 - S_MO**2)
    # note: uses the diagonal site energies e_i of A and e_j of B for each element (i,j)
    #apparently finds a symmetric lowdin orthogonalisation (?)

    #  where do A's and B's frontier orbitals end up among the PAIR's orbitals?
    # PA[i, k] = <phi_i^A | psi_k^AB>, so PA[i, k]**2 is the fraction of A's orbital i found in pair orbital k.
    # Summed over ALL pair orbitals k this should be 1 (the pair orbitals form a complete set).
    wA_H = PA[hA, :]**2          # A's HOMO
    wB_H = PB[hB, :]**2          # B's HOMO
    wA_L = PA[hA + 1, :]**2      # A's LUMO
    wB_L = PB[hB + 1, :]**2      # B's LUMO

    # Window of pair orbitals to store: HOMO-4 ... LUMO+4
    W = 5
    idxW = np.arange(hAB - W + 1, hAB + W + 1)

    # Site energies inside the pair (diagonal of JAA, JBB), in eV
    siteA = np.diag(JAA)
    siteB = np.diag(JBB)

    # Keep only the frontier block
    block = np.ix_(idxA, idxB)
    return {
        "labels": np.array(make_labels(K)),
        "E_A":    A.moenergies[0][idxA],
        "E_B":    B.moenergies[0][idxB],
        "E_AB":   E[idxAB],
        "J_eff":  1000 * Jeff[block],             # meV
        "S_MO":   S_MO[block],
        "nbasis": np.array([A.nbasis, B.nbasis, AB.nbasis]),
        "nmo":    np.array([A.nmo, B.nmo, AB.nmo]),
        "pair_window":  np.arange(-W + 1, W + 1),       # pair orbital positions: -4..0 = H-4..HOMO, 1..5 = LUMO..L+4
        "wA_HOMO":      wA_H[idxW],
        "wB_HOMO":      wB_H[idxW],
        "wA_LUMO":      wA_L[idxW],
        "wB_LUMO":      wB_L[idxW],
        "wA_HOMO_total": wA_H.sum(),
        "wB_HOMO_total": wB_H.sum(),
        "wA_LUMO_total": wA_L.sum(),
        "wB_LUMO_total": wB_L.sum(),
        "site_HOMO":    np.array([siteA[hA], siteB[hB]]),          # eV
        "site_LUMO":    np.array([siteA[hA + 1], siteB[hB + 1]]),  # eV
    }


if __name__ == "__main__":
    root = sys.argv[1]
    for dimer in sys.argv[2:]:
        result = diagnose(root, dimer)
        out_file = os.path.join(root, dimer, "J_diagnosis.npz")
        #.npz is numpys format for storing several named arrays in one file
        np.savez(out_file, **result)
        print(f"Saved {out_file}", flush=True)
