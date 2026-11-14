"""
molecule_dict.py

Molecule-specific data for trimming side chains to methyl groups
and reading GROMACS .gro files.

All atom indices are 0-based (GROMACS atom number - 1).

For each side chain:
    backbone_atom : backbone atom the chain is attached to
    methyl_C      : first chain carbon (becomes the CH3 carbon)
    H             : the two hydrogens already on methyl_C;
                    H[0] is used as the reference hydrogen in fix_hydrogens
    cap_C         : second chain carbon; kept, relabelled as H, and moved
                    to become the third methyl hydrogen
"""




MOLECULES = {
    "Y6": {
        "resname": "Y6",          
        "backbone_stop": 81,      # atoms 0..80 are kept
        "chains": [
            {"backbone_atom": 21, "methyl_C": 34, "H": (75, 76), "cap_C": 156},
            {"backbone_atom": 24, "methyl_C": 33, "H": (73, 74), "cap_C": 125},
            {"backbone_atom": 14, "methyl_C": 32, "H": (71, 72), "cap_C": 96},
            {"backbone_atom": 11, "methyl_C": 31, "H": (69, 70), "cap_C": 81},
        ],
    },
    "L8BO": {
        "resname": "L8B",      
        "backbone_stop": 81,      
        "chains": [
            {"backbone_atom": 21, "methyl_C": 34, "H": (75, 76), "cap_C": 83},
            {"backbone_atom": 24, "methyl_C": 33, "H": (73, 74), "cap_C": 112},
            {"backbone_atom": 14, "methyl_C": 32, "H": (71, 72), "cap_C": 82},
            {"backbone_atom": 11, "methyl_C": 31, "H": (69, 70), "cap_C": 81},
        ],
    },
}


def get_trim_params(which):
    """Return the lists the trimming script expects, in its original format."""
    mol = MOLECULES[which]
    chains = mol["chains"]

    # Refuse to run with incomplete data
    missing = [] if mol["backbone_stop"] is not None else ["backbone_stop"]
    missing += [f"chains[{i}].{k}" for i, ch in enumerate(chains)
                for k, v in ch.items() if v is None]
    if missing:
        raise ValueError(f"{which} is missing values in molecule_dict.py: {missing}")

    return {
        "backbone_stop":    mol["backbone_stop"],
        "carbon_to_keep":   [ch["cap_C"] for ch in chains],
        "atoms_pairs":      [[ch["backbone_atom"], ch["methyl_C"]] for ch in chains],
        "ref_hydrogen":     [ch["H"][0] for ch in chains],
        "hydrogen_triples": [[*ch["H"], ch["cap_C"]] for ch in chains],
    }