def no_limit_pairs(df, molecules):
    init,fin,fin_lattice,dist_arr = [],[],[],[]
    for i in range(len(molecules)):
        mol1 = df[(df['Mol'] == molecules[i]) & (df['Lattice'] == 0)].copy()      
        for j in range(i+1,len(molecules)):
            
            mol_j_df = (df[df['Mol'] == molecules[j]])
            lattice_nums = mol_j_df["Lattice"].unique()
            seps = []
            for lattice in lattice_nums:
                mol_check = mol_j_df[mol_j_df["Lattice"] == lattice]
                #ie only keeps entries of mol_j_df where the lattice number is equal to the current lattice number
                #picks out all the atoms
                #sep = separation
                sep = calculate_min_distance_between_C_atoms(mol1,mol_check)
                seps.append(sep)

            closest_lattice = lattice_nums[np.argmin(seps)]
            #picks up the lattice number corresponding to the minimum separation

            mol2 = mol_j_df[mol_j_df['Lattice'] == closest_lattice].copy()


            init.append(molecules[i])
            fin.append(molecules[j])
            fin_lattice.append(closest_lattice)

            dist = min(seps)
            dist_arr.append(np.round(dist,6))
    pairs = pd.DataFrame({'init' : init, 'fin' : fin,
                        'fin_lattice' : fin_lattice,
                        'dist' : dist_arr})
    return pairs

def calculate_min_distance_between_C_atoms(molecule1, molecule2):
    C_atoms_mol1 = molecule1[molecule1['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    C_atoms_mol2 = molecule2[molecule2['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    pairwise_distances = np.linalg.norm(C_atoms_mol1[:, None] - C_atoms_mol2, axis=-1)
    min_distance = np.min(pairwise_distances)
    return min_distance

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

    #old, slower method which concatinates files