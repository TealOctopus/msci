

def write_xyz_file(df, mol, Latt):
    df1 = df[['Symbol','X','Y','Z']]
    fname = "MoleculesPBC/Mol{:.0f}_Lattice{:.0f}.xyz".format(mol, Latt)
    with open(fname, "w") as f:
        for _, row in df1.iterrows():
            symbol, x, y, z = row
            f.write(f"{symbol}\t{x:.7f}\t{y:.7f}\t{z:.7f}\n")
        f.write('\n')
        

    

if __name__ == '__main__':
    t = int(sys.argv[1])
    datafile = sys.argv[2]
    fname = "Y6_300K_{}ps.gro".format(t)
    file = os.path.join(datafile,fname)
    b_vectors = get_box_vectors(file)
    
    MOO_cutoff = 5
    
    backbone_stop = 81  
    #These will be different for L8BO and Y6
    #NB, have subtracted 1 versus GROMACS numbering as python indexes from 0
    carbon_to_keep = [156,125,96,81]
    atoms_pairs = [[21,34], [24,33], [14,32], [11,31]]
    ref_hydrogen = [75,73,71,69]
    hydrogen_triples = [[75,76,156],[73,74,125],[71,72,96],[69,70,81]]
    
    atoms_to_keep = np.arange(0, backbone_stop)
    atoms_to_keep = np.append(atoms_to_keep,carbon_to_keep)
    
    df = gro_to_atoms(file)
    count = 0
    for m in df['Mol'].unique():
        df1 = df[df['Mol']==m].iloc[atoms_to_keep].copy()
        df1.index = atoms_to_keep
        for j,c in enumerate(carbon_to_keep):
            df1.loc[c,'Symbol'] = "H"
            df1.loc[c,'weight'] = 1
            df1 = fix_hydrogens(df1,atoms_pairs[j],
                                    ref_hydrogen[j],hydrogen_triples[j])
        if count == 0:
            molecules = df1
        else:
            molecules = pd.concat([molecules,df1])
        count += 1
    bigmolecules, pairs = define_pairs_PBC(molecules, b_vectors)
            
    init = []
    fin = []
    seps = []
    fin_lattices = []
    for i, pair in pairs.iterrows():
        if pair.dist <= MOO_cutoff: 
            fin_lattice = pair.fin_lattice          
            molecule1 = bigmolecules[(bigmolecules['Mol'] == pair.init) & (bigmolecules['Lattice'] == 0)]
            write_xyz_file(molecule1,pair.init,0)
            molecule2 = bigmolecules[(bigmolecules['Mol'] == pair.fin) & (bigmolecules['Lattice'] == fin_lattice)]
            write_xyz_file(molecule2,pair.fin,fin_lattice)       
            fin_lattices.append(fin_lattice)
            seps.append(pair.dist)
            init.append(pair.init)
            fin.append(pair.fin)
    
    np.savetxt('init.file', np.asarray(init), fmt='%.0f')
    np.savetxt('fin.file', np.asarray(fin), fmt='%.0f')
    np.savetxt('dist.txt', np.asarray(seps), fmt='%.8f')
    np.savetxt('finlattice.file', np.asarray(fin_lattices), fmt='%.0f')







