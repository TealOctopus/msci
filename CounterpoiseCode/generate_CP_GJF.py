"""
  Purpose:  Generates Gaussian input files (.gjf/.com) from geometry files (.xyz) of 
            two molecules for a Counterpoise calculation.
            This is a rewrite of the bash scripts in https://github.com/Jenny-Nelson-Group/counterpoise-J

"""
import os
from inspect import cleandoc  # to ignore indents when printing to file

def write_CP_com(geom_molA, geom_molB, job_name, dir_out,
                 nprocshared=8,mem = '10GB',
                 method = 'b3lyp/6-31g*',extra_word=''):
    """generate gjf / com files for counterpoise calculation."""

    # Make directories
    try:
        os.makedirs(dir_out + job_name, exist_ok=True)
        dir = dir_out + job_name + '/'
    except OSError:
        print(f"Directories job_files_Ts/{job_name} can not be created")

    # Molecule A with ghost atoms of molecule B
    with open(dir + job_name + '_part1.com', "w") as text_file:
        print(cleandoc(
            f'''
            %chk=part1.chk
            %nprocshared={nprocshared}
            %mem={mem}
            #p SCF(Tight,Conver=10) Integral(Grid=UltraFine) IOp(6/7=3) {method} nosymm {extra_word}

            CP calculation part 1 - mol A with ghost atoms of mol B

            0 1'''), file=text_file)
        print(geom_molA.to_string(header=False, index=False), file=text_file)
        ghost_molB = geom_molB.copy() # Turn molecule B into a ghost
        ghost_molB['atom'] = ghost_molB['atom'] + '-Bq' # Add -Bq to each atom, marking it as a ghost
        print(ghost_molB.to_string(header=False, index=False), file=text_file, end='\n\n')

    # Molecule B with ghost atoms of molecule A
    with open(dir + job_name + '_part2.com', "w") as text_file:
        print(cleandoc(
            f'''
            %chk=part2.chk
            %nprocshared={nprocshared}
            %mem={mem}
            #p SCF(Tight,Conver=10) Integral(Grid=UltraFine) IOp(6/7=3) {method} nosymm {extra_word}

            CP calculation part 2 - mol B with ghost atoms of mol A

            0 1'''), file=text_file)
        ghost_molA = geom_molA.copy() # Turn molecule A into a ghost
        ghost_molA['atom'] = ghost_molA['atom'] + '-Bq' # Add -Bq to each atom, marking it as a ghost
        print(ghost_molA.to_string(header=False, index=False), file=text_file)
        print(geom_molB.to_string(header=False, index=False), file=text_file, end='\n\n')

    # Molecule AB - the pair/dimer
    with open(dir + job_name + '_pair.com', "w") as text_file:
        print(cleandoc(
            f'''
            %chk=dimer.chk
            %nprocshared={nprocshared}
            %mem={mem}
            #p SCF(Tight,Conver=10) Integral(Grid=UltraFine) IOp(6/7=3) {method} nosymm {extra_word}

            CP calculation part 3 - mol AB 

            0 1'''), file=text_file)
        print(geom_molA.to_string(header=False, index=False), file=text_file)
        print(geom_molB.to_string(header=False, index=False), file=text_file, end='\n\n') # add extra line break
        # The molecular pair / dimer calc has a Link1 (Gaussian restart) to separately perform IOp(3/33=1)
        print(cleandoc(
            f'''

            --Link1--
            %chk=dimer.chk
            %nprocshared={nprocshared}
            %mem={mem}
            #p geom(allcheck) guess(read,only) IOp(3/33=1) {method} nosymm {extra_word}

            '''), file=text_file)
            

    
