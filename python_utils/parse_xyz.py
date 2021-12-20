import re
import numpy as np


PERIODIC_TABLE = ["Bq","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo","X"]

def parse_xyz_ensemble(ens_file):
    with open(ens_file) as f:
        ensemble_data = f.readlines()
        n_atoms = int(ensemble_data[0])
        n_molecules = 0

        for line in ensemble_data:
            if re.search(f'^[\s\t]*{n_atoms}', line):
                n_molecules += 1

        molecular_ens = dict()
        for mol_id in range(n_molecules):
            molecule_data = dict()

            mol_end = ((n_atoms + 1) * (mol_id + 1)) + mol_id
            mol_start = mol_end - n_atoms - 1  
            #print(f"Run: {mol_id}, Start: {mol_start+1}, End; {mol_end+1}")

            ensemble_slice = ensemble_data[mol_start:mol_end+1]

            cartesian = ensemble_slice[2:]
            cartesian =  [i.strip() for i in cartesian]
            elements = list()
            coords = list()

            for i in cartesian:
                elements.append(i.split()[0])
                coords.append(i.split()[1:])

            atomic_numbers = [PERIODIC_TABLE.index(element) for element in elements]

            molecule_data['n_atoms'] = ensemble_slice[0].strip()
            molecule_data['header'] = ensemble_slice[1].strip()
            molecule_data['elements'] = np.array(elements)
            molecule_data['atomic_numbers'] = np.array(atomic_numbers)
            molecule_data['coordinates'] = np.array(coords)

            molecular_ens[mol_id] = molecule_data
    
    return molecular_ens

def write_xyz_file(elements,coordinates, header='', filename='molecule.xyz'):
    n_atoms = len(elements)
    with open(filename, 'a') as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{header}\n")
        for i in range(n_atoms):
            element = elements[i]
            coords = '\t'.join(coordinates[i])
            f.write(f'{element}\t{coords}\n')
