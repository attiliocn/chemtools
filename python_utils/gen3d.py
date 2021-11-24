import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

sys.stdout.reconfigure(line_buffering=True)

# arguments parsing #
#####################
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file')

args = parser.parse_args()

def generate_3D_coords(smiles_string):
    molecule = Chem.MolFromSmiles(smiles_string)
    molecule_with_H = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule_with_H)
    return molecule_with_H

with open(args.file) as file:
    file.readline() # ignore header in the first line
    dataset = file.readlines()

for entry in dataset:
    id = entry.split(',')[0]
    smiles_string = entry.split(',')[1].rstrip('\n')
    dotted = False

    if "." in smiles_string:
        dotted = True

    print(f"Processing molecule ID {id}")
    print(f"Entry ID: {id}\nEntry SMILES: {smiles_string}\nDotted: {dotted}")

    if dotted:
        dot_string = smiles_string.split('.')[1]
        smiles_string = smiles_string.split('.')[0]
        print(f"DOT string: {dot_string}")

        dot_molecule = generate_3D_coords(dot_string)
        Chem.MolToXYZFile(dot_molecule, f"cat_{id}_dot.xyz")

    molecule = generate_3D_coords(smiles_string)
    Chem.MolToXYZFile(molecule, f"cat_{id}.xyz")

    print(f"End of processing for molecule ID {id}\n")
