from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from modules.xyzutils import build_xyz_file

def convert_coordinates_to_mols(elements, coordinates):
    xyz = build_xyz_file(elements, coordinates)
    mol = Chem.rdmolfiles.MolFromXYZBlock(xyz)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol
