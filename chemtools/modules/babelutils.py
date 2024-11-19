
from openbabel import openbabel
from modules import xyzutils

def convert_coordinates_to_mols(elements, coordinates, header=''):
    xyz = xyzutils.build_xyz_file(elements, coordinates, header)
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("xyz")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz)
    return mol

def convert_ensemble_to_mols(ensemble):
    ens_elements = [_['elements'] for _ in ensemble.values()]
    ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
    ens_header = [_['header'] for _ in ensemble.values()]

    mols = []
    for args in zip(ens_elements,ens_coordinates,ens_header):
        mol = convert_coordinates_to_mols(*args)
        mols.append(mol)
    
    return mols
