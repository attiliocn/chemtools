import os
import numpy as np
import json

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

with open(os.path.join(__location__,'periodictable.json')) as f:
    periodic_table = json.load(f)

element_to_Z = periodic_table['elements']
Z_to_element = {v:k for k,v in element_to_Z.items()}

def read_xyz_file(filepath):
    with open(filepath) as f:
        f.readline()
        f.readline()
        data = [i.strip().split() for i in f.readlines()]
        elements = [i[0] for i in data]
        coordinates = [i[1:] for i in data]
    elements = np.array(elements)
    coordinates = np.array(coordinates, dtype=float)
    return elements, coordinates

def build_xyz_file(elements, coordinates, header=''):
    xyz_content = []
    xyz_content.append(f"{str(len(coordinates))}\n")
    xyz_content.append(f'{header}\n')

    elements = elements.reshape(-1,1)
    
    for el,co in zip(elements,coordinates):
        xyz_content.append("{: >3} {: >15.10f} {: >15.10f} {: >15.10f}\n".format(el[0], *co))
    xyz_content = ''.join(xyz_content)
    return xyz_content

def read_xyz_ensemble(filepath):
    with open(filepath) as f:
        ensemble_data = [l.strip() for l in f.readlines()]
        
    num_atoms = int(ensemble_data[0])
    num_molecules = int(len(ensemble_data) / (num_atoms + 2))

    parsed_ensemble = dict()
    for mol_id in range(num_molecules):
        molecule_data = dict()
        start_line = (num_atoms + 2) * (mol_id)
        end_line = (num_atoms + 2) * (mol_id + 1)
        ensemble_slice = ensemble_data[start_line:end_line]

        cartesian = ensemble_slice[2:]
        elements = list()
        coordinates = list()

        for i in cartesian:
            elements.append(i.split()[0])
            coordinates.append(i.split()[1:])

        atomic_numbers = [element_to_Z[element] for element in elements]

        molecule_data['num_atoms'] = int(ensemble_slice[0])
        molecule_data['header'] = ensemble_slice[1]
        molecule_data['elements'] = np.array(elements)
        molecule_data['atomic_numbers'] = np.array(atomic_numbers)
        molecule_data['coordinates'] = np.array(coordinates, dtype='float')
        molecule_data['stringContent'] = build_xyz_file(
            molecule_data['elements'], 
            molecule_data['coordinates'], 
            molecule_data['header']
        )
        parsed_ensemble[mol_id] = molecule_data

    return parsed_ensemble