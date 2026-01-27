#!/usr/bin/env python3

import re

import pandas as pd
from scipy.constants import Avogadro, calorie, physical_constants

regex_patterns = {
    'bond energy': re.compile(r"Bond Energy +(\-?[0-9]+\.[0-9]+)"),
    'orbital energy': re.compile(r"Orbital Energy +(\-?[0-9]+\.[0-9]+)"),
    'electrostatic': re.compile(r"Electrostatic Energy +(\-?[0-9]+\.[0-9]+)"),
    'pauli': re.compile(r"Pauli Energy +(\-?[0-9]+\.[0-9]+)"),
    'delta E0(XC)': re.compile(r"Delta E\^0\(XC\) +(\-?[0-9]+\.[0-9]+)"),
    'delta dispersion': re.compile(r"Delta Dispersion +(\-?[0-9]+\.[0-9]+)"),
    'delta cpcm dieletric': re.compile(r"Delta CPCM Dielectric +(\-?[0-9]+\.[0-9]+)"),
    'smd cds': re.compile(r"SMD CDS free energy correction energy : +(\-?[0-9]+\.[0-9]+)")
}

HARTREE_TO_KCALMOL = physical_constants["Hartree energy"][0] * Avogadro / (1000 * calorie)
print(
f'''
All energies(*) are extracted from the output file in Hartree (Eh) and converted 
to kcal/mol using the conversion factor 1 Eh = {HARTREE_TO_KCALMOL} kcal/mol

(*) Except for the SMD CDS free energy correction -- this is extracted in kcal/mol'''
)

def extract_orca_energy_decomposition_analysis(filepath: str) -> dict:
    with open(filepath) as f:
        parsed_eda_energies = dict()
        for line in f:
            for name, pattern in regex_patterns.items():
                if name in parsed_eda_energies.keys():
                    continue
                match = pattern.search(line)
                if match:
                    parsed_eda_energies[name] = float(match.group(1))
    return parsed_eda_energies

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Extract the Energy Decomposition Analysis table from the ORCA output file'
    )
    parser.add_argument('files', nargs='+', help='ORCA >=6.1.0 output files')
    args = parser.parse_args()

    eda_parsed_files = {}
    for file in args.files:
        eda_parsed_files[file] = extract_orca_energy_decomposition_analysis(file)

    eda_parsed_table = pd.DataFrame(eda_parsed_files).T
    eda_parsed_table.index.name = 'filename'
    eda_parsed_table *= HARTREE_TO_KCALMOL

    if 'smd cds' in eda_parsed_table.columns:
        eda_parsed_table['smd cds'] /= HARTREE_TO_KCALMOL

    eda_parsed_table.to_csv('eda_energies.csv')