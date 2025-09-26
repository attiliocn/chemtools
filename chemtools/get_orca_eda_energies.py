#!/usr/bin/env python3

import re
from scipy.constants import physical_constants, Avogadro, calorie

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

    with open('eda_energies.csv', mode='w') as f:
        f.write(f'filename,{','.join(regex_patterns.keys())}\n')
        for k,v in eda_parsed_files.items():
            sorted_v = {s:v[s]*HARTREE_TO_KCALMOL for s in regex_patterns.keys()}
            sorted_v['smd cds'] = sorted_v['smd cds'] / HARTREE_TO_KCALMOL
            f.write(f'./{k},{','.join(str(s) for s in sorted_v.values())}\n')