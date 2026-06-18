#!/usr/bin/env python3

import re
import scipy.constants as const

# --- Unit conversions (scipy.constants) -------------------------------------
# 1 hartree in J/mol  (Hartree energy per molecule * Avogadro)
HARTREE_PER_JMOL = const.physical_constants['Hartree energy'][0] * const.N_A
# 1 cal in J
CAL_TO_J = const.calorie  # 4.184


def cal_per_mol_to_hartree(value_cal_per_mol: float) -> float:
    """cal/mol -> hartree (au)."""
    return value_cal_per_mol * CAL_TO_J / HARTREE_PER_JMOL


# --- NWChem regex patterns ---------------------------------------------------
# Energies that NWChem prints in atomic units are captured directly in au.
# Thermochemistry corrections are printed as 'X kcal/mol  ( Y au)'; we grab Y.
# Entropies are printed in cal/mol-K; we capture and convert later.

regex_patterns = {
    # 1. Total DFT energy (au). May appear multiple times -> keep LAST occurrence.
    'Total DFT energy': re.compile(
        r'Total DFT energy\s*=\s*(-?[0-9]+\.[0-9]+)'),

    # 2. Temperature (K)
    'Temperature': re.compile(
        r'Temperature\s*=\s*([0-9]+\.?[0-9]*)K'),

    # 3. Zero-Point correction to Energy (au, from the parenthetical)
    'Zero-Point correction to Energy': re.compile(
        r'Zero-Point correction to Energy\s*=.*?\(\s*(-?[0-9]+\.[0-9]+)\s*au\)'),

    # 4. Thermal correction to Energy (au)
    'Thermal correction to Energy': re.compile(
        r'Thermal correction to Energy\s*=.*?\(\s*(-?[0-9]+\.[0-9]+)\s*au\)'),

    # 5. Thermal correction to Enthalpy (au)
    'Thermal correction to Enthalpy': re.compile(
        r'Thermal correction to Enthalpy\s*=.*?\(\s*(-?[0-9]+\.[0-9]+)\s*au\)'),

    # 6. Total Entropy (cal/mol-K)
    'Total Entropy': re.compile(
        r'Total Entropy\s*=\s*(-?[0-9]+\.[0-9]+)\s*cal/mol-K'),

    # 7. Translational entropy (cal/mol-K)
    'Translational entropy': re.compile(
        r'-\s*Translational\s*=\s*(-?[0-9]+\.[0-9]+)\s*cal/mol-K'),

    # 8. Rotational entropy (cal/mol-K)
    'Rotational entropy': re.compile(
        r'-\s*Rotational\s*=\s*(-?[0-9]+\.[0-9]+)\s*cal/mol-K'),

    # 9. Vibrational entropy (cal/mol-K)
    'Vibrational entropy': re.compile(
        r'-\s*Vibrational\s*=\s*(-?[0-9]+\.[0-9]+)\s*cal/mol-K'),
}

# NWChem prints Translational/Rotational/Vibrational TWICE: once under
# "Total Entropy" and once under "Cv (constant volume heat capacity)".
# We want the entropy block (the first occurrence of each), so for these three
# keys we keep the FIRST match and ignore later ones.
ENTROPY_FIRST_ONLY = {'Translational entropy', 'Rotational entropy', 'Vibrational entropy'}


def extract_nwchem_thermo(filepath: str) -> dict:
    """Parse one NWChem output file into a dict of energies/quantities (au, K, cal/mol-K)."""
    parsed = dict()
    with open(filepath) as f:
        for line in f:
            for name, pattern in regex_patterns.items():
                match = pattern.search(line)
                if match:
                    value = float(match.group(1))
                    # Total DFT energy: keep the LAST occurrence -> always overwrite.
                    # Entropy components: keep the FIRST (the entropy block, not Cv).
                    if name in ENTROPY_FIRST_ONLY and name in parsed:
                        continue
                    parsed[name] = value

    # --- 10. Correction to Gibbs Energy (calculated) ------------------------
    # G_corr = H_corr - T * S
    #   H_corr in hartree, S in cal/mol-K, T in K  ->  T*S in cal/mol -> hartree
    if all(k in parsed for k in
           ('Thermal correction to Enthalpy', 'Temperature', 'Total Entropy')):
        H_corr = parsed['Thermal correction to Enthalpy']      # au
        T = parsed['Temperature']                              # K
        S = parsed['Total Entropy']                            # cal/mol-K
        TS_au = cal_per_mol_to_hartree(T * S)                  # au
        parsed['Correction to Gibbs Energy'] = H_corr - TS_au

    return parsed


# --- Output column layout ----------------------------------------------------
# Reference numbers from the request:
#   1 = Total DFT energy
#   2 = Temperature
#   3 = Zero-Point correction to Energy
#   4 = Thermal correction to Energy
#   5 = Thermal correction to Enthalpy
#   6 = Total Entropy
#   7 = Translational entropy
#   8 = Rotational entropy
#   9 = Vibrational entropy
#  10 = Correction to Gibbs Energy (calculated)
#
# Requested final order: 1, (1+3), (1+5), (1+10), 3, 5, 10, 2, 6, 7, 8, 9
KEY = {
    1: 'Total DFT energy',
    2: 'Temperature',
    3: 'Zero-Point correction to Energy',
    4: 'Thermal correction to Energy',
    5: 'Thermal correction to Enthalpy',
    6: 'Total Entropy',
    7: 'Translational entropy',
    8: 'Rotational entropy',
    9: 'Vibrational entropy',
    10: 'Correction to Gibbs Energy',
}

# (header label, function mapping a parsed dict -> value or None)
def _get(d, n):
    return d.get(KEY[n])

def _sum(d, *ns):
    vals = [d.get(KEY[n]) for n in ns]
    return sum(vals) if all(v is not None for v in vals) else None

OUTPUT_COLUMNS = [
    ('Total DFT energy (au)',                 lambda d: _get(d, 1)),
    ('DFT + ZPE (au)',                        lambda d: _sum(d, 1, 3)),
    ('DFT + H_corr -> Enthalpy (au)',         lambda d: _sum(d, 1, 5)),
    ('DFT + G_corr -> Gibbs (au)',            lambda d: _sum(d, 1, 10)),
    ('Zero-Point correction (au)',            lambda d: _get(d, 3)),
    ('Thermal correction to Enthalpy (au)',   lambda d: _get(d, 5)),
    ('Correction to Gibbs Energy (au)',       lambda d: _get(d, 10)),
    ('Temperature (K)',                       lambda d: _get(d, 2)),
    ('Total Entropy (cal/mol-K)',             lambda d: _get(d, 6)),
    ('Translational entropy (cal/mol-K)',     lambda d: _get(d, 7)),
    ('Rotational entropy (cal/mol-K)',        lambda d: _get(d, 8)),
    ('Vibrational entropy (cal/mol-K)',       lambda d: _get(d, 9)),
]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Parse NWChem SPE/thermochemistry output to a CSV of energies (au).'
    )
    parser.add_argument('files', nargs='+', help='NWChem output files')
    parser.add_argument('-o', '--output', default='nwchem_thermo.csv',
                        help='Output CSV filename (default: nwchem_thermo.csv)')
    args = parser.parse_args()

    parsed_files = {}
    for file in args.files:
        parsed_files[file] = extract_nwchem_thermo(file)

    with open(args.output, mode='w') as f:
        headers = ['filename'] + [label for label, _ in OUTPUT_COLUMNS]
        f.write(','.join(headers) + '\n')
        for fname, data in parsed_files.items():
            row = [f'./{fname}']
            for _, func in OUTPUT_COLUMNS:
                val = func(data)
                row.append('' if val is None else str(val))
            f.write(','.join(row) + '\n')