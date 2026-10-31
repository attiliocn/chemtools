#!/usr/bin/env python3

import re
import scipy.constants as const

# --- Unit conversions (scipy.constants) -------------------------------------
# 1 hartree in kcal/mol  (Hartree energy per molecule * Avogadro / 1000 cal)
HARTREE_TO_KCALMOL = (const.physical_constants['Hartree energy'][0] * const.N_A
                      / (1000 * const.calorie))


# --- ORCA 6 EDA-NOCV regex patterns ------------------------------------------
# The "Energy Decomposition Analysis" table prints each term as
#     '<Term>   <Hartree>   <kcal/mol>'
# We capture the Hartree column (more digits) and convert to kcal/mol later.
regex_patterns = {
    # 1. Bond (interaction) energy
    'Bond Energy': re.compile(
        r'^\s*Bond Energy\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 2. Orbital (relaxation) energy = Sum_k DE_k of the NOCV table
    'Orbital Energy': re.compile(
        r'^\s*Orbital Energy\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 3. Electrostatic energy
    'Electrostatic Energy': re.compile(
        r'^\s*Electrostatic Energy\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 4. Pauli repulsion
    'Pauli Energy': re.compile(
        r'^\s*Pauli Energy\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 5. Exchange-correlation contribution of the promolecule
    'Delta E^0(XC)': re.compile(
        r'^\s*Delta E\^0\(XC\)\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 6. Dispersion (only with D3/D4)
    'Delta Dispersion': re.compile(
        r'^\s*Delta Dispersion\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 7. CPCM dielectric term (only with implicit solvation)
    'Delta CPCM Dielectric': re.compile(
        r'^\s*Delta CPCM Dielectric\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),

    # 8. SMD non-electrostatic term (only with SMD)
    'Delta SMD CDS correction': re.compile(
        r'^\s*Delta SMD CDS correction\s+(-?[0-9]+\.[0-9]+)\s+-?[0-9]+\.[0-9]+'),
}

# NOCV table rows: 'neg. eigenvalue  pos. eigenvalue  DE_k  DT_k  DV_k'
# (eigenvalues in e, energies already in kcal/mol)
NOCV_ROW = re.compile(
    r'^\s*(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s+'
    r'(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s*$')
# 'Sum_k:  DE  DT  DV' (kcal/mol)
NOCV_SUM = re.compile(
    r'Sum_k:\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)')

NOCV_COLUMNS = ['neg. eigenvalue (e)', 'pos. eigenvalue (e)',
                'DE_k (kcal/mol)', 'DT_k (kcal/mol)', 'DV_k (kcal/mol)']


def extract_orca_eda_nocv(filepath: str) -> dict:
    """Parse one ORCA 6 EDA-NOCV output.

    Returns {'eda': {term: Hartree}, 'nocv': [row tuples], 'nocv_sum': (DE, DT, DV)}.
    """
    eda = dict()
    nocv = []
    nocv_sum = None
    in_nocv = False
    with open(filepath) as f:
        for line in f:
            # EDA terms: anchored '<Term> <Eh> <kcal/mol>' lines are unique to the
            # EDA table, so a plain line scan is safe (keep LAST occurrence).
            for name, pattern in regex_patterns.items():
                match = pattern.search(line)
                if match:
                    eda[name] = float(match.group(1))

            # --- NOCV table: starts at the 'NOCV analysis' header and ends at Sum_k
            if line.strip() == 'NOCV analysis':
                in_nocv = True
                nocv = []  # keep only the LAST table if there is more than one
                continue
            if in_nocv:
                match = NOCV_SUM.search(line)
                if match:
                    nocv_sum = tuple(float(x) for x in match.groups())
                    in_nocv = False
                    continue
                match = NOCV_ROW.match(line)
                if match:
                    nocv.append(tuple(float(x) for x in match.groups()))

    return {'eda': eda, 'nocv': nocv, 'nocv_sum': nocv_sum}


# --- Output column layout ----------------------------------------------------
# Summary CSV (one row per file), all energies in kcal/mol:
#   EDA terms, sum of terms (should reproduce the Bond Energy),
#   NOCV Sum_k (DE should reproduce the Orbital Energy) and the
#   first N_NOCV_PAIRS DE_k / eigenvalues.
EDA_TERMS = list(regex_patterns.keys())


def _kcal(d, name):
    v = d['eda'].get(name)
    return None if v is None else v * HARTREE_TO_KCALMOL


def _sum_of_terms(d):
    terms = [v for k, v in d['eda'].items() if k != 'Bond Energy']
    return sum(terms) * HARTREE_TO_KCALMOL if terms else None


def _nocv_sum(d, i):
    return None if d['nocv_sum'] is None else d['nocv_sum'][i]


def build_output_columns(n_pairs: int):
    columns = [(f'{name} (kcal/mol)', lambda d, n=name: _kcal(d, n))
               for name in EDA_TERMS]
    columns += [
        ('Sum of terms (kcal/mol)',   _sum_of_terms),
        ('NOCV Sum_k DE (kcal/mol)',  lambda d: _nocv_sum(d, 0)),
        ('NOCV Sum_k DT (kcal/mol)',  lambda d: _nocv_sum(d, 1)),
        ('NOCV Sum_k DV (kcal/mol)',  lambda d: _nocv_sum(d, 2)),
    ]
    for k in range(n_pairs):
        columns += [
            (f'NOCV {k + 1} eigenvalue (e)',
             lambda d, k=k: d['nocv'][k][1] if k < len(d['nocv']) else None),
            (f'NOCV {k + 1} DE_k (kcal/mol)',
             lambda d, k=k: d['nocv'][k][2] if k < len(d['nocv']) else None),
        ]
    return columns


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Parse ORCA 6 EDA-NOCV output to CSVs (EDA summary + NOCV pairs, kcal/mol).'
    )
    parser.add_argument('files', nargs='+', help='ORCA EDA-NOCV output files')
    parser.add_argument('-o', '--output', default='orca_eda.csv',
                        help='Summary CSV filename (default: orca_eda.csv)')
    parser.add_argument('--nocv-output', default='orca_nocv.csv',
                        help='Full NOCV pairs CSV filename (default: orca_nocv.csv)')
    parser.add_argument('-n', '--n-pairs', type=int, default=5,
                        help='Number of leading NOCV pairs in the summary CSV (default: 5)')
    args = parser.parse_args()

    parsed_files = {}
    for file in args.files:
        parsed_files[file] = extract_orca_eda_nocv(file)

    # --- Summary: one row per file
    output_columns = build_output_columns(args.n_pairs)
    with open(args.output, mode='w') as f:
        headers = ['filename'] + [label for label, _ in output_columns]
        f.write(','.join(headers) + '\n')
        for fname, data in parsed_files.items():
            row = [f'./{fname}']
            for _, func in output_columns:
                val = func(data)
                row.append('' if val is None else str(val))
            f.write(','.join(row) + '\n')

    # --- Full NOCV table: one row per (file, pair)
    with open(args.nocv_output, mode='w') as f:
        f.write(','.join(['filename', 'k'] + NOCV_COLUMNS) + '\n')
        for fname, data in parsed_files.items():
            for k, values in enumerate(data['nocv'], start=1):
                f.write(','.join([f'./{fname}', str(k)] + [str(v) for v in values]) + '\n')
