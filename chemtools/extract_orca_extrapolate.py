#!/usr/bin/env python3

import re

# ORCA 504 Patterns
regex_patterns_504 = {
    'scf tzvpp': re.compile(r'SCF energy with basis def2-TZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'scf qzvpp': re.compile(r'SCF energy with basis def2-QZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'scf cbs': re.compile(r'Extrapolated CBS SCF energy \(3/4\) : +(\-?[0-9]+\.[0-9]+)'),
    'mdci tzvpp': re.compile(r'MDCI energy with basis def2-TZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'mdci qzvpp': re.compile(r'MDCI energy with basis def2-QZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'mdci cbs': re.compile(r'Extrapolated CBS correlation energy \(3/4\) : +(\-?[0-9]+\.[0-9]+)'),
    'E(el) cbs': re.compile(r'Estimated CBS total energy \(3/4\) : +(\-?[0-9]+\.[0-9]+)')
}

# ORCA 504 Patterns for Direct Extraction (no extrapolation)
regex_patterns_504_direct = {
    'Singles Norm': re.compile(r'Singles Norm \<S\|S\>\*\*1\/2 +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'T1 diagnostic': re.compile(r'T1 diagnostic +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'E(0)': re.compile(r'E\(0\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'E(CORR)(strong-pairs)': re.compile(r'E\(CORR\)\(strong-pairs\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'E(CORR)(weak-pairs)': re.compile(r'E\(CORR\)\(weak-pairs\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'Triples Correction (T)': re.compile(r'Triples Correction \(T\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'Final correlation energy': re.compile(r'Final correlation energy +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'E(CCSD)': re.compile(r'E\(CCSD\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
    'E(CCSD(T))': re.compile(r'E\(CCSD\(T\)\) +\.\.\. +(\-?[0-9]+\.[0-9]+)'),
}

# ORCA 611 PATTERNS
regex_patterns_611 = {
    'scf tzvpp': re.compile(r'SCF energy with basis +def2-TZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'scf qzvpp': re.compile(r'SCF energy with basis +def2-QZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'scf cbs': re.compile(r'Extrapolated CBS SCF energy +\(3/4\) +: +(\-?[0-9]+\.[0-9]+)'),
    'mdci tzvpp': re.compile(r'Correlation energy with basis +def2-TZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'mdci qzvpp': re.compile(r'Correlation energy with basis +def2-QZVPP: +(\-?[0-9]+\.[0-9]+)'),
    'mdci cbs': re.compile(r'Extrapolated CBS Correlation energy +\(3/4\) +: +(\-?[0-9]+\.[0-9]+)'),
    'E(el) cbs': re.compile(r'Extrapolated CBS Total energy +\(3/4\) +: +(\-?[0-9]+\.[0-9]+)')
}

def extract_orca_extrapolate(filepath: str) -> dict:
    with open(filepath) as f:
        parsed_energies = dict()
        for line in f:
            for name, pattern in regex_patterns.items():
                match = pattern.search(line)
                if match:
                    parsed_energies[name] = float(match.group(1))
    return parsed_energies

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Extract the Extrapolate (CBS) table from the ORCA output file'
    )
    parser.add_argument('files', nargs='+', help='ORCA output files')
    parser.add_argument('--version', choices=['504', '611'], default='611', help='ORCA version (default: 611)')
    parser.add_argument('--mode', choices=['extrapolate', 'direct'], default='extrapolate', help='Extraction mode: extrapolate (default) or direct (no extrapolation, extract the raw energy values from a DLPNO-CCSD(T) calculation)')
    args = parser.parse_args()

    if args.mode == 'direct':
        if args.version != '504':
            raise ValueError('Direct extraction mode is only supported for ORCA 5.0.4')
        regex_patterns = regex_patterns_504_direct
    else:
        if args.version == '504':
            regex_patterns = regex_patterns_504
        elif args.version == '611':        
            regex_patterns = regex_patterns_611

    cbs_parsed_files = {}
    for file in args.files:
        cbs_parsed_files[file] = extract_orca_extrapolate(file)

    with open('cbs_energies.csv', mode='w') as f:
        f.write(f'filename,{','.join(regex_patterns.keys())}\n')
        for k,v in cbs_parsed_files.items():
            sorted_v = {s:v[s] for s in regex_patterns.keys()}
            f.write(f'./{k},{','.join(str(s) for s in sorted_v.values())}\n')