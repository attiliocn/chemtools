#!/usr/bin/env python3

###############################################################################
# Extract the optimised geometry from NWChem output files and write XYZ files.
#
# Priority logic:
#   1. If "Optimization converged" is found, use the last coordinate block
#      that appears before that line (the converged geometry).
#   2. Otherwise fall back to the last coordinate block in the file and warn.
#
# Coordinate blocks are identified by "Output coordinates in angstroms".
# NWChem also prints AU coordinates; the Bohr→Angstrom conversion factor is
# taken from scipy.constants so no magic numbers appear in the code.
###############################################################################

import re
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import scipy.constants as const

###############################################################################
# Unit conversion
###############################################################################
BOHR_TO_ANGSTROM = const.physical_constants["Bohr radius"][0] / const.angstrom

###############################################################################
# Regex patterns
###############################################################################
COORD_ANGSTROM_RE = re.compile(r"Output coordinates in angstroms")
COORD_AU_RE       = re.compile(r"Output coordinates in a\.u\.")
CONVERGED_FLAG    = "Optimization converged"

# Atom data line:  No.  Tag  Charge  X  Y  Z
ATOM_RE = re.compile(
    r"^\s+\d+\s+([A-Za-z]+)\s+[\d.]+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)


###############################################################################
# Parsing
###############################################################################
def parse_file(lines):
    # Returns:
    #   blocks        – list of (start_line_idx, in_angstrom, [(el, x, y, z)])
    #   converged_line – line index of first "Optimization converged", or None

    blocks = []
    converged_line = None

    i = 0
    while i < len(lines):
        line = lines[i]

        if converged_line is None and CONVERGED_FLAG in line:
            converged_line = i

        is_angstrom = bool(COORD_ANGSTROM_RE.search(line))
        is_au       = bool(COORD_AU_RE.search(line))

        if is_angstrom or is_au:
            block_start = i
            i += 4  # skip blank line, "No. Tag Charge X Y Z", and "---- ..."
            atoms = []
            while i < len(lines):
                m = ATOM_RE.match(lines[i])
                if m:
                    atoms.append((m.group(1), m.group(2), m.group(3), m.group(4)))
                    i += 1
                else:
                    break
            if atoms:
                blocks.append((block_start, is_angstrom, atoms))
            continue  # don't double-increment

        i += 1

    return blocks, converged_line


###############################################################################
# Geometry selection
###############################################################################
def select_geometry(blocks, converged_line):
    # Returns (atoms, in_angstrom, converged_flag)

    if not blocks:
        return None, True, False

    if converged_line is not None:
        candidates = [b for b in blocks if b[0] < converged_line]
        if candidates:
            start, in_ang, atoms = candidates[-1]
            return atoms, in_ang, True

    start, in_ang, atoms = blocks[-1]
    return atoms, in_ang, False


###############################################################################
# Main conversion routine
###############################################################################
def convert(input_file):
    path = Path(input_file)
    try:
        lines = path.read_text().splitlines()
    except OSError as e:
        print(f"Cannot read {input_file}: {e}", file=sys.stderr)
        return

    blocks, converged_line = parse_file(lines)

    if not blocks:
        print(f"No coordinate blocks found in {input_file}", file=sys.stderr)
        return

    atoms, in_angstrom, converged = select_geometry(blocks, converged_line)

    if not converged:
        print(
            f"Warning: {input_file} – optimization not converged; "
            f"using last available geometry (step {len(blocks)}).",
            file=sys.stderr,
        )

    status  = "Optimized geometry" if converged else f"Last geometry (step {len(blocks)}, not converged)"
    comment = f"{input_file} | {status}"

    out_path = path.with_suffix(".xyz")
    with out_path.open("w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for elem, x, y, z in atoms:
            xf, yf, zf = float(x), float(y), float(z)
            if not in_angstrom:
                xf *= BOHR_TO_ANGSTROM
                yf *= BOHR_TO_ANGSTROM
                zf *= BOHR_TO_ANGSTROM
            f.write(f"  {elem:<4s}  {xf:14.8f}  {yf:14.8f}  {zf:14.8f}\n")

    print(f"Written: {out_path}")


###############################################################################
# Entry point
###############################################################################
if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        print("Usage: outnwchem.py <file1.out> [file2.out ...]")
        sys.exit(0)

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(convert, f): f for f in files}
        for future in as_completed(futures):
            exc = future.exception()
            if exc:
                print(f"Error processing {futures[future]}: {exc}", file=sys.stderr)