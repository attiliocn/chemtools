#!/usr/bin/env python3

import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def convert(input_file):
    path = Path(input_file)
    lines = path.read_text().splitlines()

    angstroem_idx = None
    au_idx = None

    for i, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            angstroem_idx = i
        if "CARTESIAN COORDINATES (A.U.)" in line:
            au_idx = i

    if angstroem_idx is None or au_idx is None:
        return

    coord_lines = lines[angstroem_idx + 2 : au_idx - 2]

    out_path = path.with_suffix(".xyz")
    with out_path.open("w") as f:
        f.write(f"{len(coord_lines)}\n")
        f.write(f"{input_file}\n")
        f.write("\n".join(coord_lines) + "\n")


if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        sys.exit(0)
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(convert, f): f for f in files}
        for future in as_completed(futures):
            exc = future.exception()
            if exc:
                print(f"Error processing {futures[future]}: {exc}", file=sys.stderr)
