#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sort and concatenate .xyz conformation files based on an energy-sorted list."
    )
    parser.add_argument(
        "energy_file",
        type=Path,
        help="Text file with filenames in desired order (e.g. energy-sort.txt)",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("."),
        help="Directory where the .xyz files are located (default: current dir)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory to write output files (default: current dir)",
    )
    parser.add_argument(
        "--suffix",
        default="confs",
        help="Suffix appended to basename for output file (default: 'confs')",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Silently skip missing files instead of warning",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    filenames = [
        l.strip()
        for l in args.energy_file.read_text().splitlines()
        if l.strip()
    ]

    # Group filenames by basename (strip trailing _<digits>.xyz)
    groups = defaultdict(list)
    for fname in filenames:
        base = Path(fname).stem.rsplit("_", 1)[0]  # e.g. "bitet5_confs_102" -> "bitet5_confs"
        groups[base].append(fname)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for base, fnames in groups.items():
        outfile = args.output_dir / f"{base}_{args.suffix}.xyz"
        written = 0
        with outfile.open("w") as out:
            for fname in fnames:
                p = args.input_dir / fname
                if p.exists():
                    out.write(p.read_text())
                    written += 1
                elif not args.skip_missing:
                    print(f"  WARNING: {fname} not found, skipping")
        print(f"Created: {outfile}  ({written}/{len(fnames)} files written)")


if __name__ == "__main__":
    main()