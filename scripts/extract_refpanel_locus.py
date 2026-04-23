#!/usr/bin/env python3
import argparse
import json
import os
import subprocess


def _mkdir_for(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def _read_loci(loci_path: str):
    loci = {}
    with open(loci_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].lower() in {"locus", "loci", "id"}:
                continue
            if len(parts) < 4:
                continue
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue
            loci[parts[0]] = (str(parts[1]).replace("chr", ""), start, end)
    return loci


def _count_bim_rows(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path, "r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract locus-specific PLINK reference panel files")
    parser.add_argument("--ref-prefix", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--locus", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--diag-json", required=True)
    parser.add_argument("--done-file", required=True)
    args = parser.parse_args()

    _mkdir_for(args.out_prefix + ".bim")
    _mkdir_for(args.diag_json)
    _mkdir_for(args.done_file)

    loci = _read_loci(args.loci)
    if args.locus not in loci:
        raise ValueError(f"Locus '{args.locus}' not found in loci file: {args.loci}")

    chr_, start, end = loci[args.locus]

    make_bed_cmd = [
        "plink2",
        "--bfile",
        args.ref_prefix,
        "--chr",
        chr_,
        "--from-bp",
        str(start),
        "--to-bp",
        str(end),
        "--make-bed",
        "--out",
        args.out_prefix,
    ]
    subprocess.run(make_bed_cmd, check=True)

    freq_cmd = [
        "plink2",
        "--bfile",
        args.out_prefix,
        "--freq",
        "--out",
        args.out_prefix,
    ]
    subprocess.run(freq_cmd, check=True)

    n_variants = _count_bim_rows(args.out_prefix + ".bim")
    diag = {
        "locus": args.locus,
        "requested_window": {"chr": chr_, "start": start, "end": end},
        "ref_prefix_in": args.ref_prefix,
        "ref_prefix_out": args.out_prefix,
        "n_variants_extracted": int(n_variants),
        "files": {
            "bed": args.out_prefix + ".bed",
            "bim": args.out_prefix + ".bim",
            "fam": args.out_prefix + ".fam",
            "afreq": args.out_prefix + ".afreq",
        },
    }
    with open(args.diag_json, "w", encoding="utf-8") as handle:
        json.dump(diag, handle, indent=2)

    with open(args.done_file, "w", encoding="utf-8") as handle:
        handle.write("ok\n")


if __name__ == "__main__":
    main()
