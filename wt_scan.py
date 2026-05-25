#!/usr/bin/env python3
"""
ASPB-IE Phase 1: WT per-residue decomposition scan
===================================================

Runs PB decomposition on a wild-type complex to identify per-residue
energy contributions. Outputs ranked hotspot candidates for subsequent
alanine scanning (ASPB-IE Phase 2).

Usage
-----
    python wt_scan.py \\
        --topol topol.top --tpr md.tpr --traj md.xtc \\
        --index index.ndx --cg "1 13" --cr ref.pdb \\
        [--prefix WT_scan] [--top-n 20] [--threshold -1.0] \\
        [--mpi 5]

Workflow
--------
Phase 1 (this script):
    PB decomposition on WT complex → per-residue energy ranking →
    hotspot candidates → auto-generates ala.in for Phase 2.

Phase 2 (vel_cache_runner.py or similar):
    Alanine scanning on hotspot candidates with normal_cache=1.

Outputs
-------
    {prefix}_per_residue.csv    Full per-residue energy breakdown
    {prefix}_hotspots.csv       Ranked hotspot candidates (TOTAL energy)
    {prefix}_hotspots.ala       Ready-to-use alanine scanning input
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from datetime import datetime


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_MPI      = 5
DEFAULT_TOP_N    = 20
DEFAULT_THRESHOLD = -1.0   # kcal/mol; residues with TOTAL < threshold are hotspots
DEFAULT_PREFIX   = "WT_scan"

DECOMP_INPUT_TEMPLATE = """\
WT per-residue decomposition scan (ASPB-IE Phase 1)
&general
  sys_name="WT_decomp_scan"
  startframe={startframe}
{endframe_line}  verbose=1
/
&pb
  istrng=0.15
  inp=1
  radiopt=0
/
&decomp
  idecomp=1
  csv_format=1
/
"""


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="ASPB-IE Phase 1: WT per-residue decomposition scan")
    p.add_argument("--topol",  required=True, help="GROMACS topology file")
    p.add_argument("--tpr",    required=True, help="GROMACS TPR file")
    p.add_argument("--traj",   required=True, help="GROMACS trajectory (XTC/TRR)")
    p.add_argument("--index",  required=True, help="GROMACS index file")
    p.add_argument("--cg",     required=True, help="Complex groups, e.g. '1 13'")
    p.add_argument("--cr",     required=True, help="Reference PDB file")
    p.add_argument("--startframe", type=int, default=1, help="First frame")
    p.add_argument("--endframe",   type=int, default=0, help="Last frame (0 = all)")
    p.add_argument("--prefix",     default=DEFAULT_PREFIX, help="Output prefix")
    p.add_argument("--top-n",      type=int, default=DEFAULT_TOP_N,
                    help=f"Number of top hotspot candidates (default: {DEFAULT_TOP_N})")
    p.add_argument("--threshold",  type=float, default=DEFAULT_THRESHOLD,
                    help=f"Energy threshold in kcal/mol (default: {DEFAULT_THRESHOLD})")
    p.add_argument("--mpi",  type=int, default=DEFAULT_MPI, help="MPI processes")
    p.add_argument("--workdir", default=".", help="Working directory")
    p.add_argument("--no-run", action="store_true",
                   help="Only generate input file, do not run gmx_MMPBSA")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Input file generation
# ---------------------------------------------------------------------------
def make_input(args) -> str:
    """Write the PB + decomp input file. Returns path."""
    endframe_line = ""
    if args.endframe and args.endframe > 0:
        endframe_line = f"endframe={args.endframe},\n"
    text = DECOMP_INPUT_TEMPLATE.format(
        startframe=args.startframe,
        endframe_line=endframe_line,
    )
    path = Path(args.workdir) / f"{args.prefix}_decomp.in"
    path.write_text(text)
    print(f"[INPUT] {path}")
    return str(path)


# ---------------------------------------------------------------------------
# Run gmx_MMPBSA
# ---------------------------------------------------------------------------
def run_decomp(args, inp_path: str):
    """Run gmx_MMPBSA in standard PB + decomp mode."""
    out_dat = Path(args.workdir) / f"{args.prefix}_binding.dat"
    out_csv = Path(args.workdir) / f"{args.prefix}_binding.csv"
    out_do  = Path(args.workdir) / f"{args.prefix}_decomp.dat"

    cmd = (
        f"mpirun -np {args.mpi} gmx_MMPBSA "
        f"-i {inp_path} "
        f"-cp {args.topol} -cs {args.tpr} -ct {args.traj} "
        f"-ci {args.index} -cg {args.cg} -cr {args.cr} "
        f"-o {out_dat} -eo {out_csv} "
        f"-do {out_do} -nogui -O"
    )

    print(f"[RUN] {cmd}")
    sys.stdout.flush()
    ret = subprocess.run(cmd, shell=True, cwd=args.workdir)
    if ret.returncode != 0:
        print(f"[ERROR] gmx_MMPBSA returned {ret.returncode}")
        sys.exit(1)
    return str(out_do)


# ---------------------------------------------------------------------------
# Parse decomp CSV output
# ---------------------------------------------------------------------------
def parse_decomp_csv(decomp_path: str) -> list[dict]:
    """
    Parse the CSV-format decomp output from gmx_MMPBSA.
    
    The CSV format produced by DecompOut.summary('csv') looks like:
        Complex:
        Residue,Internal,,,van der Waals,,,,Electrostatic,,,,Polar Solvation,,,,Non-Polar Solv.,,,,TOTAL,,
        ,Avg.,Std. Dev.,Std. Err. of Mean,...
        LEU A 1,-0.123,0.456,...

    Returns list of dicts with residue info and TOTAL energy.
    """
    residues = []
    # Pattern: R:CHAIN:RESNAME:RESNUM  (e.g. "R:A:ALA:17")
    res_pattern = re.compile(r'^R:([A-Z]):([A-Z]{3}):(\d+)')
    # Also try the table format: "RESNAME CHAIN RESNUM" (e.g. "LEU A 1")
    table_pattern = re.compile(r'^([A-Z]{3})\s+([A-Z]?)\s+(\d+)')

    with open(decomp_path, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or not row[0]:
                continue
            # Skip header rows
            if row[0].startswith('Complex') or row[0].startswith('Residue'):
                continue
            if row[0].startswith('Receptor') or row[0].startswith('Ligand'):
                continue
            if 'idecomp' in row[0] or 'Energy' in row[0] or 'DELTAS' in row[0]:
                continue
            # Skip empty rows
            row[0] = row[0].strip()
            if not row[0]:
                continue

            # Try the R:CHAIN:RESNAME:RESNUM format first
            m = res_pattern.match(row[0])
            if m:
                resname = m.group(2)
                chain   = m.group(1)
                resnum  = int(m.group(3))
            else:
                # Fall back to table format: "RESNAME CHAIN RESNUM"
                m = table_pattern.match(row[0])
                if not m:
                    continue
                resname = m.group(1)
                chain   = m.group(2) if m.group(2) else 'A'
                resnum  = int(m.group(3))
            restitle = f"{resname} {chain}:{resnum}"

            # Extract TOTAL average (column 16, 0-indexed)
            try:
                tot_avg = float(row[16]) if len(row) > 16 and row[16] else 0.0
                tot_std = float(row[17]) if len(row) > 17 and row[17] else 0.0
            except (ValueError, IndexError):
                # Skip rows that don't have valid numeric data
                continue

            residues.append({
                'resname':  resname,
                'chain':    chain,
                'resnum':   resnum,
                'title':    restitle,
                'TOTAL_avg': tot_avg,
                'TOTAL_std': tot_std,
            })

    return residues


# ---------------------------------------------------------------------------
# Rank and filter hotspots
# ---------------------------------------------------------------------------
def rank_hotspots(residues: list[dict], top_n: int, threshold: float):
    """
    Rank residues by TOTAL energy (most negative = strongest contributor).
    Returns (all_sorted, hotspot_candidates).
    """
    # Sort by TOTAL energy (most negative first)
    sorted_res = sorted(residues, key=lambda r: r['TOTAL_avg'])

    # Hotspots: those below threshold OR top N if threshold yields few
    by_energy = [r for r in sorted_res if r['TOTAL_avg'] < threshold]
    if len(by_energy) < 1:
        by_energy = sorted_res[:top_n]

    # Keep top N
    hotspots = by_energy[:top_n]
    return sorted_res, hotspots


# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
def write_per_residue_csv(residues: list[dict], path: str):
    """Write full per-residue breakdown."""
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Residue', 'Chain', 'ResNum', 'TOTAL_avg_kcal_mol', 'TOTAL_std_kcal_mol'])
        for r in residues:
            w.writerow([r['resname'], r['chain'], r['resnum'],
                        f"{r['TOTAL_avg']:.3f}", f"{r['TOTAL_std']:.3f}"])
    print(f"[OUT] {path}  ({len(residues)} residues)")


def write_hotspots_csv(hotspots: list[dict], path: str):
    """Write ranked hotspot candidates."""
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Rank', 'Residue', 'Chain', 'ResNum', 'TOTAL_avg_kcal_mol', 'TOTAL_std_kcal_mol'])
        for i, r in enumerate(hotspots, 1):
            w.writerow([i, r['resname'], r['chain'], r['resnum'],
                        f"{r['TOTAL_avg']:.3f}", f"{r['TOTAL_std']:.3f}"])
    print(f"[OUT] {path}  ({len(hotspots)} hotspots)")


def write_ala_template(hotspots: list[dict], prefix: str, workdir: str):
    """Write an ala.in template for Phase 2 alanine scanning.

    The template follows the format expected by vel_cache_runner.py:
    only &general and &pb sections (no &alanine_scanning block —
    the runner adds it dynamically for each mutation).
    """
    lines = [
        "# ASPB-IE Phase 2: alanine scanning on WT-predicted hotspots",
        "# Generated by wt_scan.py on " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "#",
        "# Hotspot residues (ranked by PB decomposition TOTAL energy):",
    ]
    for i, r in enumerate(hotspots, 1):
        lines.append(f"#   {i:>3}. {r['title']:>12s}   TOTAL = {r['TOTAL_avg']:>8.3f} kcal/mol")

    lines.extend([
        "#",
        "# Usage: run with vel_cache_runner.py which appends &alanine_scanning",
        "#        block with normal_cache=1 for each mutation site.",
        "#",
        "# To run a single mutation manually:",
        "#   cp this file ala.in",
        "#   then edit vel_cache_runner.py SITES = [" + ",".join(str(r['resnum']) for r in hotspots[:3]) + ("..." if len(hotspots) > 3 else "") + "]",
        "#   and run: python vel_cache_runner.py",
        "",
        "&general",
        "  sys_name='ASPB-IE'",
        "  startframe=1",
        "  endframe=1000",
        "  verbose=1",
        "/",
        "&pb",
        "  istrng=0.15",
        "  inp=1",
        "  radiopt=0",
        "/",
    ])

    path = Path(workdir) / f"{prefix}_hotspots.ala"
    path.write_text('\n'.join(lines) + '\n')
    print(f"[OUT] {path}  ({len(hotspots)} hotspot sites)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  ASPB-IE Phase 1: WT per-residue decomposition scan")
    print("=" * 60)

    # 1. Build input
    inp_path = make_input(args)

    if args.no_run:
        print("[SKIP] --no-run set; input file generated. Exiting.")
        return

    # 2. Run gmx_MMPBSA decomp
    decomp_path = run_decomp(args, inp_path)
    print(f"[DECOMP] {decomp_path}")

    # 3. Parse results
    residues = parse_decomp_csv(decomp_path)
    if not residues:
        print("[ERROR] No residue decomposition data found in output.")
        print("        Check that &decomp with idecomp=1, csv_format=1 was used.")
        sys.exit(1)
    print(f"[PARSE]  Found {len(residues)} residue entries.")

    # 4. Rank hotspots
    sorted_res, hotspots = rank_hotspots(residues, args.top_n, args.threshold)
    print(f"[HOTSPOTS] Top {len(hotspots)} candidates (threshold={args.threshold} kcal/mol):")
    for i, r in enumerate(hotspots, 1):
        print(f"          {i:>3}. {r['title']:>12s}   TOTAL = {r['TOTAL_avg']:>8.3f} ± {r['TOTAL_std']:>6.3f} kcal/mol")

    # 5. Write outputs
    write_per_residue_csv(sorted_res, Path(workdir) / f"{args.prefix}_per_residue.csv")
    write_hotspots_csv(hotspots, Path(workdir) / f"{args.prefix}_hotspots.csv")
    write_ala_template(hotspots, args.prefix, str(workdir))

    print("-" * 60)
    print("  ASPB-IE Phase 1 complete.")
    print(f"  Next: use {args.prefix}_hotspots.ala for alanine scanning Phase 2")
    print(f"        with normal_cache=1 (already set in template).")
    print("=" * 60)


if __name__ == '__main__':
    main()
