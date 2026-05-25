[![Python](https://img.shields.io/badge/Python-v3.x-blue)]()
[![Pypi](https://img.shields.io/pypi/v/gmx-MMPBSA)](https://pypi.org/project/gmx-MMPBSA/)
[![Downloads](https://pepy.tech/badge/gmx-mmpbsa)](https://pepy.tech/project/gmx-mmpbsa)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.1c00645-blue)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645)

[![Help forum](https://img.shields.io/badge/Docs-Getting_Started-blue)](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/getting-started/)
[![Help forum](https://img.shields.io/badge/Help_forum-Google_group-blue)](https://groups.google.com/g/gmx_mmpbsa)
[![Issue tracking](https://img.shields.io/badge/Issue_tracking-GitHub-blue)](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues)

[![Support](https://img.shields.io/badge/Support-JetBrains-brightgreen)](https://www.jetbrains.com/?from=gmx_MMPBSA)
[![Support](https://img.shields.io/badge/Support-Sourcery-orange)](https://sourcery.ai/invite/gndRrjlo)


# gmx_MMPBSA-ASPB-IE

**GROMACS-native alanine scanning with WT cache.**

ASPB-IE is an alanine scanning workflow designed for **GROMACS users**. All inputs are standard GROMACS
format (`topol.top`, `md.tpr`, `md.xtc`, `index.ndx`). The tool converts GROMACS topologies to Amber
internally via `cpptraj`/`tleap`, runs PB/IE evaluation with **sander (CPU only)**, and outputs per-mutation ΔΔG.

> **You must have an existing MD trajectory (`.xtc`/`.trr`) before running.** The workflow does not
> perform MD — it post-processes a completed trajectory for free energy analysis.

---

## Why this fork?

| Problem | How ASPB-IE solves it |
|---------|----------------------|
| **Repeated WT calculation** — N mutations on same trajectory recompute WT N times | `normal_cache=1`: WT sander runs **once**, cached in `.wt_cache_snapshot_backup/` |
| **File overwrite** — gmx_MMPBSA deletes/regenerates WT intermediates per mutation | `CopyCalc` + `parse_output_files` restore missing WT mdout from cache backup |
| **Cache silently invalid** — frame range change gives wrong ΔΔG | Fingerprint stores `startframe`/`endframe` + file mtime; mismatch triggers full recompute |
| **No hotspot discovery** — need a predefined residue list for alanine scanning | `wt_scan.py`: PB per-residue decomposition on WT trajectory predicts hotspots (Phase 1) |

---

## Cite us

<a href="https://www.scimagojr.com/journalsearch.php?q=5100155074&amp;tip=sid&amp;exact=no" title="SCImago Journal 
&amp; Country Rank"><img border="0" align="right" src="https://www.scimagojr.com/journal_img.php?id=5100155074" 
alt="SCImago Journal &amp; Country Rank"  /></a>

Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. **gmx_MMPBSA: A New Tool to Perform 
End-State Free Energy Calculations with GROMACS**. _Journal of Chemical Theory and Computation_, 2021 17 (10), 6281-6291
https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645

Please also consider citing MMPBSA.py's paper:

Bill R. Miller, T. Dwight McGee, Jason M. Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg. **MMPBSA.py: 
An Efficient Program for End-State Free Energy Calculations**. _Journal of Chemical Theory and Computation_, 2012 8 
(9), 3314-3321. https://pubs.acs.org/doi/10.1021/ct300418h

Please, visit [Cite gmx_MMPBSA](cite_us.md) page for more information on how to cite gmx_MMPBSA and the programs/methods implemented in it.

---------------------------------------

Authors:
- [Mario Sergio Valdés Tresanco](https://www.researchgate.net/profile/Mario-Valdes-Tresanco-2), PhD. _University of Medellin, Colombia_
- [Mario Ernesto Valdés Tresanco](https://www.researchgate.net/profile/Mario-Valdes-Tresanco), PhD Student. _University of Calgary, Canada._
- [Pedro Alberto Valiente](https://www.researchgate.net/profile/Pedro-Valiente), PhD. _University of Toronto, Canada_
- [Ernesto Moreno Frías](https://www.researchgate.net/profile/Ernesto-Moreno-Frias), PhD. _University of Medellin, Colombia_

---------------------------------------

Acknowledgements:
- First of all, to Amber and GROMACS developers. Without their incredible and hard work, gmx_MMPBSA would not exist.
- Jason Swails (Amber developer and [ParmEd](https://github.com/ParmEd/ParmEd) principal developer) for his continuous support on ParmEd issues.
- Dr. Hymavathi Veeravarapu for helping with the [introductory video](https://www.youtube.com/watch?v=_2mYeffqFIo) for gmx_MMPBSA.
- To the Open Source license of the [JetBrains](https://www.jetbrains.com) programs.
- To the [Sourcery](https://sourcery.ai/invite/gndRrjlo) team for supporting us with the [Pro version](https://sourcery.ai/pro/).
- To all researchers who help improve gmx_MMPBSA with comments, feedback, and bug reports.

---------------------------------------

# Fork modifications: Normal Cache Feature

This fork adds **WT (wild-type) calculation caching** to gmx_MMPBSA for alanine scanning, avoiding redundant re-computation of the normal trajectory PB energy across different mutation runs.

## Feature: `normal_cache`

A new `&alanine_scanning` namelist parameter:

```fortran
normal_cache = 1,    # skip normal WT sander, parse results from disk
```

When enabled, `load_calc_list()` skips setting up normal (complex/receptor/ligand) sander calculations. The `mutant_only=0` path is still followed for the parsing phase, so ΔΔG (mutant − normal) output is preserved.

## Supporting modifications

Two additional patches protect the **parse-output** phase when `normal_cache=1` is active.  
(gmx_MMPBSA's `file_setup()` clears scratch files before each run, including previous WT mdout files.)

### 1. CopyCalc fallback (`calculation.py`)

`CopyCalc.run()` copies normal receptor/ligand results to mutant when those topologies are unchanged. If the source normal file has been cleaned during `file_setup`, it now falls back to reading from a `.wt_cache_snapshot_backup/` directory in the working directory before failing.

### 2. parse_output_files restore (`main.py`)

`parse_output_files()` now checks if `normal_cache=1` is active and, before attempting to read normal WT mdout files, restores them from `.wt_cache_snapshot_backup/` if they are missing.

## How to use

1. **First run** (generates WT cache):
   ```bash
   mpirun -np 10 gmx_MMPBSA -i ala.in -cp topol.top ...
   ```

2. **Subsequent mutations** (reuse cached WT):
   ```fortran
   &alanine_scanning
   mutant='ALA', mutant_res='A:126', cas_intdiel=1, normal_cache=1
   /
   ```

## How WT cache prevents overwrite-driven recomputation

gmx_MMPBSA's internal `file_setup` deletes and recreates WT intermediate files (`*_pb.mdout.*`,
`COM.prmtop`, etc.) for **each mutation run**. Without caching, every mutation on the same trajectory
repeats the full WT sander calculation — even though the WT system never changes.

### Three-layer overwrite protection

```
Layer 1 — Backup (runner level)
  vel_cache_runner.py / v5_mmpbsa_cal.py
    ├─ backup_cache(): copies WT mdout → .wt_cache_snapshot_backup/ after first FULL run
    └─ restore_cache(): restores backup before each CACHE run

Layer 2 — CopyCalc fallback (calculation.py)
  When copy(orig_name, final_name) fails because gmx_MMPBSA deleted the source:
    → reads from .wt_cache_snapshot_backup/ instead

Layer 3 — Parse restore (main.py)
  When parse_output_files() can't find WT mdout for parsing:
    → restores missing files from .wt_cache_snapshot_backup/
```

### Fingerprint validation

Each cached run stores a `.gmxmmpbsa_cache_info` fingerprint:

```json
{
  "traj": [1234567890, 52428800],
  "tpr": [1234567890, 1048576],
  "startframe": 1,
  "endframe": 1000,
  "mpi_size": 5
}
```

If any field changes (e.g. you modify `ala.in` to skip frames), the cache **automatically invalidates**
and the next run falls back to FULL. This ensures cached WT mdout always matches the current trajectory.

### Usage

Input files: `topol.top` + `md.tpr` + `md.xtc` + `index.ndx` + `ref.pdb` + `ala.in`

```bash
# Phase 1: WT decomposition scan (optional — only if you need hotspot discovery)
python wt_scan.py \
    --topol topol.top --tpr md.tpr --traj md.xtc \
    --index index.ndx --cg "1 13" --cr ref.pdb \
    --prefix MY_COMPLEX --top-n 20

# Phase 2: alanine scanning with cache
# First mutation does FULL (WT + mutant), subsequent mutations do CACHE (mutant only)
python vel_cache_runner.py
```
```

---

## ASPB-IE: Alanine Scanning Poisson Boltzmann — Interaction Entropy

> **GROMACS-native workflow.** Inputs are `topol.top`, `md.tpr`, `md.xtc`, `index.ndx`.
> All PB calculations run on **CPU via sander** (no GPU required).
> Requires an **existing MD trajectory** — the pipeline post-processes, it does not simulate.

**ASPB-IE** is a two-stage end-state free energy workflow for identifying and quantifying hotspot residues in protein–ligand complexes when no prior mutation data is available.

### Motivation

Traditional alanine scanning requires a predefined list of candidate residues. When the protein is large or no hotspot information exists, scanning every residue is computationally prohibitive, while random selection risks missing key contributors.

ASPB-IE solves this by using the wild-type (WT) trajectory itself to guide mutation selection.

### Two-Stage Workflow

```
Phase 1                          Phase 2
──────────                       ──────────
WT complex  ─→ PB decomposition  →  Ranked       ─→ Alanine scanning
(trajectory)    per-residue          hotspot           on hotspots
                energy              candidates        (normal_cache=1)
                                                    → ΔΔG per mutation
```

#### Phase 1 — WT Decomposition Scan

Run PB decomposition on the **wild-type** complex to obtain per-residue binding energy contributions:

```bash
python wt_scan.py \
    --topol topol.top --tpr md.tpr --traj md.xtc \
    --index index.ndx --cg "1 13" --cr ref.pdb \
    --prefix MY_COMPLEX --top-n 30 --threshold -1.0
```

This produces:

| File | Content |
|---|---|
| `MY_COMPLEX_per_residue.csv` | Full per-residue PB decomposition (all residues, ranked by TOTAL) |
| `MY_COMPLEX_hotspots.csv`    | Top N hotspot candidates (default: 20) |
| `MY_COMPLEX_hotspots.ala`    | Ready-to-use `ala.in` template for Phase 2 |

The script ranks residues by **TOTAL PB decomposition energy** (most negative = strongest binding contributor). Residues below `--threshold` (default -1.0 kcal/mol) or the top N are selected as candidate hotspots.

#### Phase 2 — Targeted Alanine Scanning

Use the generated `ala.in` template with the cache-aware runner:

```bash
python vel_cache_runner.py --ala MY_COMPLEX_hotspots.ala
```

The template already has `normal_cache=1` set, so the WT calculation runs only once (on the first mutation) and is reused for all subsequent mutations.

### Method Details

- **Poisson Boltzmann (PB)**: electrostatic solvation energy calculated with the PB model (`ipb=1`, `istrng=0.15`, `radiopt=0`).
- **Interaction Entropy (IE)**: configurational entropy contribution estimated via the IE method from the trajectory ensemble.
- **Decomposition**: per-residue energy decomposition (`idecomp=1`) decomposes the total PB binding energy into per-residue contributions.
- **WT caching**: `normal_cache=1` skips redundant WT sander calculations across mutations on the same trajectory (see cache documentation above).

### When to Use ASPB-IE

| Scenario | Recommendation |
|---|---|
| Large protein, no known hotspots | Full ASPB-IE: Phase 1 → Phase 2 |
| Known hotspots from literature | Skip Phase 1; directly configure `ala.in` for Phase 2 with `normal_cache=1` |
| Validation of existing hotspots | Run Phase 2 only with a targeted `ala.in` |
| Per-residue energy map desired | Phase 1 only (`--no-run` to generate input, or run and check `*_per_residue.csv`) |

### Example: Running the Full Pipeline

```bash
# Phase 1: scan WT to find hotspots
python wt_scan.py \
    --topol topol.top --tpr md.tpr --traj md.xtc \
    --index index.ndx --cg "1 13" --cr ref.pdb \
    --prefix my_protein --top-n 15

# Phase 2: alanine scan on found hotspots (using WT cache)
python vel_cache_runner.py --ala my_protein_hotspots.ala
```

---

## Install this fork

This version is **not** on PyPI/conda-forge. Install directly from source:

```bash
git clone https://github.com/lluoto/gmx_MMPBSA-ASPB-IE.git
cd gmx_MMPBSA-ASPB-IE
pip install .
```

To upgrade an existing conda environment:

```bash
conda activate gmxMMPBSA
pip install --upgrade --no-deps .
```

### Runner deployment

The `vel_cache_runner.py` and `wt_scan.py` scripts are **not part of the gmx_MMPBSA package**.
They are standalone workflow orchestrators. Deploy them to any work directory:

```bash
# Copy runner to your work directory
cp vel_cache_runner.py /path/to/project/
cp wt_scan.py /path/to/project/

# Edit these variables in vel_cache_runner.py:
#   BASE      = Path('/path/to/project')
#   SITES     = [21, 33, 68, 95]       # mutation residues
#   MPI_SIZE  = 5                       # match -np in Slurm

# Run
python vel_cache_runner.py
```

### Deployment status (em server)

| Component | Location | Status |
|-----------|----------|--------|
| `input_parser.py` (normal_cache) | site-packages/GMXMMPBSA/ | ✅ deployed, tested |
| `main.py` (skip normal + parse cache restore) | site-packages/GMXMMPBSA/ | ✅ deployed, tested |
| `calculation.py` (CopyCalc fallback) | site-packages/GMXMMPBSA/ | ✅ deployed, tested |
| `vel_cache_runner.py` | `/home/ajsali/lluoto/vel/` | ✅ deployed, job 25210991 completed |
| `v5_mmpbsa_cal.py` | `/home/ajsali/lluoto/3070_0082_rp2/` | ✅ deployed (pre-`normal_cache` era data) |
| `wt_scan.py` | `vel/` + `3070_0082_rp2/` | ✅ deployed, vel test passed |
| `mutation_list.slurm` (CPU-only) | `/home/ajsali/lluoto/` | ✅ deployed `-p cpu`, `GMXMMPBSA_HYBRID_DISABLE=1` |
