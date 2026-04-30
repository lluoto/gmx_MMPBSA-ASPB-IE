# Hybrid remote deployment artifacts

This branch stores the exact remote deployment artifacts used on em for the hybrid mutant-PB experiment.

## Remote target
- Server alias: em
- Main working directory: /home/ajsali/lluoto
- Job workdir: /home/ajsali/lluoto/3070_0082_rp2
- Installed runtime patched in place: /home/ajsali/.conda/envs/gmxMMPBSA/lib/python3.9/site-packages/MMPBSA_mods/calculation.py
- Backup created: /home/ajsali/.conda/envs/gmxMMPBSA/lib/python3.9/site-packages/MMPBSA_mods/calculation.py.BACKUP

## Hybrid runtime assumptions
- GPU minimizer: /home/ajsali/lluoto/amber22/bin/pmemd.cuda
- Required compiler runtime: gcc/10.5.0
- Required CUDA runtime: cuda/11.8.0
- PB CPU executor remains sander from the gmxMMPBSA env.

## What is uploaded here
- calculation.remote.patch: unified diff between the original and patched legacy MMPBSA_mods/calculation.py
- mutation_list.remote.slurm: the Slurm submission script patched for the remote hybrid test

## Test submission
- Submitted with: sbatch /home/ajsali/lluoto/mutation_list.slurm
- Job id: 24774426
- Observed state at upload time: PD (Priority)

## Scope note
This branch records the **actual remote legacy runtime deployment** that was patched and tested.
It does **not** claim that the same diff applies unchanged to the current upstream GMXMMPBSA/ source tree.
