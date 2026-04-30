#! /bin/bash

#SBATCH -J 0082_hybrid
#SBATCH -p gpu-c
#SBATCH -t 4-23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --output=/home/ajsali/lluoto/3070_0082_rp2/mutation_list.%j.log
#SBATCH --error=/home/ajsali/lluoto/3070_0082_rp2/mutation_list.%j.log

date

WORKDIR=/home/ajsali/lluoto/3070_0082_rp2

source /etc/profile 2>/dev/null
module purge
module load gcc/10.5.0
module load cuda/11.8.0
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

export GMXMMPBSA_HYBRID_PMEMD=/home/ajsali/lluoto/amber22/bin/pmemd.cuda
export GMXMMPBSA_HYBRID_CPPTRAJ=/home/ajsali/.conda/envs/gmxMMPBSA/bin/cpptraj

cd $WORKDIR
nvidia-smi || true
which python || true
which sander || true
python -c "import MMPBSA_mods.calculation as c; print(c.PBEnergyCalculation._DEFAULT_HYBRID_PMEMD)"
python mmpbsa_cal.py

hostname
date
