#!/bin/bash

## Yunfei Dai
## 20Sep2023

#SBATCH --time=08:00:00
#SBATCH --job-name=macs2
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/slurm_log/slurm_macs2.%N.%j.out
#SBATCH -e /work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/slurm_log/slurm_macs2.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

module load anaconda3/2021.05
source activate YD_chipseq_venv

SCRIPT=/work/geisingerlab/YD_GitProjects/geisingerlab-chipseq-almc-jb/sbatch_scripts/batch_run_macs2.py
INDIR=/work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/Mapped_BAM
OUTDIR=/work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/MACS2_output

mkdir -p ${OUTDIR}/narrowPeak ${OUTDIR}/xls ${OUTDIR}/bed

python $SCRIPT $INDIR $OUTDIR

mv ${OUTDIR}/*.narrowPeak ${OUTDIR}/narrowPeak
mv ${OUTDIR}/*.xls ${OUTDIR}/xls
mv ${OUTDIR}/*.bed ${OUTDIR}/bed

