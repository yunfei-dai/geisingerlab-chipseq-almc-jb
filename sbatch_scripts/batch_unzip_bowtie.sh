#!/bin/bash

## Yunfei Dai
## 11SEP2023

#SBATCH --array=0-9%10
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --job-name=unzip-bowtie
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH -o /work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/slurm_log/slurm_unzip-bowtie.%N.%j.out
#SBATCH -e /work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD/slurm_log/slurm_unzip-bowtie.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dai.yun@northeastern.edu

#Loading modules and venv
module load gcc/4.8.5
module load anaconda3/2021.05
module load bowtie/1.3.0
module load samtools/1.10
source activate /work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/YD_chipseq_venv

#Directories
DATA_DIR=/work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis
FASTQ_DIR=$DATA_DIR/201028-0304M_Edward_Geisinger_6095/fastq_Lane1
WORK_DIR=/work/geisingerlab/2023_Sep_ALMC_Jinna_ChIPSeq_analysis_YD
UNZIPPED=$WORK_DIR/Unzipped_fastq
CLIPPED=$WORK_DIR/Clipped_fastq
MAPPED_SAM=$WORK_DIR/Mapped_SAM
MAPPED_BAM=$WORK_DIR/Mapped_BAM
BigWig=$WORK_DIR/BigWig
REF=/work/geisingerlab/Yunfei/Ab_reference/Ab_all/Ab_all
CLIPPER=/work/geisingerlab/software/fastx_toolkit/bin/fastx_clipper
ADAPTOR=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

#Create output folders
mkdir -p $UNZIPPED
mkdir -p $CLIPPED
mkdir -p $MAPPED_SAM
mkdir -p $MAPPED_BAM
mkdir -p $BigWig

#Array of inputs
INPUT=( $(find $FASTQ_DIR -type f -name "*.gz") )
func () {
    gz_file=$1
    fname=$(basename $gz_file)
    echo "processing ${gz_file}"
    unzipped_fastq=$UNZIPPED/${fname%%.gz}
    clipped_fastq=$CLIPPED/${fname%%.fastq.gz}.clipped.fastq
    BAM_output=$MAPPED_SAM/${fname%%.fastq.gz}.bam
    BAM_sorted=$MAPPED_BAM/${fname%%.fastq.gz}.sorted.bam
    BW=$BigWig/${fname%%.fastq.gz}.bw
    gzip -cd $gz_file > $unzipped_fastq
    $CLIPPER -v -l 20 -d 0 -Q 33 -a $ADAPTOR -i $unzipped_fastq -o $clipped_fastq
    bowtie -q -m 1 -n 1 --best -y -S $REF $clipped_fastq $SAM_output | samtools view -bS - > $BAM_output
    samtools sort $BAM_output -o $BAM_sorted
    samtools index $BAM_sorted
    bamCoverage -b $BAM_sorted -o $BW --normalizeUsing CPM --binSize 20 # 16Jan2023 updated binsize and normalizing method
}

#Lauch job array
func "${INPUT[$SLURM_ARRAY_TASK_ID]}"
