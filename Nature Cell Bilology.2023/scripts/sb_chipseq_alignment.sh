#!/bin/bash

#xusheng 4-01-22 sb_chipseq_pipeline_JCui used template from rob code below 
#Jihong Cui Chipseq
#dubin 2-20-22 sb_alignToHg38 Dan Weiser ultralow

#SBATCH --partition normal
#SBATCH --job-name=alignToHg38
#SBATCH --output=di_alignToHg38.out.%j
#SBATCH --mem=15G
#SBATCH -n 4
#SBATCH --time=10:0:0

SCRIPT="sb_chipseq_alignToHg38_JCui"
DATE=`date +%Y-%m-%d`

prefixFilename=$1
TRIMMED_FASTQ_FILE_READ1=$2
TRIMMED_FASTQ_FILE_READ2=$3

echo
echo "${SCRIPT} run on ${DATE}"
echo

																	
GENOME_SELECTED="hg38"
INDEXED_GENOME="/gs/gsfs0/users/svc_wasp/wasp/metadata-new/Hsa/hg38/ncbi_no_alt_analysis_set/bwa/${GENOME_SELECTED}"
REFERENCE_FASTA="/gs/gsfs0/users/svc_wasp/wasp/metadata-new/Hsa/hg38/ncbi_no_alt_analysis_set/fasta/hg38_no_alt.fa"

echo
echo "prefixFilename = ${prefixFilename}"
echo "TRIMMED_FASTQ_FILE_READ1 = ${TRIMMED_FASTQ_FILE_READ1}"
echo "TRIMMED_FASTQ_FILE_READ2 = ${TRIMMED_FASTQ_FILE_READ2}"
echo
echo "GENOME_SELECTED = ${GENOME_SELECTED}"
echo "INDEXED_GENOME = ${INDEXED_GENOME}"
echo "REFERENCE_FASTA = ${REFERENCE_FASTA}"
echo

source ~/.bashrc; module load conda3; module list; echo; conda activate biobase; PATH=/public/apps/conda3/envs/biobase/bin:$PATH;

echo "bwa info: "
which bwa
echo
bwa   #spits out info that includes version (and more)
echo

echo
echo "starting bwa alginment"
echo

##step1 alignment with BWA mem#####################

##bwa
bwa mem -M ${INDEXED_GENOME} ${TRIMMED_FASTQ_FILE_READ1} ${TRIMMED_FASTQ_FILE_READ2} > ${prefixFilename}.bwa.${GENOME_SELECTED}.sam

echo
echo "should be all done with bwa align for ${TRIMMED_FASTQ_FILE_READ1} : ${TRIMMED_FASTQ_FILE_READ2} "
echo

echo
echo "samtools info:"
which samtools
echo 
samtools --version
echo

samtools view -bS ${prefixFilename}.bwa.${GENOME_SELECTED}.sam > ${prefixFilename}.bwa.${GENOME_SELECTED}.bam
samtools sort -T /tmp/${prefixFilename}.sorted -o ${prefixFilename}.bwa.${GENOME_SELECTED}.sorted.bam ${prefixFilename}.bwa.${GENOME_SELECTED}.bam
samtools ${prefixFilename}.bwa.${GENOME_SELECTED}.sorted.bam


echo

echo "picard MarkDuplicates info: "
which picard
echo
picard MarkDuplicates -version 
echo

echo
echo "starting MarkDuplicates : mark not remove"
echo


##mark duplicates
picard -Xmx4g MarkDuplicates \
I=${prefixFilename}.bwa.${GENOME_SELECTED}.sorted.bam \
O=${prefixFilename}.bwa.${GENOME_SELECTED}.mdup.bam \
M=${prefixFilename}_mdup_metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
TMP_DIR=/tmp

samtools sort -T /tmp/${prefixFilename}.sorted -o ${prefixFilename}.bwa.${GENOME_SELECTED}.mdup.sorted.bam ${prefixFilename}.bwa.${GENOME_SELECTED}.mdup.bam
samtools index ${prefixFilename}.bwa.${GENOME_SELECTED}.mdup.sorted.bam

echo

echo "deeptools bamCoverage info: "
which deeptools
echo
bamCoverage -version 
echo



##calculate coverage
bamCoverage -b ${prefixFilename}.bwa.${GENOME_SELECTED}.mdup.sorted.bam \
-o ${prefixFilename}.bw \
-of bigwig -bs 30 -p 6 \
--smoothLength 300 \
--extendReads 200 \
--normalizeUsing RPKM  
