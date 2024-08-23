### H3k4me1, H3k4me3 and H3K27ac, WT and MLL3 Mut MDA-MB-231 chip seq ###

### download file and unzip
unzip downloaded file from dropbox

https://www.dropbox.com/s/ly3qu1eomidm9l4/Project_12784.zip?dl=0. 


#=== Merge files from two lanes ===

echo {Input_,H3K4ME1_,H3K4ME3_,H3K27ac_}{WT,Mut}{1,2,3} >> SampleName
SampleName="Input_WT H3K4ME1_WT H3K4ME3_WT H3K27ac_WT"
SampleName={Input_,H3K4ME1_,H3K4ME3_,H3K27ac_}{WT,Mut}{1,2,3}


find . -name 'H3K4ME1_Mut1*_R1*' | xargs -P 4 cat {} > H3K4ME1_Mut1_R1.fastq.gz
find . -name 'H3K4ME1_Mut2*_R1*' | xargs -P 4 cat {} > H3K4ME1_Mut2_R1.fastq.gz
find . -name 'H3K4ME1_Mut3*_R1*' | xargs -P 4 cat {} > H3K4ME1_Mut3_R1.fastq.gz

find . -name 'H3K4ME1_Mut1*_R2*' | xargs -P 4 cat {} > H3K4ME1_Mut1_R2.fastq.gz
find . -name 'H3K4ME1_Mut2*_R2*' | xargs -P 4 cat {} > H3K4ME1_Mut2_R2.fastq.gz
find . -name 'H3K4ME1_Mut3*_R2*' | xargs -P 4 cat {} > H3K4ME1_Mut3_R2.fastq.gz


find . -name 'Input_Mut1*_R1*' | xargs -P 4 cat {} > Input_Mut1_R1.fastq.gz
find . -name 'Input_Mut2*_R1*' | xargs -P 4 cat {} > Input_Mut2_R1.fastq.gz
find . -name 'Input_Mut3*_R1*' | xargs -P 4 cat {} > Input_Mut3_R1.fastq.gz

find . -name 'Input_Mut1*_R2*' | xargs -P 4 cat {} > Input_Mut1_R2.fastq.gz
find . -name 'Input_Mut2*_R2*' | xargs -P 4 cat {} > Input_Mut2_R2.fastq.gz
find . -name 'Input_Mut3*_R2*' | xargs -P 4 cat {} > Input_Mut3_R2.fastq.gz


for sample in {Input_,H3K4ME1_,H3K4ME3_,H3K27ac_}{WT,Mut}{1,2,3}
do 
find . -name ${sample}_R1*
done


#=== Fastqc ===

fastqc -t 4 -o FastQC *.fastq.gz

#=== chipseq workflow ===

#!/bin/bash
set -e
set -u
set -o pipefail


cd /media/chizhang/Bigdata3/Chi/Mll3_SRA-download/Project_12784

conda activate encode-chip-seq-pipeline
# specify the input samples file, where the third column is the path to each sample FASTQ file
sample_info=samples.txt


REF="/media/Data/General_References/Genomes/hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"


sample_names=($(cut -f 1 "$sample_info" | uniq))

for sample in ${sample_names[@]} 
do
	
	bwa mem -M $REF Raw/${sample}_R1.fastq.gz Raw/${sample}_R2.fastq.gz > ${sample}.sam
	
	samtools view -bS ${sample}.sam > ${sample}.bam
	samtools sort -T /tmp/${sample}.sorted -o ${sample}_sorted.bam ${sample}.bam
	samtools index ${sample}_sorted.bam
	
	picard -Xmx4g MarkDuplicates \
	I=${sample}_sorted.bam \
	O=${sample}_mdup.bam \
	M=${sample}_mdup_metrics.txt \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
	TMP_DIR=/tmp

	samtools sort -T /tmp/${sample}.sorted -o ${sample}_mdup.sorted.bam ${sample}_mdup.bam
	samtools index ${sample}_mdup.sorted.bam
	bamCoverage -b ${sample}_mdup.sorted.bam \
	-o ${sample}.bw \
	-of bigwig -bs 30 -p 6 \
	--smoothLength 300 \
	--extendReads 200 \
	--normalizeUsing RPKM  
done



for sample in ${sample_names[@]} 
do
    echo $sample
done

#=== Macs2 call peaks ===
bash IDR_broad.sh H3K27ac
bash IDR_narrow.sh H3K4ME3
bash IDR_broad.sh H3K4ME1 


#====== To profile globally ========# 
# Merge bam files. Mut2&3, WT1&2
samtools merge H3K27ac_Mut.bam H3K27ac_Mut2_mdup.sorted.bam H3K27ac_Mut3_mdup.sorted.bam
samtools merge H3K4me1_Mut.bam H3K4ME1_Mut2_mdup.sorted.bam H3K4ME1_Mut3_mdup.sorted.bam
samtools merge H3K4me3_Mut.bam H3K4ME3_Mut2_mdup.sorted.bam H3K4ME3_Mut3_mdup.sorted.bam

samtools merge H3K27ac_WT.bam H3K27ac_WT1_mdup.sorted.bam H3K27ac_WT2_mdup.sorted.bam
samtools merge H3K4me1_WT.bam H3K4ME1_WT1_mdup.sorted.bam H3K4ME1_WT2_mdup.sorted.bam 
samtools merge H3K4me3_WT.bam H3K4ME3_WT1_mdup.sorted.bam H3K4ME3_WT2_mdup.sorted.bam 

# sort and index
SampleName="H3K4me3 H3K4me1 H3K27ac"

for sample in $SampleName
do

	samtools sort -T /tmp/${sample}_WT.sorted -o ${sample}_WT.sorted.bam ${sample}_WT.bam
	samtools index ${sample}_WT.sorted.bam
	samtools sort -T /tmp/${sample}_Mut.sorted -o ${sample}_Mut.sorted.bam ${sample}_Mut.bam
	samtools index ${sample}_Mut.sorted.bam
	
done

# convert merged bam to bigwig

for sample in $SampleName
do
	bamCoverage -b ${sample}_WT.sorted.bam \
	-o ${sample}_WT.bw \
	-of bigwig -bs 30 -p 6 \
	--smoothLength 300 \
	--extendReads 150 \
	--normalizeUsing RPKM  

	bamCoverage -b ${sample}_Mut.sorted.bam \
	-o ${sample}_Mut.bw \
	-of bigwig -bs 30 -p 6 \
	--smoothLength 300 \
	--extendReads 150 \
	--normalizeUsing RPKM  
done

### ROSE for H3K27ac

cd /media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/ROSE

BamDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/Bam_bw/
PeakDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/macs/
WorkDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/ROSE

SampleName="H3K27ac_WT1 H3K27ac_WT2 H3K27ac_Mut2 H3K27ac_Mut3"
for sample in $SampleName
do
	awk '{print $1"\t"$4"\t"".""\t"$2"\t"$3"\t"".""\t"".""\t"".""\t"$4}' $PeakDir/${sample}_peaks.broadPeak > ${sample}_peaks.gff
done

cd /media/b/Bigdata3/Chi/Mll3_SRA-download/ROSE-main

python ROSE_main.py -g hg38 -i $WorkDir/H3K27ac_WT1_peaks.gff -r $BamDir/Replicates/H3K27ac_WT1_mdup.sorted.bam \
-c $BamDir/Replicates/Input_WT1_mdup.sorted.bam -o H3K27ac_WT1 -s 12500 -t 2000
python ROSE_geneMapper.py -i H3K27ac_WT1/H3K27ac_WT1_peaks_AllEnhancers.table.txt -g hg38 -o H3K27ac_WT1

RepNames="Mut2 WT2 Mut3"
for id in $RepNames
do
# python ROSE_main.py -g hg38 -i $WorkDir/H3K27ac_${id}_peaks.gff -r $BamDir/Replicates/H3K27ac_${id}_mdup.sorted.bam \
# -c $BamDir/Replicates/Input_${id}_mdup.sorted.bam -o H3K27ac_${id} -s 12500 -t 2000
python ROSE_geneMapper.py -i H3K27ac_${id}/H3K27ac_${id}_peaks_AllEnhancers.table.txt -g hg38 -o H3K27ac_${id}
done



