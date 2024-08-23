#!/bin/bash 

macsDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/macs
outputDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/Bam_bw/Replicates
tmpDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/tmp

#Peak calling on pseudoreplicates
echo "Calling peaks for replicate1 "
macs2 callpeak -t ${outputDir}/$1_WT1*.bam -c ${outputDir}/Input_WT1*.bam -f BAM -g hs -n $1_WT1 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir 

echo "Calling peaks for replicate2"
macs2 callpeak -t ${outputDir}/$1_WT2*.bam -c ${outputDir}/Input_WT2*.bam -f BAM -g hs -n $1_WT2 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir   

echo "Calling peaks for replicate3"
macs2 callpeak -t ${outputDir}/$1_WT3*.bam -c ${outputDir}/Input_WT3*.bam -f BAM -g hs -n $1_WT3 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir  


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/$1_WT1_peaks.broadPeak > $macsDir/$1_WT1_sorted.broadPeak
sort -k8,8nr $macsDir/$1_WT2_peaks.broadPeak > $macsDir/$1_WT2_sorted.broadPeak
sort -k8,8nr $macsDir/$1_WT3_peaks.broadPeak > $macsDir/$1_WT3_sorted.broadPeak

###### Mut
#Peak calling on pseudoreplicates
macs2 callpeak -t ${outputDir}/$1_Mut1*.bam -c ${outputDir}/Input_Mut1*.bam -f BAM -g hs -n $1_Mut1 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir 

macs2 callpeak -t ${outputDir}/$1_Mut2*.bam -c ${outputDir}/Input_Mut2*.bam -f BAM -g hs -n $1_Mut2 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir   

macs2 callpeak -t ${outputDir}/$1_Mut3*.bam -c ${outputDir}/Input_Mut3*.bam -f BAM -g hs -n $1_Mut3 --broad --broad-cutoff 0.1 -p 1e-3 --outdir $macsDir  


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/$1_Mut1_peaks.broadPeak > $macsDir/$1_Mut1_sorted.broadPeak
sort -k8,8nr $macsDir/$1_Mut2_peaks.broadPeak > $macsDir/$1_Mut2_sorted.broadPeak
sort -k8,8nr $macsDir/$1_Mut3_peaks.broadPeak > $macsDir/$1_Mut3_sorted.broadPeak
