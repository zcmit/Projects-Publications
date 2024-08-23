#!/bin/bash 

macsDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/macs
outputDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/Bam_bw/Replicates
tmpDir=/media/b/Bigdata3/Chi/Mll3_SRA-download/Project_12784/tmp

#Peak calling on pseudoreplicates
echo "Calling peaks for replicate1 "
macs2 callpeak -t ${outputDir}/$1_WT1*.bam -c ${outputDir}/Input_WT1*.bam -f BAM -g hs -n $1_WT1 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir 

echo "Calling peaks for replicate2"
macs2 callpeak -t ${outputDir}/$1_WT2*.bam -c ${outputDir}/Input_WT2*.bam -f BAM -g hs -n $1_WT2 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir   

echo "Calling peaks for replicate3"
macs2 callpeak -t ${outputDir}/$1_WT3*.bam -c ${outputDir}/Input_WT3*.bam -f BAM -g hs -n $1_WT3 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir  


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/$1_WT1_peaks.narrowPeak > $macsDir/$1_WT1_sorted.narrowPeak
sort -k8,8nr $macsDir/$1_WT2_peaks.narrowPeak > $macsDir/$1_WT2_sorted.narrowPeak
sort -k8,8nr $macsDir/$1_WT3_peaks.narrowPeak > $macsDir/$1_WT3_sorted.narrowPeak

###### Mut
#Peak calling on pseudoreplicates
echo "Calling peaks for replicate1 "
macs2 callpeak -t ${outputDir}/$1_Mut1*.bam -c ${outputDir}/Input_Mut1*.bam -f BAM -g hs -n $1_Mut1 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir 

echo "Calling peaks for replicate2"
macs2 callpeak -t ${outputDir}/$1_Mut2*.bam -c ${outputDir}/Input_Mut2*.bam -f BAM -g hs -n $1_Mut2 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir   

echo "Calling peaks for replicate3"
macs2 callpeak -t ${outputDir}/$1_Mut3*.bam -c ${outputDir}/Input_Mut3*.bam -f BAM -g hs -n $1_Mut3 --nomodel --shift 37 --extsize 73 -p 1e-3 --outdir $macsDir  


#Sort peak by -log10(p-value)
echo "Sorting peaks..."
sort -k8,8nr $macsDir/$1_Mut1_peaks.narrowPeak > $macsDir/$1_Mut1_sorted.narrowPeak
sort -k8,8nr $macsDir/$1_Mut2_peaks.narrowPeak > $macsDir/$1_Mut2_sorted.narrowPeak
sort -k8,8nr $macsDir/$1_Mut3_peaks.narrowPeak > $macsDir/$1_Mut3_sorted.narrowPeak


# Remove the tmp directory
rm -r $tmpDir/*