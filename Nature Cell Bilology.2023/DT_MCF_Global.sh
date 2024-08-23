#!/usr/bin/env bash

bw_WT=$1
bw_KO=$2
name=$3
computeMatrix reference-point \
       --referencePoint center \
       -b 3000 -a 3000 \
       -R bed_elements_atac/MLL3_ATAC_promoter_sort.bed bed_elements_atac/MLL3_ATAC_intron_sort.bed bed_elements_atac/MLL3_ATAC_intergenic_sort.bed \
       -S $bw_WT $bw_KO \
       -p 4 --sortRegions keep \
       -o ${name}_elements.gz 
       
plotHeatmap -m ${name}_elements.gz \
--colorMap OrRd \
-out ${name}_elements.pdf

plotProfile -m ${name}_elements.gz \
              -out ${name}_elements_profile.pdf \
              --colors blue red \
              --perGroup   
