#!/usr/bin/env bash

bw_WT=$1
bw_KO=$2
name=$3

computeMatrix scale-regions -S $bw_WT $bw_KO \
-R ROSE_bed/MDA231_SE.bed \
--beforeRegionStartLength 500 \
--regionBodyLength 5000 \
--afterRegionStartLength 500 \
--skipZeros -p 4 \
-o MDA231_${name}_SE.gz

plotProfile -m MDA231_${name}_SE.gz \
              -out MDA231_${name}_SE.pdf \
              --plotTitle "MDA231 Super Enhancers" \
              --perGroup

plotHeatmap -m MDA231_${name}_SE.gz \
--colorMap RdBu_r \
--plotTitle "MDA231 Super Enhancers" \
-out MDA231_${name}_SE_heatmap.pdf 
