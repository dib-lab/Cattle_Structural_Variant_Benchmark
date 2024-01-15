multiIntersectBed -i \
../mappility/lowmappabilityall.bed.gz \
 ../gc_content/ARS_UCD1.2_l100_gclt25orgt65_slop50.bed \
../low_complexity/ARS-UCD1.2_AllTandemRepeatsandHomopolymers_slop5.bed.gz \
../segment_duplication/ARS-UCD.1.2.dup.bed \
| sortBed -faidx ../ARS-UCD1.2_Btau5.0.1Y.fa.fai -i stdin \
| mergeBed -i stdin \
| bgzip -c > ARS-UCD1.2_alldifficultregions.bed.gz

subtractBed \
-a ../ARS-UCD1.2_Btau5.0.1Y.genome.bed \
-b ARS-UCD1.2_alldifficultregions.bed.gz \
| bgzip -c > ARS-UCD1.2_notinalldifficultregions.bed.gz
