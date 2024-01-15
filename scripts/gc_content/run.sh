/home/mshokrof/sv/cattle-sv/datasets/references/gc_content

parallel --gnu -j 16 'seqtk gc -f 0.{} -l 100 ../ARS-UCD1.2_Btau5.0.1Y.fa  > l100_gc{}.bed' ::: 55 60 65 70 75 80 85

seqtk gc -wf 0.7 -l 100 ../ARS-UCD1.2_Btau5.0.1Y.fa  > l100_gc30.bed

seqtk gc -wf 0.75 -l 100 ../ARS-UCD1.2_Btau5.0.1Y.fa  > l100_gc25.bed

seqtk gc -wf 0.8 -l 100 ../ARS-UCD1.2_Btau5.0.1Y.fa  > l100_gc20.bed

seqtk gc -wf 0.85 -l 100 ../ARS-UCD1.2_Btau5.0.1Y.fa  > l100_gc15.bed


 ls *bed |parallel --gnu 'cat {} | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/'  > ARS_UCD1.2_{}'

 parallel --gnu "slopBed -i ARS_UCD1.2_l100_gc{}.bed -g ../ARS-UCD1.2_Btau5.0.1Y.genome -b 50 | awk '\$3>0' | mergeBed -i stdin > ARS_UCD1.2_l100_gc{}_slop50.bed" ::: 15 20 25 30 55 60 65 70 75 80 85





subtractBed -a  ARS_UCD1.2_l100_gc20_slop50.bed -b  ARS_UCD1.2_l100_gc15_slop50.bed >  ARS_UCD1.2_l100_gc15to20_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc25_slop50.bed -b  ARS_UCD1.2_l100_gc20_slop50.bed >  ARS_UCD1.2_l100_gc20to25_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc30_slop50.bed -b  ARS_UCD1.2_l100_gc25_slop50.bed >  ARS_UCD1.2_l100_gc25to30_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc80_slop50.bed -b  ARS_UCD1.2_l100_gc85_slop50.bed >  ARS_UCD1.2_l100_gc80to85_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc75_slop50.bed -b  ARS_UCD1.2_l100_gc80_slop50.bed >  ARS_UCD1.2_l100_gc75to80_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc70_slop50.bed -b  ARS_UCD1.2_l100_gc75_slop50.bed >  ARS_UCD1.2_l100_gc70to75_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc65_slop50.bed -b  ARS_UCD1.2_l100_gc70_slop50.bed >  ARS_UCD1.2_l100_gc65to70_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc60_slop50.bed -b  ARS_UCD1.2_l100_gc65_slop50.bed >  ARS_UCD1.2_l100_gc60to65_slop50.bed

subtractBed -a  ARS_UCD1.2_l100_gc55_slop50.bed -b  ARS_UCD1.2_l100_gc60_slop50.bed >  ARS_UCD1.2_l100_gc55to60_slop50.bed


subtractBed -a ../ARS-UCD1.2_Btau5.0.1Y.genome.bed -b ARS_UCD1.2_l100_gc55_slop50.bed | subtractBed -a stdin -b ARS_UCD1.2_l100_gc30_slop50
.bed > ARS_UCD1.2_l100_gc30to55_slop50.bed

multiIntersectBed -i  ARS_UCD1.2_l100_gc15_slop50.bed  ARS_UCD1.2_l100_gc15to20_slop50.bed  ARS_UCD1.2_l100_gc20to25_slop50.bed  ARS_UCD1.2_l100_gc65to70_slop50.bed  ARS_UCD1.2_l100_gc70to75_slop50.bed  ARS_UCD1.2_l100_gc75to80_slop50.bed  ARS_UCD1.2_l100_gc80to85_slop50.bed  ARS_UCD1.2_l100_gc85_slop50.bed | mergeBed -i stdin >  ARS_UCD1.2_l100_gclt25orgt65_slop50.bed

multiIntersectBed -i  ARS_UCD1.2_l100_gc15_slop50.bed  ARS_UCD1.2_l100_gc15to20_slop50.bed  ARS_UCD1.2_l100_gc20to25_slop50.bed  ARS_UCD1.2_l100_gc25to30_slop50.bed  ARS_UCD1.2_l100_gc55to60_slop50.bed  ARS_UCD1.2_l100_gc60to65_slop50.bed  ARS_UCD1.2_l100_gc65to70_slop50.bed  ARS_UCD1.2_l100_gc70to75_slop50.bed  ARS_UCD1.2_l100_gc75to80_slop50.bed  ARS_UCD1.2_l100_gc80to85_slop50.bed  ARS_UCD1.2_l100_gc85_slop50.bed | mergeBed -i stdin >  ARS_UCD1.2_l100_gclt30orgt55_slop50.bed


