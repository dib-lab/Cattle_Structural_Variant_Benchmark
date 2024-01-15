/home/mshokrof/sv/cattle-sv/datasets/references/segment_duplication

minimap2 -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -t 32 ../ARS-UCD1.2_Btau5.0.1Y.fa ../ARS-UCD1.2_Btau5.0.1Y.fa   | samtools sort -m10G -o {output} - 


bedtools genomecov -bga -ibam ARS-UCD.12.bam | awk '$4 < 2' | bedtools merge > ARS-UCD.12.unique.bed

bedtools bamtobed -i ARS-UCD.12.bam | awk '($3-$2) >= 50000' | bedtools merge > ARS-UCD.12.covered.bed

bedtools intersect -a ARS-UCD.12.covered.bed -b ARS-UCD.12.unique.bed > ARS-UCD.12.callable.bed
