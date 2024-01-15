#parallel --gnu --bar -j3 'python findSimpleRegions_quad.py -p {} -d 100000 -t 100000 -q 100000 ../references/ARS-UCD1.2_Btau5.0.1Y.fa  > ./ARS-UCD1.2_SimpleRepeat_p{}.bed' ::: 3 6 11 20

#parallel --gnu --bar -j3 'python findSimpleRegions_quad.py -p 100000 -d {} -t 100000 -q 100000 ../references/ARS-UCD1.2_Btau5.0.1Y.fa  > ./ARS-UCD1.2_SimpleRepeat_d{}.bed' ::: 11 51 201

#parallel --gnu --bar -j3 'python findSimpleRegions_quad.py -p 100000 -d 100000 -t {} -q 100000 ../references/ARS-UCD1.2_Btau5.0.1Y.fa  > ./ARS-UCD1.2_SimpleRepeat_t{}.bed' ::: 15 51 201

#parallel --gnu --bar -j3 'python findSimpleRegions_quad.py -p 100000 -d 100000 -t 100000 -q {} ../references/ARS-UCD1.2_Btau5.0.1Y.fa  > ./ARS-UCD1.2_SimpleRepeat_q{}.bed' ::: 20 51 201

cat ../references/ARS-UCD1.2_Btau5.0.1Y.fa.fai | cut -f 1,2 | grep -v "NKLS" | grep -v "MT"  > ../references/ARS-UCD1.2_Btau5.0.1Y.genome
awk -v OFS='\t' {'print $1,"0",$2'} ../references/ARS-UCD1.2_Btau5.0.1Y.fa.fai| grep -v "NKLS" | grep -v "MT"   > ../references/ARS-UCD1.2_Btau5.0.1Y.genome.bed


subtractBed -a ARS-UCD1.2_SimpleRepeat_p3.bed -b ARS-UCD1.2_SimpleRepeat_p6.bed  | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_homopolymer_4to6.bed.gz

subtractBed -a ARS-UCD1.2_SimpleRepeat_p6.bed -b ARS-UCD1.2_SimpleRepeat_p11.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_homopolymer_7to11.bed.gz

sed 's/^chr//' ARS-UCD1.2_SimpleRepeat_p11.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_homopolymer_gt11.bed.gz

sed 's/^chr//' ARS-UCD1.2_SimpleRepeat_p20.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_homopolymer_gt20.bed.gz



subtractBed -a ARS-UCD1.2_SimpleRepeat_d11.bed -b ARS-UCD1.2_SimpleRepeat_d51.bed  | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_diTR_11to50.bed.gz
subtractBed -a ARS-UCD1.2_SimpleRepeat_d51.bed -b ARS-UCD1.2_SimpleRepeat_d201.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_diTR_51to200.bed.gz
cat ARS-UCD1.2_SimpleRepeat_d201.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_diTR_gt200.bed.gz

subtractBed -a ARS-UCD1.2_SimpleRepeat_t15.bed -b ARS-UCD1.2_SimpleRepeat_t51.bed  | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_triTR_15to50.bed.gz
subtractBed -a ARS-UCD1.2_SimpleRepeat_t51.bed -b ARS-UCD1.2_SimpleRepeat_t201.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_triTR_51to200.bed.gz
cat  ARS-UCD1.2_SimpleRepeat_t201.bed | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_triTR_gt200.bed.gz

subtractBed -a ARS-UCD1.2_SimpleRepeat_q20.bed -b ARS-UCD1.2_SimpleRepeat_q51.bed | grep -v "NKLS" | grep -v "MT" | sed -e 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed -e 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_quadTR_20to50.bed.gz
subtractBed -a ARS-UCD1.2_SimpleRepeat_q51.bed -b ARS-UCD1.2_SimpleRepeat_q201.bed | grep -v "NKLS" | grep -v "MT" | sed  -e 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed -e 's/^30/X/;s/^31/Y/'| bgzip -c > ARS-UCD1.2_SimpleRepeat_quadTR_51to200.bed.gz
cat ARS-UCD1.2_SimpleRepeat_q201.bed  | grep -v "NKLS" | grep -v "MT" | sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | bgzip -c > ARS-UCD1.2_SimpleRepeat_quadTR_gt200.bed.gz

slopBed -i ARS-UCD1.2_SimpleRepeat_homopolymer_4to6.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_4to6_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_homopolymer_7to11.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_7to11_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_homopolymer_gt11.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_gt11_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_homopolymer_gt20.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_gt20_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_diTR_11to50.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_11to50_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_diTR_51to200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_51to200_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_diTR_gt200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_gt200_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_triTR_15to50.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_15to50_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_triTR_51to200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_51to200_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_triTR_gt200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_gt200_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_quadTR_20to50.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_20to50_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_quadTR_51to200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_51to200_slop5_withUNITS.bed.gz
slopBed -i ARS-UCD1.2_SimpleRepeat_quadTR_gt200.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_gt200_slop5_withUNITS.bed.gz

zcat ARS-UCD1.2_SimpleRepeat_homopolymer_4to6_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_4to6_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_homopolymer_7to11_slop5_withUNITS.bed.gz | cut -f1-3  | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_7to11_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_homopolymer_gt11_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_gt11_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_homopolymer_gt20_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_homopolymer_gt20_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_diTR_11to50_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_11to50_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_diTR_51to200_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_51to200_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_diTR_gt200_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_diTR_gt200_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_triTR_15to50_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_15to50_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_triTR_51to200_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_51to200_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_triTR_gt200_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_triTR_gt200_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_quadTR_20to50_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_20to50_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_quadTR_51to200_slop5_withUNITS.bed.gz| cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_51to200_slop5.bed.gz
zcat ARS-UCD1.2_SimpleRepeat_quadTR_gt200_slop5_withUNITS.bed.gz | cut -f1-3 | bgzip -c  > ARS-UCD1.2_SimpleRepeat_quadTR_gt200_slop5.bed.gz

grep 'unit=C' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>10' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_C.bed
grep 'unit=G' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>10' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_G.bed
grep 'unit=A' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>10' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_A.bed
grep 'unit=T' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>10' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_T.bed

multiIntersectBed -i ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_C.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_G.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_A.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_T.bed | 
	cut -f1-3 | grep -v "NKLS" | grep -v "MT" |
	slopBed -i stdin -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome |
	 sed 's/^X/30/;s/^Y/31/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
	bgzip -c > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz


grep 'unit=C' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>20' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_C.bed
grep 'unit=G' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>20' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_G.bed
grep 'unit=A' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>20' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_A.bed
grep 'unit=T' ARS-UCD1.2_SimpleRepeat_p3.bed | mergeBed -i stdin -d 1 | awk '$3-$2>20' > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_T.bed

multiIntersectBed -i ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_C.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_G.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_A.bed \
	ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_T.bed |
	cut -f1-3 | grep -v "NKLS" | grep -v "MT" |
	slopBed -i stdin -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome |
	sed 's/^chr//' |
	sed 's/^X/30/;s/^Y/31/' |
	sort -k1,1n -k2,2n -k3,3n |
	sed 's/^30/X/;s/^31/Y/' |
	mergeBed -i stdin |
	bgzip -c > ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt20_slop5.bed.gz

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau9/database/rmsk.txt.gz ./

# SIMPLE REPEATS
date
zgrep Simple_repeat rmsk.txt.gz |
	awk '{ print $6 "\t" $7 "\t" $8 ; }' |
	sed 's/^chr//' |
	grep -v "NKLS" | grep -v "MT" |grep -v "NW_"|
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
	bgzip -c > ARS-UCD1.2_rmsk_Simple_repeat.bed.gz
    
gunzip -c ARS-UCD1.2_rmsk_Simple_repeat.bed.gz |
	awk '$3-$2<51' | 
	bgzip -c > ARS-UCD1.2_rmsk_Simple_repeat_lt51.bed.gz
    
gunzip -c ARS-UCD1.2_rmsk_Simple_repeat.bed.gz |
	awk '$3-$2>50 && $3-$2<201' | 
	bgzip -c > ARS-UCD1.2_rmsk_Simple_repeat_51to200.bed.gz
    
gunzip -c ARS-UCD1.2_rmsk_Simple_repeat.bed.gz |
	awk '$3-$2>200' | 
	bgzip -c > ARS-UCD1.2_rmsk_Simple_repeat_gt200.bed.gz

date
zgrep Low_complexity rmsk.txt.gz |
	awk '{ print $6 "\t" $7 "\t" $8 ; }' |
	sed 's/^chr//' |
	grep -v "NKLS" | grep -v "MT" |grep -v "NW_"|
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
    bgzip -c > ARS-UCD1.2_rmsk_Low_complexity.bed.gz
    
gunzip -c ARS-UCD1.2_rmsk_Low_complexity.bed.gz |
	awk '$3-$2<51' | 
	bgzip -c > ARS-UCD1.2_rmsk_Low_complexity_lt51.bed.gz

gunzip -c ARS-UCD1.2_rmsk_Low_complexity.bed.gz |
	awk '$3-$2>50 && $3-$2<201' | 
	bgzip -c > ARS-UCD1.2_rmsk_Low_complexity_51to200.bed.gz

gunzip -c ARS-UCD1.2_rmsk_Low_complexity.bed.gz |
	awk '$3-$2>200' | 
	bgzip -c > ARS-UCD1.2_rmsk_Low_complexity_gt200.bed.gz

date
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau9/database/simpleRepeat.txt.gz ./
zcat simpleRepeat.txt.gz |
	cut -f2-4  |
	sed 's/^chr//' |
	grep -v "NKLS" | grep -v "MT"|grep -v "NW_" |
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
    bgzip -c > ARS-UCD1.2_trf_simpleRepeat.bed.gz

gunzip -c ARS-UCD1.2_trf_simpleRepeat.bed.gz |
	awk '$3-$2<51' | 
	bgzip -c > ARS-UCD1.2_trf_simpleRepeat_lt51.bed.gz

gunzip -c ARS-UCD1.2_trf_simpleRepeat.bed.gz |
	awk '$3-$2>50 && $3-$2<201' | 
	bgzip -c > ARS-UCD1.2_trf_simpleRepeat_51to200.bed.gz

gunzip -c ARS-UCD1.2_trf_simpleRepeat.bed.gz |
	awk '$3-$2>200' | 
	bgzip -c > ARS-UCD1.2_trf_simpleRepeat_gt200.bed.gz

zcat rmsk.txt.gz \
| grep "Satellite" \
| awk '{ print $6 "\t" $7 "\t" $8 ; }' \
| grep -Ev '^chr[0-9XYM]_|^chr[0-9][0-9]_|^chrUn_' \
| sed 's/^chr//' \
| sortBed -faidx ../references/ARS-UCD1.2_Btau5.0.1Y.fa.fai -i stdin \
| mergeBed -i stdin \
| bgzip -c > ARS-UCD1.2_satellites.bed.gz

slopBed -i ARS-UCD1.2_satellites.bed.gz -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome \
| mergeBed -i stdin \
| bgzip -c  > ARS-UCD1.2_satellites_slop5.bed.gz

zcat ARS-UCD1.2_satellites_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

subtractBed -a ../references/ARS-UCD1.2_Btau5.0.1Y.genome.bed -b ARS-UCD1.2_satellites_slop5.bed.gz | bgzip -c > ARS-UCD1.2_notinsatellites_slop5.bed.gz
grep -Ev '_|^chrEBV' ARS-UCD1.2_SimpleRepeat_p6.bed | grep -v "NKLS" | grep -v "MT"|grep -v "NW_"  |
	slopBed -i stdin -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome |
	mergeBed -i stdin |
	multiIntersectBed -i stdin ARS-UCD1.2_SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz |
	sed 's/^chr//' |
	cut -f1-3 | grep -v "NKLS" | grep -v "MT"|grep -v "NW_"  |
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
	bgzip -c > ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz

subtractBed -a ../references/ARS-UCD1.2_Btau5.0.1Y.genome.bed -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | bgzip -c > ARS-UCD1.2_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz

multiIntersectBed -i ARS-UCD1.2_SimpleRepeat_d11.bed \
	ARS-UCD1.2_SimpleRepeat_t15.bed \
	ARS-UCD1.2_SimpleRepeat_q20.bed \
	ARS-UCD1.2_rmsk_Simple_repeat.bed.gz \
	ARS-UCD1.2_rmsk_Low_complexity.bed.gz \
	ARS-UCD1.2_trf_simpleRepeat.bed.gz \
	ARS-UCD1.2_satellites.bed.gz | 
	sed 's/^chr//' | 
	cut -f1-3 | grep -v "NKLS" | grep -v "M" | 
	slopBed -i stdin -b 5 -g ../references/ARS-UCD1.2_Btau5.0.1Y.genome |  
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin > ARS-UCD1.2_AllTandemRepeats_intermediate.bed
    
echo "sum of regions:"
cat ARS-UCD1.2_AllTandemRepeats_intermediate.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

awk '$3-$2<61' ARS-UCD1.2_AllTandemRepeats_intermediate.bed | 
	subtractBed -a stdin -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeats_lt51bp_slop5.bed.gz
    
echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeats_lt51bp_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
awk '$3-$2>60 && $3-$2<211' ARS-UCD1.2_AllTandemRepeats_intermediate.bed | 
	subtractBed -a stdin -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeats_51to200bp_slop5.bed.gz
    
echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeats_51to200bp_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
awk '$3-$2>210 && $3-$2<10011' ARS-UCD1.2_AllTandemRepeats_intermediate.bed | 
	subtractBed -a stdin -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeats_201to10000bp_slop5.bed.gz
    
echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeats_201to10000bp_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
awk '$3-$2>10010' ARS-UCD1.2_AllTandemRepeats_intermediate.bed | 
	subtractBed -a stdin -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeats_gt10000bp_slop5.bed.gz
    
echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeats_gt10000bp_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
awk '$3-$2>110' ARS-UCD1.2_AllTandemRepeats_intermediate.bed | 
	subtractBed -a stdin -b ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeats_gt100bp_slop5.bed.gz

echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeats_gt100bp_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
multiIntersectBed -i ARS-UCD1.2_AllTandemRepeats_lt51bp_slop5.bed.gz \
	ARS-UCD1.2_AllTandemRepeats_51to200bp_slop5.bed.gz \
	ARS-UCD1.2_AllTandemRepeats_201to10000bp_slop5.bed.gz \
	ARS-UCD1.2_AllTandemRepeats_gt10000bp_slop5.bed.gz \
	ARS-UCD1.2_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz | 
	sed 's/^chr//' | 
	cut -f1-3 | grep -v "NKLS" | grep -v "MT" | 
	sed 's/^X/30/;s/^Y/31/' | 
	sort -k1,1n -k2,2n -k3,3n | 
	sed 's/^30/X/;s/^31/Y/' | 
	mergeBed -i stdin | 
	bgzip -c > ARS-UCD1.2_AllTandemRepeatsandHomopolymers_slop5.bed.gz

subtractBed -a ../references/ARS-UCD1.2_Btau5.0.1Y.genome.bed \
-b ARS-UCD1.2_AllTandemRepeatsandHomopolymers_slop5.bed.gz | 
bgzip -c > ARS-UCD1.2_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz

echo "sum of regions:"
zcat ARS-UCD1.2_AllTandemRepeatsandHomopolymers_slop5.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
mkdir -p Union
multiIntersectBed -i \
ARS-UCD1.2_AllTandemRepeats_lt51bp_slop5.bed.gz \
ARS-UCD1.2_AllTandemRepeats_51to200bp_slop5.bed.gz \
ARS-UCD1.2_AllTandemRepeats_201to10000bp_slop5.bed.gz \
ARS-UCD1.2_AllTandemRepeats_gt10000bp_slop5.bed.gz \
| sortBed -faidx ../references/ARS-UCD1.2_Btau5.0.1Y.fa.fai -i stdin \
| mergeBed -i stdin \
| bgzip -c > Union/GRCh38_allTandemRepeats.bed.gz

subtractBed \
-a ../references/ARS-UCD1.2_Btau5.0.1Y.genome.bed -b Union/GRCh38_allTandemRepeats.bed.gz \
| bgzip -c > Union/GRCh38_notinallTandemRepeats.bed.gz

echo "sum of regions:"
zcat Union/GRCh38_allTandemRepeats.bed.gz | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
