 bcftools norm -m -any /home/mshokrof/sv/cattle-sv/datasets/pangenome/NxB/multisample-vcfs/graph-filtered.vcf.gz |python ~/TheGreatGenotyper/scripts/assign-variant-ids_and_len.py  |python ~/TheGreatGenotyper/scripts/assign-variant-len.py |bcftools view -i "INFO/VARLEN >=40"  > variants/assemblies.vcf
 jasmine file_list=jasmine_list.txt out_file=merged.vcf threads=16

bcftools reheader -h header merged.vcf |grep -v "NKLS"|grep -vP "^MT" |grep -vP "^X"|grep -vP "^Y"|bcftools sort |bgzip > merged.sorted.vcf.gz

 bcftools norm -c x  -f ../datasets/references/ARS-UCD1.2_Btau5.0.1Y.fa merged.sorted.vcf.gz  > tmp.vcf
python ~/TheGreatGenotyper/ExtendedPangenome/1KG/addSamples.py tmp.vcf tmp.2.vcf

/home/mshokrof/TheGreatGenotyper/build6/pangenie/src/TheGreatGenotyper  -a  -g  -i parents.index.txt  -j 32 -t 32 -r ../datasets/references/ARS-UCD1.2_Btau5.0.1Y.fa  -y emissions -v tmp.vcf -o - 2> gg_log_2 | bgzip > gg_vcf_2.gz


bcftools view -s SAMEA9986201,SAMEA9986199  gg_vcf_2.gz |bcftools +fill-tags - -Ov  -- -t AC |bcftools view -i  'AC>1' |grep -vP "^#" |cut -f3 >> ids

bcftools view -i '(ID~"0_") || SUPP>1' 


bcftools view -i 'ID=@ids.uniq' merged.sorted.vcf.gz | bgzip > TruthSet.GxP.vcf.gz

