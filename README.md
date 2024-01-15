![image](https://github.com/dib-lab/Cattle_Structural_Variant_Benchmark/assets/5207616/0c03f970-e531-4aaf-987f-0b15ccb7ca06)# Cattle_Structural_Variant_Benchmark
Detecting structural variants (SVs) remains a significant challenge, and benchmarking SV detection methods is further complicated by the limited availability of truth sets. Numerous efforts have been made to evaluate SV callers in humans, primarily relying on the Genome in a Bottle (GIAB) and Human Genome Structural Variation (HGSV) datasets. These evaluations have demonstrated that the performance of SV callers varies considerably depending on factors such as read coverage, variant type and size, association with repetitive elements, and overall genome quality.
Although these findings provide valuable insights into the performance of SV callers, they cannot be directly extrapolated to animal genomics without further investigation. This is due to differences in genome quality and content, as well as the fact that reference genomes are not available for all breeds. Consequently, researchers often resort to utilizing reference genomes from closely related breeds. The performance of SV callers may differ when the sample and reference genomes diverge significantly from one another.
In our study, we got inspired from GIAB workflows to establish benchmark datasets and curated a Structural Variant Benchmark suite for cattle genomes. We used publicly available data from four trio assembly projects.  

## Benchmark SV:
Please download the trust sets for structural variant (SV) calling and bed files from datasets where we have confidence in the consistency of results produced by mappers and callers. We have generated four trust sets based on four cattle trios: Gaur and Piedmontese (GxP), Nellore and Brown Swiss (NxB), Brahman and Angus (BxA), and Original Braunvieh (OxO). Additionally, ensure to download the high-confidence bed file. [ARS-UCD1.2_notinalldifficultregions.bed.gz](https://github.com/dib-lab/Cattle_Structural_Variant_Benchmark/blob/main/ARS-UCD1.2_notinalldifficultregions.bed.gz)

To run the benchmark you need truvari to be installed on the system

```
 truvari bench   -b  TruthSet.OxO.vcf.gz  -c <test vcf>  --includebed ARS-UCD1.2_notinalldifficultregions.bed.gz -o truvari_final_outputs/OxO_{1} --passonly -r 2000 -C 3000  --reference ARS-UCD1.2_Btau5.0.1Y.fa
```

## Scripts and Intermediate files
Scripts used to generate the benchmark and the high confidence regions are available at scripts/. Also, Intermediate files are uploaded to the repo under Intermediate_Files/





 
