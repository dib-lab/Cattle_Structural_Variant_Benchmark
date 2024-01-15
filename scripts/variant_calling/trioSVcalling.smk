pepfile: "project_config.yaml"
configfile: "config.yaml" 

tempFolder=config["tempFolder"]
outputFolder=config["outputFolder"]
pbsvMappingTool=config["pbsvMappingTool"]

reference=list(filter(lambda x:x.library=="ref",pep.samples))[0].file[0]
reference_repeats=list(filter(lambda x:x.library=="ref_repeats",pep.samples))[0].file[0]


bioSamples=dict([(s.sample_name,s.BioSample) for s in pep.samples])


tmp=list(filter(lambda x:x.library=="assembly",pep.samples))
assemblies=[s.sample_name for s in tmp]
assemblyVCF= outputFolder+"svim-asm/diploid_"+assemblies[0]+"-"+assemblies[1]+".vcf"

mergeInput={}

tmp=list(filter(lambda x:x.library=="hifi",pep.samples))
hifiSamples=[s.sample_name for s in tmp]
mergeInput["hifi"]= expand("{out}mapping/{sample}.{tool}.bam",out=outputFolder,sample=hifiSamples,tool=pbsvMappingTool)

tmp=list(filter(lambda x:x.library=="ont",pep.samples))
ontSamples=[s.sample_name for s in tmp]
mergeInput["ont"]= expand("{out}mapping/{sample}.minimap2.bam",out=outputFolder,sample=ontSamples)
 



tmp=list(filter(lambda x:x.library=="illumina",pep.samples))
illuminaSamples=[s.sample_name for s in tmp]
mergeInput["illumina"]= expand("{out}mapping/{sample}.bwa.bam",out=outputFolder,sample=illuminaSamples)
 

tmp=list(filter(lambda x:x.library=="geno",pep.samples))
genoSamples=[s.sample_name for s in tmp]
mergeInput["geno"]= expand("{out}mapping/{sample}.bwa.bam",out=outputFolder,sample=genoSamples)
 

def getFile(name,i):
    return list(filter(lambda x:x.sample_name==name,pep.samples))[0].file[i]


vcfs=[
        assemblyVCF,
	outputFolder+"pbsv/varints.PERCISE.vcf",
        outputFolder+"cuteSV/variants.hifi.PERCISE.vcf",
	outputFolder+"sniffles/variants.hifi.PERCISE.vcf",
#        outputFolder+"cuteSV/variants.ont.PERCISE.vcf",
#	outputFolder+"sniffles/variants.ont.PERCISE.vcf",
	outputFolder+"manta/merged.illumina.candidateSV.PERCISE.vcf",
	outputFolder+"delly/merged.illumina.PERCISE.vcf",
        outputFolder+"gridss/merged.illumina.filtered.vcf"
]
vcfsNames=[
       'svim-asm',
	"pbsv",
        "cuteSV-hifi",
	"sniffles-hifi",
#        "cuteSV-ont",
#	"sniffles-ont",
	"manta",
	"delly",
        "gridss"
]


rule all:
    input:
        outputFolder+"mapping/merged.hifi.flagstat",
	outputFolder+"mapping/merged.hifi.alfred.txt",
        outputFolder+"mapping/merged.ont.flagstat",
	outputFolder+"mapping/merged.ont.alfred.txt",
        outputFolder+"mapping/merged.illumina.flagstat",
	outputFolder+"mapping/merged.illumina.alfred.txt",
        outputFolder+"SURVIVOR/allToolsUnion.all.tsv.gz",
	vcfs,
	mergeInput["geno"],
	outputFolder+"svanalyzer/allToolsUnion.all.vcf"
	



### Assembly SV Calling
rule minimap2asm:
    input:
        ref=reference,
	assembly=lambda wildcards: getFile(f"{wildcards.sample}",0)
    output:
        bam="{outputFolder}svim-asm/{sample}.bam",
        bai="{outputFolder}svim-asm/{sample}.bam.bai"	
    conda:
        "../envs/minimap2.yaml"
    log:
        "{outputFolder}svim-asm/{sample}.mapping.log"
    threads: 32
    resources:
        mem_mb=50000,cores=32
    shell:
        """
                mkdir -p {tempFolder}$$/
		minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t {threads} {input.ref} {input.assembly} 2>>{log} | samtools sort --reference {input.ref} -m4G -O BAM -T {tempFolder}$$/ -o {output.bam} - 2>> {log}
		samtools index {output.bam}
		rm {tempFolder}$$/ -rf
	"""

rule svimasDiploid:
    input:
        bam1="{out}svim-asm/{sample}.bam", 
        bai1="{out}svim-asm/{sample}.bam.bai",
	bam2="{out}svim-asm/{sample2}.bam", 
        bai2="{out}svim-asm/{sample2}.bam.bai",
	ref=reference
    output:
        sv="{out}svim-asm/diploid_{sample}-{sample2}.vcf",
        hist="{out}svim-asm/diploid_{sample}-{sample2}.png"
    conda:
        "../envs/svim-asm.yaml"
    log:
        "{out}svim-asm/diploid_{sample}-{sample2}.log"
    threads: 1
    resources:
        mem_mb=20000,cores=1
    shell:
        """
        mkdir -p {tempFolder}$$/
	svim-asm diploid {tempFolder}$$/ {input.bam1} {input.bam2} {input.ref} &> {log}
	cp {tempFolder}$$/variants.vcf {output.sv}
	cp {tempFolder}$$/sv-lengths.png {output.hist}
	cat {tempFolder}$$/SVIM_*.log >> {log}
	rm {tempFolder}$$/ -rf
	"""



### Long Read SV Calling

rule minimap2:
    input:
        reads=lambda wildcards: getFile(f"{wildcards.sample}",0),
	ref=reference+".mmi"	
    output:
        bam=outputFolder+"mapping/{sample}.minimap2.bam",
        bai=outputFolder+"mapping/{sample}.minimap2.bam.bai"       
    conda:
        "../envs/minimap2.yaml"
    log:
        outputFolder+"mapping/{sample}.minimap2.log"
    params:
        biosample= lambda wildcards: bioSamples[f"{wildcards.sample}"]
    resources:
         mem_mb= 50000, cores =64
    threads: 64
    shell:
        """
                mkdir -p {tempFolder}$$/
		mkdir -p {outputFolder}
               # cp {input.reads} {input.ref} {tempFolder}$$/

               # inputReads=$(basename {input.reads})
              #  reference=$(basename {input.ref})

		minimap2 -ax map-ont --MD  -t {threads} {input.ref} {input.reads}   2> {log}  |samtools sort -T {tempFolder}$$/tmpBam  -O BAM  -   > {output.bam}
		samtools index {output.bam}
 	#	cp {tempFolder}$$/out.bam {output.bam} 2>> {log}
	#	cp {tempFolder}$$/out.bam.bai {output.bai} 2>> {log}
		rm {tempFolder}$$/ -r
	"""


rule pbmm2_index:
    input:
        ref=reference	
    output:
        index=reference+".mmi"	
    conda:
        "../envs/pbsv.yaml"
    log:
        reference+".pbmm2Index.log"
    threads: 1
    resources:
         mem_mb= 50000, cores =1
    shell:
        """
		pbmm2 index {input.ref}  {output.index} &> {log}

	"""

rule pbmm2:
    input:
        reads=lambda wildcards: getFile(f"{wildcards.sample}",0),
	ref=reference+".mmi"	
    output:
        bam=outputFolder+"mapping/{sample}.pbmm2.bam",
        bai=outputFolder+"mapping/{sample}.pbmm2.bam.bai"       
    conda:
        "../envs/pbsv.yaml"
    log:
        outputFolder+"mapping/{sample}.pbmm2.log"
    params:
        biosample= lambda wildcards: bioSamples[f"{wildcards.sample}"]
    resources:
         mem_mb= 50000, cores =32
    threads: 32
    shell:
        """
                mkdir -p {tempFolder}$$/
		mkdir -p {outputFolder}
                cp {input.reads} {input.ref} {tempFolder}$$/

                inputReads=$(basename {input.reads})
                reference=$(basename {input.ref})

		pbmm2 align {tempFolder}$$/$reference {tempFolder}$$/$inputReads --preset HIFI --sample {params.biosample} --rg "@RG\tID:{params.biosample}" 2> {log}  |samtools sort -T {tempFolder}$$/tmpBam  -O BAM  -   > {output.bam}
		samtools index {output.bam}
 	#	cp {tempFolder}$$/out.bam {output.bam} 2>> {log}
	#	cp {tempFolder}$$/out.bam.bai {output.bai} 2>> {log}
		rm {tempFolder}$$/ -r
	"""


## pbsv


rule discoverSVSignature:
    input:
         bam=outputFolder+"mapping/{file}"+pbsvMappingTool+".bam"
    output:
         outputFolder+"pbsv/{file}"+pbsvMappingTool+".svsig.gz"
    conda:
         "../envs/pbsv.yaml"
    log:
         outputFolder+"pbsv/{file}"+pbsvMappingTool+".log"
    threads: 1
    resources:
        mem_mb=30000, cores=1
    shell:
       	        """
		mkdir -p {outputFolder}pbsv/ 
		pbsv discover {input.bam} {output} &> {log} 
                """

rule svCall:
    input:
         samples=[expand("{prefix}pbsv/{file}.{map}.svsig.gz",prefix=outputFolder,file=hifiSamples,map=pbsvMappingTool)],
	 ref=reference
    output:
         outputFolder+"pbsv/varints.vcf"
    conda:
         "../envs/pbsv.yaml"
    log:
         outputFolder+"pbsv/log"
    threads: 32
    resources:
        mem_mb=100000, cores=32
    shell:
       	        """
		mkdir -p {outputFolder}/pbsv/
		pbsv call -j {threads} {input.ref}  {input.samples} {output} &> {log}
                """

rule cuteSVONT:
    input:
        bam = outputFolder+"mapping/merged.ont.bam",
	bai = outputFolder+"mapping/merged.ont.bam.bai",
	ref = reference
    output:
        outputFolder+"cuteSV/variants.ont.vcf"
    conda:
        "../envs/cuteSV.yaml"
    log:
        outputFolder+"cuteSV/ont.log"
    threads: 32
    resources:
        mem_mb=100000, cores=32
    shell:
        """
		mkdir -p {tempFolder}$$/

		bamName=$(basename {input.bam})
		cp {input.bam}* {tempFolder}$$/
		cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL	100 --diff_ratio_merging_DEL 0.3 -t {threads} {tempFolder}$$/$bamName {input.ref} {output} {tempFolder}$$/ &>{log}
		rm {tempFolder}$$/ -r
	"""

rule cuteSVHIFI:
    input:
        bam = outputFolder+"mapping/merged.hifi.bam",
	bai = outputFolder+"mapping/merged.hifi.bam.bai",
	ref = reference
    output:
        outputFolder+"cuteSV/variants.hifi.vcf"
    conda:
        "../envs/cuteSV.yaml"
    log:
        outputFolder+"cuteSV/hifi.log"
    threads: 32
    resources:
        mem_mb=100000, cores=32
    shell:
        """
		mkdir -p {tempFolder}$$/

		bamName=$(basename {input.bam})
		cp {input.bam}* {tempFolder}$$/
		cuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t {threads} {tempFolder}$$/$bamName {input.ref} {output} {tempFolder}$$/ &>{log}
		rm {tempFolder}$$/ -r
	"""

rule sniffles:
    input:
        bam = outputFolder+"mapping/merged.{tech}.bam",
	bai = outputFolder+"mapping/merged.{tech}.bam.bai",
	ref = reference
    output:
        outputFolder+"sniffles/variants.{tech}.vcf"
    conda:
        "../envs/sniffles.yaml"
    log:
        outputFolder+"sniffles/{tech}.log"
    threads: 16
    resources:
        mem_mb=30000, cores=16
    shell:
        """
		mkdir -p {outputFolder}/sniffles/
		mkdir -p {tempFolder}$$/

		cp {input.ref} {tempFolder}$$/
#		cp {input.bam} {tempFolder}$$/
		bamName=$(basename {input.bam})
		refName=$(basename {input.ref})

		samtools calmd -@ {threads} -b {input.bam} {tempFolder}$$/$refName  > {tempFolder}$$/$bamName.md.bam 2> {log}
		samtools index  {tempFolder}$$/$bamName.md.bam

		sniffles -t {threads} --input {tempFolder}$$/$bamName.md.bam  --vcf  {output}  &> {log}
		rm {tempFolder}$$/ -r
	"""




### Short read calling
ruleorder: merge > bwa

rule bwaIndex:
    input:
         reference
    output:
         reference+".bwt"
    log: config["outputFolder"] + "bwaIndex.log"
    threads: 1
    resources:
        mem_mb=10000, cores=1
    conda:
         "../envs/bwa.yaml"
    shell:
       	        """		
		mkdir -p {tempFolder}$$/		
		cp {input} {tempFolder}$$/
		inputName=$(basename {input})
		bwa index {tempFolder}$$/$inputName &> {log}
		outDir=$(dirname {output})
		cp {tempFolder}$$/$inputName.* $outDir
		rm {tempFolder}$$/ -r		
        """

rule bwa:
    input:
        readsForward=lambda wildcards: getFile(f"{wildcards.sample}",0),
	readsReverse=lambda wildcards: getFile(f"{wildcards.sample}",1),
	ref=reference,
	index=reference +".bwt"
    output:
        bam=outputFolder+"mapping/{sample}.bwa.bam",
        bai=outputFolder+"mapping/{sample}.bwa.bam.bai"	
    conda:
        "../envs/bwa.yaml"
    log:
        outputFolder+"{sample}.bwa.log"
    threads: 32
    resources:
        mem_mb=30000, cores=32
    shell:
        """
		mkdir -p {outputFolder}
		mkdir -p {tempFolder}$$/
		cp {input.readsForward} {input.readsReverse} {tempFolder}$$/
		inputFQ1=$(basename {input.readsForward})
		inputFQ2=$(basename {input.readsReverse})
		bwa mem -t {threads} {input.ref} {tempFolder}$$/$inputFQ1 {tempFolder}$$/$inputFQ2   2> {log} |samtools sort -T {tempFolder}$$/tmpBam -o {output.bam} -O BAM 2>> {log}
		samtools index {output.bam} 2>> {log}
		rm {tempFolder}$$/ -r
	"""

rule manta:  
    input:
        ref = reference,
        bam = outputFolder+"mapping/{file}.bam",
	bai = outputFolder+"mapping/{file}.bam.bai"
    output:
        vcf = outputFolder+"manta/{file}.candidateSV.vcf.gz",
        vcfUncompressed = outputFolder+"manta/{file}.candidateSV.vcf",	
	vcfIndex = outputFolder+"manta/{file}.candidateSV.vcf.gz.tbi",
	diploidVCF = outputFolder+"manta/{file}.diploidSV.vcf.gz",
	diploidVCFIndex = outputFolder+"manta/{file}.diploidSV.vcf.gz.tbi",
	resultsTar= outputFolder+"manta/{file}.allresults.tar.gz"
    log: outputFolder+"manta/{file}.log"
    conda:
        "../envs/manta.yaml"
    threads:
        32
    resources:
        mem_mb = 20000,cores=32, mem_gb=20
    shell:
        """
	mkdir -p {outputFolder}manta/
	mkdir -p {tempFolder}$$/
	samtools index {input.bam}
	configManta.py --bam={input.bam} --referenceFasta={input.ref} --runDir={tempFolder}$$/ &> {log}
	{tempFolder}$$/runWorkflow.py -j {threads}   -g {resources.mem_gb} &> {log}
	tar -czvf {output.resultsTar} {tempFolder}$$/
	cp {tempFolder}$$/results/variants/candidateSV.vcf.gz  {output.vcf}
	cp {tempFolder}$$/results/variants/candidateSV.vcf.gz.tbi  {output.vcfIndex}
	cp {tempFolder}$$/results/variants/diploidSV.vcf.gz  {output.diploidVCF}
	cp {tempFolder}$$/results/variants/diploidSV.vcf.gz.tbi  {output.diploidVCFIndex}
	gzip -dc {output.vcf} > {output.vcfUncompressed}
	rm -rf {tempFolder}$$/
        """

ruleorder: filterPRECISE > delly
ruleorder: filterPRECISE > manta
ruleorder: filterPRECISE > cuteSVONT
ruleorder: filterPRECISE > cuteSVHIFI
ruleorder: filterPRECISE > sniffles


rule filterPRECISE:
    input: "{prefix}.vcf"
    output: "{prefix}.PERCISE.vcf"
    log: "{prefix}.PERCISE.vcf.log"
    resources:
        mem_mb = 1000, cores=1
    shell:
        """
		grep -av  "IMPRECISE" {input} >{output} 2>{log}
        """     

rule delly:  
    input:
        ref = reference,
        bam = outputFolder+"mapping/{file}.bam",
	bai = outputFolder+"mapping/{file}.bam.bai"
    output:
        vcf = outputFolder+"delly/{file}.vcf"
    log:
        outputFolder+"delly/{file}.log"
    conda:
        "../envs/delly.yaml"
    threads:
        32
    resources:
        mem_mb = 20000, cores=32
    shell:
        """
	mkdir -p {outputFolder}delly/
	mkdir -p {tempFolder}$$/
	cp {input.bam} {input.bai} {tempFolder}$$/
	bam=$(basename {input.bam})
	bam="{tempFolder}$$/$bam"
	delly call -o {tempFolder}$$/sv.bcf -g {input.ref} $bam &> {log}
	bcftools convert -o {output.vcf} -O v {tempFolder}$$/sv.bcf   
	rm -rf {tempFolder}$$/
        """

ruleorder: filterGridss > gridss 

rule filterGridss:
    input: outputFolder+"gridss/{file}.vcf.gz"
    output: outputFolder+"gridss/{file}.filtered.vcf"
    log: outputFolder+"gridss/{file}.filtered.vcf.log"
    resources:
        mem_mb = 3000, cores=1
    conda:
        "../envs/gridss.yaml"
    shell:
        """
		gzip -dc {input} | grep -v  "LOW_QUAL"  > tmp.$$.vcf 2>{log}
		Rscript --vanilla ../addSVTypeGridss.R tmp.$$.vcf {output} 2>> {log}
        """     



rule gridss:  
    input:
        ref = reference,
        bam = outputFolder+"mapping/{file}.bam",
	bai = outputFolder+"mapping/{file}.bam.bai"
    output:
        vcf = outputFolder+"gridss/{file}.vcf.gz",
	bam = outputFolder+"gridss/{file}.assembly.bam"
    log: outputFolder+"gridss/{file}.log"
    conda:
        "../envs/gridss.yaml"
    threads:
        8
    resources:
        mem_mb = 40000, cores=8
    shell:
        """ 
	mkdir -p {outputFolder}gridss/
	mkdir -p {tempFolder}$$/
	#touch {input.bam}.bai
	../gridss/bin/gridss --reference {input.ref} --output {tempFolder}$$/out.vcf --assembly {output.bam} --threads {threads} --workingdir {tempFolder}$$/  {input.bam} --jvmheap {resources.mem_mb}M  &> {log}
	gzip -c {tempFolder}$$/out.vcf > {output.vcf} 
	rm -rf {tempFolder}$$/
        """



### Merging Results

rule check:
    input:
         vcf="{prefix}.vcf",
	 ref=reference
    output:
         "{prefix}.checked.vcf"
    conda:
         "../envs/survivor.yaml"
    log:"{prefix}.check.vcf.log"
    resources:
        mem_mb = 3000, cores=1
    shell:
       	 """
	 /home/mshokrof/sv/vcflib/build/vcfcheck -x -f {input.ref} {input.vcf} > {output} 2> {log}
         """



ruleorder: check > SurvivorUnion 

rule SurvivorUnion:
    input:
         vcfs=vcfs
    output:
         outputFolder+"SURVIVOR/allToolsUnion.all.vcf"
    params: names=vcfsNames
    conda:
         "../envs/delly.yaml"
    log: outputFolder+"SURVIVOR/allToolsUnion.all.vcf.log"
    resources:
        mem_mb =4000, cores=1
    shell:
       	 """
		echo {input.vcfs} |tr -s ' ' $'\n' > input.$$.lst
	 	../SURVIVOR-master/Debug/SURVIVOR merge input.$$.lst 2000 1 1 1 1 50 tmp.$$.vcf &> {log}
		bcftools sort -O v -o tmp2.$$.vcf tmp.$$.vcf
		echo {params.names} |tr -s ' ' $'\n' > input.$$.names
		bcftools reheader -s input.$$.names tmp2.$$.vcf > {output}
		rm input.$$.lst	 tmp2.$$.vcf tmp.$$.vcf
         """
rule svAnalyzerUnioin:
    input:
         vcfs=vcfs,
	 ref=reference
    output:
         vcf=outputFolder+"svanalyzer/allToolsUnion.all.vcf",
	 distances=outputFolder+"svanalyzer/allToolsUnion.distances"
    params: names=vcfsNames
    conda:
         "../envs/svanalyzer.yaml"
    log: outputFolder+"svanalyzer/allToolsUnion.all.vcf.log"
    resources:
        mem_mb =10000, cores=1
    shell:
       	 """
		mkdir -p {tempFolder}$$/
		ls {input.vcfs} >  {tempFolder}$$/source.lst
		grep -oP "[^/]*/[^/]*$" {tempFolder}$$/source.lst > {tempFolder}$$/tmp1 2> {log}
		cut -f1 tmp1 -d "/" | sed -e "s#^#{tempFolder}$$/#" > 	{tempFolder}$$/folders.lst 2>> {log}
		parallel --gnu 'mkdir {{}}' :::: {tempFolder}$$/folders.lst 2>> {log}
		paste {tempFolder}$$/source.lst {tempFolder}$$/folders.lst |   xargs -n2 cp
		paste -d '/' {tempFolder}$$/folders.lst <(cut -f2 tmp1 -d "/") >  {tempFolder}$$/input.lst 2>> {log}
		cp {input.ref} {tempFolder}$$/
		refName=$(basename {input.ref})
		svanalyzer merge --ref {tempFolder}$$/$refName  --fof {tempFolder}$$/input.lst --prefix {tempFolder}$$/tmp &>> {log}
		cat {tempFolder}$$/tmp.log >> {log}
		mv {tempFolder}$$/tmp.distances {output.distances}
		mv {tempFolder}$$/tmp.clustered.vcf {output.vcf}
		rm -rf  {tempFolder}$$/
         """


rule convertToMLTSV:
    input:
         vcf=outputFolder+"SURVIVOR/allToolsUnion.all.vcf",
	 manta=outputFolder+"manta/merged.illumina.candidateSV.PERCISE.vcf",
	 delly=outputFolder+"delly/merged.illumina.PERCISE.vcf",
	 gridss=outputFolder+'gridss/merged.illumina.filtered.vcf',
	 pbsv=outputFolder+"pbsv/varints.PERCISE.vcf",
         cuteSV=outputFolder+"cuteSV/variants.hifi.PERCISE.vcf",
	 sniffles=outputFolder+"sniffles/variants.hifi.PERCISE.vcf",
	 repeats= reference_repeats
    output:
         training = outputFolder+"SURVIVOR/allToolsUnion.all.tsv.gz"
    conda:
         "../envs/survivor.yaml"
    resources:
        mem_mb =4000, cores=1
    log: outputFolder+"SURVIVOR/allToolsUnion.all.tsv.log"
    shell:
       	 """
	 python ../createBenchmarkTable.py {input.vcf} svim pbsv cuteSV manta delly gridss  > tmp.$$.tsv 2> {log}
	 python ../addPbsvSpecificData.py tmp.$$.tsv {input.pbsv} > tmp2.$$.tsv 2>> {log}
	 python ../addCallerSpecificData.py tmp2.$$.tsv {input.sniffles} sniffles ../snifflesInfoFields ../snifflesSampleFields > tmp3.$$.tsv 2>> {log}
	 python ../addCallerSpecificData.py tmp3.$$.tsv {input.cuteSV} cuteSV ../cuteVSInfo ../cuteVSExtra > tmp4.$$.tsv 2>> {log}
	 python ../addCallerSpecificData.py tmp4.$$.tsv {input.manta} manta ../mantaInfo  ../mantaSample > tmp5.$$.tsv 2>> {log}
	 python ../addCallerSpecificData.py tmp5.$$.tsv {input.delly} delly ../dellyInfo  ../dellySamples > tmp6.$$.tsv 2>> {log}
	 python ../addCallerSpecificData.py tmp6.$$.tsv {input.gridss} gridss ../gridssInfo  ../gridssSample > tmp7.$$.tsv 2>> {log}
	 python ../addRepeatAnnotation.py tmp7.$$.tsv {input.repeats} > tmp8.$$.tsv
	 python ../createTrainingSet.py tmp8.$$.tsv all training  2>> {log} |gzip > {output.training}
	 rm tmp*.$$.tsv
         """





### BAM utilities

rule merge:
    input:
         samples=lambda wildcards: mergeInput[f"{wildcards.type}"],
	 ref = reference
    output:
         bam=outputFolder+"mapping/merged.{type}.bam",
         bai=outputFolder+"mapping/merged.{type}.bam.bai"
    conda:
         "../envs/pbsv.yaml"
    log:
        outputFolder+"mapping/merged.{type}.log"
    resources:
         mem_mb= 200000, cores =8
    threads: 8	 
    shell:
       	        """
		mkdir -p {tempFolder}$$/		
		cp {input.samples} {tempFolder}$$/
		samtools merge   -O BAM --reference {input.ref} --threads {threads} {output.bam} {tempFolder}$$/*bam &> {log}
		samtools index {output.bam}
		rm {tempFolder}$$/ -r		
                """




rule alfredQC:
    input:
         bam=outputFolder+"{file}.bam",
         ref=reference
    output:
         tsv=outputFolder+"{file}.alfred.qc.tsv.gz",
	 txt=outputFolder+"{file}.alfred.txt"
    conda:
         "../envs/alfred.yaml"
    resources:
         mem_mb= 10000, cores =1
    log:
        outputFolder+"{file}.alfred.qc.log"
    threads: 1	 
    shell:
       	        """
                alfred qc -i  -r {input.ref} -o {output.tsv} {input.bam} &> {log}
		zgrep "^ME" {output.tsv} 2>>{log} | datamash transpose  > {output.txt} 2>>{log}
                """

rule samtoolsStat:
    input:
         outputFolder+"{file}.bam"
    output:
         outputFolder+"{file}.flagstat"
    conda:
         "../envs/bwa.yaml"
    resources:
         mem_mb= 5000, cores =1
    log:
        outputFolder+"{file}.pbmm2.flagstat.log"
    threads: 1	 
    shell:
       	        """
		samtools flagstat {input} > {output} 2>{log}		
                """


