configfile: 'config.yaml'

#module load python/3.5.1 samtools freebayes vcftools bwa bedtools r/3.5.0 rstudio/1.1.453 sratoolkit bedops; PATH=$PATH:/nas/longleaf/home/csoeder/modules/vcflib/bin:/nas/longleaf/home/csoeder/modules/parallel/bin


sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]

sra_by_name = {}
sra_by_path = {}

for s in sample_by_name.keys():
	if "SRA" in sample_by_name[s].keys():
		sra_by_name[s] = sample_by_name[s]["SRA"]
		sra_by_path[sample_by_name[s]["path"]] = sample_by_name[s]["SRA"]

annotation_by_name = { a['name'] : a for a in config['annotations']}


rule annotation_reporter:
	input:
		annot = lambda wildcards: annotation_by_name[wildcards.annot_name]["bed_path"]
	output:
		report_out = "summaries/annotations/{annot_name}.stats"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		shell(""" mkdir -p summaries/annotations/ """)
		shell(""" rm -f {output.report_out} """)
		shell(""" cat {input.annot} | cut -f 1 | sort | uniq -c | tr -s " " | tr " " "\t" | awk '{{print"count\t"$2"\t"$1}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | wc -l | awk '{{print"count\ttotal\t"$0}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | awk '{{print$3-$2;}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "size\ttotal\t",sum; print "size\tavg\t",sum/NR; print "size\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}'  >> {output.report_out} """)


rule summon_annotation_summaries:
	input:
		refgen_reports = lambda wildcards: expand("summaries/annotations/{ref_ann}.stats", ref_ann= [a["name"] for a in config['annotations']  ] ) # annotation_by_name.keys())
	output:
		refann_summary = "summaries/reference_annotations.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		print([a["name"] for a in config['annotations'] ])
		shell(""" rm -f {output.refann_summary} """)
		for anne in [a["name"] for a in config['annotations']  ]:
			shell(""" cat summaries/annotations/{anne}.stats | awk '{{print"{anne}\t"$0}}' >> {output.refann_summary}""")





def return_filename_by_sampname(sampname):
	filenames = []
	if sample_by_name[sampname]['paired']:
		filenames.append(sample_by_name[sampname]['readsfile1'])
		filenames.append(sample_by_name[sampname]['readsfile2'])
	else:
		filenames.append(sample_by_name[sampname]['readsfile'])
	return filenames

def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = return_filename_by_sampname(sampname)
	pathsout = filesin
	return pathsout



rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "summaries/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"""
		mkdir -p summaries/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("summaries/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "summaries/reference_genomes.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"

rule summon_reads_SRA:
	output: 
		reads1='FASTQs/{path}/{prefix}_1.fastq',
		reads2='FASTQs/{path}/{prefix}_2.fastq',
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
	run:
		try:
#			sra = sra_by_path["%s/" % tuple([wildcards.path])]
			sra = sra_by_path[wildcards.path]
			shell(""" mkdir -p FASTQs/{wildcards.path}/ """)
			shell("""
				fasterq-dump  --split-3 --outdir FASTQs/{wildcards.path}/ {sra}
			""")
		except KeyError:
			raise KeyError("Sample is listed as empirical but no reads found and no SRA to download!" )




rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R0.fastq"],
		jason = "{pathprefix}/{samplename}.False.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.False.json",# --html meta/FASTP/{samplename}.html", 
		se_params = "",
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "{pathprefix}/{samplename}.True.json"
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.True.json",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"

#request FASTP reportBacks instead, process them into a unified SQL upload?
# rule clean_all_samps:
# 	input: 
# 		clean_se = [sample_by_name[nom]['path']+nom+".clean.R0.fastq" for nom in sample_by_name.keys() if not sample_by_name[nom]['paired']],
# 		clean_pe = [sample_by_name[nom]['path']+nom+".clean.R"+arr+".fastq" for nom in sample_by_name.keys() if sample_by_name[nom]['paired'] for arr in ["1","2"]],
# 	output:
# 		outflg = "allclean.flag"
# 	shell:
# 		"touch {output.outflg}"

#lambda wildcards: return_file_relpath_by_sampname(wildcards)

rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}/{samp}.{pairt}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, pairt = sample_by_name[wildcards.samplename]['paired'])
	output:
		jason_pruned = "summaries/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		cp {input.jason} summaries/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("summaries/FASTP/{samplename}.json.pruned", samplename = sampname_by_group["all"])
	output:
		summary = "summaries/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"


# #	request stats/idxstats/flagstats?  
rule bwa_mem:
	input:
		reads_in = lambda wildcards: expand("{path}/{sample}.clean.R{arr}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, arr=[ [1,2] if sample_by_name[wildcards.sample]['paired'] else [0] ][0]),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwaMEM.sort.bam",
	params:
		runmem_gb=64,
		runtime="64:00:00",
		cores=8,
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
	run:
		shell(""" mkdir -p mapped_reads""")
		print(input.reads_in[0])
		shell(""" bwa mem {input.ref_genome_file} {input.reads_in[0]} {input.reads_in[1]} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - """)
		shell("samtools index {output.bam_out}")


# #	request stats/idxstats/flagstats?  

# rule bwa_uniq:
# 	input:
# 		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam"
# 	output:
# 		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwaUniq.sort.bam"
# 	params:
# 		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 -F 0x04",
# 		uniqueness="XT:A:U.*X0:i:1.*X1:i:0",
# 		runmem_gb=16,
# 		runtime="6:00:00",
# 		cores=4,
# 	message:
# 		"filtering alignment of {wildcards.sample} to {wildcards.ref_genome} for quality and mapping uniqueness.... "	
# 	run:
# 		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
# 		#original; no dedupe
# 		#"samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome} - | samtools sort -o {output.bam_out} - "
# 		shell("samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome_file} - | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -n - | samtools fixmate -m - - | samtools sort - | samtools markdup -r - {output.bam_out}")
# 		shell("samtools index {output.bam_out}")



rule bam_reporter:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam",
		amps_in = "utils/squashed_amplicons.bed",
	output:
		report_out = "summaries/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary",
		amplicov = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam.ampliconcov",
	params:
		runmem_gb=16,
		runtime="4:00:00",
		cores=1,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -ibam {input.bam_in}  > {input.bam_in}.genomcov")
		if sample_by_name[wildcards.sample]["library"] == "amplicon":
			shell("bedtools coverage -a  {input.amps_in} -b {input.bam_in} > {input.bam_in}.ampliconcov")
		else:
			shell(" touch {input.bam_in}.ampliconcov")
#		shell("bedtools coverage -a  {input.amps_in} -b {input.bam_in} > {input.bam_in}.ampliconcov")

		#change the -max flag as needed to set 
		shell("""samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats""")
		#https://www.biostars.org/p/5165/
		#save the depth file and offload the statistics to the bam_summarizer script?? 
		shell("python3 scripts/bam_summarizer.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.genomcov -d {input.bam_in}.dpthStats -o {output.report_out} -t {wildcards.sample}")


rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("summaries/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sampname_by_group['all'], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
#		amplicovs = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam.ampliconcov", sample=list(set(sampname_by_group['backcross']).intersection(sampname_by_group['all'])), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
	output:
		full_report = "summaries/alignments.vs_{ref_genome}.{aligner}.summary",
		amplicov_report = "summaries/amplicon_regions.vs_{ref_genome}.{aligner}.coverage",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
		shell("""cat {input.bam_reports} > {output.full_report}""")
		shell(""" rm -rf {output.amplicov_report} """)
		for sample in sampname_by_group['all']:
			if sample_by_name[sample]["library"] == "amplicon":
				shell(""" cat mapped_reads/{sample}.vs_{wildcards.ref_genome}.{wildcards.aligner}.sort.bam.ampliconcov | awk '{{bases_covered+=$5; bases_total+=$6}} END {{ print "{sample}\tampliconic_coverage_breadth\t"bases_covered/bases_total}}' >> {output.amplicov_report}""")
				shell(""" samtools depth -a -b utils/squashed_amplicons.bed mapped_reads/{sample}.vs_{wildcards.ref_genome}.{wildcards.aligner}.sort.bam  | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "ampliconic_coverage_depth_avg\t",sum/NR; print "ampliconic_coverage_depth_std\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{sample}\t"$0}}' >> {output.amplicov_report}""")


 #	cat test.samtools.idxstats | sed \$d | awk '{print $1, $3/$2}' > per_contig_coverage_depth
# #	echo  $( cat test.samtools.idxstats |  sed \$d | cut -f 2 | paste -sd+ | bc) $(cat test.samtools.idxstats |  sed \$d | cut -f 3 | paste -sd+ | bc) | tr  " " "\t" | awk '{print "total\t"$1/$2}'


rule bam2Bed:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam"
	output:
		bed_out = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.bed"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
		shell(""" samtools view -hb -F 4 {input.bam_in} | bedtools bamtobed -bed12 -i  - | bedtools merge -i - > {output.bed_out} """)

rule squashed_amplicons:
	input:
		beds_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.bed", sample=set(sampname_by_group['backcross']).intersection(set(sampname_by_group['all'])), ref_genome="dps3", aligner="bwaMEM")
	output:
		bed_out = "utils/squashed_amplicons.bed",
		fai_out = "utils/squashed_amplicons.fa.fai"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
		shell(""" cat {input.beds_in} | bedtools sort -i - | bedtools merge -i - > {output.bed_out} """)
		ref_gen_path = ref_genome_by_name["dps3"]['path']
		shell("""bedtools getfasta -fi {ref_gen_path} -bed {output.bed_out} -fo  utils/squashed_amplicons.fa""")
		shell("""samtools faidx utils/squashed_amplicons.fa """)


# # rule joint_vcf_caller:
# # 	input:
# # 		bams_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sample_by_name.keys(), ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
# # 	output:
# # 		vcf_out = "variants/all_samples.vs_{ref_genome}.{aligner}.vcf"
# # 	params:
# # 		freebayes="--standard-filters",
# # 		runmem_gb=32,
# # 		runtime="96:00:00",
# # 	message:
# # 		"Jointly calling variants from all samples mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
# # 	run:
# # 		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
# # 		shell("freebayes {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}")


# ## VCF summary report? eg, shared fraction of genome covered by x% of samples at mindepth?


rule window_maker:
	output:
		windowed='utils/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	params:
		runmem_gb=8,
		runtime="5:00",
		cores=1,
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell("mkdir -p utils")
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} -i winnum | bedtools sort -i - > {output.windowed}'
		)


rule joint_vcf_caller_parallel:
	input:
		bams_in = lambda wildcards: expand("mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam", sample=sampname_by_group["all"], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner),
		windows_in = "utils/{ref_genome}_w100000_s100000.windows.bed",
		amplicon_regions = "utils/squashed_amplicons.bed",

	output:
		vcf_out = "variants/all_samples.vs_{ref_genome}.{aligner}.vcf"
	params:
		freebayes="--standard-filters",
		runmem_gb=32,
		runtime="96:00:00",
		cores=24,
	message:
		"Jointly calling variants from all samples mapped to \ {wildcards.ref_genome} \ with \ {wildcards.aligner} \ "
	run:
		shell(""" mkdir -p variants """)
		ref_genome_file=ref_genome_by_name[wildcards.ref_genome]['path']
		shell("""cat {input.windows_in}| awk '{{print$1":"$2"-"$3}}' > {input.windows_in}.rfmt""")
		shell("scripts/freebayes-parallel {input.windows_in}.rfmt {params.cores} {params.freebayes} --targets {input.amplicon_regions} -f {ref_genome_file} {input.bams_in} | vcftools --bed {input.amplicon_regions} --remove-indels --min-alleles 2 --max-alleles 2 --vcf - --recode --recode-INFO-all --stdout | uniq > {output.vcf_out} ")
		#shell("~/modules/freebayes/scripts/freebayes-parallel {input.windows_in}.rfmt {params.cores}  --standard-filters -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out} ")
		#shell("freebayes {params.freebayes} -f {ref_genome_file} {input.bams_in} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}")

rule vcf_reporter:
	input:
		vcf_in = "variants/{prefix}.vs_{ref_genome}.{aligner}.vcf"
	output:
		report_out = "summaries/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary",
		frq_out = "summaries/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary.frq"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=4,
#	message:
#		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	shell:
		"""
		cat {input.vcf_in}  | grep -v "#" | cut -f 1 | sort | uniq -c > {output.report_out}.snpsPerContig.tmp
		cat {output.report_out}.snpsPerContig.tmp | awk '{{sum+=$1}} END {{ print sum,"\ttotal"}}' | cat {output.report_out}.snpsPerContig.tmp - > {output.report_out}.snpsPerContig
		rm {output.report_out}.snpsPerContig.tmp
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --freq 
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --counts
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-indv
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --missing-site
		vcftools  --vcf {input.vcf_in} --out {output.report_out} --singletons

		ref_genome={wildcards.ref_genome}

		tail -n 1 {output.report_out}.snpsPerContig | awk '{{print "total_snp_count\t"$1}}' | sed -e 's/^/'$ref_genome'\t/g' > {output.report_out}
		"""
	#cat  all_samples.vs_droSec1.bwaUniq.summary.frq.count| cut -f 3 | tail -n +2 | sort | uniq -c
	#####	bi, tri, and quadralelic counts ^^ 
	#replace some of this with vcftools::vcf-stats ?

rule summon_VCF_analytics_base:
	input:
		bam_reports = lambda wildcards: expand("summaries/VCFs/{prefix}.vs_{ref_genome}.{aligner}.summary", prefix=wildcards.prefix, ref_genome=["dps3"], aligner="bwaMEM")
	output:
		full_report = "summaries/{prefix}.calledVariants.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"""
		prefix={wildcards.prefix}
		cat {input.bam_reports} | sed -e 's/^/'$prefix'\t/g'> {output.full_report}
		"""






# rule parse_indv:
# 	input:
# #		vcf_in = "variants/all_samples.vs_dps3.bwaMEM.vcf",
# 		vcf_in = "variants/informative_sites.vcf",
# 		bed_in = "variants/informative_sites.bed",
# 		windoze= "utils/dps3_w1000000_s50000.windows.bed",
# 	output:
# 		bed_out = "variants/{sample}.informative.bed",
# 		bins_out = "variants/{sample}.informative.binned.bed",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=1,
# 	message:
# 		"collecting all alignment metadata.... "
# 	run:
# 		shell(""" cat {input.vcf_in} | vcftools --vcf - --indv {wildcards.sample} --max-missing-count 0  --recode --recode-INFO-all  --stdout | bedtools intersect -wa -wb -b - -a {input.bed_in} | cut -f 1-8,18 | sed -e 's/0\/0:[0-9:,.-]*/-1/g' | sed -e 's/1\/1:[0-9:,.-]*/1/g'  | awk '{{print$1,$2,$3,$4,$5,$6,$7*$9,$8*$9}}' | tr " " "\t" > {output.bed_out} """)
# 		shell(""" cat {output.bed_out} | awk '{{print$1,$2,$3,$4,$5,$6,$7*$9,$8*$9}}' | tr " " "\t" | bedtools map -b - -a {input.windoze} -o mean,mean -c 7,8 > {output.bins_out} """)



# rule parse_everyone:
# 	input:
# 		bins_in = expand("variants/{sample}.informative.binned.bed", sample = set(sampname_by_group['backcross']).intersection(set(sampname_by_group['all'])) ),
# 	output:
# 		parse_flag = "utils/parse.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	message:
# 		"collecting all alignment metadata.... "
# 	run:
# 		shell(""" touch {output.parse_flag} """)


rule naiive_heterozyg:#simply take an individual, and calculate heterozygosity. 
	input:
		vcf_in = "variants/all_samples.vs_dps3.bwaMEM.vcf",
		windoze= "utils/dps3_w1000000_s50000.windows.bed",
	output:
		vcf_out = "variants/{sample}.naiive_heterozyg.vcf",
		bg_het_out = "variants/{sample}.naiive_heterozyg.bedgraph",
		bg_hom_out = "variants/{sample}.naiive_homozyg.bedgraph",
#		vcf_out = "",
#		bed_out = "variants/{sample}.hetZyggy.bed",
#		wig_out = "variants/{sample}.naiive_heterozyg.windowed.wig",
#		count_wig = "variants/{sample}.counts.windowed.wig",
#		bed_out = "variants/informative_sites.bed"
#		vcf_out = "variants/informative_sites.vcf"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
#		shell(""" cat {input.vcf_in} | vcftools --vcf - --indv {wildcards.sample} --max-missing-count 0  --recode --recode-INFO-all  --stdout | vcf2bed | awk '{{print$1,$2,$3,$6">"$7,0,"*",$11}}' | tr " " "\t" | sed -e 's/0\/0:[0-9:,.-]*/0/g' | sed -e 's/1\/1:[0-9:,.-]*/0/g' | sed -e 's/0\/1:[0-9:,.-]*/1/g' | sed -e 's/1\/0:[0-9:,.-]*/1/g' > {output.bed_out}""")
		shell(""" cat {input.vcf_in} | vcftools --vcf - --indv {wildcards.sample} --max-missing-count 0  --recode --recode-INFO-all  --stdout  > {output.vcf_out}""")
#		shell(""" cat {output.vcf_out} | vcf2bed | awk '{{print$1,$2,$3,$6">"$7,0,"*",$11}}' | tr " " "\t" | sed -e 's/0\/0:[0-9:,.-]*/0/g' | sed -e 's/1\/1:[0-9:,.-]*/0/g' | sed -e 's/0\/1:[0-9:,.-]*/1/g' | sed -e 's/1\/0:[0-9:,.-]*/1/g' """)
		shell(""" cat {output.vcf_out} |  grep -v "#" | awk '{{print$1,$2,$2+1,$4">"$5,0,"*",$10}}' | tr " " "\t"  | sed -e 's/0\/0:[0-9:,.-]*/0/g' | sed -e 's/1\/1:[0-9:,.-]*/0/g' | sed -e 's/0\/1:[0-9:,.-]*/1/g' | sed -e 's/1\/0:[0-9:,.-]*/1/g' | bedtools sort | uniq > variants/{wildcards.sample}.naiive_heterozyg.bed """)
		shell(""" echo "track type=bedGraph name='{wildcards.sample} naiive heterozygosity' description='{wildcards.sample} naiive heterozygosity' graphType=bar color=255,0,0" > {output.bg_het_out} """)
		shell(""" cat variants/{wildcards.sample}.naiive_heterozyg.bed | cut -f 1,2,3,7 >> {output.bg_het_out} """)
		shell(""" echo "track type=bedGraph name='{wildcards.sample} naiive homozygosity' description='{wildcards.sample} naiive heterozygosity' graphType=bar color=0,0,255" > {output.bg_hom_out} """)
		shell(""" cat variants/{wildcards.sample}.naiive_heterozyg.bed | cut -f 1,2,3,7 | awk '{{print$1"\t"$2"\t"$3"\t"1-$4}}' >> {output.bg_hom_out} """)

#		shell(""" bedtools map -o mean -c 7 -a {input.windoze} -b {output.bed_out} | awk '{{if($5!=".")print;}}' | awk '{{if($1==2)print;}}' | cut -f 2,5 >> {output.wig_out} """)
#		shell(""" cp utils/wig.header {output.count_wig} """)
#		shell(""" bedtools map -o count -c 7 -a {input.windoze} -b {output.bed_out} | awk '{{if($5!=".")print;}}' | awk '{{if($1==2)print;}}' | cut -f 2,5 >> {output.count_wig} """)

rule naiive_melody:
	input:
		bins_in = expand("variants/{sample}.naiive_heterozyg.bedgraph", sample = set(sampname_by_group['backcross']).intersection(set(sampname_by_group['all'])) ),
	output:
		parse_flag = "utils/naiive.flg"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
		shell(""" touch {output.parse_flag} """)





rule informative_sites:
	input:
		vcf_in = "variants/all_samples.vs_dps3.bwaMEM.vcf"
	output:
		bed_out = "variants/informative_sites.bed",
		vcf_out = "variants/informative_sites.vcf"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	run:
 		shell(""" cat {input.vcf_in} | vcftools --vcf - --indv MV2_25 --indv Flag14  --recode --recode-INFO-all  --stdout | grep  -v "#" | grep "0/0" | grep "1/1" | vcf2bed | sed -e 's/0\/0:[0-9:,.-]*/-1/g' |sed -e 's/1\/1:[0-9:,.-]*/1/g' | awk '{{if($11!=$12)print$1,$2,$2+1,"potato",0,"*",$11,$12}}' | tr " " "\t" > {output.bed_out} """)
 		shell(""" grep "#" {input.vcf_in} > {output.vcf_out} """)
 		shell(""" bedtools intersect -wa -a {input.vcf_in} -b {output.bed_out} >> {output.vcf_out} """)











