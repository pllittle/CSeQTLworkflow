# CSeQTLworkflow

## To Dos

* eQTL mapping code
* Torus enrichment
* GWAS enrichment

## Introduction and Outline

This repository contains template codes for 
runnings workflow steps including

* [Reference data](https://github.com/pllittle/CSeQTLworkflow#reference-data)
* [CommonMind Consortium data preparation](https://github.com/pllittle/CSeQTLworkflow#commonmind-consortium)
* [Blueprint data preparation](https://github.com/pllittle/CSeQTLworkflow#blueprint)
* [BAM workflow to TReC/ASReC](https://github.com/pllittle/CSeQTLworkflow#bam-workflow)
* [Deconvolution](https://github.com/pllittle/CSeQTLworkflow#deconvolution) with
	CIBERSORT and ICeDT
* eQTL mapping
* Torus and GWAS enrichment

for our [manuscript](https://www.biorxiv.org/content/10.1101/2022.03.31.486605v1) 
**Cell Type-Specific Expression Quantitative Trait Loci**.

## Reference Data

* [GTEx-pipeline for RNA-seq](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md)
* [Reference fasta](https://personal.broadinstitute.org/francois/topmed/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.tar.gz)
* [GTF](https://personal.broadinstitute.org/francois/topmed/gencode.v26.GRCh38.ERCC.genes.gtf.gz)

## GTEx Brain and Blood

Since GTEx samples are available on the NHGRI AnVIL cloud, they
were processed using a combination of Docker and WDL with details
provided [here](https://github.com/Sun-lab/gtex_AnVIL).

## CommonMind Consortium

<details>
<summary>Click to expand!</summary>

Code to obtain CMC inputs

```R
# Install package
install.packages("synapser",
	repos = c("https://sage-bionetworks.github.io/ran","http://cran.fhcrc.org"))

library(synapser)

# Login to Synapse
username = "." # your username
password = "." # your password
synLogin(username,password)

down_dir = "." # user-specified directory to download files

# Function used to download files
entity = "syn4600989" # an example SYN ID unique to a file
synGet(entity = entity,downloadLocation = down_dir)

# List of SYN_IDs, use synGet() to download
"syn4600989" 	# CMC Human sampleIDkey metadata
"syn4600985" 	# QC'd genotype bed file data
"syn4600987" 	# QC'd bim file
"syn4600989" 	# QC fam file
"syn4935690" 	# README file
"syn5600756" 	# list of outlier samples
"syn18080588"	# SampleID key
"syn3354385"	# clinical data
"syn18358379"	# RNAseq meta/QC data

# List of BAM SYN_IDs
tableId = "syn11638850" 
results = synTableQuery(smart_sprintf("select * from %s 
	where species='Human' and dataType='geneExpression'", tableId))
# results$filepath
aa = as.data.frame(results)
aa = aa[grep("bam",aa$name),]
aa = aa[grep("accepted",aa$name),] # excluding unmapped bam files
dim(aa)
saveRDS(aa,"aligned_bam_files.rds")

# Download BAM files one by one
BAM_dir = "./BAM"; dir.create(BAM_dir)
for(samp in sort(unique(aa$specimenID))){
	samp_dir = file.path(BAM_dir,samp)
	dir.create(samp_dir)
	samp_syn_ids = aa$id[which(aa$specimenID == samp)]
	for(samp_syn_id in samp_syn_ids){
		synGet(entity = samp_syn_id,down_dir = samp_dir)
	}
	cat("\n")
}

# SNP6 annotation CSV file
tmp_link = "http://www.affymetrix.com/Auth/analysis"
tmp_link = file.path(tmp_link,"downloads/na35/genotyping")
tmp_link = file.path(tmp_link,"GenomeWideSNP_6.na35.annot.csv.zip")
system(sprintf("wget %s",tmp_link))

```

</details>

## Blueprint

<details>
<summary>Click to expand!</summary>

* Download metadata and BAMs with 
	[Pyega3](https://github.com/EGA-archive/ega-download-client)
* `cred_file.json` contains user login information
* To obtain metadata,

	```Shell

	# for example
	dataset=EGAD00001002663 # genotypes
	# dataset=EGAD00001002671 # purified cell type 1
	# dataset=EGAD00001002674 # purified cell type 2
	# dataset=EGAD00001002675 # purified cell type 3

	# Download metadata code
	pyega3 -cf cred_file.json files $dataset
	```

* To obtain BAM files one by one,

	```Shell
	# num threads/cores
	nt=1

	# Example file id per BAM
	id=EGAF00001330176

	# Download code
	pyega3 -cf cred_file.json -c $nt fetch $id
	```

</details>

## BAM workflow

<details>
<summary>Click to expand!</summary>

* BamToFastq
	
	<details>
	<summary>Click to expand!</summary>
	
	For paired-end reads
	
	```Shell
	java -Xmx4g -jar picard.jar SamToFastq \
		INPUT=input.bam FASTQ=output_1.fastq.gz \
		SECOND_END_FASTQ=output_2.fastq.gz \
		UNPAIRED_FASTQ=output_unpair.fastq.gz \
		INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
	```
	
	For single-end reads

	```Shell
	java -Xmx4g -jar picard.jar SamToFastq \
		INPUT=input.bam FASTQ=output.fastq \
		INCLUDE_NON_PF_READS=true \
		VALIDATION_STRINGENCY=SILENT
	```
	
	</details>

* Build STAR index
	
	<details>
	<summary>Click to expand!</summary>

	```Shell
	fasta_fn=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
	gtf_fn=gencode.v26.GRCh38.ERCC.genes.gtf
	
	STAR --runMode genomeGenerate \
		--genomeDir ./star_index \
		--genomeFastaFiles $fasta_fn \
		--sjdbGTFfile $gtf_fn \
		--sjdbOverhang 99 --runThreadN 1
	```

	</details>

* Run STAR for read alignment
	
	<details>
	<summary>Click to expand!</summary>
	
	For paired-end reads
	
	```Shell
	STAR --runMode alignReads \
		--runThreadN 1 --genomeDir ./star_index --twopassMode Basic \
		--outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
		--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
		--outFilterType BySJout --outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 \
		--readFilesIn output_1.fastq.gz output_2.fastq.gz \
		--readFilesCommand zcat --outFileNamePrefix output_hg38 \
		--outSAMstrandField intronMotif --outFilterIntronMotifs None \
		--alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM \
		GeneCounts --outSAMtype BAM Unsorted --outSAMunmapped Within \
		--genomeLoad NoSharedMemory --chimSegmentMin 15 \
		--chimJunctionOverhangMin 15 --chimOutType WithinBAM SoftClip \
		--chimMainSegmentMultNmax 1 --outSAMattributes NH HI AS nM NM ch \
		--outSAMattrRGline ID:rg1 SM:sm1
	```

	For single-end reads
	
	```Shell
	STAR --runMode alignReads \
		--runThreadN 1 --genomeDir ./star_index --twopassMode Basic \
		--outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
		--outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
		--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
		--outFilterType BySJout --outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 \
		--readFilesIn output.fastq.gz --readFilesCommand zcat \
		--outFileNamePrefix output_hg38 --outSAMstrandField intronMotif \
		--outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes \
		--quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted \
		--outSAMunmapped Within --genomeLoad NoSharedMemory \
		--chimSegmentMin 15 --chimJunctionOverhangMin 15 \
		--chimOutType WithinBAM SoftClip --chimMainSegmentMultNmax 1 \
		--outSAMattributes NH HI AS nM NM ch \
		--outSAMattrRGline ID:rg1 SM:sm1
	```
	
	</details>

* MarkDuplicates
	
	<details>
	<summary>Click to expand!</summary>
	
	```Shell
	# Sort reads by coordinate
	samtools sort -@ 0 -o output_hg38.sortedByCoordinate.bam \
		output_hg38Aligned.out.bam
	
	# Make bam index
	samtools index -b -@ 1 output_hg38.sortedByCoordinate.bam
	
	# MarkDuplicates
	java -Xmx4g -jar picard.jar MarkDuplicates \
		INPUT=output_hg38.sortedByCoordinate.bam \
		OUTPUT=output_hg38.sortedByCoordinate.md.bam \
		M=out.marked_dup_metrics.txt ASSUME_SORT_ORDER=coordinate
	```
	
	Process post-markduplicate sorted bam
	
	```Shell
	# Count num reads in bam
	samtools view -c -@ 0 output_hg38.sortedByCoordinate.md.bam
	
	# Re-index bam
	samtools index -b -@ 1 output_hg38.sortedByCoordinate.md.bam
	```
	
	</details>

* Get TReC and ASReC
	
	Download asSeq source package
	
	<details>
	<summary>Click to expand!</summary>
	
	```Shell
	url=https://github.com/Sun-lab/asSeq/raw
	url=$url/master/asSeq_0.99.501.tar.gz
	
	wget $url
	```
	
	</details>
	
	Install R package and run `asSeq` to get unique reads and filter
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	install.packages(pkgs = "asSeq_0.99.501.tar.gz",
		type = "source",repos = NULL)
	
	PE = TRUE 
		# set TRUE for paired-end samples
		# set FALSE for single-end
	
	flag1 = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
		isSecondaryAlignment = FALSE,isDuplicate = FALSE,
		isNotPassingQualityControls = FALSE,
		isSupplementaryAlignment = FALSE,isProperPair = PE)
	
	param1 = Rsamtools::ScanBamParam(flag = flag1,
		what = "seq",mapqFilter = 255)
	
	bam_file = "output_hg38.sortedByCoordinate.md.bam"
	bam_filt_fn = "output.filtered.asSeq.bam"
	Rsamtools::filterBam(file = bam_file,
		destination = bam_filt_fn,
		param = param1)
	```
	
	</details>
	
	Create exon image file
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	gtf_fn = "gencode.v26.GRCh38.ERCC.genes.gtf.gz"
	exdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,
		format = "gtf")
	exons_list_per_gene = GenomicFeatures::exonsBy(exdb,
		by = "gene")
	
	gtf_rds_fn = "exon_by_genes_gencode_v26.rds"
	saveRDS(exons_list_per_gene,gtf_rds_fn)
	```
	
	</details>
	
	Get total read count (TReC)
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	genes = readRDS(gtf_rds_fn)
	bamfile = Rsamtools::BamFileList(bam_filt_fn,
		yieldSize = 1000000)
	se = GenomicAlignments::summarizeOverlaps(features = genes,
		reads = bamfile,mode = "Union",singleEnd = !PE,
		ignore.strand = TRUE,fragments = PE)
	ct = as.data.frame(SummarizedExperiment::assay(se))
	```
	
	</details>
	
	Filter reads by Qname
	
	<details>
	<summary>Click to expand!</summary>
	
	```Shell
	samtools sort -n -o output.filtered.asSeq.sortQ.bam \
		output.filtered.asSeq.bam
	```
	
	</details>
	
	Extract allele-specific reads, outputs hap1.bam, hap2.bam, hapN.bam
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	het_snp_fn = "<tab delimited filename of heterozygous SNPs for sample>"
		# Columns: chr, position, hap1 allele, hap2 allele
		# no header
	
	bam_filt_sortQ_fn = "output.filtered.asSeq.sortQ"
	asSeq::extractAsReads(input = sprintf("%s.bam",bam_filt_sortQ_fn),
		snpList = het_snp_fn,min.avgQ = 20,min.snpQ = 20)
	```
	
	</details>
	
	Count allele-specific read counts (ASReC)
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	se1 = GenomicAlignments::summarizeOverlaps(features = genes,
		reads = sprintf("%s_hap1.bam",bam_filt_sortQ_fn),mode = "Union",
		singleEnd = !PE,ignore.strand = TRUE,fragments = PE)
	se2 = GenomicAlignments::summarizeOverlaps(features = genes,
		reads = sprintf("%s_hap2.bam",bam_filt_sortQ_fn),mode = "Union",
		singleEnd = !PE,ignore.strand = TRUE,fragments = PE)
	seN = GenomicAlignments::summarizeOverlaps(features = genes,
		reads = sprintf("%s_hapN.bam",bam_filt_sortQ_fn),mode = "Union",
		singleEnd = !PE,ignore.strand = TRUE,fragments = PE)
	```
	
	</details>
	
	Save read counts
	
	<details>
	<summary>Click to expand!</summary>
	
	```R
	ct1 = as.data.frame(SummarizedExperiment::assay(se1))
	ct2 = as.data.frame(SummarizedExperiment::assay(se2))
	ctN = as.data.frame(SummarizedExperiment::assay(seN))
	cts = cbind(ct,ct1,ct2,ctN) # trec, hap1, hap2, hapN
	dim(cts); cts[1:2,]
	out_fn = "output.trecase.txt"
	write.table(cts,file = out_fn,quote = FALSE,
		sep = "\t", eol = "\n")
	```
	
	</details>

</details>

## Deconvolution

<details>
<summary>Click to expand!</summary>

* Signature expression derived from single cell RNAseq
	* [Middle Temporaral Gyrus](https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq)
	* [Pipeline/Workflow](https://github.com/Sun-lab/scRNAseq_pipelines/tree/master/MTG)
		to perform quality control on samples and genes, perform 
		dimension reduction, clustering, and compare cluster assignments 
		to existing cell type labeling
	* [Downstream differential expression analysis](https://github.com/Sun-lab/scRNAseq_pipelines/tree/master/_brain_cell_type)

* Comments:
	* Transcripts per million (TPM), gene count normalized 
	by exonic gene length and then all genes are normalized by 
	total normalized gene count, then multipled by 1 million
	* Cell sizes: For brain, summation of gene count normalized 
	by exonic gene length. For blood, obtained from EPIC and 
	ICeDT papers.
	
* Deconvolution Inputs
	* Bulk RNA-seq TPM (`bulk_tpm`),
	* scRNA-seq TPM (`sig_tpm`),
	* cell sizes

* Template code
	
	Input variables and objects. Make sure the genes are 
	ordered consistently between `sig_tpm` and `bulk_tpm`
	
	```R
	work_dir = "." # working directory
	setwd(work_dir)
	
	sig_tpm 	# TPM signature expression matrix
	bulk_tpm 	# TPM bulk expression matrix
	
	```
	
	CIBERSORT: [[Paper](https://www.nature.com/articles/nmeth.3337), 
		[Software](https://cibersort.stanford.edu/)]
	
	```R
	sig_fn = file.path(work_dir,"signature.txt")
	mix_fn = file.path(work_dir,"mixture.txt")
	write.table(cbind(rowname=rownames(sig_tpm),sig_tpm),
		file = sig_fn,sep = "\t",quote = FALSE,row.names = FALSE)
	write.table(cbind(rowname=rownames(bulk_tpm),bulk_tpm),
		file = mix_fn,sep = "\t",quote = FALSE,row.names = FALSE)
	
	source("CIBERSORT.R") # obtained from CIBERSORT website
	results = CIBERSORT(sig_matrix = sig_fn,mixture_file = mix_fn,
		perm = 0,QN = FALSE,absolute = FALSE,abs_method = 'sig.score',
		filename = "DECON")
	unlink(x = c(sig_fn,mix_fn))
	ciber_fn = sprintf("CIBERSORT-Results_%s.txt","DECON")
	unlink(ciber_fn)
	QQ = ncol(sig_tpm)
	
	# Extract proportion of transcripts per cell type per sample
	pp_bar_ciber = results[,seq(QQ)]
	
	# Calculate proportion of cell types per sample
	pp_hat_ciber = t(apply(pp_bar_ciber,1,function(xx){
		yy = xx / cell_sizes; yy / sum(yy)
	}))
	
	```

	ICeDT: [[Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7529339/),
		[Software](https://github.com/Sun-lab/ICeDt)]
	
	```R
	fit = ICeDT::ICeDT(Y = bulk_tpm,Z = sig_tpm,
		tumorPurity = rep(0,ncol(bulk_tpm)),refVar = NULL,
		rhoInit = NULL,maxIter_prop = 4e3,
		maxIter_PP = 4e3,rhoConverge = 1e-2)
	
	# Extract proportion of transcripts per cell type per sample
	pp_bar_icedt = t(fit$rho)[,-1]
	
	# Calculate proportion of cell type per sample
	pp_hat_icedt = t(apply(pp_bar_icedt,1,function(xx){
		yy = xx / cell_sizes; yy / sum(yy)
	}))
	
	```

	</details>

## eQTL mapping

* **BULK** mode: If running bulk analyses, set `RHO` to a 
matrix with `N` rows and 1 column and set `XX_trecPC` to 
residual TReC PCs calculated **without accounting** for cell types. 
* **Cell type-specific (CTS)** mode: If running cell type-specific 
analyses, set `RHO` to estimated cell type proportions and 
set `XX_trecPC` to **proportion-adjusted** residual TReC PCs.
* Coding
	
	Inputs
	
	```R
	devtools::install_github("pllittle/smarter")
	devtools::install_github("pllittle/CSeQTL")
	
	# Inputs
	XX_base 	# baseline covariates, continuous variables centered
	XX_genoPC # genotype PCs, centered
	XX_trecPC # residual TReC PCs, centered
	XX = cbind(Int = 1,XX_base,XX_genoPC,XX_trecPC)
	
	TREC 	# TReC vector
	SNP		# phased genotype vector
	hap2	# 2nd haplotype counts
	ASREC # total haplotype counts = hap1 + hap2
	PHASE # Indicator vector of whether or not to use haplotype counts
	RHO 	# cell type proportions matrix
	trim 	# TRUE for trimmed analysis, FALSE for untrimmed
	```
	
	eQTL mapping for one gene
	
	```R
	gs_out = CSeQTL_GS(XX = XX,TREC = TREC,SNP = SNP,hap2 = hap2,
		ASREC = ASREC,PHASE = PHASE,RHO = RHO,trim = trim,
		thres_TRIM = 20,numAS = 5,numASn = 5,numAS_het = 5,
		cistrans = 0.01,ncores = 1,show = TRUE)
	
	```
	
	Hypothesis testing output. Matrices where rows are genomic loci 
	and columns are cell types for **CTS** mode or a single column for 
	**BULK** mode.
	
	```R
	# TReC-only likelihood ratio test statistics with 1 DF
	gs_out$LRT_trec
	
	# TReC+ASReC likelihood ratio test statistics with 1 DF
	gs_out$LRT_trecase
	
	# Cis/Trans likelihood ratio test statistics with 1 DF
	gs_out$LRT_cistrans
	
	```

###

