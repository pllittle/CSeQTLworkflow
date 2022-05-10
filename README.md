# CSeQTLworkflow

This repository contains template codes for 
runnings workflow steps.

## Links

* [GTEx-pipeline for RNA-seq](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md)

## To Dos

Publish template codes for pipeline

* Where datasets obtained and how
* GTEx pipeline
* asSeq: TReC, ASReC
* Deconvolution
* eQTL mapping code
* Torus enrichment
* GWAS enrichment

## Template Codes

## CommonMind Consortium

* Code to obtain CMC inputs
	<details>
	<summary>Click to expand!</summary>

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

* Download metadata and BAMs with [Pyega3](https://github.com/EGA-archive/ega-download-client)
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

## BAM workflow

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
	
	</details>

* Get TReC and ASReC
	
	Download asSeq source package
	
	```Shell
	url=https://github.com/Sun-lab/asSeq/raw
	url=$url/master/asSeq_0.99.501.tar.gz
	
	wget $url
	```
	
	Install R package
	
	```R
	install.packages(pkgs = "asSeq_0.99.501.tar.gz",
		type = "source",repos = NULL)
	```

	
###

