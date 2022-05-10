# CSeQTLworkflow

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

## Obtaining Datasets

### CommonMind Consortium

* Code to obtain CMC inputs

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

### Blueprint

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

