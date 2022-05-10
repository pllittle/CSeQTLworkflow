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

* Synapse R code

### Blueprint

* Download metadata and BAMs with [Pyega3](https://github.com/EGA-archive/ega-download-client)
* `cred_file.json` contains user login information
* To obtain metadata,

```Shell

dataset=EGAD00001002663 # for example
# dataset=EGAD00001002671
# dataset=EGAD00001002674
# dataset=EGAD00001002675

# Download metadata code
pyega3 -cf cred_file.json files $dataset \
	| grep "EGAF" | grep -v "cram" | tr -s ' ' | sed 's/^ //g' \
	| sed 's/ /\t/g' > file_metadata
```

* To obtain bams one by one,

```Shell
# num threads/cores
nt=1

# Example file id per BAM
id=EGAF00001330176

# Download code
pyega3 -cf cred_file.json -c $nt fetch $id
```

