
# Filter and translate the Yazar genes considered for testing for each cell type

# (1) only keep the genes existing gene location files sent by Yazar et al. 
# (2) translate Yazar genes from gene names to gene ensemblID, 
#     using a processed file from hg19 gene reference


library(data.table)
library(readxl)
library(stringr)
is.unq <- function(x) length(x)==length(unique(x))


data_dir = "../../CSeQTL/Yazar_2022/"
Yazar_2022_dir   = "science.abf3041_tables_s6_to_s19/"
Yazar_2022_table = "science.abf3041_tables_s6_to_s19.xlsx"
geneloc_dir = "OneK1K_cohort_cis_eQTL_input_files/Gene_Location_Files/"
hg19_symbol_ensemblID_file = "hg19_gene_symbol_ensemblID.txt"


# Load gene location files sent by Yazar et al. 

col_names = c("geneid", "chr", "start", "end", "strand")

Yazar_gene_loc = data.frame(matrix(ncol=5, nrow=0, 
                                   dimnames=list(NULL, col_names)))

for (i in c(1:22)){
  gene_loc_i = read.csv(paste0(data_dir, geneloc_dir, 
                               "geneloc_chr", as.character(i), ".tsv"), 
                        sep = "\t", header=TRUE)
  
  Yazar_gene_loc = rbind(Yazar_gene_loc, gene_loc_i)
}

dim(Yazar_gene_loc)


# read hg19 gene symbol and ensemblID

hg19_symbol_ensembID <- read.delim(paste0(data_dir, 
                                          hg19_symbol_ensemblID_file))
dim(hg19_symbol_ensembID)
head(hg19_symbol_ensembID)

# all gene names in Yazar_gene_loc fall into hg19_symbol_ensembID
table(Yazar_gene_loc$geneid %in% hg19_symbol_ensembID$GeneSymbol)


# Filter and translate

Yazar_cell_types = c("B IN", "B Mem", "CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 NC", 
                     "CD8 ET", "CD8 S100B", "DC", "Mono C", "Mono NC", "NK R",     
                     "NK", "Plasma")

# load file for Yazar 2022 genes considered for testing
Yazar_2022_genes <- read_excel(paste0(data_dir, Yazar_2022_dir, Yazar_2022_table), 
                               sheet = "Table.S8", skip = 2)
colSums(!is.na(Yazar_2022_genes))

n_con_Yazar_raw_list = list()
n_con_Yazar_loc_list = list()
n_con_Yazar_list     = list()
gene_con_Yazar_list  = list()

id_unique_flag_v = rep(NA, length(Yazar_cell_types))

for (i in 1:length(Yazar_cell_types)){
  
  ct = Yazar_cell_types[i]
  ct_string = sub(" ", "_", ct)
  
  Yazar_ct_eGenes_path = paste0(data_dir, "Yazar_eSNP1_50kb/", 
                                "Yazar_eSNP1_50kb_", ct_string, ".csv")
  df_Yazar_ct_eGenes = read.csv(Yazar_ct_eGenes_path, header = TRUE)
  
  Yazar_ct_eGene_names = df_Yazar_ct_eGenes$Gene.ID
  
  # verify whether all kept Yazar eGenes have name belonging to the genes considered
  stopifnot(all(Yazar_ct_eGene_names %in% Yazar_2022_genes[[ct]]))
  
  # filter out genes that have name not existing in any of the gene_loc files
  ct_loc_Yazar = intersect(Yazar_2022_genes[[ct]], Yazar_gene_loc$geneid)
  # it has been verified in earlier  that all gene names in 
  # Yazar_gene_loc fall into hg19_symbol_ensembID
  # so all genes in ct_loc_Yazar also fall into hg19_symbol_ensembID 
  mat1 = match(ct_loc_Yazar, hg19_symbol_ensembID$GeneSymbol)
  ct_considered_Yazar = hg19_symbol_ensembID[mat1,]$gene_id
  
  n_con_Yazar_raw_list[[ct_string]] = sum(!is.na(Yazar_2022_genes[[ct]]))
  n_con_Yazar_loc_list[[ct_string]] = length(ct_loc_Yazar)
  n_con_Yazar_list[[ct_string]]     = length(ct_considered_Yazar)
  gene_con_Yazar_list[[ct_string]]  = ct_considered_Yazar
  
  # see whether different gene symbols have different ensembl id 
  id_unique_flag_v[i] = is.unq(ct_considered_Yazar)
  
}

as.numeric(unlist(n_con_Yazar_raw_list))
as.numeric(unlist(n_con_Yazar_loc_list))
as.numeric(unlist(n_con_Yazar_list))

id_unique_flag_v

saveRDS(gene_con_Yazar_list, 
        file = paste0(data_dir, "Yazar_gene_considered.rds"))

gc()

sessionInfo()
q(save="no")
