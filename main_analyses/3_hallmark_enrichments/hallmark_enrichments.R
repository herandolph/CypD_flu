library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(fgsea)
library(qusage)
library(data.table)
library(biomaRt)

current = getwd()
folder = "3_hallmark_enrichments"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
inputs_dir <- paste0("inputs/",folder,"/")

## outputs
results_dir <- paste0("outputs/",folder,"/")

## input genesets
hallmark <- readRDS(paste0(inputs_dir,"Mm.h.all.v7.1.entrez.rds"))
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)

## read in effect of CypD-/- genotype results + subset on only t statistics in the flu condition
df <- read.table(paste0(inputs_dir,"resultsALL_with_qvalues.txt"), header = TRUE, sep = ",")
df$genes <- rownames(df)
df <- subset(df, select = c(genes, t_stat_flu))

## convert genes in this data to entrez ids instead of mgi
entrez <- unique(getBM(attributes = c("mgi_symbol","entrezgene_id"),    
                      filters = "mgi_symbol",
                      values = df$genes,
                      mart = ensembl))
colnames(entrez)[1] <- "genes"
entrez <- entrez[complete.cases(entrez),]

## add in the entrez ids to the data
df <- join(df, entrez, by = "genes")
t_statistics <- df[!duplicated(df$genes),]
t_statistics <- df[!duplicated(df$entrezgene_id),]
t_statistics <- t_statistics[complete.cases(t_statistics),]
rownames(t_statistics) <-t_statistics$entrezgene_id; t_statistics$genes <- NULL; t_statistics$entrezgene_id <- NULL
colnames(t_statistics)[1] <- "t_stat"

## run enrichments for the hallmark geneset
run_GSEA <- function(geneset, geneset_name){

	set.seed(2020)
	## calculate gene set enrichment analysis (gsea)
	## take the t stats of condition i and rank them
	t_i <- t_statistics
	t_i$gene <- rownames(t_i)
	t_i <- data.table(t_i[,c(2,1)])
	t_i_rank <- t_i[order(t_stat), ]

	input <- setNames(t_i_rank$t_stat, t_i_rank$gene)

	fgseaRes <- fgsea(pathways = geneset, 
		                  stats = input,
		                  minSize=15,
		                  maxSize=500,
		                  nperm=100000)

	fgseaRes_csv <- data.frame(lapply(fgseaRes, as.character), stringsAsFactors = FALSE)
	write.csv(fgseaRes_csv, file = paste0(results_dir,geneset_name,".csv"))
}

run_GSEA(hallmark, "hallmark")
