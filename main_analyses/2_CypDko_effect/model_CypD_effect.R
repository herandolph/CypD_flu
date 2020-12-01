library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)
library(limma)
library(edgeR)

current = getwd()
folder = "2_CypDko_effect"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
inputs_dir <- paste0("inputs/",folder,"/")
setwd("common_functions")
source("permFDR.R")
setwd(current)

## outputs
results_dir <- paste0("outputs/",folder,"/")

## meta data
meta_data <- read.csv(paste0(inputs_dir,"meta_data.csv"), header = TRUE)
meta_data$infection_status = factor(meta_data$infection_status, levels = c("NI","flu"))
meta_data$mouse_type = factor(meta_data$mouse_type, levels = c("WT","CypDko"))
meta_data$M_aligned_scale <- scale(meta_data$M_aligned)
rownames(meta_data) <- meta_data$sample_ID

# permutations
permutations_NI <- read.table(paste0(inputs_dir,"permutation_outs_NI.txt"), header = TRUE, sep = ",")
permutations_flu <- read.table(paste0(inputs_dir,"permutation_outs_flu.txt"), header = TRUE, sep = ",")

## gene expression
GE <- read.table(paste0(inputs_dir,"Mus_musculus.mm10.81.kallisto.gene_level.lengthScaledTPM_counts_MGIsymbols.txt"), header = TRUE, sep = ",")

## subset on protein coding genes
coding_ids = read.table(paste0(inputs_dir,"Mouse.mm10.81.HUGOgene_transcriptAssoc.protein_coding.txt"), header = TRUE)
GE <- GE[which(rownames(GE) %in% coding_ids$gene_name), ]

## reorder colnames
GE <- GE[meta_data$mouse_ID]
colnames(GE) <- meta_data$sample_ID

### remove lowly expressed genes, keep median expression > 1
dge <- DGEList(counts = GE)
dge <- calcNormFactors(dge)
design <- model.matrix(~ M_aligned, data = meta_data)
v <- voom(dge, design, plot = FALSE)

tab = data.frame(genes = rownames(GE), medians=apply(v$E, 1, median), order = 1:nrow(GE))
tab = tab[order(-tab$medians), ]
tab$order_by_median = 1:nrow(tab)

## threshold at median = 1
tab = tab[order(tab$order), ]

length(which(rownames(GE) != rownames(tab)))
GE <- GE[which(tab$medians > 1), ]
## 12451 genes

## voom after removal of lowly expressed genes
dge <- DGEList(counts = GE)
dge <- calcNormFactors(dge)
vv <- voomWithQualityWeights(dge, design, plot = FALSE)

## fit the model
design = model.matrix(~ M_aligned_scale + infection_status + mouse_type:infection_status, data = meta_data)
vfit <- lmFit(vv, design)
vfit <- eBayes(vfit)

## get results
## beta for the infection effect
infection_betas_sign = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("infection_statusflu"))]); colnames(infection_betas_sign)[1] <- "infection_betas_sign"

## beta for the effect of the CypD-/- genotype within each condition
get_geno_res <- function(beta_col, name){
	betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% beta_col)]); colnames(betas)[1] <- paste0("betas_",name)
	p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% beta_col)]); colnames(p_values)[1] <- paste0("pvalues_",name)
	fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- paste0("fdrs_",name)
	out <- cbind(betas, p_values, fdrs)
	return(out)
}

NI_CypD <- get_geno_res("infection_statusNI:mouse_typeCypDko", "NI")
flu_CypD <- get_geno_res("infection_statusflu:mouse_typeCypDko", "flu")

## get t statistics
get_t_statistics <- function(beta_col, name){
	t_stats = as.data.frame(cbind(rownames(topTable(vfit, coef = beta_col, number = Inf)), topTable(vfit, coef = beta_col, number = Inf)$t))
	colnames(t_stats) <- c("genes",paste0("t_stat_",name))
	return(t_stats)
}

t_stats_NI <- get_t_statistics("infection_statusNI:mouse_typeCypDko", "NI")
t_stats_flu <- get_t_statistics("infection_statusflu:mouse_typeCypDko", "flu")

## bind all results
results <- cbind(infection_betas_sign, NI_CypD, flu_CypD)
results$genes <- rownames(results)
results <- join(results, t_stats_NI, by = "genes")
results <- join(results, t_stats_flu, by = "genes")
rownames(results) <- results$genes; results$genes <- NULL

## CALCULATE FDRS USING PERMUTED DATA
perm_fdrs = permFDR(full_data = results, full_column_id = "pvalues_NI", perm_data = permutations_NI, perm_column_ids = "all", output_name = paste0(results_dir))
fdrs_intermediate <- perm_fdrs$fdrs
colnames(fdrs_intermediate)[c(10,11)] <- paste0(colnames(fdrs_intermediate)[c(10,11)],"_NI")
perm_fdrs = permFDR(full_data = fdrs_intermediate, full_column_id = "pvalues_flu", perm_data = permutations_flu, perm_column_ids = "all", output_name = paste0(results_dir))
fdrs_final <- perm_fdrs$fdrs
colnames(fdrs_final)[c(12,13)] <- paste0(colnames(fdrs_final)[c(12,13)],"_flu")

## write results
write.table(fdrs_final, file = paste0(results_dir,"resultsALL_with_qvalues.txt"), quote = FALSE, sep = ",")

## write genelists
genes_10 <- as.character(rownames(fdrs_final[fdrs_final$Fdr_SAB_flu < 0.10 & abs(fdrs_final$betas_flu) > 0.5, ]))
genes_10_higherCypD <- as.character(rownames(fdrs_final[fdrs_final$Fdr_SAB_flu < 0.10 & fdrs_final$betas_flu > 0.5, ]))
genes_10_higherWT <- as.character(rownames(fdrs_final[fdrs_final$Fdr_SAB_flu < 0.10 & fdrs_final$betas_flu < -0.5, ]))

write.table(genes_10, paste0(results_dir,"significant_genes_fdr0.10_with_qvalues.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_10_higherCypD, paste0(results_dir,"significant_genes_fdr0.10_higherCypD_with_qvalues.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_10_higherWT, paste0(results_dir,"significant_genes_fdr0.10_higherWT_with_qvalues.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(fdrs_final), paste0(results_dir,"background_geneset.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
