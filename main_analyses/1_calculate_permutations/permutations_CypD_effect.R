library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(limma)
library(edgeR)

current = getwd()
folder = "1_calculate_permutations"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))
system(paste0("mkdir -p outputs/",folder,"/by_permutation_outs/"))

## inputs
inputs_dir <- paste0("inputs/",folder,"/")
num_permutations <- c(1:10)
set.seed(2020)

## outputs
results_dir <- paste0("outputs/",folder,"/")

## meta data
meta_data <- read.csv(paste0(inputs_dir,"meta_data.csv"), header = TRUE)
meta_data$infection_status = factor(meta_data$infection_status, levels = c("NI","flu"))
meta_data$mouse_type = factor(meta_data$mouse_type, levels = c("WT","CypDko"))
meta_data$M_aligned_scale <- scale(meta_data$M_aligned)
rownames(meta_data) <- meta_data$sample_ID

for (j in 1:length(num_permutations)){

	perm_number <- num_permutations[j]

	## read in gene expression
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

	## PERMUTE GENOTYPE WITHIN CONDITION	
	meta_data_NI <- meta_data[which(meta_data$infection_status == "NI"),]
	meta_data_flu <- meta_data[which(meta_data$infection_status == "flu"),]
	meta_data_NI$mouse_type_perm <- sample(meta_data_NI$mouse_type)
	meta_data_flu$mouse_type_perm <- sample(meta_data_flu$mouse_type)

	## add back into one dataframe
	meta_data_PERM <- rbind(meta_data_NI, meta_data_flu)
	meta_data_PERM <- meta_data_PERM[colnames(vv$E),]

	length(which(colnames(vv$E)!=rownames(meta_data_PERM)))
	length(which(rownames(vv$targets)!=rownames(meta_data_PERM)))

	## fit the model with permuted data
	design = model.matrix(~ M_aligned_scale + infection_status + mouse_type_perm:infection_status, data = meta_data_PERM)
	vfit <- lmFit(vv, design)
	vfit <- eBayes(vfit)

	## get permuted p-values
	p_values_NI = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("infection_statusNI:mouse_type_permCypDko"))]); colnames(p_values_NI)[1] <- paste0("pvalues_NI_",perm_number)
	p_values_NI$genes <- rownames(p_values_NI)
	p_values_flu = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("infection_statusflu:mouse_type_permCypDko"))]); colnames(p_values_flu)[1] <- paste0("pvalues_flu_",perm_number)
	p_values_flu$genes <- rownames(p_values_flu)

	write.table(p_values_NI, paste0(results_dir,"by_permutation_outs/permutation",perm_number,"_NI.txt"), quote = FALSE, sep = ",", row.names = FALSE)
	write.table(p_values_flu, paste0(results_dir,"by_permutation_outs/permutation",perm_number,"_flu.txt"), quote = FALSE, sep = ",", row.names = FALSE)

	if(j == 1){
		results_perm_NI <- cbind(p_values_NI)
		results_perm_flu <- cbind(p_values_flu)
	}else{
		results_perm_NI <- merge(results_perm_NI, p_values_NI, by = "genes")
		results_perm_flu <- merge(results_perm_flu, p_values_flu, by = "genes")
	}

	print(perm_number)
}

write.table(results_perm_NI, paste0(results_dir,"permutation_outs_NI.txt"), quote = FALSE, sep = ",", row.names = FALSE)
write.table(results_perm_flu, paste0(results_dir,"permutation_outs_flu.txt"), quote = FALSE, sep = ",", row.names = FALSE)
