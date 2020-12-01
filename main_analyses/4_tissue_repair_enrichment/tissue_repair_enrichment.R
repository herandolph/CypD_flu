library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)

current = getwd()
folder = "4_tissue_repair_enrichment"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
inputs_dir <- paste0("inputs/",folder,"/")

## outputs
results_dir <- paste0("outputs/",folder,"/")

## read in TiRe genes
tire <- as.data.frame(read.csv(paste0(inputs_dir,"mouse_entries_TIREdatabase.csv")))
tire <- as.character(tire$Gene.symbol)
tire <- unique(tire[tire != ""])

## read in differentially-expressed genes in WT and CypD-/- genotypes in the flu condition
higherKO <- unique(as.character(read.table(paste0(inputs_dir,"significant_genes_fdr0.10_higherCypD_with_qvalues.txt"), header = FALSE)$V1))
higherWT <- unique(as.character(read.table(paste0(inputs_dir,"significant_genes_fdr0.10_higherWT_with_qvalues.txt"), header = FALSE)$V1))
background_geneset <- unique(as.character(read.table(paste0(inputs_dir,"background_geneset.txt"), header = FALSE)$V1))

tire <- tire[tire %in% background_geneset]

enrich <- function(genelist, name, name2){

		set.seed(2020)
		## calculate the number of observed genes
		obs_num <- length(genelist)

		## calculate the percentage of genes that are in tire list
		tire_num <- length(genelist[genelist %in% tire])

		## calculate percentage
		obs_percentage <- tire_num/obs_num

		## calculate a null distribution
		npermutations <- 1000
		for(i in 1:npermutations){

			## sample number rows equal to the number of unique genes
			null_vec <- sample(background_geneset, obs_num)

			## calcualte the number these genes that are in the tire dataset
			null_num <- length(null_vec[null_vec %in% tire])

			## null percent
			null_percentage <- null_num/obs_num

			if(i == 1){
				permutation_outs <- null_percentage
			}else{
				permutation_outs <- rbind(permutation_outs, null_percentage)
			}
		}

		permutation_outs <- as.data.frame(permutation_outs)
		colnames(permutation_outs) <- "null"
		equal_or_greater <- as.numeric(table(permutation_outs$null >= obs_percentage)["TRUE"])
		permutation_outs$pval <- equal_or_greater/npermutations
		permutation_outs$population <- paste0(name2)
		permutation_outs$obs_percentage <- obs_percentage * 100
		permutation_outs$null <- permutation_outs$null * 100
		return(permutation_outs)
}

KO <- enrich(higherKO, "higher_KO", "KO")
WT <- enrich(higherWT, "higher_WT", "WT")

combine <- rbind(KO, WT)
combine$population <- factor(combine$population, levels = c("KO","WT"))
write.table(combine, paste0(results_dir,"tissue_repair_enrichment_results.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

## plot 
plot <- ggplot() +
			geom_density(data = subset(combine, population == "KO"), aes(null), color = "#228B22", fill = "#228B22", alpha = 0.1, adjust = 2) +
			geom_density(data = subset(combine, population == "WT"), aes(null), color = "#0044d0", fill = "#0044d0", alpha = 0.1, adjust = 2) +
			geom_vline(data = subset(combine, population == "KO"), aes(xintercept = obs_percentage), color = "#228B22", linetype = "dashed") +
			geom_vline(data = subset(combine, population == "WT"), aes(xintercept = obs_percentage), color = "#0044d0", linetype = "dashed") +
			scale_fill_manual(values = c("#0044d0","#228B22")) +
			theme_bw() +
			ggtitle(paste0("p KO = ", unique(subset(combine, population == "KO")$pval),"\np WT = ", unique(subset(combine, population == "WT")$pval))) +
			theme(plot.title = element_text(size = 8.5),
				  axis.title.x = element_text(size = 8.5),
				  axis.title.y = element_text(size = 8.5),
				  axis.text.x = element_text(size = 7.5),
				  axis.text.y = element_text(size = 7.5), 
				  panel.grid.minor = element_blank(),
				  panel.grid.major = element_blank(),
				  panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
				  legend.title = element_blank(),
				  legend.position = "none") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
			xlab("% genes in TiRe set")


pdf(paste0(results_dir,"TiRe_enrichments.pdf"), width = 4, height = 2.5)
print(plot)
dev.off()
