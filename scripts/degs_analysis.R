library(limma)
library(AnnotationDbi)
library(pheatmap)
library(org.Rn.eg.db)
source("plots_utils.R")
source("degs_utils.R")


exprs <- read.table("../exprs/exprs_statin.tsv", header = TRUE, sep = "\t")
pd <- read.table("../pdata/pdata_statin.tsv", header = TRUE, sep = "\t", stringsAsFactors = F)
rownames(pd) <- pd$SampleID

rownames(exprs) <- exprs$ENTREZID

outlier=TRUE
if (outlier) {
  outsample = "SDd50_17023"
  exprs <- exprs[,-which(colnames(exprs) %in% outsample)]
  pd <- pd[-which(pd$SampleID %in% outsample),]
  outname="_nooutlier"
}


# Get differentially expressed genes
pd$Group <- factor(pd$Group, levels=c("SD", "PE"))
pd$Day <- factor(pd$Day, levels=c("d21", "d50"))
GD <- paste(pd$Group, pd$Day, sep=".")
GD <- factor(GD, levels=c("SD.d21", "SD.d50", "PE.d21", "PE.d50"))
design <- model.matrix(~0+GD)
colnames(design) <- levels(GD)
fit <- lmFit(exprs[,-c(1,2)], design)
cont.matrix <- makeContrasts(
                SDd21vsd50=SD.d21-SD.d50,
                PEd21vsd50=PE.d21-PE.d50,
                d21SDvsPE=SD.d21-PE.d21,
                d50SDvsPE=SD.d50-PE.d50,
                levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get all the genes with logFC, p-values, no filtering
degs <- topTable(fit2, adjust.method="fdr", number=nrow(fit2))

coefs = c("SDd21vsd50", "PEd21vsd50", "d21SDvsPE", "d50SDvsPE")
d <- c()
i=4
degs <- topTable(fit2, adjust.method="fdr", number=nrow(fit2), coef=coefs[i])

# Filter
degs <- filterDEGS(degs, 0.05, 0.5)

# Add description
anno_name<-select(org.Rn.eg.db, rownames(degs), c("GENENAME", "SYMBOL"), keytype="ENTREZID")
degs$description <- anno_name$GENENAME
degs$SYMBOL <- anno_name$SYMBOL
degs$ENTREZID <- rownames(degs)

# Add FC and rearrange columns
degs$FC <- sign(degs$logFC)*(2^abs(degs$logFC))
degs <- degs[,c("SYMBOL", "ENTREZID","logFC", "FC", "P.Value", "adj.P.Val", "description")]

write.table(degs, paste("../degs/degs_", coefs[i],".tsv", sep=""),sep="\t", quote=F, row.names=FALSE)

d <- rbind(d, degs)
# Heatmap on significant genes
joint_degs <- unique(d$SYMBOL)

rownames(exprs) <- exprs$SYMBOL
exprs <- exprs[,-c(1,2)]
exprs.degs <- exprs[which(rownames(exprs) %in% joint_degs),]
ht_matrix <- exprs.degs

anno_col <- pd[,c("Group", "Day")]

pl <- pheatmap(ht_matrix, cluster_cols = T, cluster_rows = T, drop_levels = F,
               annotation_col = anno_col, scale = "row", treeheight_row = 40, treeheight_col = 20,
               cellheight = 4, cellwidth=17, border_color = NA,
               fontsize_row = 4, fontsize = 8, silent=T)
save_plot(paste("../plots/degs_statin_heatmap.pdf", sep=""),
          base_height=23, base_width=8, pl, ncol=1)
save_plot(paste("../plots/degs_statin_heatmap.svg", sep=""),
          base_height=23, base_width=8, pl, ncol=1)
