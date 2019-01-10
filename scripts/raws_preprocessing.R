library(oligo)
library(clariomdratrnentrezg.db)
library(arrayQualityMetrics)
library(factoextra)
source("plots_utils.R")

#install.packages("http://brainarray.mbni.med.umich.edu/bioc/src/contrib/pd.clariomdrat.rn.entrezg_23.0.0.tar.gz",
#                 repos = NULL, type = "source")
#install.packages("http://brainarray.mbni.med.umich.edu/bioc/src/contrib/clariomdratrnentrezg.db_23.0.0.tar.gz",
#                 repos = NULL, type = "source")

pd <- read.table("../pdata/pdata_statin.tsv", header = TRUE, sep = "\t")
celFiles <- paste("../raws/", pd$FileName, sep="")
rawData <- oligo::read.celfiles(filenames=celFiles, pkgname="pd.clariomdrat.rn.entrezg")

rownames(pd) <- sampleNames(rawData)
pdata <- AnnotatedDataFrame(data=pd)
phenoData(rawData) <- pdata

eset = oligo::rma(rawData)

# Extract expression matrix
exprs <- exprs(eset)
colnames(exprs) <- pd$SampleID
rownames(pd) <- pd$SampleID

# Get annotations
anno <- AnnotationDbi::select(clariomdratrnentrezg.db, rownames(exprs), c("ENTREZID", "SYMBOL", "GENENAME"))
anno <- anno[which(anno$ENTREZID!="NA"),]
n_occur <- data.frame(table(anno$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
anno <- anno[which(anno$PROBEID %in% uniques),]

# Replace expression matrix row names with gene symbol
exprs <- exprs[which(rownames(exprs) %in% anno$PROBEID),]
ind <- match(rownames(exprs), anno$PROBEID)
rownames(exprs) <- anno[ind,]$SYMBOL

e <- cbind(ENTREZID=anno[ind,]$ENTREZID, SYMBOL=anno[ind,]$SYMBOL, exprs)
write.table(e, "../exprs/exprs_statin.tsv",sep="\t", quote=F, row.names=F)

###############
exprs <- read.table("../exprs/exprs_statin.tsv", header = TRUE, sep = "\t")
pd <- read.table("../pdata/pdata_statin.tsv", header = TRUE, sep = "\t")
rownames(pd) <- pd$SampleID
rownames(exprs) <- exprs$SYMBOL
exprs <- exprs[,-c(1,2)]


outlier=TRUE
outname=""

if (outlier) {
  outsample = "SDd50_17023"
  exprs <- exprs[,-which(colnames(exprs) %in% outsample)]
  pd <- pd[-which(pd$SampleID %in% outsample),]
  outname="_nooutlier"
}

eset = ExpressionSet(assayData=as.matrix(exprs), phenoData = AnnotatedDataFrame(pd))
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/AQM_report", outname, sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Group"))

# Simple PCA
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, c("Day", "Group"))
save_plot(paste("../plots/statin_PCA", outname, ".pdf", sep=""),
          base_height=3.4, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)
save_plot(paste("../plots/statin_PCA", outname, ".svg", sep=""),
          base_height=3.4, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)

# PCA percentage of variance
pl <- fviz_eig(pca)
save_plot(paste("../plots/statin_PCA_var_percentage", outname, ".pdf", sep=""),
          base_height=3, base_width=6, pl, ncol=1)
save_plot(paste("../plots/statin_PCA_var_percentage", outname, ".svg", sep=""),
          base_height=3, base_width=6, pl, ncol=1)

# PCA loadings
pl <- ggloadings(pca, num_pca=1)
save_plot(paste("../plots/statin_PCA_loadings_1", outname, ".pdf", sep=""),
          base_height=3, base_width=8, pl, ncol=1)
save_plot(paste("../plots/statin_PCA_loadings_1", outname, ".svg", sep=""),
          base_height=3, base_width=8, pl, ncol=1)
pl <- ggloadings(pca, num_pca=2)
save_plot(paste("../plots/statin_PCA_loadings_2", outname, ".pdf", sep=""),
          base_height=3, base_width=8, pl, ncol=1)
save_plot(paste("../plots/statin_PCA_loadings_2", outname, ".svg", sep=""),
          base_height=3, base_width=8, pl, ncol=1)

# PCA quality of representation
pl <- fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE)
save_plot(paste("../plots/statin_PCA_cos", outname, ".pdf", sep=""),
          base_height=5, pl, ncol=1)
save_plot(paste("../plots/statin_PCA_cos", outname, ".svg", sep=""),
          base_height=5, pl, ncol=1)

# PCA loadings on 2D
pl <- fviz_pca_biplot(pca,
                      geom=c("point", "text"),
                      label="all",
                      col.var = "contrib",
                      fill.ind = factor(pd$Group),
                      #palette = indCol,
                      col.ind = "black",
                      pointshape=21, pointsize = 3,
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE,     # Avoid text overlapping
                      select.var = list(contrib=5),
                      invisible = "quali"
) + labs(fill = "Group")
save_plot(paste("../plots/statin_PCA_loadings", outname, ".pdf", sep=""),
          base_height=6, pl, ncol=1)
save_plot(paste("../plots/statin_PCA_loadings", outname, ".svg", sep=""),
          base_height=6, pl, ncol=1)
