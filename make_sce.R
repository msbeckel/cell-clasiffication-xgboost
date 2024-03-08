#############################################################################################
#Aim: Make sce object from neuronal-cell-type-classification and make preprocess analysis
#Author: Maximiliano S. Beckel
#Creación: 2023/20/20
#Última modificación: 2023/28/20
#############################################################################################

library("edgeR")
library("ggplot2")
library("matrixStats")
library(data.table)
library(ggpubr)
library(stringr)
library(SingleCellExperiment)
library(scuttle)
library(robustbase)
library(ggplot2)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(Rmagic)
library(ALRA)
library(SAVER)
#library(DescTools)
memory.limit(60000)
options(bitmapType='cairo')
setwd("/home/maxibeckel/maestria_datos/neuronal-cell-type-classification/")

# Loads in data
intron_matrix <- fread("./datasets/human_MTG_2018-06-14_intron-matrix.csv")
exon_matrix <- fread("./datasets/human_MTG_2018-06-14_exon-matrix.csv")
full_matrix <- rbindlist(list(intron_matrix, exon_matrix))[, lapply(.SD, sum, na.rm = TRUE), by = V1]
#full_matrix <- intron_matrix + exon_matrix

# ------------------------------- External CSVs ------------------------------- #

# Get cell labels from external csv
# NOTE: This code was for only the highest cell subtype
# labels <- read.csv("final_two_sample_columns.csv")
# cell.types <- c("Exc", "Inh", "Oli", "Ast", "OPC", "End", "Mic")

# Create cluster columns
samples <- fread("./datasets/human_MTG_2018-06-14_samples-columns.csv")
labels <- samples[,.(sample_name, cluster)]
labels <- labels[cluster != "no class",]
labels <- cbind(labels, str_split_fixed(labels$cluster, " ", n=Inf))
colnames(labels) <- c("sample_name", "cluster", "higher", "layer", "intermediate", "granular")

cell.types <- unique(labels$cluster)

# Gets list of column indices for a particular cell type to subset
for(type in cell.types) {
  t(assign(paste0(type, ".celltypes"), labels[which(labels$clu == type),1]))
}

# Read in gene.rows file
gene.rows <- fread("./datasets/human_MTG_2018-06-14_genes-rows.csv")

# Replace row names of full matrix with genes
#rownames(full_matrix) <- gene.rows[,gene]

# ------------------------------- Build SCE ------------------------------- #

sce <- SingleCellExperiment(assays  = list(counts = as.matrix(full_matrix[,-1])), 
                            rowData = gene.rows, 
                            colData = samples)
#elimino celulas no asignadas a un cluster
sce <- sce[,sce$cluster != "no class"]

# Identifying the mitochondrial transcripts in our SingleCellExperiment.
is.mito <- rowData(sce)$chromosome == "MT"

df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

stats <- cbind(log10(df$sum), log10(df$detected),
               df$subsets_Mito_percent)

colData(sce) <- cbind(colData(sce), df)

#robustbase
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)
colData(sce)[, "multi.outlier"] <- as.logical(multi.outlier)

if(FALSE){

  p1 <- plotColData(sce, x="donor", y="sum", colour_by="multi.outlier") + 
    scale_y_log10() + ggtitle("Cuentas totales")+theme(axis.title.y=element_blank())
  p2 <- plotColData(sce, x="donor", y="detected", colour_by="multi.outlier") + 
    scale_y_log10() + ggtitle("Features detectados")+theme(axis.title.y=element_blank())
  p3 <- plotColData(sce, x="donor", y="subsets_Mito_percent", colour_by="multi.outlier") + 
    scale_y_log10() + ggtitle("% Mito")+theme(axis.title.y=element_blank())
  
  ggarrange(p1, p2, p3, ncol=1, nrow=3, common.legend = TRUE, legend="right")
}

sce <- sce[,!multi.outlier]

if(!dir.exists(paste0(getwd(), "/datasets/scimpute_results"))){
  sce_csv <- assay(sce)
  rownames(sce_csv) <- rowData(sce)[, "gene"]
  write.csv(sce_csv, file = paste0(getwd(), "/datasets/whole_sce.csv"), row.names=TRUE, col.names=TRUE)
  
  scimpute(count_path = paste0(getwd(), "/datasets/whole_sce.csv"), infile = "csv", outfile = "csv", out_dir = paste0(getwd(), "/datasets/scimpute_results/"), ncores = 30, Kcluster = 75)
}
  
  
  
# Numero de células por tipo celular
ggplot(as.data.frame(colData(sce)), aes(x = factor(cluster, levels = rev(names(sort(table(cluster))))), fill = cluster)) +
  geom_bar(stat = "count") +
  labs(title = "",
       x = "Tipo celular",
       y = "Frecuencia")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        legend.position = "none")

# Scaling and log-transforming ----
set.seed(100)
clust.norm <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust.norm, min.mean=0.1)
sce <- logNormCounts(sce)
assayNames(sce)

#Plot
if(FALSE){
  markers <- unique(strsplit2(colData(sce)[, "cluster"], split = " ")[, 3])
  tmp_sce <- assay(sce, "logcounts")
  rownames(tmp_sce) <- rowData(sce)[, "gene"]
  #me quedo con los markers
  tmp_sce <- tmp_sce[which(rowData(sce)[, "gene"] %in% markers), ]
  
  #
  tmp_sce <- data.table(t(tmp_sce))
  tmp_sce$cell_type <- colData(sce)[, "cluster"]
  
  tmp_sce           <- aggregate(. ~ cell_type, data=tmp_sce, FUN=median)
  rownames(tmp_sce) <- tmp_sce$cell_type; tmp_sce <- tmp_sce[,-1]
  
  #Heatmap----
  library(dendextend)
  library(circlize)
  col_fun = colorRamp2(seq(min(tmp_sce), max(tmp_sce), length = 3),  
                       c("ghostwhite", "thistle1", "thistle4"))
  row_dend = as.dendrogram(hclust(dist(tmp_sce)))
  
  #Annotations
  row_ha = rowAnnotation(N = anno_barplot(as.numeric(table(colData(sce)[, "cluster"]))), width = unit(1.5, "cm"))
  
  draw(Heatmap(tmp_sce,
          name = "Expresión",
          column_title = "Genes Marcadores", 
          row_title = "Tipos celulares", 
          column_title_side = "bottom",
          show_column_dend = FALSE,
          rect_gp = gpar(col = "gray90", lwd = 1),
          row_names_side = "right",
          column_names_side = "top",
          row_names_gp = gpar(fontsize = 8),
          row_dend_side = "right",
          cluster_rows = color_branches(row_dend, k = 5),
          #row_km = 5,
          right_annotation = row_ha,
          col = col_fun), heatmap_legend_side = "right")


}
# Feature selection ----
dec <- modelGeneVar(sce)

# Visualizing the fit:
fit.sce <- metadata(dec)
plot(fit.sce$mean, fit.sce$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)

#
rowSubset(sce) <- getTopHVGs(dec, prop=0.1) # stored in the default 'subset'.
rowSubset(sce, "HVGs.more") <- getTopHVGs(dec, prop=0.2)
rowSubset(sce, "HVGs.less") <- getTopHVGs(dec, prop=0.3)
rowSubset(sce, "HVGs.0.5") <- getTopHVGs(dec, prop=0.5)
rowSubset(sce, "HVGs.0.8") <- getTopHVGs(dec, prop=0.8)
colnames(rowData(sce))

#Make cluster aggregation

#level 1
tmp <- strsplit2(sce$cluster, " ")
sce$cluster_agg <- apply(tmp,1, FUN = function(x){paste(x[1], x[2], collapse=" ")})

#level 1
sce$cluster_agg2 <- strsplit2(sce$cluster_agg, "-")[,1]

#level 1
sce$cluster_agg3 <- strsplit2(sce$cluster_agg, " ")[,1]

#Dimensionality reduction----

#PCA
set.seed(100) # See below.
sce <- fixedPCA(sce, subset.row=TRUE)

#Variance explicada por componentes PCA
percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Varianza Explicado (%)")

plotReducedDim(sce, dimred="PCA", colour_by="cluster_agg3")

plotReducedDim(sce, dimred="PCA", ncomponents=3,
               colour_by="cluster_agg")

#TSNE
set.seed(100)
sce <- runTSNE(sce, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce, dimred="TSNE",
                       colour_by="cluster_agg3") + ggtitle("perplexity = 5")

set.seed(100)
sce <- runTSNE(sce, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce, dimred="TSNE",
                        colour_by="cluster_agg3") + ggtitle("perplexity = 20")

set.seed(100)
sce <- runTSNE(sce, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce, dimred="TSNE", 
                        colour_by="cluster_agg3") + ggtitle("perplexity = 80")

gridExtra::grid.arrange(out5, out20, out80, ncol=3)

#UMAP
set.seed(100)
sce <- runUMAP(sce, dimred="PCA")
plotReducedDim(sce, dimred="UMAP", colour_by="cluster_agg")


#Clustering----

#analizo la pureza de los clusters
pure.cluster <- neighborPurity(reducedDim(sce, "PCA"), sce$cluster)

pure.data <- as.data.frame(pure.cluster)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- sce$cluster

ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
  ggbeeswarm::geom_quasirandom(method="smiley")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        legend.position = "none")
#Within-cluster sum of squares
rmsd <- clusterRMSD(reducedDim(sce, "PCA"), sce$cluster)
barplot(rmsd, ylab="RMSD", xlab="Cluster", ylim = c(0,100))

#clustering jerarquico para tener el dendrograma
# Setting the seed due to the randomness of k-means.
set.seed(100)
khclust.info <- clusterCells(sce, use.dimred="PCA",
                             BLUSPARAM=TwoStepParam(
                               first=KmeansParam(centers=1000),
                               second=HclustParam(method="ward.D2", cut.dynamic=TRUE,
                                                  cut.param=list(deepSplit=3)) # for higher resolution.
                             ),
                             full=TRUE
)
table(khclust.info$clusters)

plotTSNE(sce, colour_by=I(khclust.info$clusters), 
         text_by=I(khclust.info$clusters))

#
k.stats <- khclust.info$objects$first
tree.sce <- khclust.info$objects$second$hclust

m <- match(as.integer(tree.sce$labels), k.stats$cluster)
final.clusters <- khclust.info$clusters[m]

# TODO: expose scater color palette for easier re-use,
# given that the default colors start getting recycled.
dend <- as.dendrogram(tree.sce, hang=0.1)
#labels_colors(dend) <- as.integer(final.clusters)[order.dendrogram(dend)]

plot(dend)

saveRDS(sce, file = "./datasets/sce.rds")

#Imputation methods----
sce <- readRDS("./datasets/sce.rds")

#magic
impMAGIC            <- magic(t(logcounts(sce)))
impMAGIC            <- t(as.matrix(impMAGIC$result));rownames(impMAGIC) <- NULL 
assay(sce, "magic") <- impMagic

#alra
k_choice <- choose_k(logcounts(sce))

library(ggplot2)
library(gridExtra)
df <- data.frame(x=1:100,y=k_choice$d)
g1<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k)   + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle('Singular values')
df <- data.frame(x=2:100,y=diff(k_choice$d))[3:99,]
g2<-ggplot(df,aes(x=x,y=y),) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k+1)   + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
grid.arrange(g1,g2,nrow=1)

alra_norm_completed <- alra(logcounts(sce),k=k_choice$k)
assay(sce, "alra") <- alra_norm_completed[[3]]


#SAVER
sce.saver <- saver(assay(sce, "counts"), ncores = 30)
rownames(sce.saver$estimate) <- NULL
assay(sce, "saver") <- sce.saver$estimate

#
saveRDS(sce, file = "./datasets/sce.rds")
