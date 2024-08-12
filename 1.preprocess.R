# Preparation of single-cell data sets from public data sets
# raw counts and original sample information and annotation were retained
# data were normalised using the standard Seurat analysis pipeline

library(Seurat)
library(dplyr)

project_dir <- "/path/to/project/"


########################################################
################# data set 1 ###########################
# Decoding myofibroblast origins in human kidney fibrosis. Nature. 2021
# https://zenodo.org/record/4059315 
# dataset name: Human_CD10negative.tar.gz

path_deposit <- paste0(project_dir, "dataDeposit/")
path_output <- paste0(project_dir, "results/")

dat <- Matrix::readMM("/inputPath/CD10negative/kidneyMap_UMI_counts.mtx"))
rowDat <- read.table("/inputPath/dataset/CD10negative//kidneyMap_UMI_counts_rowData.txt",sep=",",header=TRUE, stringsAsFactors = FALSE)
colDat <- read.table("/inputPath/dataset/CD10negative/kidneyMap_UMI_counts_colData.txt", sep=",",header=TRUE, stringsAsFactors = FALSE)

# Genes
dat1 <- dat[!duplicated(rowDat$Gene.Symbol),]
rowDat1 <- rowDat[!duplicated(rowDat$Gene.Symbol),]
rownames(dat1) <- rowDat1$Gene.Symbol
rownames(rowDat1) <- rowDat1$Gene.Symbol
rm(dat, rowDat)

# Cells
colnames(dat1) <- paste0("cell",1:ncol(dat1))
rownames(colDat) <- paste0("cell",1:ncol(dat1))

summary(as.data.frame(colDat))

sce <- SingleCellExperiment(assays=list("counts"=as.data.frame(dat1)),
                colData=colDat,
                rowData=rowDat1)
sce = scran::computeSumFactors(sce, 
                   sizes = seq(10, 200, 20), 
                   clusters = sce$Annotation.Level.2, 
                   positive = TRUE)
sce <- logNormCounts(sce)


logcounts(sce) <- as(logcounts(sce), "dgCMatrix")
#note that logcounts has to be transferred, otherwise has below error message 
#Error in checkSlotAssignment(object, name, value) : 
#  assignment of an object of class “dgTMatrix” is not valid for slot ‘data’ in an object of class “Assay”; is(value, "AnyMatrix") is not TRUE
logcounts(sce) <- as(logcounts(sce), "dgCMatrix")

cd10neg.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")


cd10neg.seurat[["percent.mt"]] <- PercentageFeatureSet(cd10neg.seurat, pattern = "^MT-")
VlnPlot(cd10neg.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cd10neg.seurat <- NormalizeData(cd10neg.seurat, normalization.method = "LogNormalize", scale.factor = 10000)

cd10neg.seurat <- FindVariableFeatures(cd10neg.seurat, selection.method = "vst", nfeatures = 5000)
cd10neg.seurat <- ScaleData(cd10neg.seurat)

cd10neg.seurat <- RunPCA(cd10neg.seurat, features = VariableFeatures(object = cd10neg.seurat))

ElbowPlot(cd10neg.seurat)
  cd10neg.seurat <- RunUMAP(cd10neg.seurat, dims = 1:38)



########################################################
################# data set 2 ###########################
# data downloaded from GSE211785 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785


obj <- readRDS(paste0(project_dir, "GSE211785_EXPORT_scRNAseq_snRNAseq_snATACseq_counts.rds"))
meta.data <- read.csv(paste0(project_dir,"GSE211785_scRNA-seq_snRNA-seq_snATAC-seq_metadata.txt"))
rownames(meta.data) <- meta.data$X
colnames(obj) <- meta.data$X
scRNA_meta <- meta.data %>% filter(tech == "SC_RNA") %>% select(X)

seu <- CreateSeuratObject(counts=obj, meta.data = meta.data)

#umap_info <- read.table(paste0(project_dir,"GSE211785_EXPORT_scRNA-seq_snRNA-seq_snATAC-seq_umap.txt"))
#umap_sc <- umap_info[grep(pattern = "SC_", x = rownames(umap_info)),]
########################################################

seu[["CellName"]] <- seu$X
sc <- subset(seu,  subset = CellName %in%  scRNA_meta$X)
rm(seu,obj)

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
sc<- PercentageFeatureSet(sc, pattern = "^Hb[abq]", col.name = "percent.hb")


sc <- CellCycleScoring(sc,
                g2m.features = cc.genes.updated.2019$g2m.genes,
                s.features = cc.genes.updated.2019$s.genes)

sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst")
sc <- ScaleData(sc)
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
sc <- RunUMAP(sc, dims = 1:20)

saveRDS(sc,file = paste0(project_dir, "scRNA.RDS"))

########################################################
# record on the original cell type names
tt <- data.frame(names(table(sc$Cluster_Idents)))
write.table(tt, file = paste0(project_dir, "celltype_name.txt"), quote = F, col.names = NA)
########################################################

p <- DimPlot(sc, group.by = "Cluster_Idents")



