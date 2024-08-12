#integration of two kidney data sets that were preprocessed by 1.preprocess.R

library(Seurat)
library(dplyr)

project_dir <- "/path/to/project/"
# dataset 1 
sub.seurat <- readRDS(file = paste0(project_dir, "subset.seurat.RDS"))
# dataset 2 
sc <- readRDS(file = paste0(project_dir, "sc.RDS"))


# remove certain cell types
#updated gene list in 2024 june
genelist.reorder2 <- c("JAK1", "N4BP2L2", "LUC7L3", "COL4A1","COL1A2", "COL6A3", "COL8A1", "FLRT2", "GXYLT2", "LONRF1", "CDH6", "SLCO3A1", "IGFBP3",  "MAFB", "NTNG1", "BMP2", "PAX8-AS1", "RYR3", "HLX", "THSD7A", "PTPN13", "SLC6A13", "FGF9", "TFCP2L1", "PKHD1", "GRB14", "CCSER1", "PROM1", "LY86", "RB1", "RNASET2", "FLI1", "CD53", "HLA-DRB1", "HLA-DQA1", "PRKX",  "SAMHD1", "CXCL8", "CHST11", "NMT2", "L3MBTL3", "SORCS1", "LINC00113", "ETV1", "RRM2", "ANKRD30B")



############### Integration of two datasets ###############
integrate.list <- list(sub.seurat, sc)
features <- SelectIntegrationFeatures(object.list = integrate.list)
integrate.list <- lapply(X = integrate.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = integrate.list, reference = c(1, 2), reduction = "rpca", dims = 1:30)
saveRDS(anchors, file=paste0(project_dir, "anchor.RDS"))

seu.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, dims = 1:30)




###################### re label integreated ######################
# merge the relevent biological terms (cell types)
library(dplyr)

metadata.integrated <- seu.integrated@meta.data
df <- cbind(metadata.integrated$Annotation.Level.2, metadata.integrated$Cluster_Idents)
result <- data.frame(combined = c(t(df)))
result <- result[!is.na(result$combined), , drop = FALSE]
seu.integrated$anno <- result$combined

DimPlot(seu.integrated, group.by = "anno")
tt <- result %>%
  mutate(across(where(is.character), 
                ~case_when(
                  . == "Arterioral Endothelium" ~ "Endothelium",
                  . == "B_Naive" ~ "B Cells",
                  . == "B_memory" ~ "B Cells",
                  . == "Baso/Mast" ~ "Basophils",
                  . == "CD14_Mono" ~ "Monocytes",
                  . == "CD16_Mono" ~ "Monocytes",
                  . == "C_TAL" ~ "TAL",
                  . == "M_TAL" ~ "TAL",
                  . == "Thick Ascending Limb" ~ "TAL",
                  . == "CD4T" ~ "T Cells",
                  . == "CD8T" ~ "T Cells",
                  . == "cDC" ~ "Dendritic Cells",
                  . == "pDC" ~ "Dendritic Cells",
                  . == "DCT1" ~ "DCT",
                  . == "DCT2" ~ "DCT",
                  . == "Des-Thin_Limb" ~ "LOH",
                  . == "Ascending_Thin_LOH" ~ "LOH",
                  . == "Mac" ~ "Macrophages",
                  . == "MyoFib/VSMC " ~ "Myofibroblast",
                  . == "Natural Killer Cells" ~ "NK",
                  . == "Collecting Duct Principal Cells" ~ "PC",
                  . == "Connecting Tubule" ~ "CNT",
                  . == "Podo" ~ "Podocytes",    
                  . == "Plasma_Cells" ~ "Plasma Cells",
                  . == "Fibroblast_1" ~ "Fibroblast",
                  . == "Fibroblast_2" ~ "Fibroblast",
                  . == "Distal Convoluted Tubule" ~ "DCT",
                  . == "IC_A" ~ "Intercalated Cells",
                  . == "IC_B" ~ "Intercalated Cells",
                  . == "Proximal Tubule" ~ "Proximal Tubule",
                  . == "PT_S1" ~ "Proximal Tubule",
                  . == "PT_S2" ~ "Proximal Tubule",
                  . == "PT_S3" ~ "Proximal Tubule",
                  . == "iPT" ~ "Proximal Tubule",
                  . == "Endo_GC" ~ "Endothelium",
                  . == "Endo_Peritubular" ~ "Endothelium",
                  . == "Venular Endothelium" ~ "Endothelium",
                  . == "Fibroblast" ~ "Myofibroblast",
                  . == "MyoFib/VSMC" ~ "Myofibroblast",
                  . == "Macula_Densa" ~ "Macula Densa",
                  TRUE ~ .
                )))

    #remove RBC, Neural_Cells, Prolif_Lym, Plasma_Cells

seu.integrated2 <- subset(seu.integrated, subset = anno2 %in% c("RBC", "Neural_Cells", "Prolif_Lym", "Plasma_Cells"), invert = T)


seu.integrated$anno2 <- tt$combined
DimPlot(seu.integrated, group.by = "anno2", label = T)
saveRDS(seu.integrated, file = paste0(project_dir, "seu.integrated.RDS"))

###################### re label BioRXiv dataset ######################
sc <- readRDS(file = paste0(project_dir, "scRNA.RDS"))

t1 <- sc@meta.data$Cluster_Idents
Idents(sc) <- "Cluster_Idents"
sc2 <- subset(sc, Cluster_Idents %in% c("Neural_Cells", "Plasma_Cells", "Prolif_Lym", "RBC", "Endo_Lymphatic", "Endo_Lymphatic", "Endo_Peritubular"), invert = TRUE)

# need to list the full options
cluster_names <- c(
  "B_memory" = "B_Cells",
  "B_Naive" = "B_Cells",
  "Baso/Mast" = "Baso/Mast",
  "CD14_Mono" = "Monocytes",
  "CD16_Mono" = "Monocytes",
  "CD4T" = "T_Cells",
  "CD8T" = "T_Cells",
  "CNT" = "CNT",
  "DCT1" = "DCT",
  "DCT2" = "DCT",
  "Des-Thin_Limb" = "LoH",
  "Ascending_Thin_LOH" = "LoH",
  "Endo_GC" = "Endothelium",
  "Fibroblast_1" = "Fibroblast",
  "Fibroblast_2" = "Fibroblast",
  "Mes" = "Messangial",
  "MyoFib/VSMC" = "myofibroblast",
  "Neutrophil" = "Neutrophil",
  "NK" = "NK",
  "Podo" = "Podocytes",
  "PT_S1" = "PT",
  "PT_S2" = "PT",
  "PT_S3" = "PT",
  "iPT" = "PT",
  "GS_Stromal" = "GS",
  "IC_A" = "IC",
  "IC_B" = "IC",
  "PC" = "PC",
  "Mac" = "Macrophage",
  "Macula_Densa" = "Macula_Densa",
  "cDC" = "DC",
  "pDC" = "DC",
  "C_TAL" = "TAL",
  "M_TAL" = "TAL",
  "PEC" = "PEC"
)
sc2$cluster_name <- cluster_names[as.character(Idents(sc2))]

# If you want to set this as the new identity
Idents(sc2) <- sc2$cluster_name
p2 <- DotPlot(sc2, features = genelist.reorder2, group.by = "cluster_name") + theme(axis.title =element_blank())+ coord_flip() + RotatedAxis()


sc2$group2 <- factor(sc2$group2, levels = c("Control", "Disease"))
Idents(sc2) <- "cluster_name"
# Get all unique cell types
cell_types <- unique(sc2$cluster_name)

# Create a list to store results
marker_lists <- list()

# Iterate through each cell type
for (cell_type in cell_types) {
  print(cell_type)
  # Find markers between health and disease for this cell type
  markers <- FindMarkers(sc2,
                         ident.1 = "Disease",
                         ident.2 = "Control",
                         group.by = "group2",
                         subset.ident = cell_type,
                         min.pct = 0.1,
                         logfc.threshold = 0.01)
  
  # Store the results
  marker_lists[[cell_type]] <- markers
}

saveRDS(marker_lists, file = paste0(path_deposit, "DE_markers.RDS"))

combined.de.markers <- NULL
for(i in names(marker_lists)){
  tt <- marker_lists[[i]] 
  tt$gene <- rownames(tt)
  tt$celltpe <- i
  tt1 <- tt[rownames(tt) %in% genelist.reorder2, ]
  combined.de.markers <- rbind(combined.de.markers, tt1)
}
write.table(combined.de.markers,file = paste0(project_dir,"dif_exp_test.txt" ) , sep = "\t", quote = F, col.names = NA)

# dotplot by two condiitons 
p <- DotPlot(sc2, 
        features = genelist.reorder2, 
        cols = c("blue", "red"),
        split.by = "group2") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("") + coord_flip() + RotatedAxis()
write.table(p$data, file = paste0(project_dir, "dotplot_data.txt"), sep = "\t", quote = F, col.names = NA)
pdf_and_png(p, paste0(project_dir, "dotplot_by_condition"), width = 15, h_to_w_ratio = 0.6, png = T)

#dotplot overall 

p2 <- DotPlot(sc2, 
             features = genelist.reorder2, 
          ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("") + coord_flip() + RotatedAxis()
write.table(p2$data, file = paste0(project_dir, "dotplot_data.txt"), sep = "\t", quote = F, col.names = NA)
pdf_and_png(p2, paste0(project_dir, "dotplot_markers"), width = 8, h_to_w_ratio = 1, png = T)

##################################################
library(plot1cell)

png(filename =  paste0(path_output, gene, ".png"), width = 4, height = 6,units = 'in', res = 300)
complex_dotplot_multiple(seu_obj = sc2, features = genelist.reorder2, group = "group2", celltypes = cell_types)
dev.off() 

png(filename =  paste0(path_output_genes, gene, "overall.png"), width = 3.5, height = 6,units = 'in', res = 300)
complex_dotplot_single(seu_obj = total_new_final, feature = gene) # default ident "cell_label_detail"
dev.off() 
