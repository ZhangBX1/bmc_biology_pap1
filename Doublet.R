wt_object=NormalizeData(wt_object)
wt_object <- ScaleData(wt_object , verbose = FALSE)
wt_object =FindVariableFeatures(wt_object , selection.method = "vst", nfeatures = 2000)
wt_object  <- RunPCA(wt_object, npcs = 40)
wt_object  <- FindNeighbors(wt_object , reduction = "pca", dims = 1:40)
wt_object <- RunUMAP(wt_object, reduction = "pca", dims = 1:40)
ElbowPlot(wt_object , ndims = 50)

sweep.res <- paramSweep(wt_object, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
optimal_pk <- as.numeric(as.character(optimal_pk))
nExp_poi <- round(0.1 * ncol(wt_object))

wt_object <- doubletFinder(wt_object, 
                        PCs = 1:40, 
                        pN = 0.25, 
                        pK = optimal_pk, 
                        nExp = nExp_poi, 
                        reuse.pANN = FALSE, 
                        sct = FALSE)

colnames(wt_object@meta.data)

df_columns <- grep("DF.classifications", colnames(wt_object@meta.data), value = TRUE)
print(df_columns)

DimPlot(wt_object, 
        reduction = "umap",
        group.by = df_columns,  
        cols = c("Singlet" = "grey", "Doublet" = "red")) +
        ggtitle("Doublets")

table(wt_object@meta.data[, df_columns])

VlnPlot(wt_object, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        group.by = df_columns,
        pt.size = 0)

FeatureScatter(wt_object, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = df_columns)

wt_object_filtered <- subset(wt_object, 
                           cells = rownames(wt_object@meta.data)[wt_object@meta.data[, df_columns] == "Singlet"])

mutant_object=NormalizeData(mutant_object)
mutant_object <- ScaleData(mutant_object , verbose = FALSE)
mutant_object =FindVariableFeatures(mutant_object , selection.method = "vst", nfeatures = 2000)
mutant_object  <- RunPCA(mutant_object, npcs = 40)
mutant_object  <- FindNeighbors(mutant_object , reduction = "pca", dims = 1:40)
mutant_object <- RunUMAP(mutant_object, reduction = "pca", dims = 1:40)
ElbowPlot(mutant_object , ndims = 50)

sweep.res <- paramSweep(mutant_object, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
optimal_pk <- as.numeric(as.character(optimal_pk))

nExp_poi <- round(0.1 * ncol(mutant_object))

mutant_object <- doubletFinder(mutant_object, 
                        PCs = 1:40, 
                        pN = 0.25, 
                        pK = optimal_pk, 
                        nExp = nExp_poi, 
                        reuse.pANN = FALSE, 
                        sct = FALSE)


colnames(mutant_object@meta.data)
df_columns <- grep("DF.classifications", colnames(mutant_object@meta.data), value = TRUE)
print(df_columns)

DimPlot(mutant_object, 
        reduction = "umap",
        group.by = df_columns,  
        cols = c("Singlet" = "grey", "Doublet" = "red")) +
        ggtitle("Doublets")

table(mutant_object@meta.data[, df_columns])

VlnPlot(mutant_object, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        group.by = df_columns,
        pt.size = 0)

FeatureScatter(mutant_object, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               group.by = df_columns)

mutant_object_filtered <- subset(mutant_object, 
                           cells = rownames(mutant_object@meta.data)[mutant_object@meta.data[, df_columns] == "Singlet"])


object.list <- list(wt_object_filtered , mutant_object_filtered )
