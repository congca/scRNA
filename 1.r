################################################
#      GEO   2021                              #
#      Cong Cao                                #
################################################

memory.limit()
memory.limit(102400)
### download table from GEO and generate
#an expression set, get annotations, create data matrix and target file

## set working directory and load libraries for processing
dir()
library(parallel)
cl.cores <- detectCores()
cl.cores
#Creates a set of copies of R running in parallel 
#and communicating over sockets.
cl <- makeCluster(cl.cores)
#stopCluster(cl)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("GEOquery")
library(GEOquery)
###########################################################################
## Series GSE134691
library(Seurat)  
library(dplyr)
library(patchwork)
new_counts <- read.table(file="./GSM3946323_SIA_0_Seurat_count_matrix.txt",row.names = 1)
head(new_counts)
dim(new_counts)#[1] 11414  7542
mydata <- CreateSeuratObject(counts = new_counts, min.cells = 3, project = "mydata_scRNAseq")
mydata#11413 features across 7542 samples within 1 assay 

pbmc.data <- new_counts
at_least_one <- apply(pbmc.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
hist(colSums(pbmc.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
##3.2 ??????pbmc???????????????Seurat??????
pbmc <- mydata
head(pbmc$RNA@data[,1:5])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
head(pbmc@meta.data, 5)
#????????????????????????, UMI??????, ?????????????????????
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# ??????
# 
# ?????????????????????????????????2500?????????200?????????
# ????????????????????????????????????????????????>5%

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hist(colSums(pbmc$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

##4.2 ???????????????

#????????????????????????????????????LogNormalize, ??????????????????????????????????????????10000,??????log?????????;???????????????pbmc[["RNA"]]@data???
#????????????,???????????????????????????
 
hist(colSums(pbmc$RNA@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#After normalization, the total expression of each cell
hist(colSums(pbmc$RNA@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")  

##4.3 Change gene identification

#To identify genes with highly variable expression between cells, 
#follow-up studies need to focus on these genes.
#Seurat's built-in FindVariableFeatures() function first calculates the mean and variance
#of each gene, and directly simulates the relationship. 2000 genes are returned by default.

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 10????????????????????????????????????
top10 <- head(VariableFeatures(pbmc), 10) #head(pbmc$RNA@var.features,10)
# "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"   "GNG11"  "S100A8"

# ???????????????????????????,?????????????????????
plot1 <- VariableFeaturePlot(pbmc)
# ???????????????????????????,?????????10?????????
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

##4.4 ????????????


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")



pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmcTsne <- RunTSNE(pbmc,dims.use=1:10,do.fast=TRUE)
TSNEPlot(pbmcTsne)

#???????????????????????????????????????????????????

 

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
 

# PC_ 1 
# Positive:  Cd74, H2-Aa, H2-Ab1, H2-Eb1, H2-DMa 
# Negative:  Fxyd2, Vsig4, Aqp1, Ctsk, Lyve1 
# PC_ 2 
# Positive:  Lyz2, Apoe, Wfdc17, Vsig4, Ctsd 
# Negative:  Cd74, H2-Aa, H2-Eb1, H2-Ab1, Stmn1 
# PC_ 3 
# Positive:  Stmn1, Birc5, Ube2c, Tubb5, Cdk1 
# Negative:  Pf4, Ccl8, Retnla, Cd209f, Cd209d 
# PC_ 4 
# Positive:  Lyz2, Ctsd, Vsig4, Sparc, H2-Aa 
# Negative:  Pf4, Malat1, Fxyd2, F13a1, Aqp1 
# PC_ 5 
# Positive:  Mmp9, Acp5, Ctsk, Atp6v0d2, Mt3 
# Negative:  Pf4, Tppp3, Lyve1, Aqp1, F13a1 
#Visualization of gene sets that have a greater impact on each principal component
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")



DimPlot(pbmc, reduction = "pca",split.by = 'ident')




DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)



pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)




pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


# Number of nodes: 5901
# Number of edges: 171107
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.7437
#    Number of communities: 7
#    Elapsed time: 0 seconds
 
#??????????????????????????????
head(Idents(pbmc), 5)
 


# V2 V5 V6 V7 V8 
# 1  2  0  1  3 
# Levels: 0 1 2 3 4 5 6

# install UMAP: reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


#??????????????????xiba?????????
DimPlot(pbmc, reduction = "umap",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "umap"),id = 'ident')





saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

#######################   KEGG ?????? ################################
# ???????????????????????? (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# 
# p_val avg_log2FC pct.1 pct.2
# Cd74   9.685467e-188  -1.941918 0.600 0.824
# H2-Ab1 7.660187e-143  -1.810298 0.381 0.655
# H2-Aa  6.808151e-138  -1.780533 0.364 0.642
# Vsig4  7.752818e-136   1.005536 0.837 0.538
# H2-Eb1 1.490916e-129  -1.832933 0.257 0.559
# p_val_adj
# Cd74   1.105402e-183
# H2-Ab1 8.742571e-139
# H2-Aa  7.770143e-134
# Vsig4  8.848291e-132
# H2-Eb1 1.701582e-125
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#            p_val        avg_log2FC pct.1 pct.2
# Mmp9      0.000000e+00   4.405342 0.747 0.009
# Acp5      0.000000e+00   5.019188 0.928 0.046
# Ctsk     1.955350e-283   5.006415 0.945 0.102
# Atp6v0d2 2.332366e-154   2.307371 0.321 0.006
# Atp6v1b2  5.341348e-36   1.518753 0.329 0.079

# p_val_adj
# Mmp9      0.000000e+00
# Acp5      0.000000e+00
# Ctsk     2.231641e-279
# Atp6v0d2 2.661929e-150
# Atp6v1b2  6.096080e-32


# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

saveRDS(pbmc, "./pbmc_5k_v3.rds")

 #As shown in the figure, after processing, the cells can be divided into 6 cell subgroups, clusters 0-5.
 
pbmc<- readRDS("./pbmc_5k_v3.rds")

# ??????GSEA,???????????????????????????
 library(devtools)
library(RPresto)
library(Rcpp)
# pbmc.genes <- wilcoxauc(pbmc, 'seurat_clusters')
# head(pbmc.genes)
# 
# 
# # ??????????????????cluster???????????????
# dplyr::count(pbmc.genes, group)
# 
# ?????????????????????????????????,??????????????????????????????????????????????????????Broad Institute?????????MsigDB[3] ,?????????:
# 
# ??????fgsea?????????????????????
# install.packages("msigdbr")
library(msigdbr)
# install.packages("fgsea")
library(fgsea)
library(ggplot2)

msigdbr_show_species()#???????????????????????????????????????

# [1] "Bos taurus"               "Caenorhabditis elegans"   "Canis lupus familiaris"   "Danio rerio"              "Drosophila melanogaster"  "Gallus gallus"           
# [7] "Homo sapiens"             "Mus musculus"             "Rattus norvegicus"        "Saccharomyces cerevisiae" "Sus scrofa"              
#  


m_df<- msigdbr(species = "Mus musculus", category = "C7") 
head(m_df)

# # A tibble: 6 x 17
# gs_cat gs_subcat gs_name entrez_gene gene_symbol human_entrez_ge~ human_gene_symb~ gs_id gs_pmid gs_geoid gs_exact_source gs_url gs_description species_name species_common_~
#   <chr>  <chr>     <chr>         <int> <chr>                  <int> <chr>            <chr> <chr>   <chr>    <chr>           <chr>  <chr>          <chr>        <chr>           
#   1 C7     ""        GOLDRA~       11305 Abca2                     20 ABCA2            M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# 2 C7     ""        GOLDRA~       27416 Abcc5                  10057 ABCC5            M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# 3 C7     ""        GOLDRA~       68644 Abhd14a                25864 ABHD14A          M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# 4 C7     ""        GOLDRA~       11364 Acadm                     34 ACADM            M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# 5 C7     ""        GOLDRA~       11433 Acp5                      54 ACP5             M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# 6 C7     ""        GOLDRA~       66659 Acp6                   51205 ACP6             M3044 164927~ ""       GSE1000002_158~ ""     Genes down-re~ Mus musculus house mouse     
# # ... with 2 more variables: ortholog_sources <chr>, num_ortholog_sources <dbl>

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets
#fgsea()????????????????????????????????????????????????,  ???????????????,???????????????,AUC?????????
devtools::install_github("immunogenomics/presto")
library('RPresto')
pbmc.genes <- wilcoxauc(pbmc, 'seurat_clusters')
# Naive CD4+ T cells
pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)      #??????????????????


# ???cluster0?????????????????????IL7R???CCR7,????????????CD4 + T?????????????????????
# ?????????fgsea???feature???auc???
cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)
ranks<- deframe(cluster0.genes)
head(ranks)


fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

#????????????
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()

# ??????top20????????????
ggplot(fgseaResTidy %>% filter(padj < 0.008) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####???7.5??????????????????


#GSEA style plot
plotEnrichment(fgsea_sets[["GSE10325_CD4_TCELL_VS_MYELOID_UP"]],
               ranks) + labs(title="GSE10325 CD4 TCELL VS MYELOID UP")
############################ Another way do KEGG #######################################
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(ggplot2)

for( j in 0:12)
{
  cluster.markers <- FindMarkers(object = pbmc, ident.1 =j, logfc.threshold = 0.25, test.use = "bimod", only.pos = TRUE)
  cluster<- row.names.data.frame(cluster.markers)
  cluster=bitr(cluster,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
  cluster.go<-enrichGO(gene=cluster[,"ENTREZID"],keyType = "ENTREZID",OrgDb=org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = TRUE)
  assign(paste0("cluster",j,".go"),cluster.go)
  pdf(file = paste0("cluster",j,"go.pdf"),,width=20,height=10)
  barplot(cluster.go,showCategory=50)
  dev.off()
  cluster.kegg<-enrichKEGG(gene = cluster[,"ENTREZID"],organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  assign(paste0("cluster",j,".kegg"),cluster.kegg)
  pdf(file = paste0("cluster",j,"kegg.pdf"),,width=20,height=10)
  dotplot(cluster.kegg,showCategory=50)
  dev.off()
  write.csv(x=cluster.markers,file=paste0("cluster",j,".csv"))
}







for( j in 0:12)
{
  cluster.markers <- FindMarkers(object = RA.integrated, ident.1 =j, logfc.threshold = 0.25, test.use = "bimod", only.pos = TRUE)
  cluster<- row.names.data.frame(cluster.markers)
  cluster=bitr(cluster,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
  cluster.go<-enrichGO(gene=cluster[,"ENTREZID"],keyType = "ENTREZID",OrgDb=org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = TRUE)
  assign(paste0("cluster",j,".go"),cluster.go)
  pdf(file = paste0("cluster",j,"go.pdf"),,width=20,height=10)
  barplot(cluster.go,showCategory=50)
  dev.off()
  cluster.kegg<-enrichKEGG(gene = cluster[,"ENTREZID"],organism = 'hsa', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  assign(paste0("cluster",j,".kegg"),cluster.kegg)
  pdf(file = paste0("cluster",j,"kegg.pdf"),,width=20,height=10)
  dotplot(cluster.kegg,showCategory=50)
  dev.off()
  write.csv(x=cluster.markers,file=paste0("cluster",j,".csv"))
}


VlnPlot(pbmc, features = c("Mmp9", "Ctsk"))
FeaturePlot(pbmc, features = c("Mmp9", "Acp5", "Ctsk", "Atp6v0d2", "Atp6v1b2"))
new.cluster.ids <- c("Mmp9", "Acp5", "Ctsk", "Atp6v0d2", "Atp6v1b2")
# 
# 

