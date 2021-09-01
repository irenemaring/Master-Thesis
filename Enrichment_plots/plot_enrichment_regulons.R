# SCRIPT FOR ENRICHMENT REGULONS

# Reexplore GSEA
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("viridis")
# input -------------------------------------------------------------------
# Regulons promote 50 boots
regulons <- readRDS(paste0(getwd(),'/ordered_filered_regulons_50.rds'))
# Annotation
OrgDb <- org.Hs.eg.db

# Gene nomenclature
genelist <- lapply(regulons, function(x) mapIds(OrgDb, column="ENTREZID",
                                                keytype="SYMBOL",
                                                multiVals="first",keys=x$targets))


# enrichment --------------------------------------------------------------

# Compare cluster
# BP
BP <- compareCluster(geneCluster = genelist,
                     fun = enrichGO,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     OrgDb=OrgDb,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     maxGSSize = 2000, 
                     readable=TRUE)

pdf("./doptplot_GOBP_regulons.pdf", height=15, width=20)

clusterProfiler::dotplot(clusterProfiler::simplify(BP),title="GO-BP ORA enrichment analysis of regulons")+scale_color_viridis()
dev.off()
