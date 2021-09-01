# SCRIPT ENRICHMENT REWIRED TFs
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("viridis")
#load file
drivers_fisher <- rownames(read.delim('./promote/sm_3500_graphs_fisherpvals_newV.txt'))

# Gene nomenclature
gene_entrez <- mapIds(OrgDb,
                      keys= drivers_fisher,
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")

# enrichGO: GO over-representation test (Hypergeometric)
biol <- enrichGO(gene = gene_entrez,
                 keyType = "ENTREZID",
                 OrgDb = OrgDb,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 readable = TRUE)

dotplot(pairwise_termsim(biol), x = "Count" ,showCategory=30,title = "GO enrichment analysis")+scale_color_viridis()

