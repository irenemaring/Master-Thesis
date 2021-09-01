              ### Comparison SGMM and TraRe plots script ###

# Load packages -----------------------------------------------------------
library(ggplot2)
library(cowplot)
library(VennDiagram)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)
library(hues)

# Load data ---------------------------------------------------------------
lo1 <- readRDS("./promote/sm_3500_linker_output.rds") #TraRe GRN output
lo2 <- readRDS("./SparseGMM/promote_10b_SGMM_lambda325_2nd.rds") #SGMM GRN output
# Number of regulators and target genes
n_regs1 <- sapply(lo1$modules$VBSR,function(x) length(x$regulators))
n_regs2 <- sapply(lo2$modules$VBSR,function(x) length(x$regulators))
n_targs1 <- sapply(lo1$modules$VBSR,function(x) length(x$target_genes))
n_targs2 <- sapply(lo2$modules$VBSR,function(x) length(x$target_genes))
n_mods1 <- sapply(lo1$raw_results$VBSR$bootstrapResults, function(x)x$NrModules)
n_mods2 <- sapply(lo2$raw_results$SGMM$bootstrapResults, function(x)x$NrModules)

# Histograms ---------------------------------------------------------------
par(mfrow(c(2,1)))
myhistogram <- function(data, main="", sub="", mybreak=NULL, num=FALSE, xlab = ""){
  if (!is.null(mybreak)){
    h <- hist(data,main=main, sub=sub,breaks=mybreak,xlab = xlab)
  }else{h <- hist(data,main=main, sub=sub, xlab = xlab)}
  if (num) {
    text(h$mids,h$counts,labels=h$counts,adj=c(0.5, -0.5))}
  return(h)
}

h1 <- myhistogram(n_regs1,main='Histogram of number of regulators in GRNs',sub= "TraRe", num=T)
h2 <- myhistogram(n_regs2,main='Histogram of number of regulators in GRNs',sub= "SGMM", num=T,mybreak = median(regs2))
h3 <- myhistogram(n_targs1,main='Histogram of number of targets in GRNs',sub= "TraRe", num=F,mybreak = 30)
h4 <- myhistogram(n_targs2,main='Histogram of number of targets in GRNs',sub= "SGMM", num=F,mybreak = 30)

par(mfrow=c(1,2))
h1 <- myhistogram(n_regs1,sub= "TraRe", num=T)
h2 <- myhistogram(n_regs2,sub= "SGMM", num=T,mybreak = max(n_regs2))
title("Histogram of number of regulators in GRNs", line = -1, outer = TRUE)

h3 <- myhistogram(n_targs1,sub= "TraRe", num=F,mybreak = 30)
h4 <- myhistogram(n_targs2,sub= "SGMM", num=F,mybreak = 30)
title("Histogram of number of targets in GRNs", line = -1, outer = TRUE)
par(mfrow=c(1,1))

myhistogram(n_mods1)
myhistogram(n_mods2)
barplot(n_mods1,ylim = c(0,100))
barplot(n_mods2,ylim = c(0,100))

t1 <- summary(n_mods1)
t1 <- c(t1,sd=round(sd(n_mods1),2))
t2 <- summary(n_mods2)
t2 <- c(t2,sd=round(sd(n_mods2),2))
tsum <- rbind(t1,t2)
rownames(tsum) <- c("TraRe", "SGMM")
pdf("./summary_table_nmods.pdf", height=2, width=6)
plot.new()
title(main="Number of modules",sub="Descriptive statistics in 10 bootstraps",outer=F,adj =0.5)
gridExtra::grid.table(tsum)
dev.off()

plot.new()+title(main="Number of modules",sub="Descriptive statistics in 10 bootstraps",outer=F,adj =0.5)+gridExtra::grid.table(tsum)

# Violin plots --------------------------------------------------------------

C1 <- data.frame(value=n_regs1, variable = "TraRe")
C2 <- data.frame(value=n_regs2, variable = "SGMM")
data1 <- rbind(C1, C2)
C3 <- data.frame(value=n_targs1, variable = "TraRe")
C4 <- data.frame(value=n_targs2, variable = "SGMM")
data2 <- rbind(C3, C4)
# ggplot(data, aes(x = variable, y = value, fill=variable))+
#   geom_boxplot(outlier.shape = NA,)+
#   geom_jitter( size=0.4, alpha=0.5,) +
#   ggtitle("A boxplot with jitter") +
#   theme_bw()+
#   xlab("")

myviolplot <-function(d){
  ggplot(d, aes(x = variable, y = value,fill=variable)) +
  geom_violin() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=12, hjust = 0.5))+
    xlab("")
    
}

vp1 <- myviolplot(data1)+ylab("Counts")+
ggtitle("Regulators")

vp2 <- myviolplot(data2)+ylab("Counts")+
ggtitle("Targets")


pdf(file = "./Comparison_genes_sGMM_TraRe_violinplot.pdf",
    width = 5, height = 3,pointsize = 1, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "rgb")
plot_grid(vp1,vp2,
          labels= c("1", "2"),
          scale = 0.95,
          hjust=0.5, vjust=0.5,
          label_x=0.05, label_y=0.1)+
  ggtitle("Violin chart of gene counts distribution in GRNs")+
  theme(plot.title = element_text(size=15, hjust = 0.5,face = "bold"))
dev.off()
vv <- plot_grid(vp1,vp2,
          labels= c("1", "2"),
          scale = 0.95,
          hjust=0.5, vjust=0.5,
          label_x=0.05, label_y=0.1)+
  ggtitle("Violin chart of gene counts distribution in GRNs")+
  theme(plot.title = element_text(size=15, hjust = 0.5,face = "bold"))

ttplot <- plot.new()+title(main="Number of modules",sub="Descriptive statistics in 10 bootstraps",outer=F,adj =0.5)+gridExtra::grid.table(tsum)
plot_grid(ttplot,
          vv,nrow = 2,ncol = 1)
# VennDiagrams ------------------------------------------------------------

# Targets sigmods 10b promote sm3500
sigmodl10b <- read.delim('./promote/sigmodules_VBSR_1_10b.txt', header = FALSE)
sigmodl10b <- as.numeric(sapply(sigmodl10b,gsub,pattern = "VBSR.1.mod.", replacement = ""))
cluster_prom <- c(547,183,470,517,623,690,70,156,261,791,396)

trare_reg_module <- unique(unlist(sapply(cluster_prom,function(x) c(lo1$modules$VBSR[[x]]$target_genes,lo1$modules$VBSR[[x]]$regulators))))

regs1 <- unique(unlist(sapply(cluster_prom,function(x) lo1$modules$VBSR[[x]]$regulators)))

# Targets sigmods SGMM
path <- './rewirings/rewiring_promote_10b_lambda325_SGMM_2nd_0.05/'
sigmodl_sgmm <- read.delim(paste0(path,'sigmodules_VBSR_1.txt'), header = FALSE)
sigmodl_sgmm <- as.numeric(sapply(sigmodl_sgmm,gsub,pattern = "VBSR.1.mod.", replacement = ""))
cluster1 <- c(65,404,718,278,502,630,140,304,560)
cluster2 <- c(468,517,364,15,177,525,634,696,141,398)
#cluster3 <- c(63,105,240,704, 264, 43,406, 513,42,673,347,614,118,274,520,297,528)
cluster3 <- c(264,43,406, 513)

sgmm_reg_module1 <- unique(unlist(sapply(cluster1,function(x) c(lo2$modules$VBSR[[x]]$target_genes,lo2$modules$VBSR[[x]]$regulators))))
sgmm_reg_module2 <- unique(unlist(sapply(cluster2,function(x) c(lo2$modules$VBSR[[x]]$target_genes,lo2$modules$VBSR[[x]]$regulators))))
sgmm_reg_module3 <- unique(unlist(sapply(cluster3,function(x) c(lo2$modules$VBSR[[x]]$target_genes,lo2$modules$VBSR[[x]]$regulators))))

regs_gmm1 <- unique(unlist(sapply(cluster1,function(x) lo2$modules$VBSR[[x]]$regulators)))
regs_gmm2 <- unique(unlist(sapply(cluster2,function(x) lo2$modules$VBSR[[x]]$regulators)))


intersect(regs1,regs_gmm2)
setdiff(regs1,regs_gmm2)
setdiff(regs_gmm2,regs1)

intersect(regs_gmm2,intersect(regs_prom,regs_c22))
# Plot

myvenD2<-function(sgmm,TR,name="",myCol2){ VennDiagram::venn.diagram(
  #venn
  x=list(SGMM=sgmm,TraRe=TR),filename = paste0(name,".png"),  
  height = 480 , 
  width = 480 , 
  resolution = 300, 
  # Circles
  lwd = 1,
  lty = 'blank',
  fill=myCol2,
  fontface=c("plain","bold","plain"),
  fontfamily = "sans",
  cex=0.6,
  rotation.degree=180,
  # out
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.pos = c(35, -35),
  cat.dist = c(0.055, 0.055),
  output=TRUE)}

myvenD3<-function(cl1,cl2,cl3,name="",myCol3){
  VennDiagram::venn.diagram(
  #venn
  x=list(cluster1=cl1,cluster2=cl2,cluster3=cl3),
  filename = paste0(name,".png"),  
  height = 480 , width = 480 , resolution = 300, 
  # Circles
  lwd = 1,
  lty = 'blank',
  fill=myCol3,
  # Numbers
  fontface=c("plain","bold","plain","bold","bold","bold", "plain"),
  fontfamily = "sans",
  cex=0.6,
  #rotation.degree=180,
  # out
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.pos = c(-30, 30, 180),
  cat.dist = c(0.09, 0.09, 0.085),
  #cat.fontfamily = "sans",
  #label.col = c("white", "darkblue", "white", "coral1", 
                #"black", "lightcyan4", "white"),
  output=TRUE)
  }

#Plot
myCol2 <- c("#F8766D","#619CFF")
scales::show_col(myCol2)
myCol3 <- brewer.pal(3,"Reds")
scales::show_col(myCol3)

myvenD2(sgmm_reg_module2,trare_reg_module,"vennD_TraReSGMM_coincident_cluster", myCol2)
myvenD2(sgmm_reg_module1,trare_reg_module,"vennD_TraReSGMM_cluster1", myCol2)
myvenD2(sgmm_reg_module3,trare_reg_module,"vennD_TraReSGMM_cluster3", myCol2)

myvenD2(Com22,trare_reg_module,"pruega_C22_prom", myCol2)


myvenD3(sgmm_reg_module1,sgmm_reg_module2,sgmm_reg_module3, "SGMM_clusters", myCol3)

# statistics --------------------------------------------------------------


# Function create contingency table
create_conting_tbl <- function(universe_size, gr1,gr2){
contig_tbl <- as.table(matrix(c(length(intersect(gr1, gr2)),
                                length(setdiff(gr1, gr2)),
                                length(setdiff(gr2,gr1)),
                                universe_size - length(gr1) - length(gr2) + length(intersect(gr1, gr2))),
                              ncol = 2, 
                              byrow = TRUE))
}
universe_size<- length(rownames(preparerew$datasets[[1]]$norm_expr_mat_keep))
# For TraRe and SGMM
contig_tbl <- create_conting_tbl(gr1=sgmm_reg_module2,gr2=trare_reg_module, universe_size = universe_size)
res <- stats::fisher.test(contig_tbl, alternative = "g")

mygrid <- as.data.frame(t(combn(x=c("sgmm_reg_module1","sgmm_reg_module2", "sgmm_reg_module3"),m = 2)))
mygrid[4,] <- c("sgmm_reg_module1", "trare_reg_module")
mygrid[5,] <- c("sgmm_reg_module2", "trare_reg_module")
mygrid[6,] <- c("sgmm_reg_module3", "trare_reg_module")
# For SGMM clusters
for (i in seq(dim(mygrid)[1])){
  data <- split(as.character(mygrid[i,c(1,2)]),f=factor(mygrid[i,c(1,2)]))
  contig_tbl <- create_conting_tbl(gr1=get(data[[1]]),gr2=get(data[[2]]), universe_size = universe_size)
  res<- stats::fisher.test(contig_tbl, alternative = "g")
  mygrid[i,3] <- as.character(round(res$p.value,3))
}
colnames(mygrid)[3] <- "pvalue"
mygrid[5,3] <- "<2.2e-16"
pdf("./stats_rewiring_clusters.pdf", height=3, width=6)
plot.new()
title(main="Regulatory modules overlap test",sub ="Hypergeometric test",outer=F,adj =0.5)
gridExtra::grid.table(mygrid, rows=NULL)
dev.off()


# GO ----------------------------------------------------------------------

library("org.Hs.eg.db")
library(clusterProfiler)
OrgDb <- org.Hs.eg.db 
gene_ensembl <- na.omit(AnnotationDbi::mapIds(OrgDb,
                                              keys= sgmm_reg_module2,
                                              column="ENSEMBL",
                                              keytype="SYMBOL",
                                              multiVals="first"))
biol <- enrichGO(gene = gene_ensembl,
                 keyType = "ENSEMBL",
                 OrgDb = OrgDb,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 readable = TRUE)
pdf("./goBP_sgmm_cluster2.pdf", height=12, width=6)
barplot(biol, showCategory=30,title = "GO enrichment analysis")
dev.off()
goplot(biol)
