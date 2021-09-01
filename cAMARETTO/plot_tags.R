# PLOT ENRICHMENT TAGS from cAMARETTO enrichments
# Read, tag and normalize the cAmaretto's .csv files
library(dplyr)

#promote
# comms_imp from cAMARETOOhtml_margaret
#comm_imp <- c(8,13,22,40)
comm_imp <- c(1,7,8,9,13,21,22,40,48)

first_df<- lapply(seq_along(comm_imp), function(x){
  print(comm_imp[x])
  file <- read.csv(paste0('./promote/csv_comm_rewired/enrichment_amareto_com',comm_imp[x],'.csv'))
  # Create Community col
  mutate(file,Community= factor(comm_imp[x])) %>% 
    # Cut first word and lowercase of geneset name adn description
    mutate(Gene.Set.Name=stringr::word(gsub(tolower(Gene.Set.Name), pattern = "_",replacement =  " "), start=2,end=-1))
})

#SU2C
# comms_su2c <- c(4,6,8,9,21,44)
# 
# first_df<- lapply(seq_along(comms_su2c), function(x){
#   print(comms_su2c[x])
#   file <- read.csv(paste0('./su2c/csv_comm_rewired/enrichment_amareto_com',comms_su2c[x],'.csv'))
#   # Create Community col
#   mutate(file,Community= factor(comms_su2c[x])) %>% 
#     # Cut first word and lowercase of geneset name adn description
#     mutate(Gene.Set.Name=stringr::word(gsub(tolower(Gene.Set.Name), pattern = "_",replacement =  " "), start=2,end=-1))
# })
###
first_df <- do.call(rbind,first_df)
str(filtered_df)

###--####

# Define tags:
Hemato <- c("heme", "hemato","blood","erythro","bone marrow", "myeloid", "megakaryocyte")

Liver <- c("liver","hepato","bile acid", "glycogen","gluconeogenesis","xenobiotic metabolism", "drug metabolism", "liver tissue", "glucuronidation", "cytochrome p450", "glutathione")

Neuro <- c("Alzheimer","Parkinson","Huntington","neural","neuro", "brain","astrocyte", "glial","cerebral cortex", "prefrontal cortex", "temporal cortex", 
           "hypothalamus", "synapsis", "synaptic", "cerebellum")

Immune <- c("immune","neutrophil", "leukocyte", "macrophage", "lympho", "cytokine", "antigen response", " il2 ", " il6 ",
            "antigen processing","natural killer", "mhc class","dendritic cell", "NF-kB", "monocyte", 
            " t cell"," T-cell", " b cell"," B cell", " b-cells", "Complement cascade", "Interleukins", "class I MHC")

Cancer <- c("cancer","carcino", "onco", "neoplasm", "tumor", "adenoma", "melanoma","sarcoma", "leukemia", "myeloma",
            "prostate cancer", "breast cancer","lung cancer", "colorrectal cancer","lymphoma","mesothelioma",
            "metast", "epitelial-mesenchymal","ECM", "matrisome", "neuroblastoma", "glioblastoma", "glioma", "medulloblastoma")

Mitochon <- c("mitochondria", "energy", "oxidative", "phosphorilation","oxidation", "respiratory electron transport", "atp")

Protein <- c("protein synthesis","protein modification", "transcription", "ribosome", "housekeeping",
             "Peptide chain elongation", "metabolism of proteins","metabolism of amino acids", "mtorc1", "proteasome")

Cellcycle <- c("cell division","TP53", " mTOR " ,"Cellular Proliferation",
               "G1 to S Transition","PI3K/AKT/mTOR pathway" ," RB ", "E2F","apoptosis", "mapk", "rb1", "cell cycle", "cell-cycle" , "p53", "mitotic")

GPCR <- c(" GPCR ", " G protein ", " rhodopsin ", " G alpha ", " G-protein ")

Embryo <- c("embryo", "developmental biology", " stem ", "Wnt")

Muscle <- c("skeletal", "smooth", "myogen", "myocyte", "heart", "muscle", "amyotrophic", "myofibers")

Micro <- c("myco", "bacteria", "virus", "infection", "interferon signaling", "virulence", "salmonella")

Adipo <- c("adipo", "fatty acid", "triacylglycerol","ketone", "lipol","lipid metabolism", "lipids")

Sexorg <- c("mammary", "prostate", "breast", "estrogen", "androgen", "oocyte","testis",
            "ovulation", "fertility", "progesterone", "myometri", "spermatogonial", "spermatogenesis " )
Skin <- c(" skin ", "melanoma","epidermis", "epidermal", "keratinocyte")
###--####
# Function create tags
create_tags <-function(tag_list, text){
  
  # For each tag family
  # For each term
  #grepl in Gene.Set.Name + Gene.Set.Description (text)
  #check if any of the tags of the family is present to assign the tag to the genset
  
  sapply(tag_list, function(x){
    
    tag_true <- sapply(x, function(y){
      
      if ("GPCR"%in%x){
        grepl(pattern=y,text, ignore.case=F)
      }else{
        grepl(pattern=y,text, ignore.case=T)
      }
      
    })
    
    idx_tag <- apply(tag_true,1,any)
    
  })
}

tag_list <- list(Hemato=Hemato, Liver=Liver,
                 Neuro=Neuro, Immune=Immune,
                 Cancer=Cancer, Mitochon=Mitochon,
                 Cellcycle=Cellcycle, Protein=Protein, 
                 GPCR=GPCR, Embryo=Embryo, Muscle=Muscle, Micro=Micro, Adipo=Adipo, Sexorg=Sexorg, Skin=Skin)

text <- paste(first_df$Gene.Set.Name, first_df$Gene.Set.Description)

tags <- create_tags(tag_list, text)

rownames(tags) <- first_df$Community
head(tags)
##
sum(apply(tags,1,function(x) sum(x))==0) 
length(unique(text[which(apply(tags,1,function(x) sum(x))==0)]))

# % without tags 0.30

length(unique(text[which(apply(tags,1,function(x) sum(x))==0)]))/length(unique(text))
##
# Take away non-tagged
idx <- which(rowSums(tags)!=0) 
filt_tags <- tags[idx,]

tots_tags <- colSums(filt_tags) #total to normalize

rownames(filt_tags) <- first_df$Community[idx]

filtered_df <- first_df[idx,]

# Filter significant genesets
idx2 <- which(filtered_df$X..Genes.in.Gene.Set > 100 &
                filtered_df$X..Genes.in.Overlap > 20 &
                filtered_df$P.value<0.001 &
                filtered_df$FDR.Q.value<0.001)

refiltered_df <- filtered_df[idx2,c(1,2,6,7,9)]
refilt_tags <- filt_tags[idx2,]

# Generate table for barplot
build_multiple_barplot_df <- function(bin_df){
  
  comm <- rownames(bin_df)
  tissues <- colnames(bin_df)
  
  df <- matrix(0,1,2)
  
  for (row in seq(nrow(bin_df))){
    
    pos <- which(bin_df[row,])
    
    for (x in as.numeric(pos)){
      df <- rbind(df,c(comm[row],tissues[x]))
    }
    
  }
  
  df <- df[-1,]
  rownames(df) <- df[,1]
  colnames(df) <- c('Comm','Tissue')
  
  return(as.data.frame(df))
  
}

mul_filt_tags <- build_multiple_barplot_df(refilt_tags)

table(mul_filt_tags$Tissue)


df_table <- table(mul_filt_tags$Tissue, mul_filt_tags$Comm)
# Reorder Comms
s_df_table <- df_table[,order(as.numeric(colnames(df_table)))]

# Norm per number of tags
#norm_intra <- t(s_df_table/tots_tags[sort(names(tots_tags))])

# Take Cancer out
not_cancer_df_table <- t(s_df_table[which(rownames(s_df_table)!="Cancer"),])
# Norm per total community

comm_count <- apply(not_cancer_df_table,1,sum)


norm_intra_within <- not_cancer_df_table/comm_count

colorset <- hues::iwanthue(14)
idx3 <- rownames(norm_intra_within)%in%mycoms


#plot export pdf
pdf("./plot_tags_coms_promote.pdf", height=6, width=15)
#pdf("./plot_tags_coms_SU2C.pdf", height=6, width=15)

barplot(t(norm_intra_within), 
        col=colorset,
        beside=TRUE, 
        ylab="Counts", xlab="Regulatory Module", ylim = c(0,1))
legend("topright", legend=as.character(levels(as.factor(colnames(norm_intra_within)))),
       fill=colorset,cex=1.5,y.intersp = 0.8,x.intersp = 0.4,ncol = 2)

dev.off()
