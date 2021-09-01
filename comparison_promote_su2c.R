# SCRIPT COMPARISON COMMUNITIES PROMOTE-SU2C
#Here we are going to check if the communities from promote 21k 
#is enriched in the communities from the common genes of su2c-promote.

library("CommunityAMARETTO")
library("foreach")
library(igraph)
library(dplyr)
source('./scripts/cAMARETTO_modified_functions_margaret.R')

# Function generate communities
generate_comms <- function(cAMARETTOresults,jaccard=0.2,m=80*10){
  
  #Module Network
  cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork_TraRe(prom_su2c_results,
                                                   pvalue = 0.05/m/m,
                                                   inter = 5,
                                                   jaccard=j,
                                                   edge_method="overlap",
                                                   plot_network = FALSE)
  
  #Identify communities
  cAMARETTOnetworkC<-cAMARETTO_IdentifyCom2(cAMARETTOnetworkM,
                                            filterComm = FALSE,
                                            ratioCommSize = 0.01,
                                            MinRuns = 2,
                                            ratioRunSize = 0.1,
                                            ratioEdgesInOut = 0.5,
                                            plot_network=F)
  
  return(cAMARETTOnetworkC)
  
}

# load data PROMOTE
prom_21k_results_comms <- readRDS('./promote/network_promote_21k_50b.rds')
linkeroutput_prom_21k <- readRDS('./promote/sm_3500_linker_output_50bootstraps.rds')
#load data SU2C
prom_su2c_results <- readRDS('./su2c/13k_cAMARETTOresults.rds')
prom_su2c_results_comms <- generate_comms(prom_su2c_results)
linkeroutput_su2c <- readRDS(paste0(getwd(),'/su2c/su2c_50bootstraps.rds'))


# Function adj matrix -----------------------------------------------------

generate_adjmatrix <- function(df_1,df_2, dfl_1,dfl_2, lo_1,lo_2){
  
  get_targs_regs_comm <- function(imp_boot_mod,lo){
    
    regs_targs <- lapply(imp_boot_mod,function(x){
      
      x_num <- as.numeric(x)
      regprog <- lo$raw_results$VBSR$bootstrapResults[[x_num[1]]]$RegulatoryPrograms[x_num[2],]
      regs <- names(regprog[regprog!=0])
      
      modmember <- lo$raw_results$VBSR$bootstrapResults[[x_num[1]]]$ModuleMembership
      targs <- rownames(modmember)[modmember%in%x_num[2]]
      
      return(list(regs=regs,targs=targs))
      
    })
    
    tot_regs <- unique(unlist(sapply(regs_targs,function(x) x$regs)))
    tot_targs <- unique(unlist(sapply(regs_targs,function(x) x$targs)))
    
    return(list(regs=tot_regs,targs=tot_targs))
  }
  
  calculate_hiptest_jacc <- function(group1,group2,universe_size=13e3){
    
    group1 <- unique(group1)
    group2 <- unique(group2)
    
    
    # building contigency table.
    contig_tbl <- as.table(matrix(c(length(intersect(group1, group2)), 
                                    length(setdiff(group1, group2)), 
                                    length(setdiff(group2, group1)), 
                                    universe_size - length(group2) - length(group1) + length(intersect(group1, group2))), 
                                  ncol = 2, 
                                  byrow = TRUE))
    
    res <- stats::fisher.test(contig_tbl, alternative = "g")
    pvals <- res$p.value
    
    jacc <- length(intersect(group1,group2))/length(union(group1,group2))
    
    #return(jacc)
    return(pvals)
    #return(list(pvals=pvals,jacc=jacc))
    
  }
  
  adjmatrix <- c()
  message('Generating reference')
  #preprocess dfl_2
  imp_comms <- df_2$community_list[dfl_2]
  
  for (x in imp_comms){
    
    
    imp_mods <- x[sapply(x,function(x) substr(x,1,1)=='S')]
    imp_boot_mod <- lapply(strsplit(imp_mods,'\\|'),function(x) stringr::str_extract(x,'[0-9]+'))
    ref_targs_regs <- get_targs_regs_comm(imp_boot_mod,lo_2)
    
    #append regs and targs
    regs <- unique(ref_targs_regs$regs)
    targs <- unique(ref_targs_regs$targs)
    genes_ref <- c(regs,targs)
    
    for (comm in dfl_1){
      
      comm_i_info <- df_1$community_list[[comm]]
      imp_boot_mod <- lapply(strsplit(comm_i_info,'\\|'),function(x) stringr::str_extract(x,'[0-9]+'))
      ref_targs_regs <- get_targs_regs_comm(imp_boot_mod,lo_1)
      
      regs <- unique(ref_targs_regs$regs)
      targs <- unique(ref_targs_regs$targs)
      genes <- c(regs,targs)
      
      adjmatrix <- c(adjmatrix,calculate_hiptest_jacc(genes_ref,genes))
      
    }
    
  }
  
  adjmat <-matrix(adjmatrix,length(dfl_1),length(dfl_2))
  rownames(adjmat) <- paste0('P_',dfl_1)
  colnames(adjmat) <- paste0('S_',dfl_2)
  return (adjmat)
  
  
  
}

## make comparison
adjmatt_metrics <- generate_adjmatrix(df_1=prom_21k_results_comms, df_2=prom_su2c_results_comms,
                                      dfl_1= c(8,13,22,40),dfl_2= c(4,6,8,9,21,44),
                                      lo_1= linkeroutput_prom_21k, lo_2=linkeroutput_su2c)

#plot heatmap
library(pheatmap)
library(viridis)
addidx <- adjmatt_metrics
addidx[addidx>0.1] <- 0.1

viridis::cividis(n = 10)

pheatmap::pheatmap(addidx,color =viridis::cividis(n = 20),cluster_rows = c(F,F,F,F), cluster_cols = c(F,F,F,F,F,F))
