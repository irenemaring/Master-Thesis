## SCRIPT - SGMM ## 

#TRANSFORM OUTPUT FROM SGMM TO WORKING INPUT FOR TRARE_REWIRING

# for 10 botrstraps
object_rew <- list()

for (i in seq(1,10)){
  # Import results from SparseGMM
  weights <- read.csv(paste0("./SparseGMM/myoutput_lambda325_10b/weights_",i,"_promote_100.csv"),header = T,row.names = 1)
  lab <- read.csv(paste0("./SparseGMM/myoutput_lambda325_10b/labels_",i,"_promote_100.csv"),header = T,row.names = 1)
  lab <- lab+1 # because of python indexing
  
  idx <- sort(unique(lab$X0)) #Number of unique modules
  
  # For rewiring we need the modules
  boots_results <- lapply(seq(1,length(idx)), function(x){
    targs <- rownames(lab)[lab$X0==idx[x]]
    #message("targs", x)
    cl <- weights[,idx[x],drop=F]
    #cl_ord <- cl[order(abs(cl),decreasing = T),,drop=F]
    #cl_res <- scale(cl_ord,scale = T,center = T)
    #regs <- rownames(cl_res)[abs(cl_res)>1.96] #Regs needed for Regulatoryprograms
    cl_res <- cl[which(cl!=0),1, drop=F]
    regs <- rownames(cl_res)
    #message("regs", x)
    outc <- list(target_genes=targs, regulators=regs)
  })
  regs <- sapply(boots_results, function(x)x$regulators)
  
  if (i==1){
    # for bootrstrap 1
    #SPGMM
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$weights <- weights
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$labels <- lab
    # FOR AMARETTO:
    #NrModules
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$NrModules <- length(idx)
    
    # Regulatory Programs (matrix)
    # Initialise with 0 
    data <- matrix(nrow = length(idx),
                   ncol = length(rownames(weights)),
                   dimnames = list(paste0("Module_",seq_along(idx)),rownames(weights)),
                   data = 0)
    # Add weights at the corresponding modules in bootstrap 1
    for(j in seq(idx)){
      for (k in regs[[j]]){
        data[j,k] <- t(weights)[idx[j],k]
      }}
    # Insert in obect
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$RegulatoryPrograms <- data
    #AllRegulators
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$AllRegulators <- rownames(weights)
    
    #AllGenes
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$AllGenes <- rownames(lab)
    
    #ModuleMembership
    object_rew$raw_results$SGMM$bootstrapResults[[i]]$ModuleMembership <- as.matrix(lab)
    for (j in seq(idx)){
      object_rew$raw_results$SGMM$bootstrapResults[[i]]$ModuleMembership[which(lab==idx[j]),]=j
    }
    colnames(object_rew$raw_results$SGMM$bootstrapResults[[i]]$ModuleMembership) <- "ModuleNr"
    #TrainingStats (blank)
    
  }else{
    modmem <- as.matrix(lab)
    for (j in seq(idx)){
      modmem[which(lab==idx[j]),]=j
    }
    colnames(modmem) <- "ModuleNr"
    # Regulatory programs
    data <- matrix(nrow = length(idx),
                   ncol = length(rownames(weights)),
                   dimnames = list(paste0("Module_",seq_along(idx)),rownames(weights)),
                   data = 0)
    # Add weights at the corresponding modules in bootstrap 1
    for(j in seq(idx)){
      for (k in regs[[j]]){
        data[j,k] <- t(weights)[idx[j],k]
      }}
    actualiz <- list(list(weights=weights,
                          labels=lab,
                          NrModules=length(idx),
                          RegulatoryPrograms= data,
                          AllRegulators=rownames(weights),
                          AllGenes = rownames(lab),
                          ModuleMembership = modmem))
    
    object_rew$raw_results$SGMM$bootstrapResults <- c(object_rew$raw_results$SGMM$bootstrapResults,actualiz)
  }
  
  # Add to object the modules
  object_rew$modules$VBSR <- c(object_rew$modules$VBSR, boots_results)
}

saveRDS(object_rew, "./SparseGMM/promote_10b_SGMM_lambda325_2nd.rds")


