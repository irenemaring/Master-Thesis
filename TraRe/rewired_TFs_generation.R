# SCRIPT REWIRING TFs

source('./scripts/format_graphs_to_modules_def.R')
#read the geneinfo
gene_info_p <-'./promote/promote_v1.gene_info.txt'

#Read the linker output using 50 bootstraps
linker_output_p <- './promote/sm_3500_linker_output_50bootstraps.rds'

#read the exp matrix
lognorm_est_counts_p <- './promote/mrm_sm_3500.txt'

#read the clinic file
phenotype_p <- './promote/promote_clinical_ours_fixed.txt'


drivers_gene_level <- function(linker_output_p, 
                               lognorm_est_counts_p,
                               gene_info_p,
                               phenotype_p,
                               final_sig_th=0.05,
                               include_cliques=FALSE,
                               ImpTH=0.05,
                               corrTH=0.7){
  
  #Initialize
  #geneinfo
  gi <- read.delim(gene_info_p)[,c('uniq_isos','regulator')]
  rownames(gi) <- gi[,1]
  #graphtomodules
  linker_object <- graph_to_modules(readRDS(linker_output_p),gi)
  linker_object <- graph_to_modules(linkeroutput,gi)
  #bootdist
  bootdist <- unlist(sapply(linker_object[[1]][[1]],function(x) x$bootstrap_idx))
  #get the expression matrix
  exp_mrm <- read.delim(lognorm_est_counts_p)
  
  #Do fast rewiring
  fast_rew <- function(ObjectList){
    
    #initialize common parameters
    regulator_info_col_name<-ObjectList$regulator_info_col_name
    phenotype_class_vals<-ObjectList$phenotype_class_vals
    phenotype_class_vals_label<-ObjectList$phenotype_class_vals_label
    outdir<-ObjectList$outdir
    orig_test_perms<-ObjectList$orig_test_perms
    retest_thresh<-ObjectList$retest_thresh
    retest_perms<-ObjectList$retest_perms
    logfile <- ObjectList$logfile
    
    #set.seed(1)
    dqrng::dqset.seed(1)
    
    #we create rundata and combine the modules of both parsers.
    #Lets generate the dupla (dataset's method - dataset's number)
    
    duplas <- unlist(lapply(seq_along(ObjectList$datasets),
                            function(i) paste(names(ObjectList$datasets[[i]]$rundata$modules),i)))
    
    #Initialize c_allstats and statsnames variables
    c_allstats <- c() #This is for the combined heatmap
    c_module_membership_list <- c() #This is for the combined heatmap
    
    statsnames <- c("module-method", "module-index", "orig-pval",
                    "revised-pvalue", "num-targets", "num-regulators",
                    "regulator-names","target-names", "num-samples", "num-genes",
                    "num-class1", "num-class2")
    
    for (dupla in duplas) {
      
      #Initialize rewired module's hash table
      module_membership_list <- hash::hash()
      
      #Initialize allstats array
      allstats <- NULL
      
      #For instance: 'VBSR X' to 'VBSR' and 'X'
      modmeth_i <- unlist(strsplit(dupla,' '))
      
      modmeth <- modmeth_i[1]
      i <- as.numeric(modmeth_i[2])
      
      #Output to the user which dupla we are working with
      message(modmeth,' ',i)
      
      rundata<-ObjectList$'datasets'[[i]]$rundata
      norm_expr_mat_keep<-ObjectList$'datasets'[[i]]$norm_expr_mat_keep
      keepsamps<-ObjectList$'datasets'[[i]]$keepsamps
      keeplabels<-ObjectList$'datasets'[[i]]$keeplabels
      class_counts<-ObjectList$'datasets'[[i]]$class_counts
      final_signif_thresh<-ObjectList$'datasets'[[i]]$final_signif_thresh
      responder<-ObjectList$'datasets'[[i]]$responder
      gene_info_df_keep<-ObjectList$'datasets'[[i]]$gene_info_df_keep
      name2idx<-ObjectList$'datasets'[[i]]$name2idx
      
      
      # This will register nr of cores/threads, keep this here
      # so the user can decide how many cores based on
      # their hardware.
      
      parallClass <- BiocParallel::bpparam()
      parallClass$workers <- ObjectList$NrCores
      
      GenerateStats <- function(mymod){
        
        signify <- NULL
        modregs <- unique(rundata$modules[[modmeth]][[mymod]]$regulators)
        modtargs <- unique(rundata$modules[[modmeth]][[mymod]]$target_genes)
        regnames <- paste(collapse = ", ", modregs)
        targnames <- paste(collapse = ", ", modtargs)
        keepfeats <- unique(c(modregs, modtargs))
        modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps])
        modmeth_i_c <- paste(modmeth_i,collapse=' ')
        
        orig_pval <- TraRe::rewiring_test(modmat, keeplabels + 1, perm = orig_test_perms)
        new_pval <- orig_pval
        stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3), signif(new_pval, 3),
                   length(modtargs), length(modregs), regnames,targnames, dim(modmat),
                   class_counts)
        if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
          #methods::show(paste(c("ModNum and NumGenes", mymod, length(keepfeats))))
          result <- TraRe::rewiring_test_pair_detail(modmat, keeplabels + 1,perm = retest_perms)
          new_pval <- result$pval
          stats <- c(modmeth_i_c, mymod, signif(orig_pval, 3),
                     signif(new_pval, 3), length(modtargs),
                     length(modregs), regnames,targnames, dim(modmat), class_counts)
          
          if (new_pval <= final_signif_thresh | new_pval == 1) {
            # save as list
            modname <- paste0(modmeth,'.',i, ".mod.", mymod)
            #module_membership_list[[modname]] <- keepfeats
            signify <- list(modname,keepfeats)
          }
          
        }
        return(list(stats,signify))
        
      }
      foreach_allstats <-BiocParallel::bplapply(seq_along(rundata$modules[[modmeth]]),
                                                GenerateStats, BPPARAM = parallClass)
      
      for (elements in foreach_allstats){
        
        #now we recover first allstats matrix
        foreach_stats <- elements[[1]]
        allstats<-rbind(allstats,foreach_stats)
        
        #and then update the module_membership dictionary
        hashtable <- elements[[2]]
        
        if (!is.null(hashtable)){
          
          module_membership_list[[hashtable[[1]]]] <- hashtable[[2]]
          
        }
      }
      
    }
    
    return(module_membership_list)
    
  }
  
  #generate dataframe with driver scores
  score_driv <- function(linkeroutput,driver,rew_list_names,bootdist,showmess=FALSE){
    
    #split due to cliques
    drivers <- unlist(strsplit(driver,'\\|'))
    
    oddsratio <- function(contig_tbl){
      #+0.5 using the Haldane Anscombe correction
      #(C/D) / (A/B)
      oddsr <- ((contig_tbl[2,1]+0.5) / (contig_tbl[2,2]+0.5)) / ((contig_tbl[1,1]+0.5) / (contig_tbl[1,2]+0.5))
      if (showmess){
        message(oddsr)
      }
      return(oddsr)
    }
    
    #go through every bootstrap
    pvals <- sapply(seq(50),function(x) {
      
      #message('Bootstrap: ', x)
      
      #All the modules within a bootstrap
      pos <- which(bootdist==x)
      
      #A (rewired modules)
      A_modnames <- intersect(pos,rew_list_names)-min(pos)+1
      #print(length(A_modnames))
      
      #Ac (non rewired modules)
      Ac_modnames <- setdiff(pos,rew_list_names)-min(pos)+1
      #print(length(Ac_modnames))
      
      bootstrap_universe <- sapply(pos,function(y) sum(drivers%in%linkeroutput[[1]][[1]][[y]]$regulators))
      
      A <- bootstrap_universe[A_modnames]
      Ac <- bootstrap_universe[Ac_modnames]
      
      A_0 <- sum(A==0)
      A_1 <- sum(A!=0)
      Ac_0 <- sum(Ac==0)
      Ac_1 <- sum(Ac!=0)
      
      # message(' A_0: ', A_0,
      #         ' A_1: ', A_1,
      #         ' Ac_0: ', Ac_0,
      #         ' Ac_1: ', Ac_1)
      
      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = TRUE))
      if (showmess){
        message('Every bootstrap')
        print(contig_tbl)
      }
      
      #show(contig_tbl)
      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative = 'l')
      
      #Pvalue!
      c(res$p.value, oddsratio(contig_tbl))
      
    })
    
    #evaluate all bootstraps at once
    pvals_all <- function(){
      
      #A (rewired modules)
      A_modnames <- rew_list_names
      #print(length(A_modnames))
      
      #Ac (non rewired modules)
      Ac_modnames <- setdiff(seq_along(bootdist),rew_list_names)
      #print(length(Ac_modnames))
      
      bootstrap_universe <- sapply(seq_along(bootdist),
                                   function(y) sum(drivers%in%linkeroutput[[1]][[1]][[y]]$regulators))
      
      A <- bootstrap_universe[A_modnames]
      Ac <- bootstrap_universe[Ac_modnames]
      
      A_0 <- sum(A==0)
      A_1 <- sum(A!=0)
      Ac_0 <- sum(Ac==0)
      Ac_1 <- sum(Ac!=0)
      
      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = TRUE))
      
      if (showmess){
        message('All bootstraps')
        print(contig_tbl)
      }
      
      #show(contig_tbl)
      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative='l')
      
      #Pvalue!
      c(res$p.value, oddsratio(contig_tbl))
      
      
    }
    
    #evaluate every 5 bootstraps
    pvals_b_5 <- sapply(seq(50/5),function(x) {
      
      #message('Bootstrap: ', x)
      
      #All the modules within a bootstrap
      pos <- which(bootdist%in%seq(x*5-4,x*5))
      
      #A (rewired modules)
      A_modnames <- intersect(pos,rew_list_names)-min(pos)+1
      #print(length(A_modnames))
      
      #Ac (non rewired modules)
      Ac_modnames <- setdiff(pos,rew_list_names)-min(pos)+1
      #print(length(Ac_modnames))
      
      bootstrap_universe <- sapply(pos,function(y) sum(drivers%in%linkeroutput[[1]][[1]][[y]]$regulators))
      
      A <- bootstrap_universe[A_modnames]
      Ac <- bootstrap_universe[Ac_modnames]
      
      A_0 <- sum(A==0)
      A_1 <- sum(A!=0)
      Ac_0 <- sum(Ac==0)
      Ac_1 <- sum(Ac!=0)
      
      # message(' A_0: ', A_0,
      #         ' A_1: ', A_1,
      #         ' Ac_0: ', Ac_0,
      #         ' Ac_1: ', Ac_1)
      
      
      
      contig_tbl <- as.table(matrix(c(A_0,Ac_0,A_1,Ac_1),
                                    ncol = 2, byrow = TRUE))
      
      if (showmess){
        message('Every 5 bootstraps')
        print(contig_tbl)
      }
      
      #Null hypothesis tail!
      res <- stats::fisher.test(contig_tbl,alternative='l')
      
      #Pvalue and odds ratio #(c/d) / (a/b)
      c(res$p.value, oddsratio(contig_tbl))
      
    })
    
    fishermeth <- function(pvals_v,showmess=FALSE){
      
      if (showmess){
        print(pvals_v)
      }
      
      #degrees of freedom
      df <- 2*length(pvals_v)
      
      #approximation by fisher
      newvar <- (-2) * sum(log(pvals_v))
      
      #Why upper tail??
      pchisq(newvar,df,lower.tail = FALSE)
      
    }
    
    
    #Reformat
    pvals <- c(fishermeth(pvals[1,]),
               mean(pvals[2,]))
    pvals_b_5 <- c(fishermeth(pvals_b_5[1,]),
                   mean(pvals_b_5[2,]))
    
    #return(min(pvals))
    results_l <- list(single_b = pvals,
                      all_b = pvals_all(),
                      every5_b = pvals_b_5)
    #return as a dataframe
    return(t(data.frame(do.call(c,results_l))))
    
  }
  
  #Include cliques (Without duplicities)
  cliques_w_duplicities <- function(regs_exp,regsofinterest,th){
    
    corr_m <- abs(cor(t(regs_exp)))
    
    corr_m <- (corr_m>=th) + 0
    
    g <- igraph::graph_from_adjacency_matrix(corr_m,mode = "undirected", weighted = TRUE)
    cliques_g <- igraph::max_cliques(g)
    
    # Add cliques
    interest_idx <- sapply(seq_along(cliques_g),
                           function(x) if (any(names(cliques_g[[x]])%in%regsofinterest)) return (x) else return(0))
    interest_idx <- interest_idx[interest_idx!=0]
    
    list_interest <- lapply(cliques_g[interest_idx],names)
    
    list_interest_c <- sapply(list_interest,function(x) paste0(x,collapse='|'))
    
    
    return(list_interest_c)
    
  }
  
  if (!file.exists(paste0(getwd(),'/Module 3500/50bootstraps/50b_graphs_rew_names.txt'))){
    
    message('Preparing rewiring...')
    
    #Prepare rewiring creating object
    preparedrewiring <- TraRe::preparerewiring(linker_output_p = linker_output_p,
                                               lognorm_est_counts_p = lognorm_est_counts_p,
                                               gene_info_p = gene_info_p,
                                               phenotype_p = phenotype_p,
                                               final_signif_thresh = final_sig_th)
    
    #Do the fast rewiring
    rew_list <- fast_rew(preparedrewiring)
    rew_list_names <- sub('[a-zA-Z]+.[0-9].mod.','',ls(rew_list))
    
    #write rew_list_names
    write.table(rew_list_names,file='./promote/50b_graphs_rew_names.txt',
                quote=F,sep='\t',row.names = FALSE,col.names=FALSE)
    
    
  }
  
  
  #read rew_list_names 
  rew_list_names <- read.delim('./promote/50b_graphs_rew_names.txt')[,1]

  #Intersect genes
  regs_mrm <- intersect(rownames(exp_mrm),gi$uniq_isos[gi[,2]==1])
  
  #confirm regulator genes
  rew_driv <- unique(unlist(sapply(rew_list_names,function(x) linker_object[[1]][[1]][[x]]$regulators)))
  
  #Build score matrix
  score_matrix <- t(sapply(rew_driv,function(x) score_driv(linker_object,x,rew_list_names,bootdist)))
  
  #Add colnames
  colnames(score_matrix) <- c('1B - FS','1B - OR',
                    'AllB - FS','AllB - OR',
                    '5B - FS','5B - OR')
  
  #List refinement
  imp <- as.numeric(which(apply(score_matrix,1,function(x) any(x[c(1,3,5)]<ImpTH))))
  
  #filter by refinement
  score_matrix_refined <- score_matrix[imp[order(score_matrix[imp,3])],]
  
  if (include_cliques){
    
    cliques_sm_3500 <- cliques_w_duplicities(exp_mrm[regs_mrm,],
                                             regsofinterest = rownames(score_matrix_refined),
                                             th = corrTH)
    
    cliques_sm_3500 <- union(rownames(score_matrix_refined),cliques_sm_3500)
    message('Drivers db augmented from ',nrow(score_matrix_refined), ' to ',length(cliques_sm_3500))
    
    #Define our final score matrix
    final_score <- t(sapply(cliques_sm_3500,function(x){
      #Select the first match (best score as they are sortered by significance)
      fmatch <- which(rownames(score_matrix_refined)%in%unlist(strsplit(x,'\\|')))[1]
      return(score_matrix_refined[fmatch,])
    }))
    #reorder again
    final_score <- final_score[order(final_score[,3]),]
    
  }else{
    final_score <- score_matrix_refined
  }
  
  return(final_score)
  
  
}

#write.table(final_score,'./promote/sm_3500_graphs_fisherpvals_newV.txt')
# impgenes <- drivers_gene_level(linker_output_p ,
#                                lognorm_est_counts_p,
#                                gene_info_p,
#                                phenotype_p,
#                                include_cliques=TRUE)
