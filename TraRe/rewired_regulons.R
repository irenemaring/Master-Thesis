# Script para Jesus

library(TraRe)
library(igraph)

source(paste0(getwd(),'/scripts/format_graphs_to_modules_def.R'))
source(paste0(getwd(),'/scripts/rewiring_single_module.R'))

## INPUT FILE'S path:

# linker_output_p <- './promote/sm_3500_linker_output_50bootstraps.rds'
# 
# exp_p <- './promote/mrm_sm_3500.txt'
# 
# gene_info_p <-'./promote/promote_v1.gene_info.txt'
# 
# clinic_p <- './promote/promote_v1_clinical_ours_improved.txt'
# 
# rewired50boots_p <- './50b_graphs_rew_names.txt'

### SUMMARY of FUNCTIONS:

## function get_regulons: 
        # Takes as input linkeroutput's path, geneinfo's path and possibly 
        # the number of rewired modules' path (if wanted to get the regulons from only those modules).
        # Retrieves an object with the regulon subgraphs of the modules of linkeroutput and the boostrap idx
## function rewiring_regulon:
        # Takes as input the output of get_regulons, expression matrix's path, linkeroutput's path,
        # geneinfo's path, phenotype info's path and the significant threshold of the rewiring method.
        # Retreives an object that includes the regulon information and the pvalues obtained in the rewiring metho
## function fishermethod:
      # calculates the Fisher's pvalue (used in the following functions)
## targets_genelevel:
      # Takes as input the output of rewiring_regulon and retrieves a data.frame with the info related to 
      # each target (global multiplicity, mult. per regulons, associated regulators and pvals(Fisher).
## function merge_regulons:
      # Takes as input the output of rewirin_regulon and retrieves an object containing the info for 
      # the unified regulon (merged regulons withe the same driver).
## function filter_regulons: takes as input the output of merge_regulons and retrieves an object 
      # with a data.frame for each regulon with filtered targets by pvalue and its corresponding info



# get_regulons ---------------------------------------------------------------

get_regulons <- function(linker_output_p,
                         gene_info_p,
                         rewired50boots_p=NULL,
                         rewired = TRUE){
  
  
  #Read the linkeroutput
  
  linkeroutput <- readRDS(linker_output_p)
  
  #Read the geneinfo
  geneinfo <- read.delim(gene_info_p)[,c('uniq_isos','regulator')]
  rownames(geneinfo) <- geneinfo[,1]
  
  #Convert graph to modules
  linkeroutput_g <- graph_to_modules(linkeroutput,geneinfo)
  
  #Retrieve the conversion from old to original
  oldnumber <- linkeroutput_g$modules$VBSR_original
  
  if (rewired){
    #Read the rewired modules from the 50 bootstraps
    #rewired50boots <- read.delim(rewired50boots_p,header = F)[,1]
    
    #changed bc Jesus gave me the object of equivalences
    equiv_mods <- readRDS(rewired50boots_p)
    rewired50boots <- unlist(lapply(equiv_mods, function(x) x$orig))
    
    #Sort the index of the rewired modules
    modulesold <- sort(oldnumber[rewired50boots])
    
  }else{modulesold <- oldnumber}
  
  
  # Rewired graphs
  graph_objs <- lapply(modulesold,function(x)linkeroutput$graphs$VBSR$VBSR[[x]])
  
  #Get the bootstraps of the rewired graphs
  boost_rewired_mods <- sapply(modulesold,function(x) linkeroutput$modules$VBSR[[x]]$bootstrap_idx)
  
  #Get the regulons from all graphs 
  regulons_all_modules <- lapply(seq_along(graph_objs),function(idx){
                                                      g <- graph_objs[[idx]]
                                                      
                                                      # Regulators in module
                                                      regs_mod<- igraph::V(g)[type]$name
                                                      
                                                      out <- lapply(regs_mod,function(reg){
                                                        v <- igraph::V(g)[.nei(reg)]
                                                        regulons <- igraph::induced_subgraph(g,c(names(v),reg))
                                                      })
                                                      
                                                      return(list(regulons=out,boosts_idx=boost_rewired_mods[idx]))
                                                    })
  #Change regulons names
  names(regulons_all_modules) <- paste0("module_",modulesold)
  total_regulons <- sum(unlist(lapply(regulons_all_modules, function(x) length(x$regulons))))
  
  if (rewired) {
  message(paste(total_regulons,"regulons extracted from",length(rewired50boots), "rewired modules"))
    }else{message(paste(length(regulons_all_modules),"regulons extracted from all modules"))}

  return(regulons_all_modules)
}
# use
regulons_all_modules <- get_regulons(linker_output_p,gene_info_p,rewired50boots_p=rewired50boots_p,rewired = TRUE)

# rewiring_regulon ---------------------------------------------------------------

rewiring_regulon <- function(regulons_all_modules,
                             linker_output_p,
                             exp_p,
                             gene_info_p,
                             clinic_p,
                             final_signif_thresh=0.05){
  # Load rewiring object
  prepared <- TraRe::preparerewiring(linker_output_p = linker_output_p,
                                     lognorm_est_counts_p = exp_p,
                                     gene_info_p = gene_info_p,
                                     phenotype_p = clinic_p,
                                     final_signif_thresh = final_signif_thresh)
  
  #Unlist the regulons
  graphs_regulons <- lapply(regulons_all_modules, function(x) x$regulons)
  unlisted_regulons <- unlist(graphs_regulons,recursive = F,use.names = F)
  
  # Rewiring of regulon (slow)
  set.seed(1)
  message("Calculating the rewiring of the regulons")
  #Compute pvals for regulons on every rewired module
  pvals <- lapply(seq_along(unlisted_regulons),function(r){
    rewiring_single_module(prepared,unlisted_regulons[[r]])
  })
  pvals_vector <- unlist(pvals)
  
  # Calculate padj
  message('Calculating padjust values')
  padj <- p.adjust(pvals_vector,"BH")
  
  # Track of module number and bootstrap
    # Compute n? regulons/module:
  nr_regulons_per_module<- sapply(regulons_all_modules,function(x)length(x$regulons))
    # Create a vector with the number of regulons/modules
  module_nr_vector <- rep(names(nr_regulons_per_module),nr_regulons_per_module)
    # Extract the bootstrap of the regulon
  boost_rewired_mods <- as.numeric(sapply(regulons_all_modules,function(x) x$boosts_idx))
    # Create a vector with the number of bootstrap of the regulon
  boots_idx_vector <- rep(boost_rewired_mods,nr_regulons_per_module)
  
  
  # Compile everything in  list
  finals_adj <- lapply(seq_along(unlisted_regulons),function(r){
    
    list(driver=igraph::V(unlisted_regulons[[r]])[type]$name,
         
         targets=igraph::V(unlisted_regulons[[r]])[!type]$name,
         
         p_value=pvals_vector[r],
         
         padj=padj[r],
         
         module_nr=module_nr_vector[r],
         
         boots_idx=boots_idx_vector[r],
         
         G=unlisted_regulons[[r]])
  })
  return(finals_adj)
}

#use
finals_adj <- rewiring_regulon(regulons_all_modules,linker_output_p,exp_p,gene_info_p,clinic_p)
#saveRDS(finals_adj,  "./promote/regulons_finals_adj2.rds")
# fishermeth --------------------------------------------------------------

fishermeth <- function(pvals_v,showmess=FALSE,method='NA'){
  if (showmess){
    print(pvals_v)
  }
  if (method!='NA'){
    pvals_v <- p.adjust(pvals_v,method=method)
  }
  #degrees of freedom
  df <- 2*length(pvals_v)
  #approximation by fisher
  newvar <- (-2) * sum(log(pvals_v))
  #Why upper tail??
  pchisq(newvar,df,lower.tail = FALSE)
}


# list of TARGETS at the GENE level ---------------------------------------

targets_genelevel <- function(finals_adj){

# Get vector to iterate over
alltargs <- unique(unlist(sapply(finals_adj,function(x) x$targets)))
# Get list of targets in regulon to iterate over
targslist <- sapply(finals_adj,function(x) x$targets)

#Compute for each target within rewired regulons
#the corresponding stats
message("Calculating target's multiplicity across modules and pvalues")
data <- as.data.frame(t(as.data.frame(sapply(alltargs,function(x){
  
  bool_arr <- sapply(targslist,function(y) x%in%y) #Bool wether a target is in that element of the targslist
  
  pvals <- sapply(which(bool_arr),function(z){
    finals_adj[[z]]$p_value
  }) # if the target is in the list (bool_arr = TRUE) extract the padj of the regulon
  
  regs <- sort(unlist(sapply(which(bool_arr),function(k){
    finals_adj[[k]]$driver
  })))
  
  pvals_raw <- fishermeth(pvals) #calculate the p_value of Fisher method
  
  #Multiplicity within regulons
  m <- sum(bool_arr)
  #Multiplicity of regulators
  mult_regs  <- table(regs)
  n_occur <- paste(as.character(mult_regs),collapse=', ')
  
  regs <- paste(names(mult_regs), collapse='. ')
 
  
  return(list(mult_total=m, mult_per_reg=n_occur,regs=regs,pvals=pvals_raw))
  
}))))

}


# Use ---------------------------------------------------------------------

data <- targets_genelevel(finals_adj)

# TARGETS FILTERED by pval
filt_targ <- data[which(data$pvals<0.05),] 
head(filt_targ)
# Order by multiplicity
m <-filt_targ[order(as.data.frame(filt_targ$mult_total),decreasing = TRUE),]
head(rownames(m),10)
head(m[,c(1,4)])

# Order by pval
pfish <- filt_targ[order(as.data.frame(filt_targ$pvals),decreasing = FALSE),]
head(rownames(pfish),10)
head(pfish[,c(1,4)])


# Regulons ----------------------------------------------------------------

merge_regulons <- function(finals_adj){
  message("Calculating unique regulons")
  #((repeat))
  targs <- sapply(finals_adj,function(x) x$targets)  
  # should be equal to targslist: identical(targslist,targs)
  padjs <- sapply(finals_adj,function(x) x$p_value)    
  # should be equal to padj: identical(padj,padjs)
  
  #Calculate drivers within regulons
  uni_drivers <- unique(sapply(finals_adj,function(x){x$driver}))
  
  #Evaluate each driver that appear in regulons within rewired drivers
  regulons_uniq <- lapply(uni_drivers, function(x){
    
    # Idx where there is a regulon for that driver
    bool_driv<- sapply(finals_adj, function(y) x%in%y$driver)
    
    iter_targs <- targs[bool_driv]
    iter_padjs <- padjs[bool_driv]
    
    # List of targets of this driver
    targs_vec <- unique(unlist(iter_targs))
    
    # For each target get the pvals if there is an edge
    info_edge <- sapply(targs_vec,function(j){
      bool <- sapply(iter_targs, function(k){j%in%k})
      m <- sum(bool) #multiplicity
      pval_targ <- iter_padjs[bool]
      pfish <- fishermeth(unlist(pval_targ))
      list(pfish_targ=pfish, multiplicity=m)
    })
    
    p_value <- as.numeric(info_edge[1,])
    multiplicity <- as.numeric(info_edge[2,])
    
    list(driver=x,target=targs_vec, p_value=p_value,multiplicity=multiplicity)
  })
  
  #Give name to the drivers
  
  names(regulons_uniq) <- uni_drivers
  return(regulons_uniq)
}

# use
regulons_uniq <- merge_regulons(finals_adj)

# function filter ----------------------------------------------------------
filter_regulons <- function(regulons_uniq){
  # Boolean vector to filter the regulons <0.05
  bool <- sapply(regulons_uniq, function(x) x$p_value<0.05)
  
  # Save the names of the regulons
  uni_drivers <- names(regulons_uniq)
  
  
  filtered_regulons <- lapply(seq_along(regulons_uniq),function(x){
    filt_t <- regulons_uniq[[x]]$target[bool[[x]]]
    filt_p <- regulons_uniq[[x]]$p_value[bool[[x]]]
    filt_m <- regulons_uniq[[x]]$multiplicity[bool[[x]]]
    list(driver=regulons_uniq[[x]]$driver,targets=filt_t, p_value=filt_p,multiplicity=filt_m)
  })
  #Give name to the regulons
  names(filtered_regulons) <- uni_drivers
  
  # Take out empty regulons
  bool_empty_regulons <- sapply(seq_along(filtered_regulons),function(x){
    length(filtered_regulons[[x]]$targets)>1
  })
  
  filtered_regulons <- filtered_regulons[bool_empty_regulons]
  return(filtered_regulons)
}

#use
 filtered_regulons <- filter_regulons(regulons_uniq)

# Transform to data frame and sort by multiplicity
order_filtr_regulons <- lapply(filtered_regulons, function(x){
  y <- as.data.frame(x)
  y[order(y$multiplicity,decreasing = TRUE),]
})
#saveRDS(order_filtr_regulons, "./promote/order_filt_regulons_50.rds")

# igraphs

# order_filtr_regulons_graphs <-  lapply(order_filtr_regulons,function(x){
#   g <- igraph::graph_from_data_frame(x[,1:2],directed = F)
#   igraph::set_edge_attr(g, "weight",index = E(g), x$multiplicity)
# })

