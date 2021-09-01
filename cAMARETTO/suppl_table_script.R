# Script DATA for SUPPLEMENTARY TABLE

# 1. [[rewiring_equivalency]]: indexing object to transform the number of a
#    rewired module across all bootstraps (1-4088) to the number of the module
#    within a bootstrap .
# 2. [[rewcom]]: vector with the index number of the Communities in cAmaretto results
#    object that contain rewired modules.
# 3. [[comm_summ]]: object with the info shown in the supplementary table 2 (Community summary table 10b) 

# 1. Get index objects -------------------------------------------------------


# GRN inferred method output
linker_output_p <- './promote/sm_3500_linker_output_50bootstraps.rds'

linkeroutput <- readRDS(linker_output_p)

## Get rewiring equivalency of modules in bootstrap and total modules taking into account the lost graphs
# Script graphs to modules
source("./scripts/format_graphs_to_modules_def.R")

# Rewired modules in 50 bootstraps
rewired50boots <- read.delim('promote/50b_graphs_rew_names.txt',header = F)$V1

geneinfo <- read.delim('promote/promote_v1.gene_info.txt')[,c('uniq_isos','regulator')]
rownames(geneinfo) <- geneinfo[,1]

#Get graphs
linkeroutput_g <- graph_to_modules(linkeroutput,geneinfo)

#Get number of rewired module within bootstrap

oldnumber <- linkeroutput_g$VBSR$VBSR_original_index
rewiredmodulesold <- oldnumber[rewired50boots]

# Bootstrap of rewired modules
boost_rewired_mods <- sapply(rewiredmodulesold,function(x) linkeroutput$modules$VBSR[[x]]$bootstrap_idx)

# Number of modules in each bootstrap
nrmod <- sapply(linkeroutput$raw_results$VBSR$bootstrapResults,function(x)x$NrModules)

accum <- c(0)
rewiring_equivalency <- list()

for (i in seq(1,length(nrmod))){
  idx <- which(boost_rewired_mods==i)
  idx2 <- accum-rewiredmodulesold[idx]
  accum <- c(accum+nrmod[i])
  rewiring_equivalency[[i]] <- list(n_mod=rewiredmodulesold[idx],n_mod_bootstr=abs(idx2))
}
# Conversion object
rewiring_equivalency

# cAmaretto's output
network02_p <- './network_promote_21k_50b.rds'
network02 <- readRDS(network02_p)

check_rewiring_Amaretto <- function(Amaretto_communities,Jaccard,rew_equivalency,Nboot=10){
  message('Jaccard index set to ',Jaccard)
  amaretto_length <- sapply(Amaretto_communities$community_list,length)
  #More than 5 modules
  message('Communities with more than 5 modules ',sum(amaretto_length>5))
  message('% of that modules ', sum(amaretto_length>5)/length(amaretto_length))
  #Define current dataset
  current_dataset <- c(rep(1,Nboot/2),rep(2,Nboot/2))
  #Rewired modules
  rewired_modules_comms <- sapply(seq_along(Amaretto_communities$community_list),function(i){
    x <- Amaretto_communities$community_list[[i]]
    #split by |
    split_x <- strsplit(x,'\\|')
    run_mods_numbers <- lapply(split_x,function(x) as.numeric(gsub('[^0-9]','',x)))
    rew_mods_in_community <- sum(sapply(run_mods_numbers,function(x){
      x[2]%in%rew_equivalency[[x[1]]]$'n_mod_bootstr'
    }))
    return(rew_mods_in_community)
  })
  
  message('% of modules that have rewiring modules ',sum(rewired_modules_comms!=0)/sum(rewired_modules_comms))
  return(list(rew_mods = rewired_modules_comms,
              main= paste0('Jaccard index of ',Jaccard)))
}

nr_modules_comms <- check_rewiring_Amaretto(Amaretto_communities=network02,Jaccard=0.2,rew_equivalency=rewiring_equivalency)
barplot(nr_modules_comms$rew_mods,)

rewcom <- which(nr_modules_comms$rew_mods!=0)

# Known rewired Communities
#rewcom <- c(8,9,13,21,22,24,25,40,46,48,64)
# 3. summary table object ----------------------------------------------------

# Rewiring equiv for 10 first bootstraps
rew_mods_10b <- rewiring_equivalency[1:10]

# Number of modules in bootstrap (for the 10 firs bootstraps)
nrmod10 <- nrmod[1:10]

#Module result of 10 first bootstraps
linkeroutput10 <- linkeroutput$modules$VBSR[1:sum(nrmod10)]


# Object containing the info for supplementary table and community summary
#rewcom
mycoms <- c(1,7,8,9,13,21,22,40,48)
comm_summ <- lapply(mycoms, function(comm){
  
  string <- network02$community_list[[comm]]
  string_split <- strsplit(string,'|', fixed = T)
  string_list<- lapply(string_split,function(x)as.numeric(gsub('[^0-9]','',x)))
  
  string_df <- as.data.frame(do.call(rbind,string_list))
  colnames(string_df) <- c("Bootstrap", "Module")
  data <- string_df[order(string_df$Bootstrap),]
  rownames(data) <- NULL
  
  # get the regs and targs
  targets_uniq <- c()
  for (j in seq(1,length(data$Bootstrap))){
    boot <- data$Bootstrap[j]
    rew <-  data$Module[j]
    # For 10 boot
    mod <- sum(nrmod10[1:boot-1])+data[j,2]
    
    targets_mod <- length(linkeroutput10[[mod]]$target_genes)
    regulators_mod <- linkeroutput10[[mod]]$regulators
    data$n_targs[j] <- targets_mod
    data$regs[j] <- paste(regulators_mod, collapse = ' ')
    
    data$rewired[j] <- rew%in%rew_mods_10b[[boot]]$n_mod_bootstr
    targets_uniq <- c(targets_uniq,linkeroutput10[[mod]]$target_genes)
  }
  
  targets_uniq <- length(unique(targets_uniq))
  #message(targets_uniq)
  regs_uniq <- sort(table(unlist(stringr::str_split(paste(data$regs, collapse = ' '), " "))), decreasing = T)
  return(list(data=data, n_targs=targets_uniq, regs_u=regs_uniq,n_rewired=sum(data$rewired)))
})

names(comm_summ) <- paste0("Com", mycoms) #rewcom