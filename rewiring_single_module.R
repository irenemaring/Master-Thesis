rewiring_single_module <- function(ObjectList,module,mymod=1){
  if (inherits(module,'igraph')){
    genes <- names(igraph::V(module))
  }
  #We are working with just 1
  i<-1
  #Needed info
  regulator_info_col_name<-ObjectList$regulator_info_col_name
  phenotype_class_vals<-ObjectList$phenotype_class_vals
  phenotype_class_vals_label<-ObjectList$phenotype_class_vals_label
  outdir<-ObjectList$outdir
  orig_test_perms<-ObjectList$orig_test_perms
  retest_thresh<-ObjectList$retest_thresh
  retest_perms<-ObjectList$retest_perms
  logfile <- ObjectList$logfile
  rundata<-ObjectList$'datasets'[[i]]$rundata
  norm_expr_mat_keep<-ObjectList$'datasets'[[i]]$norm_expr_mat_keep
  keepsamps<-ObjectList$'datasets'[[i]]$keepsamps
  keeplabels<-ObjectList$'datasets'[[i]]$keeplabels
  class_counts<-ObjectList$'datasets'[[i]]$class_counts
  final_signif_thresh<-ObjectList$'datasets'[[i]]$final_signif_thresh
  responder<-ObjectList$'datasets'[[i]]$responder
  gene_info_df_keep<-ObjectList$'datasets'[[i]]$gene_info_df_keep
  name2idx<-ObjectList$'datasets'[[i]]$name2idx
  #Now rewiring test!
  signify <- NULL
  modregs <- intersect(ObjectList$datasets[[i]]$allregs,genes)
  modtargs <- intersect(ObjectList$datasets[[i]]$alltargs,genes)
  regnames <- paste(collapse = ", ", modregs)
  targnames <- paste(collapse = ", ", modtargs)
  keepfeats <- unique(c(modregs, modtargs))
  modmat <- t(norm_expr_mat_keep[keepfeats, keepsamps])
  orig_pval <- TraRe::rewiring_test(modmat, keeplabels + 1, perm = orig_test_perms)
  new_pval <- orig_pval
  # stats <- c( mymod, signif(orig_pval, 3), signif(new_pval, 3),
  #             length(modtargs), length(modregs), regnames,targnames, dim(modmat),
  #             class_counts)
  if (orig_pval < retest_thresh | orig_pval == 1 | mymod %% 300 == 0) {
    #methods::show(paste(c("ModNum and NumGenes", mymod, length(keepfeats))))
    result <- TraRe::rewiring_test_pair_detail(modmat, keeplabels + 1,perm = retest_perms)
    new_pval <- result$pval}
  #   stats <- c( mymod, signif(orig_pval, 3),
  #               signif(new_pval, 3), length(modtargs),
  #               length(modregs), regnames,targnames, dim(modmat), class_counts)
  #   if (new_pval <= final_signif_thresh | new_pval == 1) {
  #     # save as list
  #     modname <- mymod
  #     message('TRUE')
  #     #module_membership_list[[modname]] <- keepfeats
  #     signify <- list(modname,keepfeats)
  #   }
  # }
  return(new_pval)
}