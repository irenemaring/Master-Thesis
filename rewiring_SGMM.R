## rewiring
path <- getwd()

# Load files
linker_output_p <- "./SparseGMM/promote_10b_SGMM_lambda325_2nd.rds"
lo <- readRDS(linker_output_p)
geneinfo_p <- paste0(path,'/promote/promote_v1.gene_info.txt')
exp_p <- paste0(path,'/promote/mrm_sm_3500.txt')
clinic_p <- paste0(path,'/promote/promote_clinical_ours_fixed.txt')


preparerew <- TraRe::preparerewiring(name = "rewiring_promote_10b_lambda325_SGMM_2nd",
                                     linker_output_p = linker_output_p,
                                     lognorm_est_counts_p = exp_p,
                                     gene_info_p = geneinfo_p,
                                     phenotype_p = clinic_p,
                                     final_signif_thresh = 0.05,
                                     outdir = "rewirings/",
                                     nrcores = 20)

TraRe::runrewiring(preparerew)

