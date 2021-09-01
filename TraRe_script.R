# SCRIPT TO RUN ALL STEPS OF TraRe

library(TraRe)

# INPUT -------------------------------------------------------------------

path <- paste0(getwd(),'/promote')


# PROMOTE
# Expression matrix
exp_p <- paste0(path,'/mrm_sm_3500.txt')
lognorm_est_counts <- read.table(exp_p)

# Gene info: drivers and targets
gene_info_p <- paste0(path,'/promote_v1.gene_info.txt')
gene_info <- read.delim(gene_info_p)
regulator_name <- gene_info[which(gene_info$regulator==1),1]
target_name <- gene_info[which(gene_info$regulator==0),1]

regulator_filtered_idx<-which(rownames(lognorm_est_counts)%in%regulator_name)
target_filtered_idx <- which(rownames(lognorm_est_counts)%in%target_name)

# Clinical info (Responder - non-responders)
clinic_p <- paste0(path,'/promote_clinical_ours_fixed.txt')


# GRN inference method ----------------------------------------------------

linkeroutput <- LINKER_run(lognorm_est_counts =  as.matrix(lognorm_est_counts),
                            target_filtered_idx = target_filtered_idx,
                            regulator_filtered_idx = regulator_filtered_idx,
                            link_mode = "VBSR", #phase 1
                            graph_mode = "VBSR", #phase 2
                            NrModules = 100, #default
                            Nr_bootstraps = 10, #default
                            NrCores = 32)
#saveRDS(linkeroutput,file=paste0(path,"sm_3500_linker_output.rds"))
# Rewiring ----------------------------------------------------------------

prew_output <- preparerewiring(name = "rewiring_test_promote",
                               linker_output_p = paste0(path,'sm_3500_linker_output.rds'),
                               lognorm_est_counts_p = exp_p,
                               gene_info_p = gene_info_p,
                               phenotype_p = clinic_p,
                               outdir = paste0(path,'/rewiring_tests/'), 
                               final_signif_thresh = 0.05,
                               nrcores = 20)
runrewiring(prew_output)


# excels -------------------------------------------------------------------


drivers <- lognorm_est_counts[rownames(lognorm_est_counts)%in%regulator_name,]
excel_generation(gpath = paste0(path,'/rewiring_tests/rewiring_test_promote_0.05/supermodule_1/refinedsumm.rds'),
                 wpath = paste0(path,'/rewiring_tests/rewiring_test_promote_0.05/'),
                 dataset = drivers)

