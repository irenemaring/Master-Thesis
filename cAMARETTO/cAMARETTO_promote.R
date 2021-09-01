# cAMARETTO SCRIPT
# Example for PROMOTE

library("CommunityAMARETTO")
library("foreach")
library(igraph)
library(dplyr)
source("./scripts/cAMARETTO_modified_functions_margaret.R")

#load data
P<-readRDS("./promote/sm_3500_linker_output_50bootstraps.rds") #TraRe GRN result
Res<-P$raw_results$VBSR$bootstrapResults[1:10]
names(Res)<-paste0("Run_",seq(length(Res)),sep = '')

output_directory_cAMARETTO = "./promote/cAMARETTO_report_TraRe/"
dir.create(output_directory_cAMARETTO)

genesets_database_reference <- "H_C2_genesets.gmt"
download.file(url="https://www.broadinstitute.org/~npochet/NotebookExample/ExampleData/H_C2_genesets.gmt", destfile=genesets_database_reference) 

# Run clustering
cAMARETTOresults<-cAMARETTO_Results_TraRe(AMARETTOresults_all = Res,
                                          gmt_filelist = NULL,
                                          NrCores = 30,
                                          drivers = FALSE)

saveRDS(cAMARETTOresults, paste0('./promote/prom_21k_50b_first_cAMARETTOresults.rds'))
#

#paramenters
m=80*10
j=0.2
cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork_TraRe(cAMARETTOresults,
                                                 pvalue = 0.05/m/m,
                                                 inter = 5,
                                                 jaccard=j,
                                                 edge_method="overlap",
                                                 plot_network = F)

cAMARETTOnetworkC<-cAMARETTO_IdentifyCom2(cAMARETTOnetworkM,
                                          filterComm = FALSE,
                                          ratioCommSize = 0.01,
                                          MinRuns = 2,
                                          ratioRunSize = 0.1,
                                          ratioEdgesInOut = 0.5)

saveRDS(cAMARETTOnetworkC,"./promote/cAMARETTO_report_TraRe/network_promote_21k_50b.rds")
# Run html - enrichment communities ---------------------------------------


cAMARETTO_HTMLreport(cAMARETTOresults=cAMARETTOresults,
                     cAMARETTOnetworkM=cAMARETTOnetworkM,
                     cAMARETTOnetworkC=cAMARETTOnetworkC,
                     PhenotypeTablesList = NULL,
                     output_address = output_directory_cAMARETTO,
                     HTMLsAMARETTOlist = NULL,
                     CopyAMARETTOReport = FALSE,
                     hyper_geo_reference = genesets_database_reference,
                     hyper_geo_reference_gp = NULL,
                     hyper_geo_reference_cp = NULL,
                     driverGSEA = TRUE,
                     NrCores=20)

