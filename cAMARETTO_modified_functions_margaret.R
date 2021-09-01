#' @title cAMARETTO_Results
#' 
#' To initiate the community AMARETTO this results functions performs
#'  hyper geometric tests between all modules and gene sets.
#'
#' @param AMARETTOresults_all A list of multiple
#'  AMARETTO_Run outputs. The names are run names.
#' @param NrCores Nr of Cores that can be used to
#'  calculated the results.
#' @param output_dir A directory that stores gmt files
#'  for the runs
#' @param gmt_filelist A list with gmt files that are added into
#'  the communities. The names of the list are used in the networks.
#'   NULL if no list is added.
#' @param drivers Boolean that defines if only targets or drivers
#'  and targets are used to calculate the HGT.
#'
#' @return a list with cAMARETTOresults data objects from multiple runs
#' @importFrom gtools combinations
#' @importFrom dplyr arrange group_by left_join mutate select summarise 
#' rename  filter everything pull distinct case_when
#' @importFrom stats p.adjust phyper
#' @importFrom devtools install_github
#' @examples 
#' try(
#' cAMARETTOresults <- cAMARETTO_Results(AMARETTOresults_all,
#'   gmt_filelist=list(ImmuneSignature = Cibersortgmt),
#'   NrCores = 4 ,
#'   output_dir = "./")
#' )
#' @export
cAMARETTO_Results_TraRe <- function(AMARETTOresults_all,
                              NrCores=1,
                              output_dir="./",
                              gmt_filelist=NULL,
                              drivers = FALSE){
  
  #test if names are matching
  RunName1<-Geneset1<-RunName2<-Geneset2<-NULL
  runnames <- names(AMARETTOresults_all)
  if (!length(unique(runnames)) == length(runnames)){
    stop("The run names are not unique. Give unique names.")
  }
  dir.create(file.path(output_dir, "gmt_files"), recursive = FALSE,
             showWarnings = FALSE)
  # for each file a gmt for the modules
  create_gmt_filelist<-c()
  
  for (run in runnames){
    gmt_file <- file.path(output_dir,
                          "gmt_files",paste0(run, "_modules.gmt"))
    GmtFromModules(AMARETTOresults_all[[run]], gmt_file, run,
                   Drivers = drivers)
    create_gmt_filelist <- c(create_gmt_filelist,gmt_file)
  }
  names(create_gmt_filelist)<-runnames
  
  # add extra gmt files to compare with
  given_gmt_filelist <- c()
  if (!is.null(gmt_filelist)){
    for (gmt_file in gmt_filelist){
      gmt_file_path <- file.path(gmt_file)
      given_gmt_filelist <- c(given_gmt_filelist,gmt_file_path)
    }
  }
  names(given_gmt_filelist)<-names(gmt_filelist)
  all_gmt_files_list <- c(create_gmt_filelist,given_gmt_filelist)
  
  if (
    ! length(unique(names(all_gmt_files_list))) == 
    length(names(all_gmt_files_list))){
    stop("There is overlap between the gmt file names and run names")
  }
  
  if (length(all_gmt_files_list) < 2){
    stop("There are none or only one group given,
         community AMARETTO needs at least two groups.")
  }
  # compare gmts pairwise between runs
  all_run_combinations <- as.data.frame(gtools::combinations(
    n=length(all_gmt_files_list),
    r=2,
    v=names(all_gmt_files_list),
    repeats.allowed=FALSE))
  output_hgt_allcombinations <- apply(all_run_combinations, 1, function(x) {
    gmt_run1 <- all_gmt_files_list[x["V1"]]
    gmt_run2 <- all_gmt_files_list[x["V2"]]
    output_hgt_combination <- HyperGTestGeneEnrichment_JACCARD(gmt_run1,
                                                       gmt_run2,
                                                       runname1=as.character(x["V1"]),
                                                       runname2=as.character(x["V2"]),
                                                       NrCores=NrCores)
    return(output_hgt_combination)
  }
  )
  genelists<-lapply(all_gmt_files_list, function(x){
    genelist<-readGMT(x)
    names(genelist)<-gsub(" ","_",names(genelist))
    return(genelist)
  }
  )
  
  output_hgt_allcombinations <- do.call(rbind, output_hgt_allcombinations)
  output_hgt_allcombinations$padj <- stats::p.adjust(
    output_hgt_allcombinations$p_value,
    method="BH")
  output_hgt_allcombinations <- output_hgt_allcombinations %>% 
    dplyr::mutate(p_value=dplyr::case_when(
      Geneset1 == Geneset2~NA_real_, TRUE~p_value))
  output_hgt_allcombinations <- output_hgt_allcombinations %>%
    dplyr::mutate(Geneset1=ifelse(RunName1%in%names(given_gmt_filelist),
                                  paste0(RunName1,"|",gsub(" ","_",Geneset1))
                                  ,Geneset1),
                  Geneset2=ifelse(RunName2%in%names(given_gmt_filelist),
                                  paste0(RunName2,"|",gsub(" ","_",Geneset2)),
                                  Geneset2))
  
  # Extract relationship between genes and modules for all
  #AMARETTo files and the given gmt files.
  all_genes_modules_df<-Extract_Genes_Modules_All(AMARETTOresults_all,
                                                  gmt_filelist)
  
  return(list(runnames=runnames,
              gmtnames=names(given_gmt_filelist),
              hgt_modules=output_hgt_allcombinations,
              genelists = genelists,
              all_genes_modules_df=all_genes_modules_df,
              NrCores=NrCores))
  }



#' @title HyperGTestGeneEnrichment_JACCARD
#'
#' @param gmtfile1 A gmtfilename that you want to compare
#' @param gmtfile2 A second gmtfile to compare with.
#' @param runname1 name of the first dataset.
#' @param runname2 name of the second dataset.
#' @param NrCores  Number of cores for parallel computing. 
#' @param ref.numb.genes The reference number of genes.
#' @importFrom doParallel registerDoParallel 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats p.adjust phyper
#' 
#' @return Creates resultfile with p-values and padj when
#'  comparing two gmt files with a hyper geometric test.
#' @examples 
#' try(
#' HyperGTestGeneEnrichment(gmtfile1, gmtfile2,
#' runname1, runname2,NrCores, ref.numb.genes=45956)
#' )
#' @export
HyperGTestGeneEnrichment_JACCARD<-function(gmtfile1, gmtfile2,
                                   runname1, runname2,
                                   NrCores, ref.numb.genes=45956){
  i<-j<-NULL
  gmtfile1<-readGMT(gmtfile1) # our gmt_file_output_from Amaretto
  gmtfile2<-readGMT(gmtfile2)  # the hallmarks_and_co2...
  ###########################  Parallelizing :
  cluster <- parallel::makeCluster(c(rep("localhost", NrCores)),
                                   type = "SOCK")
  doParallel::registerDoParallel(cluster,cores=NrCores)
  resultloop<-foreach::foreach(j=seq_len(length(gmtfile2)),
                               .combine='rbind') %do% {
                                 #print(j)
                                 foreach::foreach(i=seq_len(length(gmtfile1)),.combine='rbind') %dopar% {
                                   #print(i)
                                   l<-length(gmtfile1[[i]])
                                   k<-sum(gmtfile1[[i]] %in% gmtfile2[[j]])
                                   u<-length(union(gmtfile1[[i]],gmtfile2[[j]]))
                                   m<-ref.numb.genes
                                   n<-length(gmtfile2[[j]])
                                   p1<-stats::phyper(k-1,l,m-l,n,lower.tail=FALSE)
                                   overlapping.genes<-gmtfile1[[i]][gmtfile1[[i]] %in% gmtfile2[[j]]]
                                   overlapping.genes<-paste(overlapping.genes,collapse = ', ')
                                   c(RunName1=runname1,
                                     RunName2=runname2,
                                     Geneset1=names(gmtfile1[i]),
                                     Geneset2=names(gmtfile2[j]),
                                     p_value=p1,n_Overlapping=k,
                                     Overlapping_genes=overlapping.genes,
                                     jaccard_idx=k/u,
                                     union=u)
                                 }
                               }
  parallel::stopCluster(cluster)
  resultloop<-as.data.frame(resultloop,stringsAsFactors=FALSE)
  resultloop$p_value<-as.numeric(resultloop$p_value)
  resultloop$n_Overlapping<-as.numeric((resultloop$n_Overlapping))
  resultloop$jaccard_idx<-as.numeric((resultloop$jaccard_idx))
  resultloop$union<-as.numeric((resultloop$union))
  resultloop[,"padj"]<-stats::p.adjust(resultloop[,"p_value"],method='BH')
  return(resultloop)
}

#' @title cAMARETTO_ModuleNetwork
#' 
#' Creates a module network.
#'
#' @param cAMARETTOresults cAMARETTOresults object
#' @param pvalue pvalue cut-off for each gene set
#' @param inter minimal overlap between two gene sets
#' @param color_list An optional list with colors
#' @param edge_method Define edge weights based 
#' p-values or overlap between the gene sets.
#' @param plot_network If TRUE, plots the Community Network at the end.
#'
#' @return a list with the module network, layout for the network,
#'  used p-value, used overlap en colors
#' @importFrom utils globalVariables
#' @importFrom randomcoloR randomColor
#' @importFrom igraph as_data_frame degree E graph_from_data_frame 
#' layout_with_fr V graph.data.frame norm_coords edge.betweenness.community
#' @importFrom dplyr filter mutate select
#' 
#' @examples 
#' try(
#' cAMARETTOnetworkM<-cAMARETTO_ModuleNetwork(cAMARETTOresults,0.10,5)
#' )
#' @export
cAMARETTO_ModuleNetwork_TraRe<-function(cAMARETTOresults,
                                  pvalue = 0.05,
                                  inter = 5,
                                  jaccard = 0.1,
                                  color_list = NULL,
                                  edge_method = "pvalue",
                                  plot_network = TRUE){
  
  ############################################################
  #if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
  ############################################################
  padj<-n_Overlapping<-modulenames<-NULL
  RunName1<-RunName2<-NULL
  output_hgt_allcombinations_filtered <- cAMARETTOresults$hgt_modules %>% 
    dplyr::filter(padj<=pvalue & n_Overlapping>=inter & jaccard_idx>=jaccard)
  node_information <- as.data.frame(unique(
    c(output_hgt_allcombinations_filtered$Geneset1,
      output_hgt_allcombinations_filtered$Geneset2)
  ))
  colnames(node_information) <- c("modulenames")
  node_information <- node_information %>%
    dplyr::mutate(run=sub("\\|.*$","",modulenames))
  module_network <- graph_from_data_frame(
    d=output_hgt_allcombinations_filtered%>%
      dplyr::select(-RunName1,-RunName2),
    vertices=node_information, directed=FALSE)
  if (is.null(color_list)){
    color_list <- randomcoloR::randomColor(length(c(
      cAMARETTOresults$runnames,cAMARETTOresults$gmtnames)),
      luminosity="dark")
    names(color_list) <- c(cAMARETTOresults$runnames,
                           cAMARETTOresults$gmtnames)
  } else {
    c(cAMARETTOresults$runnames,
      cAMARETTOresults$gmtnames) %in% names(color_list)
  }
  igraph::V(module_network)$color <- color_list[igraph::V(module_network)$run]
  igraph::V(module_network)$size <- 2*sqrt(igraph::degree(module_network,
                                                          mode="all"))
  if (edge_method=="pvalue"){
    igraph::E(module_network)$width <-
      -(log10(igraph::E(module_network)$p_value))*0.2
  } else if (edge_method=="overlap"){
    igraph::E(module_network)$width <- igraph::E(
      module_network)$n_Overlapping/16
  } else {
    stop("The edge method is not properly defined.")
  }
  layoutMN <- layout_with_fr(module_network)
  if(plot_network){
    plot(module_network, vertex.frame.color=NA,
         layout=layoutMN, vertex.label=NA, main="Module network",
         edge.color="gray80")
    legendMN <- legend(x = -5.5,
                       y = -1.1+0.05*length(color_list),
                       legend = names(color_list),
                       col = color_list,
                       pch=19,
                       bty="n",
                       ncol=ceiling(length(color_list)/5))
    legendMN
  }
  return(list(module_network=module_network,
              layoutMN=layoutMN, pvalue=pvalue,
              inter=inter, colMN=color_list))
}


#' Title HyperGeoEnrichmentTest
#'
#' @param ResultsGRN GRN results output
#' @param hyper_geo_reference GMT file with gene sets to compare with.
#' @param driverGSEA if TRUE, module drivers will also be included in the hypergeometric test.
#' @param NrCores Number of cores for parallel processing. 
#'
#' @return Hyper-Geometric Enrichment Test table
#' @export
#'
#' @examples HyperGeoEnrichmentTest(AMARETTOresults=NULL, hyper_geo_reference, driverGSEA=TRUE, MSIGDB=TRUE, NrCores=4)
HyperGeoEnrichmentTest_TraRe<-function(ResultsGRN, hyper_geo_reference, driverGSEA=TRUE, NrCores=4){
  output_hgt_all<-NULL
  for(i in 1:length(hyper_geo_reference)){
    if (is.null(ResultsGRN)){
      return(1)
    }
    CommunityAMARETTO::GmtFromModules(ResultsGRN, driverGSEA)
    output_hgt <- CommunityAMARETTO::HyperGTestGeneEnrichment(hyper_geo_reference[i], "./Modules_genes.gmt", NrCores)
    load("./MsigdbMapping.rda")#utils::data(MsigdbMapping)
    MsigdbMapping<-MsigdbMapping%>%dplyr::mutate(url=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',geneset,'">',gsub("_"," ",geneset),'</a>'))
    output_hgt<-output_hgt%>%dplyr::left_join(MsigdbMapping,by=c("Geneset"="geneset"))%>%
      dplyr::mutate(description=ifelse(is.na(description),Geneset,description))%>%
      dplyr::mutate(Geneset=ifelse(is.na(url),Geneset,url))%>%dplyr::rename("Description"="description")%>%dplyr::select(-url)
    cat("The hyper geometric test results are calculated.\n")
    output_hgt_all<-rbind(output_hgt_all,output_hgt)
  }
  return(output_hgt_all)
}


# plot communities --------------------------------------------------------

cAMARETTO_IdentifyCom2<- function (cAMARETTOnetworkM, color_list = NULL, filterComm = TRUE, 
                                   ratioCommSize = 0.01, MinRuns = 2, ratioRunSize = 0.1, ratioEdgesInOut = 0.5, 
                                   plot_network = TRUE) 
{
  from <- to <- Community <- numEdgesInComm <- NULL
  totalNumEdges <- numEdgesNotInComm <- numTotalEdgesInCommunity <- NULL
  numTotalEdgesNotInCommunity <- numDatasetsPerCommunity <- CommSize <- NULL
  NewComNumber <- CommsizeFrac <- fractDatasetsSize <- fractEdgesInVsOut <- NULL
  . <- nodeName <- NULL
  comm <- igraph::edge.betweenness.community(cAMARETTOnetworkM$module_network, 
                                             directed = FALSE, merges = TRUE, modularity = TRUE, membership = TRUE)
  message("There are ", length(unique(comm$membership)), 
          " different communities detected using weighted edges.")
  names(comm$membership) <- igraph::V(cAMARETTOnetworkM$module_network)$name
  membership <- as.data.frame(cbind(c(seq_len(length(comm$membership))), 
                                    comm$membership))
  colnames(membership) <- c("nodeID", "Community")
  numCommunitiesOrig <- length(unique(membership[, "Community"]))
  membership <- tibble::rownames_to_column(membership, "nodeName") %>% 
    dplyr::mutate(run = sub("|Module_.*$", "", 
                            nodeName))
  Edges_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, 
                                          what = "edges")
  Nodes_Mnetwork <- igraph::as_data_frame(cAMARETTOnetworkM$module_network, 
                                          what = "vertices")
  for (m in seq_len(nrow(membership))) {
    commNum <- membership[m, "Community"]
    Id <- membership[m, "nodeName"]
    edgeMatrixVector <- unlist(Edges_Mnetwork %>% dplyr::filter(from == 
                                                                  Id | to == Id) %>% dplyr::select(from, to))
    edgeMatrixVector <- edgeMatrixVector[-which(edgeMatrixVector == 
                                                  Id)]
    membership[m, "totalNumEdges"] <- length(edgeMatrixVector)
    membership[m, "numEdgesInComm"] <- nrow(membership[match(edgeMatrixVector, 
                                                             membership[, "nodeName"]), ] %>% dplyr::filter(Community == 
                                                                                                              commNum))
  }
  membership <- membership %>% dplyr::mutate(fractEdgesInOut = numEdgesInComm/totalNumEdges, 
                                             numEdgesNotInComm = totalNumEdges - numEdgesInComm)
  commEdgeInfo <- membership %>% dplyr::group_by(Community) %>% 
    dplyr::summarise(numTotalEdgesInCommunity = sum(numEdgesInComm)/2, 
                     numTotalEdgesNotInCommunity = sum(numEdgesNotInComm), 
                     fractEdgesInVsOut = numTotalEdgesInCommunity/(numTotalEdgesNotInCommunity + 
                                                                     numTotalEdgesInCommunity), numDatasetsPerCommunity = length(unique(run)), 
                     CommSize = n(), fractDatasetsSize = numDatasetsPerCommunity/CommSize, 
                     CommsizeFrac = CommSize/nrow(Nodes_Mnetwork))
  commEdgeInfo <- commEdgeInfo %>% dplyr::arrange(-CommSize) %>% 
    dplyr::mutate(NewComNumber = row_number())
  suppressMessages(membership <- dplyr::left_join(membership, 
                                                  commEdgeInfo %>% dplyr::select(Community, NewComNumber)))
  if (filterComm == TRUE) {
    KeepCommEdgeInfo <- commEdgeInfo %>% dplyr::filter(CommsizeFrac >= 
                                                         ratioCommSize & numDatasetsPerCommunity >= MinRuns & 
                                                         fractDatasetsSize >= ratioRunSize & fractEdgesInVsOut >= 
                                                         ratioEdgesInOut)
    message("There are ", nrow(commEdgeInfo) - nrow(KeepCommEdgeInfo), 
            " communities to remove.")
  }
  else {
    KeepCommEdgeInfo <- commEdgeInfo %>% dplyr::filter(numDatasetsPerCommunity >= 
                                                         MinRuns)
  }
  Nodes_Mnetwork <- dplyr::left_join(Nodes_Mnetwork, membership %>% 
                                       select(-run), by = c(name = "nodeName"))
  CommGraph <- igraph::graph.data.frame(Edges_Mnetwork, directed = FALSE, 
                                        vertices = data.frame(Nodes_Mnetwork))
  graph.degrees <- igraph::degree(CommGraph)
  igraph::V(CommGraph)$size <- 2 * sqrt(graph.degrees)
  community_list_df <- membership %>% dplyr::filter(Community %in% 
                                                      KeepCommEdgeInfo$Community) %>% dplyr::select(nodeName, 
                                                                                                    Community) %>% split(.$Community)
  community_list <- lapply(community_list_df, function(x) unlist(x$nodeName))
  names(community_list) <- names(community_list_df)
  if (is.null(color_list)) {
    color_list <- randomcoloR::randomColor(length(community_list), 
                                           luminosity = "light")
    names(color_list) <- names(community_list)
  }
  else {
    length(color_list) >= length(community_list)
  }
  if (plot_network) {
    layout_1 <- cAMARETTOnetworkM$layoutMN
    layout_1 <- igraph::norm_coords(layout_1, ymin = -1, 
                                    ymax = 1, xmin = -1, xmax = 1)
    plot(CommGraph, rescale = FALSE, layout = layout_1*0.5, 
         vertex.color = as.character(Nodes_Mnetwork$color),
          
         vertex.label = NA, mark.groups = community_list,
         #mark.border=color_list,
         mark.border="black",
         mark.col = color_list,
        # main = "Rewired Modules Community network",
         main = "", #Sparse Gaussian Mixture Models
         sub = paste("Jaccard index", j))
    legendMN <- legend(x = -1.85, y = 0.5,# + 0.05 * length(cAMARETTOnetworkM$colMN),
                       legend = names(cAMARETTOnetworkM$colMN), col = cAMARETTOnetworkM$colMN,
                       pch = 19, bty = "n", ncol = ceiling(length(cAMARETTOnetworkM$colMN)/5))
    legendMN
    # legendCOM <- legend(x = 0.85, y = 1, legend = names(color_list), 
    #                     col = color_list, pch = 19, bty = "n",
    #                     cex = max(0.9,1/(1 + 0.02 * length(color_list))),
    #                     ncol = ceiling(length(color_list)/15))
    legendCOM <- legend(x =1, y = 0.5, legend = names(color_list),
                        col = color_list, pch = 19, bty = "n",
                        cex = max(0.9,1/(1 + 0.02 * length(color_list))),
                        ncol = 4)
    legendCOM
  }
  return(list(CommGraph = CommGraph, community_list = community_list, 
              commEdgeInfo = commEdgeInfo, color_list = color_list))
}


