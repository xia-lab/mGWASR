##################################################
## R script for mGWAS
## Description: GO/Pathway ORA 
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformNetEnrichment
#' @export 
PerformNetEnrichment <- function(file.nm, fun.type, IDs){
    # prepare query
    ora.vec <- NULL;
     # net is tf/mir/drug, they already in entrez
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- as.character(ora.vec);
    res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec);
    return(res);
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param ora.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[qs]{qread}}
#'  \code{\link[rjson]{toJSON}}
#' @rdname PerformEnrichAnalysis
#' @export 
#' @importFrom qs qread
#' @importFrom rjson toJSON
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec){
  org.code <- "hsa";
  if(fun.type %in% c("keggc", "smpdb")){
    .load.enrich.compound.lib(org.code, fun.type);
  }else{
    .load.enrich.lib(org.code, fun.type);
  }
  
  # prepare query
  ora.nms <- names(ora.vec);
  ora <<- ora.vec;
  
  current.geneset <- qs::qread("current.geneset.qs");
  current.universe <- unique(unlist(current.geneset));

  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes 
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  
  if(fun.type %in% c("keggc", "smpdb")){
    hits.query <- lapply(current.geneset, function(x){x[names(x) %in% ora.vec]});
    hits.query<- lapply(hits.query, function(x){names(x)});
  }else{
    hits.query <- lapply(current.geneset, function(x) { ora.nms[ora.vec%in%unlist(x)];});
  }
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat) == 0){
    AddErrMsg("No hits found for your query!");
    return(0);
  }
  
  # order by p value
  ord.inx<-order(res.mat[,4]);
  res.mat <- signif(res.mat[ord.inx,],3);
  hits.query <- hits.query[ord.inx];
  imp.inx <- res.mat[,4] <= 0.05;
  if(sum(imp.inx) < 10){ # too little left, give the top ones
    topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
    res.mat <- res.mat[1:topn,];
    hits.query <- hits.query[1:topn];
  }else{
    res.mat <- res.mat[imp.inx,];
    hits.query <- hits.query[imp.inx];
    if(sum(imp.inx) > 120){
      # now, clean up result, synchronize with hit.query
      res.mat <- res.mat[1:120,];
      hits.query <- hits.query[1:120];
    }
  }
  
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  AddMsg("Functional enrichment analysis was completed");
  
  # write json
  fun.anot = hits.query; 
  fun.pval = resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  hit.num = resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  if(fun.type %in% c("keggc", "smpdb")){
    fun.ids <- as.vector(current.setids[which(current.setids %in% names(fun.anot))]); 
    names(fun.anot) = as.vector(names(current.setids[which(current.setids %in% names(fun.anot))]));
  }else{
    fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  }
  
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    hit.num = hit.num
  );
  json.mat <- rjson::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- resTable[,5];
  hit.num <<- resTable[,4];
  csv.nm <- paste(file.nm, ".csv", sep="");
  fast.write.csv(resTable, file=csv.nm, row.names=F);
  return(1);
}

# these are met set libraries
.load.enrich.compound.lib<-function(org.code, fun.type){

  # prepare lib
  is.go <- FALSE;
  mSetObj <- .get.mSet();
  if(tolower(fun.type) == 'smpdb'){ 
    nm <- mSetObj$org;
    sub.dir <- paste0("smpdb");
  }else if(tolower(fun.type) == 'keggc'){ 
    nm <- "hsa";
    sub.dir <- paste0("kegg/metpa");
  }
  
  lib.nm <- paste0(nm, ".qs");
  
  my.lib <- .get.my.lib(lib.nm, sub.dir);
  
  current.setlink <- my.lib$link;
  current.mset <- my.lib$mset.list;
  set.ids<- my.lib$path.ids;
  
  current.setlink <<- current.setlink;
  current.setids <<- set.ids;
  qs::qsave(current.mset, "current.geneset.qs");

}

# these are geneset libraries
.load.enrich.lib<-function(org.code, fun.type){

  # prepare lib
  is.go <- FALSE;
  if(tolower(fun.type) == 'kegg'){ 
    nm <- "kegg";
  }else if(tolower(fun.type) == 'reactome'){ 
    nm <- "reactome";
  }else if(tolower(fun.type) == 'motif'){ 
    nm <- "motif_set";
  }else if(tolower(fun.type) == 'vset'){ 
    nm <- "disgenet_snp_disease";
  }else if(tolower(fun.type) == 'disgenet'){ 
    nm <- "disgenet";
  }else if(tolower(fun.type) == 'orphanet'){ 
    nm <- "orphanet";
  }else if(tolower(fun.type) == 'drugmatrix'){ 
    nm <- "drugmatrix";
  }else if(tolower(fun.type) == 'dsigdb'){ 
    nm <- "dsigdb";
  }else{ # GO
    is.go <- TRUE;
    nm <- paste0("go_", tolower(fun.type));
  }
  lib.nm <- paste0(nm, ".qs");
  sub.dir <- paste0("genesets/", org.code);
  if(fun.type=="vset"){
    sub.dir <- paste0("variantsets/");
  }
  my.lib <- .get.my.lib(lib.nm, sub.dir);
  if(is.go){ # fix some issue in go lib
    if(is.null(names(my.lib))){
      names(my.lib) <- c("link", "term", "sets");
    }
  }
  
  current.setlink <- my.lib$link;
  current.mset <- my.lib$sets;
  set.ids<- names(current.mset); 
  names(set.ids) <- names(current.mset) <- my.lib$term;
  
  current.setlink <<- current.setlink;
  current.setids <<- set.ids;

  qs::qsave(current.mset, "current.geneset.qs");
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param symbol.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doSymbol2EntrezMapping
#' @export 
doSymbol2EntrezMapping <- function(symbol.vec){
  db.path <- paste("../../data/hsa/entrez.rds", sep="");
  if(!.on.public.web){
    nmdb <- basename(db.path);
    download.file(db.path, destfile = nmdb, method="libcurl", mode = "wb");
    db.path <- nmdb;
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(symbol.vec, db.map[, "symbol"]);
  entrezs <- db.map[hit.inx, "gene_id"];

  # if not gene symbol, use id by itself
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- symbol.vec[na.inx];
  mode(entrezs) <- "character";
  return(entrezs);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ldclumpOpt PARAM_DESCRIPTION
#' @param ldProxies PARAM_DESCRIPTION
#' @param ldThresh PARAM_DESCRIPTION
#' @param pldSNPs PARAM_DESCRIPTION
#' @param mafThresh PARAM_DESCRIPTION
#' @param harmonizeOpt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[TwoSampleMR]{clump_data}}, \code{\link[TwoSampleMR]{extract_outcome_data}}, \code{\link[TwoSampleMR]{harmonise_data}}, \code{\link[TwoSampleMR]{mr}}, \code{\link[TwoSampleMR]{mr_heterogeneity}}, \code{\link[TwoSampleMR]{mr_pleiotropy_test}}, \code{\link[TwoSampleMR]{mr_singlesnp}}, \code{\link[TwoSampleMR]{mr_leaveoneout}}
#' @rdname PerformMRAnalysis
#' @export 
#' @importFrom TwoSampleMR clump_data extract_outcome_data harmonise_data mr mr_heterogeneity mr_pleiotropy_test mr_singlesnp mr_leaveoneout
PerformMRAnalysis <- function(ldclumpOpt, ldProxies, ldThresh, pldSNPs, mafThresh, harmonizeOpt){
  #ldclumpOpt<<-ldclumpOpt;
  #ldProxies<<-ldProxies;
  #ldThresh<<-ldThresh;
  #pldSNPs<<-pldSNPs;
  #mafThresh<<-mafThresh;
  #harmonizeOpt<<-harmonizeOpt;
  #save.image("PerformMRAnalysis.RData")
  mSetObj <- .get.mSet(mSetObj);
  
  # get instruments
  exposure.dat <- mSetObj$dataSet$exposure;
  # missing sample size, 
  exposure.dat <- exposure.dat[,c("P-value", "Chr", "SE","Beta","BP","HMDB","SNP","A1","A2","Common Name")]
  colnames(exposure.dat) <- c("pval.exposure","chr.exposure","se.exposure","beta.exposure","pos.exposure","id.exposure","SNP","effect_allele.exposure","other_allele.exposure","exposure")
  exposure.snp <- mSetObj$dataSet$exposure$SNP;
  outcome.id <- mSetObj$dataSet$outcome$id;
  # LD clumping
  if(ldclumpOpt!="no_ldclump"){
    exposure.dat <- TwoSampleMR::clump_data(exposure.dat);
  }
  # get effects of instruments on outcome
  outcome.dat <- TwoSampleMR::extract_outcome_data(snps=exposure.snp, outcomes = outcome.id, proxies = as.logical(ldProxies),
                                                   rsq = as.numeric(ldThresh), palindromes=as.numeric(as.logical(pldSNPs)), maf_threshold=as.numeric(mafThresh))
  fast.write.csv(outcome.dat, file="mr_outcome_data.csv", row.names=FALSE);
  # harmonise the exposure and outcome data
  # need to have allele frequency to harmonize
  exposure.dat$eaf.exposure <- rep("NA",nrow(exposure.dat))
  dat <- TwoSampleMR::harmonise_data(exposure.dat, outcome.dat, action = as.numeric(harmonizeOpt));
  # perform mr
  method.type <- mSetObj$dataSet$methodType;
  mr.res <- TwoSampleMR::mr(dat, method_list = method.type);
  rownames(mr.res) <- mr.res$method;
  #Analysing 'HMDB0000042' on 'ebi-a-GCST007799'
  # Heterogeneity tests
  mr_heterogeneity.res <- TwoSampleMR::mr_heterogeneity(dat);
  rownames(mr_heterogeneity.res) <- mr_heterogeneity.res$method;
  fast.write.csv(mr_heterogeneity.res, file="mr_heterogeneity_results.csv", row.names=FALSE);
  mSetObj$dataSet$mr.hetero_mat <- round(data.matrix(mr_heterogeneity.res[6:8]),3) #"Q"           "Q_df"        "Q_pval"
  
  # Test for directional horizontal pleiotropy
  mr_pleiotropy_test.res <- TwoSampleMR::mr_pleiotropy_test(dat);
  fast.write.csv(mr_pleiotropy_test.res, file="mr_pleiotropy_results.csv", row.names=FALSE);
  mSetObj$dataSet$mr.pleio_mat <- signif(data.matrix(mr_pleiotropy_test.res[5:7]),2)
  mSetObj$dataSet$mr_results <- mr.res;
  fast.write.csv(mr.res, file="mr_results.csv", row.names=FALSE);
  mSetObj$dataSet$mr_dat <- dat;
  mSetObj$dataSet$mr.res_mat <- signif(data.matrix(mr.res[6:9]), 3) #"nsnp","b","se","pval" 

  res_single <- TwoSampleMR::mr_singlesnp(dat);
  mSetObj$dataSet$mr_res_single <- res_single;
  res_loo <- TwoSampleMR::mr_leaveoneout(dat);
  mSetObj$dataSet$mr_res_loo <- res_loo;

  .set.mSet(mSetObj);
  if(.on.public.web){
    return(1);
  }else{
    return(current.msg);
  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMRRes.rowNames
#' @export 
GetMRRes.rowNames<-function(mSetObj=NA){
  #save.image("GetMRRes.rowNames.RData")
  mSetObj <- .get.mSet(mSetObj);
  nms <- rownames(mSetObj$dataSet$mr.res_mat);
  if(is.null(nms)){
    return("NA");
  }
  return (nms);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMRRes.mat
#' @export 
GetMRRes.mat<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$mr.res_mat);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetHeteroRes.rowNames
#' @export 
GetHeteroRes.rowNames<-function(mSetObj=NA){
  #save.image("GetHeteroRes.rowNames.RData")
  mSetObj <- .get.mSet(mSetObj);
  nms <- rownames(mSetObj$dataSet$mr.hetero_mat);
  if(is.null(nms)){
    return("NA");
  }
  return (nms);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetHeteroRes.mat
#' @export 
GetHeteroRes.mat<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$mr.hetero_mat);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetPleioRes.mat
#' @export 
GetPleioRes.mat<-function(mSetObj=NA){
  mSetObj <- .get.mSet(mSetObj);
  return(mSetObj$dataSet$mr.pleio_mat);
}
