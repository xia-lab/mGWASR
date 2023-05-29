#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @param imgName PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param width PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Cairo]{Cairo}}
#' @rdname PlotScatter
#' @export 
#' @importFrom Cairo Cairo
PlotScatter<-function(mSetObj=NA, imgName, format="png", dpi=72, width=NA){
  imgName <<-imgName;
  #save.image("PlotScatter.RData")
  
  
  mSetObj <- .get.mSet(mSetObj);
  
  # get mr_results and mr_dat
  mr.res <- mSetObj$dataSet$mr_results;
  mr.dat <- mSetObj$dataSet$mr_dat;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 12;
  }else if(width == 0){
    w <- 7;
  }else{
    w <-width;
  }
  h <- 9;
  
  #record img
  mSetObj$imgSet$mr.scatter <- imgName
  mSetObj$imgSet$current.img <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  .mr_scatterPlot(mr.res, mr.dat);
  .set.mSet(mSetObj);
  if(.on.public.web){
    return(1);
  }
}

.mr_scatterPlot <- function(mr_results, dat)
{
  # mr_results<<-mr_results;
  # dat<<-dat;
  # save.image("mr_scatterPlot.RData")
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- TwoSampleMR::mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- TwoSampleMR::mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    p <- ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position="right", legend.direction="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
    print(p);
    dev.off();
  })
  #mrres
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @param imgName PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param width PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Cairo]{Cairo}}
#' @rdname PlotForest
#' @export 
#' @importFrom Cairo Cairo
PlotForest<-function(mSetObj=NA, imgName, format="png", dpi=72, width=NA){
  #imgName <<-imgName;
  #save.image("PlotForest.RData")

  mSetObj <- .get.mSet(mSetObj);
  
  # get mr_results and mr_dat
  mr.res_single <- mSetObj$dataSet$mr_res_single;

  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 7;
  }else{
    w <-width;
  }
  h <- w;
  
  #record img
  mSetObj$imgSet$mr.scatter <- imgName
  mSetObj$imgSet$current.img <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  .mr_forestPlot(mr.res_single);
  .set.mSet(mSetObj);
  if(.on.public.web){
    return(1);
  }
}

.mr_forestPlot <- function(singlesnp_results, exponentiate=FALSE)
{
  #singlesnp_results<<-singlesnp_results;
  #save.image("mr_forestPlot.RData")
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
    levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
    am <- grep("All", d$SNP, value=TRUE)
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 0.01
    d$tot[d$SNP %in% am] <- 1
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[! d$SNP %in% am]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c(am, "", nom))
    
    xint <- 0
    if(exponentiate)
    {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }
    
    p <- ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
      # ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      ggplot2::scale_colour_manual(values=c("black", "red")) +
      ggplot2::scale_size_manual(values=c(0.3, 1)) +
      # xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
      ggplot2::theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=8), 
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=8)) +
      ggplot2::labs(y="", x=paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
    
    print(p);
    dev.off();
  })
  #res
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @param imgName PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param width PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Cairo]{Cairo}}
#' @rdname PlotLeaveOneOut
#' @export 
#' @importFrom Cairo Cairo
PlotLeaveOneOut<-function(mSetObj=NA, imgName, format="png", dpi=72, width=NA){
  #imgName <<-imgName;
  #save.image("PlotLeaveOneOut.RData")
  
  mSetObj <- .get.mSet(mSetObj);
  
  mr.res_loo <- mSetObj$dataSet$mr_res_loo;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 7;
  }else{
    w <-width;
  }
  h <- w;
  
  #record img
  mSetObj$imgSet$mr.scatter <- imgName
  mSetObj$imgSet$current.img <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  .mr_looPlot(mr.res_loo);
  .set.mSet(mSetObj);
  if(.on.public.web){
    return(1);
  }
}

.mr_looPlot <- function(leaveoneout_results)
{
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    # Need to have at least 3 SNPs because IVW etc methods can't be performed with fewer than 2 SNPs
    if(sum(!grepl("All", d$SNP)) < 3) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 1
    d$tot[d$SNP != "All"] <- 0.01
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[d$SNP != "All"]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c("All", "", nom))
    
    p <- ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=0, linetype="dotted") +
      # ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      ggplot2::scale_colour_manual(values=c("black", "red")) +
      ggplot2::scale_size_manual(values=c(0.3, 1)) +
      # xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
      ggplot2::theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=8), 
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=8)) +
      ggplot2::labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
    print(p);
    dev.off();
    
  })
  #res
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mSetObj PARAM_DESCRIPTION, Default: NA
#' @param imgName PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param width PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[Cairo]{Cairo}}
#' @rdname PlotFunnel
#' @export 
#' @importFrom Cairo Cairo
PlotFunnel<-function(mSetObj=NA, imgName, format="png", dpi=72, width=NA){
  #imgName <<-imgName;
  #save.image("PlotFunnel.RData")
  
  mSetObj <- .get.mSet(mSetObj);
  
  mr.res_single <- mSetObj$dataSet$mr_res_single;
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 12;
  }else if(width == 0){
    w <- 7;
  }else{
    w <-width;
  }
  h <- 9;
  
  #record img
  mSetObj$imgSet$mr.scatter <- imgName
  mSetObj$imgSet$current.img <- imgName;
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  .mr_funnelPlot(mr.res_single);
  .set.mSet(mSetObj);
  if(.on.public.web){
    return(1);
  }
}

.mr_funnelPlot <- function(singlesnp_results)
{
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value=TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    p <- ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(data=subset(d, SNP %in% am), ggplot2::aes(xintercept=b, colour = SNP)) +
      # ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                              "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                              "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(y=expression(1/SE[IV]), x=expression(beta[IV]), colour="MR Method") +
      ggplot2::theme(legend.position="right", legend.direction="vertical")
    print(p);
    dev.off();
  })
  #res
}
