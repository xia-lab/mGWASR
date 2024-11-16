#' Batch download VCF files from ieu openGWAS database
#'
#' @param IDlist The GWAS ID in the ieu database you want to download
#' @param path The path where the file is saved
#' @param timeout Download timeout setting
#'
#' @return
#' @export
#'
#' @examples download_ieuvcf(c('bbj-a-86','bbj-a-94','bbj-a-99'), './', timeout = 600)
download_ieuvcf <- function(IDlist, path, timeout=600) {
  options(timeout = timeout)
  if (!dir.exists(path)) {
    dir.create(path)
  }
  download_status <- data.frame(gwasID = IDlist, status = "Failed", stringsAsFactors = FALSE)
  failed_downloads <- c()  # 用于存储下载失败的文件ID

  for (i in 1:length(IDlist)) {
    gwasID <- IDlist[i]
    local_file <- paste(path, gwasID, '.vcf.gz', sep = "")
    local_tbi_file <- paste(path, gwasID, '.vcf.gz.tbi', sep = "")
    remote_file <- paste("https://gwas.mrcieu.ac.uk/files/", gwasID, '/', gwasID, '.vcf.gz', sep = '')
    remote_tbi_file <- paste("https://gwas.mrcieu.ac.uk/files/", gwasID, '/', gwasID, '.vcf.gz.tbi', sep = '')

    tryCatch({
      # 检查文件是否已经存在
      if (file.exists(local_file)) {
        local_file_size <- file.info(local_file)$size
        remote_file_info <- RCurl::getURL(remote_file, nobody = 1, header = 1)
      } else {
        download.file(remote_file, local_file, mode = "wb")
      }

      # 如果文件大小不匹配，重新下载
      if (!grepl(local_file_size, remote_file_info)) {
        message(paste(gwasID, "file size mismatch. Redownloading..."))
        download.file(remote_file, local_file, mode = "wb")
        download.file(remote_tbi_file, local_tbi_file, mode = "wb")
      } else {
        message(paste(gwasID, "already downloaded and file size matches. Skipping."))
      }

      # 如果.tbi文件不存在，则下载
      if (!file.exists(local_tbi_file)) {
        download.file(remote_tbi_file, local_tbi_file, mode = "wb")
      }

      # 下载成功时更新状态
      download_status$status[i] <- "Success"

    }, error = function(e) {
      # 捕获错误并将 gwasID 添加到失败列表中，同时更新状态
      message(paste("Failed to download", gwasID, ":", e$message))
      failed_downloads <<- c(failed_downloads, gwasID)
    })
  }

  # 打印所有下载失败的文件
  if (length(failed_downloads) > 0) {
    message("The following files failed to download:")
    print(failed_downloads)
  } else {
    message("All files downloaded successfully.")
  }

  # 返回下载状态数据框
  assign("download_status", download_status, envir = .GlobalEnv)
}

#' Title
#'
#' @param gwas1
#' @param gwas2
#' @param rsid_or_pos
#' @param region
#' @param plot
#' @param population
#'
#' @return
#' @export
#'
#' @examples
colocal <- function(gwas1,gwas2,rsid_or_pos,region,plot=FALSE,population){
  if (grepl("^rs[0-9]+$", rsid_or_pos)) {
    snppos <- variants_rsid(rsid_or_pos, opengwas_jwt = get_opengwas_jwt())
    pos_region = paste0(snppos$chr, ":", snppos$pos - region, "-", snppos$pos + region)
  } else if (grepl("^[0-9XYM]+:[0-9]+$", rsid_or_pos)) {
    pos_region = paste0(sub(":.*", "", rsid_or_pos), ":", sub(".*:", "", rsid_or_pos) - region, "-", sub(".*:", "", rsid_or_pos) + region)
  } else {
    stop('Please input a valid rsID or SNP position.')
  }
  if (file.exists(gwas1) & file.exists(gwas2)) {
    out <- gwasglue::gwasvcf_to_coloc(gwas1, gwas2, pos_region);
    print("Perform co-localization of vcf files")
  } else {
    out <<- gwasglue::ieugwasr_to_coloc(id1=gwas1, id2=gwas2, chrompos=pos_region)
    print("Perform co-localization of ieu GWAS ID")
  }
  res <- coloc::coloc.abf(out[[1]], out[[2]])
  if (plot==TRUE) {
    print(coloc_plot(out,gwas1name = basename(gwas1),gwas2name = basename(gwas2),population = population))
  }
  return(res)
}

#共定位画图
#' Title
#'
#' @param out
#' @param gwas1name
#' @param gwas2name
#'
#' @return
#' @export
#'
#' @examples
coloc_plot <- function(out=out,gwas1name='gwas1',gwas2name='gwas2',population = 'EUR'){
  gwas1_fn <- data.frame(rsid = as.character(out[[1]]$snp),
                        pval = as.numeric(out[[1]]$pvalues),
                        stringsAsFactors = FALSE)
  gwas2_fn <- data.frame(rsid = as.character(out[[2]]$snp),
                        pval = as.numeric(out[[2]]$pvalues),
                        stringsAsFactors = FALSE)
  locuscomparer::locuscompare(in_fn1 = gwas1_fn, in_fn2 = gwas2_fn,
                              title1 = gwas1name, title2 = gwas2name)

}

#vcf finemap
#' Title
#'
#' @param region
#' @param vcf
#' @param bfile
#' @param plink_bin
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
vcf_finemap <- function (region, vcf, bfile, plink_bin = genetics.binaRies::get_plink_binary(),
          threads = 1)
{
  message("Extracting data from vcf")
  ext <- gwasvcf::query_gwas(vcf = vcf, chrompos = region)
  out <- parallel::mclapply(unique(region), function(i) {
    message(i)
    m <- list()
    temp <- gwasvcf::query_gwas(vcf = ext, chrompos = i)
    m[["ld"]] <- ieugwasr::ld_matrix(names(temp), bfile = bfile,
                                     plink_bin = plink_bin, with_alleles = FALSE) %>%
      gwasglue:::greedy_remove()
    tib <- gwasvcf::vcf_to_tibble(temp)
    m[["z"]] <- tib[match(rownames(m[["ld"]]),tib$rsid),] %>%
      dplyr::mutate(z = ES/SE) %>% dplyr::select(snp = rsid,
                                                 zscore = z)
    m[["n"]] <- tib[["SS"]]
    return(m)  # 直接返回结果m
  }, mc.cores = threads)
  class(out) <- "FinemaprList"
  return(out)
  }

#ieufinemap
#' Title
#'
#' @param region
#' @param id
#' @param bfile
#' @param plink_bin
#'
#' @return
#' @export
#'
#' @examples
ieu_finemap <- function (region, id, bfile, plink_bin = genetics.binaRies::get_plink_binary())
{
  id <- unique(id)
  message("Getting rsids in region")
  as <- ieugwasr::associations(region, id, proxies = 0)
  rsid_avail <- unique(as$rsid)
  message("Calculating LD for ", length(rsid_avail), " variants")
  ld <- ieugwasr::ld_matrix(rsid_avail, bfile=bfile, plink_bin=plink_bin, with_alleles = FALSE) %>%
    gwasglue:::greedy_remove()
  rsid_avail <- rownames(ld)
  message("Data available for ", length(rsid_avail), " variants")
  as <- subset(as, rsid %in% rsid_avail)
  out <- list()
  for (i in 1:length(unique(id))) {
    dat <- list()
    x <- as[as[["rsid"]] %in% rsid_avail & as[["id"]] ==
              id[i], ]
    dat[["z"]] <- dplyr::tibble(snp = x[["rsid"]], zscore = x[["beta"]]/x[["se"]])
    index <- match(x[["rsid"]], rsid_avail)
    dat[["ld"]] <- ld[index, index]
    stopifnot(all(x[["rsid"]] == rownames(dat[["ld"]])))
    n <- x[["n"]]
    if (all(is.na(n))) {
      g <- ieugwasr::gwasinfo(id[i])
      n <- g[["sample_size"]]
    }
    dat[["n"]] <- n
    out[[id[i]]] <- dat
  }
  class(out) <- "FinemaprList"
  return(out)
}

#findcausalSNP: multiple region, single GWAS
#' Title
#'
#' @param regionlist
#' @param bfile
#' @param id_or_vcf 
#' @param threshold 
#'
#' @return
#' @export
#'
#' @examples
findcausalSNP <- function(regionlist,id_or_vcf,threshold = 0.8,bfile){
  finemap_df <- data.frame(matrix(ncol = 1, nrow = length(regionlist)))
  rownames(finemap_df) <- regionlist
  colnames(finemap_df) <- id_or_vcf
  # 先检查 id_or_vcf 是文件还是 ID
  vcf_exists <- file.exists(id_or_vcf)
  for (i in 1:length(regionlist)) {
    # 尝试执行代码，并捕获错误
    dat <- tryCatch({
      if (vcf_exists) {
        # 如果存在 VCF 文件，使用 vcf_finemap
        vcf_finemap(region = regionlist[i], vcf = id_or_vcf, bfile = bfile)
      } else {
        # 如果不存在 VCF 文件，使用 ieu_finemap
        ieu_finemap(region = regionlist[i], id = id_or_vcf, bfile = bfile)
      }
    }, error = function(e) {
      NULL  # 返回 NULL 以继续循环
    })
    # 如果 dat 有效，则进行 susie_rss 分析
    fitted_rss <- tryCatch({
      susieR::susie_rss(
        dat[[1]]$z$zscore,
        dat[[1]]$ld,
        mean(dat[[1]]$n)
      )
    }, error = function(e) {
      NULL  # 返回 NULL 以继续循环
    })
    merged_data <- tryCatch(unlist(strsplit(summary(fitted_rss)$cs$variable, ",")), error = function(e) NULL)
    new_finemap_snp <- tryCatch(names(which(fitted_rss$pip>threshold))[which(fitted_rss$pip>threshold) %in% merged_data], error = function(e) NULL)
    # 将找到的SNP作为单元格填入数据框
    if (length(new_finemap_snp) > 0) {
      finemap_df[i, 1] <- paste(new_finemap_snp, collapse = ", ")
    } else {
      finemap_df[i, 1] <- NA  # 若未找到SNP，填入NA
    }
    gc()
  }
  return(finemap_df)
}

#batch_coloc: single pheno1, multiple pheno2
#' Title
#'
#' @param pheno1
#' @param pheno2
#' @param region
#'
#' @return
#' @export
#'
#' @examples
batch_coloc <- function(pheno1, pheno2, rsid_or_pos, region, output_dir) {
  for (snppos in rsid_or_pos) {
  csvname <- paste0(basename(pheno1),'_',snppos, '.csv')
  csvfile <- paste0(output_dir, '/',csvname)
  if (all(grepl("^rs[0-9]+$", rsid_or_pos))) {
    snppos <- variants_rsid(rsid_or_pos, opengwas_jwt = get_opengwas_jwt())
    pos_region = paste0(snppos$chr, ":", snppos$pos - region, "-", snppos$pos + region)
  } else if (all(grepl("^[0-9XYM]+:[0-9]+$", rsid_or_pos))) {
    pos_region = paste0(sub(":.*", "", rsid_or_pos), ":", sub(".*:", "", rsid_or_pos) - region, "-", sub(".*:", "", rsid_or_pos) + region)
  } else {
    stop('Please input a valid rsID or SNP position.')
  }
  vcf_exists <- c(file.exists(pheno1), sapply(pheno2, file.exists))

  res_list <- parallel::mclapply(1:length(pheno2), function(i) {
    tryCatch({
    if (all(vcf_exists)) {
      out <- gwasglue::gwasvcf_to_coloc(vcf1 = pheno1, vcf2 = pheno2[i], chrompos = pos_region);
    } else if (!any(vcf_exists)) {
      out <- gwasglue::ieugwasr_to_coloc(id1 = pheno1, id2 = pheno2[i], chrompos = pos_region)
    }
    res <- coloc::coloc.abf(out[[1]], out[[2]])
    }, error = function(e) NA)

    if (isTRUE(is.na(res))) {
      new_res <- data.frame(pheno1 = pheno1,
                            pheno2 = pheno2[i],
                            PP.H0.abf = NA,
                            PP.H1.abf = NA,
                            PP.H2.abf = NA,
                            PP.H3.abf = NA,
                            PP.H4.abf = NA)
    } else {
      new_res <- data.frame(pheno1 = pheno1,
                            pheno2 = pheno2[i],
                            as.data.frame(t(res$summary)))
    }
    new_res
  }, mc.cores = detectCores() - 1)

  all_res <- bind_rows(res_list)
  # 保存结果到CSV文件
  write.csv(all_res, csvfile, row.names = FALSE)
  print(paste0("Saved results for ", csvname))
  }
}

#vcf_ldsc文件转换
#' Title
#'
#' @param vcf_file_path
#' @param saveRDS
#' @param output_dir
#' @param num_cores
#' @param plan
#'
#' @return
#' @export
#'
#' @examples
vcf2ldsc <- function(vcf_file_path,saveRDS=FALSE,output_dir, num_cores = 2, plan='multicore'){
  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # 设置并行计划
  plan(plan, workers = num_cores)
  # 获取vcf.gz文件路径
  vcf_files <- list.files(vcf_file_path, pattern = "\\.vcf\\.gz$", full.names = TRUE)
  # 创建进度条
  progressr::handlers(global = TRUE)
  with_progress({
    p <- progressor(steps = length(vcf_files))

    # 并行处理每个 VCF 文件的核心函数
    process_vcf <- function(vcf) {
      tryCatch({
        # 读取 VCF 文件并转换为 GRanges 对象，然后转为 tibble
        data <- readVcf(vcf) %>%
          vcf_to_granges() %>%
          dplyr::as_tibble()
        # 使用 mutate 一次性生成 Z 列
        data <- data %>%
          dplyr::mutate(Z = ES / SE) %>%
          dplyr::select(SNP = ID, N = SS, Z, A1 = ALT, A2 = REF)  # 直接选择所需的列并重命名
        name <- tools::file_path_sans_ext(basename(vcf))
        # 初始化 datalist，每个文件一个单独的列表
        datalist <- list()
        datalist[[name]] <- data

        if (saveRDS == TRUE) {
          # 保存 RDS 文件
          saveRDS(datalist, file = file.path(output_dir, paste0(name, ".rds")))
          result <- name  # 返回文件名
          rm(data)  # 删除中间数据
          gc()  # 进行垃圾回收
        } else {
          result <- datalist  # 返回子列表
        }
        # 更新进度条
        p(sprintf("Processed file: %s", basename(vcf)))
        return(result)
      }, error = function(e) {
        warning(sprintf("Error processing file: %s", basename(vcf)))
        warning("Error message:", e$message)
        return(NULL)
      })
    }

    # 使用 future.apply 包并行处理
    processed_files <- future_lapply(vcf_files, process_vcf)
  })

  processed_files <- Filter(Negate(is.null), processed_files)

  if (saveRDS == FALSE) {
    # 将所有子列表合并成一个大的列表
    result_list <- list()
    for (res in processed_files) {
      result_list <- c(result_list, res)
    }
    return(result_list)  # 返回最终的合并列表
  }

  return(processed_files)  # 返回成功处理的文件名
}



#ldsc
#' Title
#'
#' @param data
#' @param ancestry
#'
#' @return
#' @export
#'
#' @examples
uni_ldsc <- function(data,ancestry="EUR"){
  require('ldscr')
  ldscdata <-  vcf2ldsc(data)
  h2_res <- ldsc_h2(munged_sumstats = ldscdata, ancestry = ancestry)
}

#' Title
#'
#' @param data
#' @param ancestry
#'
#' @return
#' @export
#'
#' @examples
mul_ldsc <- function(data, ancestry="EUR"){
  datalist <- list()
  for (i in 1:length(data)) {
    ldscdata <-  vcf2ldsc(data[i])
    name <- tools::file_path_sans_ext(basename(data[i]))
    datalist[[name]] <- ldscdata
    print(paste0(name,' has been readed'))
  }
  rg_res <- ldsc_rg(munged_sumstats = datalist,ancestry = ancestry)
}

#' Title
#'
#' @param path
#' @param ancestry
#' @param para_plan
#' @param num_cores
#'
#' @return
#' @export
#'
#' @examples
batch_ldsc_h2 <- function(path, ancestry = "EUR", para_plan='multicore',num_cores = 3){
  # 设置并行计划
  plan(para_plan, workers = num_cores)
  progressr::handlers(global = TRUE)

  datalist <- list.files(path, full.names = TRUE)
  process_h2 <- function(data){
    tryCatch({
      name <- tools::file_path_sans_ext(basename(data))
      h2_res <- ldsc_h2(munged_sumstats = readRDS(data)[[1]], ancestry = ancestry)
      h2_res$trait <- name
      return(h2_res)
    }, error = function(e) {
      h2_res <- data.frame(
        mean_chisq = NA, lambda_gc = NA, intercept = NA, intercept_se = NA,
        ratio = NA, ratio_se = NA, h2_observed = NA, h2_observed_se = NA,
        h2_Z = NA, h2_p = NA, trait = name
      )
      return(h2_res)
    })
  }
  with_progress({
    p <- progressor(steps = length(datalist))
    results_list <- future_lapply(datalist, function(file) {
      res <- process_h2(file)
      p()  # 每次处理完一个文件后，更新进度条
      return(res)
    })
    })
  results_df <- bind_rows(results_list) %>% dplyr::select(trait, everything())
  return(results_df)
}

#' Title
#'
#' @param path
#' @param output_dir
#' @param ancestry
#' @param para_plan
#' @param num_cores
#' @param log_file
#'
#' @return
#' @export
#'
#' @examples
pair_ldsc <- function(path, output_dir, ancestry = "EUR", para_plan='multicore',num_cores = 3, log_file = "ldsc.log") {
  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # 设置并行计划
  plan(para_plan, workers = num_cores)
  # 获取所有两两配对
  data <- list.files(path, full.names = TRUE)
  pairs <- combn(length(data), 2, simplify = FALSE)

  # 获取指定路径下的所有文件名
  existing_files <- tools::file_path_sans_ext(basename(list.files(output_dir, full.names = TRUE)))
  existing_hash <- hash::hash()  # 创建一个空的哈希表
  # 将文件名存入哈希表
  for (file in existing_files) {
    existing_hash[[file]] <- TRUE  # 键为文件名，值为 TRUE
  }

  # 创建进度追踪器
  progressr::handlers(global = TRUE)
  with_progress({
    p <- progressor(steps = length(pairs))  # 进度追踪器
  # 定义处理每对的函数
  process_pair <- function(pair) {
    result <- NULL
    tryCatch({
      name1 <- tools::file_path_sans_ext(basename(data[pair[1]]))
      name2 <- tools::file_path_sans_ext(basename(data[pair[2]]))
      pair_name <- paste0(name1, '_', name2)
      file_name <- file.path(output_dir, paste0(pair_name, '.csv'))

      # 跳过已存在的配对
      if (!is.null(existing_hash[[pair_name]])) {
        p()  # 更新进度条
        return(NULL)
      }

      # 直接在函数调用时读取并传递文件
      rg_res <- ldsc_rg(munged_sumstats = c(readRDS(data[pair[1]]), readRDS(data[pair[2]])),
                        ancestry = ancestry)
      rg_res <- rg_res$rg
      # 将结果写入文件
      write.csv(rg_res, file = file_name)
      # 更新进度条
      p()
    }, error = function(e) {
      rg_res <- data.frame(
        trait1 = name1,
        trait2 = name2,
        rg = NA_real_,
        rg_se = NA_real_,
        rg_p = NA_real_
      )
      # 保存包含 NA 的 tibble 到 CSV 文件
      write.csv(rg_res, file = file_name)
    })
    return(file_name)
  }

  # 并行计算
  results_list <- future_lapply(pairs, process_pair)
  })
  return(results_list)  # 返回结果列表
}


#' Combine and Filter CSV Files
#'
#' This function reads and combines CSV files from a specified directory or a single CSV file. It allows filtering of the combined data based on specified traits and additional conditions.
#'
#' @param path A character string specifying the path to a directory containing CSV files, or a single CSV file path. 
#' @param traits A character vector of traits for filtering. Only rows where 'trait1' and 'trait2' match the specified traits will be retained.
#' @param ... conditions for filtering, including trait1, trait2, rg, rg_se, rg_p
#'
#' @return
#' @export
#'
#' @examples
combine_filter <- function(path, traits = NULL, ...) {
  require(dplyr)
  # 获取所有 CSV 文件路径
  if (file.exists(path) && grepl("\\.csv$", path)) {
    # 如果 path 是 CSV 文件路径，直接读取
    all_data <- read.csv(path, header = TRUE)
  } else {
    # 否则，获取路径下的所有 CSV 文件
    csv_files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
    # 读取 CSV 文件，并跳过空文件或有格式错误的文件
    all_data <- do.call(rbind, lapply(csv_files, function(file) {
      tryCatch({
        # 尝试读取文件
        read.csv(file, header = TRUE, row.names = 1)
      }, error = function(e) {
        message("跳过空文件或无法读取的文件: ", file)
        NULL  # 返回 NULL 以便在 rbind 时跳过该文件
      })
    }))
    # 如果没有有效数据，返回空数据框
    if (is.null(all_data)) {
      message("未读取到任何有效数据。")
      return(data.frame())  # 返回空数据框
    }
  }

  # 根据 traits 参数筛选
  if (!is.null(traits)) {
    all_data <- all_data %>%
      filter(trait1 %in% traits & trait2 %in% traits)
  }
  # 获取传入的条件
  conditions <- rlang::exprs(...)
  # 使用 dplyr 的 filter 根据条件筛选
  filtered_data <- all_data %>% filter(!!!conditions)
  return(filtered_data)
}


#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
cor_plot <- function(data){
  # 准备相关性矩阵和显著性水平矩阵
  traits <- unique(c(data$trait1, data$trait2))
  n <- length(traits)

  # 初始化空矩阵
  rg_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(traits, traits))
  display_matrix <- matrix("", nrow = n, ncol = n, dimnames = list(traits, traits))

  # 填充相关性矩阵和显示矩阵
  for (i in 1:nrow(data)) {
    rg_matrix[data$trait1[i], data$trait2[i]] <- data$rg[i]
    rg_matrix[data$trait2[i], data$trait1[i]] <- data$rg[i]

    # 下三角显示相关性值
    display_matrix[data$trait1[i], data$trait2[i]] <- sprintf("%.2f", data$rg[i])

    # 上三角显示显著性水平
    p_value <- data$rg_p[i]
    if (is.na(p_value)) {
      display_matrix[data$trait2[i], data$trait1[i]] <- ""
    } else if (p_value <= 0.001) {
      display_matrix[data$trait2[i], data$trait1[i]] <- "***"
    } else if (p_value <= 0.01) {
      display_matrix[data$trait2[i], data$trait1[i]] <- "**"
    } else if (p_value <= 0.05) {
      display_matrix[data$trait2[i], data$trait1[i]] <- "*"
    } else {
      display_matrix[data$trait2[i], data$trait1[i]] <- ""
    }
  }

  # 为对角线填充 1
  diag(rg_matrix) <- 1
  diag(display_matrix) <- ""

  # 绘制热图
  pheatmap::pheatmap(rg_matrix,
           display_numbers = display_matrix,   # 自定义显示矩阵
           color = colorRampPalette(c("#4575b4", "white", "#d73027"))(50), # Nature风格配色
           cluster_rows = FALSE,          # 禁止行聚类
           cluster_cols = FALSE,          # 禁止列聚类
           na_col = "white",              # NA 值用白色填充
           main = "Genetic Correlation Heatmap (rg)",  # 热图标题
           number_color = "black",        # 显示数字的颜色
           fontsize_number = 10)          # 数字字体大小
}
