##### PRECAST Benchmark Helper Func. #####

PRECAST_test <- function(seu_ls, k=7, gene_num=2000, sel_gene_method = "SPARK-X"){
  tic <- Sys.time()
  PRECAST_obj <- my_CreatePRECASTObject(seu_ls, project = "SpatialLIBD", gene.number = gene_num, selectGenesMethod = sel_gene_method,
                                        premin.spots = 20, premin.features = 20, postmin.spots = 1,
                                        postmin.features = 10, verbose = F)

  PRECAST_obj <- suppressMessages(AddAdjList(PRECAST_obj, platform = "Visium"))
  PRECAST_obj <- AddParSetting(PRECAST_obj, Sigma_equal = FALSE, verbose = F, int.model = NULL)

  PRECAST_obj <- PRECAST(PRECAST_obj, K = k)
  PRECAST_obj <- SelectModel(PRECAST_obj)
  seuInt <- IntegrateSpaData(PRECAST_obj, species = "Human")
  toc <- Sys.time()
  message(toc - tic)
  seuInt
}

draw_PRECAST <- function(seuInt_obj, k=7){
  if(k <= 20){
    cols_cluster <- chooseColors(palettes_name = "Classic 20", n_colors = k, plot_colors = TRUE)
  }else{
    cols_cluster <- chooseColors(palettes_name = "hue n", n_colors = k, plot_colors = TRUE)
  }

  pList <- SpaPlot(seuInt_obj, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = FALSE,
                   nrow.legend = 7)
  drawFigs(pList, layout.dim = c(2, 2), common.legend = TRUE, legend.position = "right", align = "hv")
}

# Record results
rkd_res <- function(bench_df, tic, toc, mem_usage, method_name,
                    seu_ls=NULL, label_col=NULL, result_col=NULL,
                    label_vec=NULL, result_vec=NULL, sample_vec = NULL){
  idx <- which(bench_df$method == method_name)
  bench_df$tic[idx] <- tic
  bench_df$toc[idx] <- toc
  bench_df$time[idx] <- toc - tic
  bench_df$mem[idx] <- mem_usage$Peak_RAM_Used_MiB
  if(!is.null(seu_ls) & !is.null(label_col) & !is.null(result_col)){
    ari_vec <- sapply(seu_ls,function(seu_obj) mclust::adjustedRandIndex(seu_obj@meta.data[[label_col]], seu_obj@meta.data[[result_col]]))
    nmi_vec <- sapply(seu_ls,
                      function(seu_obj){
                        v1 <- seu_obj@meta.data[[label_col]]
                        v2 <- seu_obj@meta.data[[result_col]]
                        aricode::NMI(v1[!is.na(v1)&!is.na(v2)],v2[!is.na(v1)&!is.na(v2)])
                      })
  }else if(!is.null(sample_vec) & !is.null(label_vec) & !is.null(result_vec)){
    sample_vec <- as.factor(sample_vec)
    ari_vec <- sapply(levels(sample_vec),function(i) mclust::adjustedRandIndex(label_vec[sample_vec==i], result_vec[sample_vec==i]))
    nmi_vec <- sapply(levels(sample_vec),
                      function(i){
                        v1 <- label_vec[sample_vec==i]
                        v2 <- result_vec[sample_vec==i]
                        aricode::NMI(v1[!is.na(v1)&!is.na(v2)],v2[!is.na(v1)&!is.na(v2)])
                      })
  }else stop("No clustering results and ground truth labels for ARI/NMI calculation")

  bench_df$ari[idx] <- ari_vec %>% mean()
  bench_df$ari_vec[idx] <- ari_vec %>% paste(collapse = ",")
  bench_df$nmi[idx] <- nmi_vec %>% mean()
  bench_df$nmi_vec[idx] <- nmi_vec %>% paste(collapse = ",")
  bench_df
}


draw_slide_bench <- function(x,y,label_vec,
                             col_name = NULL, sample_id="",
                             pal = NULL,
                             flip = F){
  df_ <- data.frame(x=x,y=y,label=label_vec)
  if(!is.null(col_name))colnames(df_)[3] <- col_name
  else title = "label"
  g <- ggplot()
  g <- g + geom_point(mapping = aes(x = df_[["x"]], y = df_[["y"]], color=df_[[col_name]]),
                      data = df_) +
    labs(x="x",y="y",title = paste(sample_id, col_name))
  if(is.numeric(df_[[col_name]])){
    require(viridis)
    g <- g + scale_colour_gradientn(colours = viridis::viridis(n=10))
  }else{
    if(is.null(pal)) g <- g + scale_color_discrete()
    else g <- g + scale_colour_manual(values = pal)
  }
  g <- g + theme_classic()
  g <- g + guides(color=guide_legend(title=col_name))
  if(flip) g + coord_flip()
  else g
}


# Fix some bugs in PRECAST

my_CreatePRECASTObject <- function (seu_ls, project = "PRECAST", gene.number = 800,
                selectGenesMethod = "SPARK-X", numCores_sparkx = 1, customGenelist = NULL,
                premin.spots = 20, premin.features = 20, postmin.spots = 15,
                postmin.features = 15, rawData.preserve = FALSE, verbose = TRUE)
{
  if (!inherits(seu_ls, "list"))
    stop("CreatePRECASTObject: check the argument: seu_ls! it must be a list.")
  flag <- sapply(seu_ls, function(x) !inherits(x, "Seurat"))
  if (any(flag))
    stop("CreatePRECASTObject: check the argument: seu_ls! Each component of seu_ls must be a Seurat object.")
  exist_spatial_coods <- function(seu) {
    flag_spatial <- all(c("row", "col") %in% colnames(seu@meta.data))
    return(flag_spatial)
  }
  flag_spa <- sapply(seu_ls, exist_spatial_coods)
  if (any(!flag_spa))
    stop("CreatePRECASTObject: check the argument: seu_ls! Each Seurat object in seu_ls must include  the spatial coordinates saved in the meta.data, named 'row' and 'col'!")
  if (numCores_sparkx < 0)
    stop("CreatePRECASTObject: check the argument: numCores_sparkx! It must be a positive integer.")
  if (!is.null(customGenelist) && (!is.character(customGenelist)))
    stop("CreatePRECASTObject: check the argument: customGenelist! It must be NULL or a character vector.")
  object <- new(Class = "PRECASTObj", seuList = seu_ls, seulist = NULL,
                AdjList = NULL, parameterList = list(), resList = NULL,
                project = project)
  seu_ls <- object@seuList
  if (verbose)
    message("Filter spots and features from Raw count data...")
  tstart <- Sys.time()
  seu_ls <- lapply(seu_ls, PRECAST:::filter_spot, premin.features)
  seu_ls <- pbapply::pblapply(seu_ls, PRECAST:::filter_gene, premin.spots)
  message(" \n ")
  PRECAST:::.logDiffTime(sprintf(paste0("%s Filtering step for raw count data finished!"),
                       "*****"), t1 = tstart, verbose = verbose)
  if (verbose)
    message("Select the variable genes for each data batch...")
  tstart <- Sys.time()
  if (is.null(customGenelist)) {
    if (tolower(selectGenesMethod) == "spark-x") {
      seu_ls <- pbapply::pblapply(seu_ls, PRECAST:::.findSVGs,
                                   nfeatures = gene.number, num_core = numCores_sparkx,
                                   verbose = verbose)
      spaFeatureList <- lapply(seu_ls, PRECAST:::.topSVGs, ntop = gene.number)
    }
    else if (tolower(selectGenesMethod) == "hvgs") {
      seu_ls <- pbapply::pblapply(seu_ls, FindVariableFeatures,
                                   nfeatures = gene.number, verbose = verbose)
      getHVGs <- function(seu) {
        assay <- DefaultAssay(seu)
        if (!inherits(seu[[assay]], "Assay5")) {
          vf <- seu[[assay]]@var.features
        }
        else {
          vf_dat <- seu[[assay]]@meta.data[, c("vf_vst_counts_rank",
                                               "var.features")]
          vf_dat <- vf_dat[complete.cases(vf_dat), ]
          vf <- vf_dat$var.features[order(vf_dat$vf_vst_counts_rank)]
        }
        return(vf)
      }
      spaFeatureList <- lapply(seu_ls, getHVGs)
    }
    else {
      stop("CreatePRECASTObject: check the argument: selectGenesMethod! It only support 'SPARK-X' and 'HVGs' to select genes now. You can provide self-selected genes using customGenelist argument.")
    }
    spaFeatureList <- lapply(spaFeatureList, function(x) x[!is.na(x)])
    if (any(sapply(spaFeatureList, length) < gene.number)) {
      gene.number_old <- gene.number
      gene.number <- min(sapply(spaFeatureList, length))
      warning(paste0("Number of genes in one of sample is less than ",
                     gene.number_old, ", so set minimum number of variable genes as gene.number=",
                     gene.number))
    }
    if (verbose)
      message("Select common top variable genes  for multiple samples...")
    genelist <- my_selectIntFeatures(seu_ls, spaFeatureList = spaFeatureList,
                                  IntFeatures = gene.number)
  }
  else {
    geneNames <- Reduce(intersect, (lapply(seu_ls, row.names)))
    if (any(!(customGenelist %in% geneNames)))
      message("CreatePRECASTObject: remove genes:", paste0(setdiff(customGenelist,
                                                                   geneNames), "  "), "with low count reads in seu_ls.")
    genelist <- intersect(customGenelist, geneNames)
  }
  PRECAST:::.logDiffTime(sprintf(paste0("%s Gene selection finished!"),
                       "*****"), t1 = tstart, verbose = verbose)
  seu_ls <- lapply(seu_ls, function(x) x[genelist, ])
  if (verbose)
    message("Filter spots and features from SVGs(HVGs) count data...")
  tstart <- Sys.time()
  seu_ls <- lapply(seu_ls, PRECAST:::filter_spot, postmin.features)
  seu_ls <- pbapply::pblapply(seu_ls, PRECAST:::filter_gene, postmin.spots)
  seu_ls <- lapply(seu_ls, NormalizeData, verbose = verbose)
  object@seulist <- seu_ls
  PRECAST:::.logDiffTime(sprintf(paste0("%s Filtering step for count data with variable genes finished!"),
                       "*****"), t1 = tstart, verbose = verbose)
  if (!rawData.preserve) {
    object@seuList <- NULL
  }
  return(object)
}




my_selectIntFeatures <- function (seu_ls, spaFeatureList, IntFeatures = 2000)
{
  if (length(seu_ls) != length(spaFeatureList))
    stop("The length of suelist and spaFeatureList must be equal!")
  if (length(seu_ls) == 1) {
    if (length(spaFeatureList[[1]]) >= IntFeatures) {
      genelist <- spaFeatureList[[1]][1:IntFeatures]
    }
    else {
      genelist <- spaFeatureList[[1]]
      warning("The IntFeature is larger than the  number of elements in FeatureList!")
    }
    return(genelist)
  }
  if (any(sapply(spaFeatureList, length) < IntFeatures))
    stop("Feature list exists number of features less than IntFeatures!")
  geneUnion <- unique(unlist(lapply(spaFeatureList, function(x) x[1:IntFeatures])))
  gene_delete <- unique(unlist(lapply(seu_ls, function(x) setdiff(geneUnion,
                                                                   row.names(x)))))
  geneUnion <- setdiff(geneUnion, gene_delete)
  genes_zeroVar <- unique(unlist(lapply(seu_ls, function(x) {
    assay <- DefaultAssay(x)
    count_matrix <- GetAssayData(x, assay = assay, slot = "counts")
    geneUnion[Matrix::rowSums(count_matrix[geneUnion, ]) ==
                0]
  })))
  gene_Var <- setdiff(geneUnion, genes_zeroVar)
  nsample <- length(seu_ls)
  numVec <- rep(0, length(gene_Var))
  rankMat <- matrix(NA, length(gene_Var), nsample)
  row.names(rankMat) <- gene_Var
  for (i in 1:length(gene_Var)) {
    for (j in 1:nsample) {
      if (is.element(gene_Var[i], spaFeatureList[[j]])) {
        numVec[i] <- numVec[i] + 1
        rank1 <- which(spaFeatureList[[j]] == gene_Var[i])
        rankMat[i, j] <- rank1
      }
    }
  }
  cutNum <- sort(numVec, decreasing = T)[min(IntFeatures,
                                             length(numVec))]
  if (max(numVec) > cutNum) {
    genelist1 <- gene_Var[numVec > cutNum]
  }else {
    genelist1 <- NULL
  }
  num_rest_genes <- min(IntFeatures, length(numVec)) - length(genelist1)
  gene2 <- gene_Var[numVec == cutNum]
  # Add case that length(gene2) <=1
  if(length(gene2) <= 1){
    genelist <- c(genelist1, gene2)
  }else{
    rankMat2 <- rankMat[gene2, ]
    rowMedian <- function(xmat, na.rm = TRUE) {
      apply(xmat, 1, median, na.rm = na.rm)
    }
    genes1rank <- gene2[order(rowMedian(rankMat2, na.rm = T))[1:num_rest_genes]]
    genelist <- c(genelist1, genes1rank)
  }

  return(genelist)
}





