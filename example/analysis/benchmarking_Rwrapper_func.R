require(reticulate)
require(anndata)
require(stringr)
require(aricode)
require(reshape)
require(SpatialExperiment)
require(scater)
require(harmony)
require(tidyr)
require(dplyr)
require(Seurat)
require(Banksy)
require(data.table)


# Function for creating anndata objects for
# each sample in spe object and saving as h5ad files.
create_annFiles <- function(spe, annots_label, sample_label,
                            fileSavepath, genePanel = TRUE) {

  # ####  Options explanation
  # spe = spe object containing gene expressions of
  # cells from all samples
  # annots_label = name of colData(spe) column with domain annotations
  # sample_label = name of colData(spe) column with sample names
  # fileSavepath = path where anndata files should be saved
  # genePanel = TRUE when data is from a gene panel, and not
  #   whole transcriptome (genePanel = FALSE)

  spe = logNormCounts(spe)
  if (genePanel == TRUE) {
    spe = runPCA(spe)
    reducedDim(spe, "PCA_harmony") <- harmony::RunHarmony(reducedDim(spe, "PCA"),
                                                          colData(spe), 'sample')
  } else { # for data from whole transcriptome identify top 2000 HVGs
    seu = as.Seurat(spe)
    seu = FindVariableFeatures(seu, nfeatures = 2000)
    spe = spe[VariableFeatures(seu), ]
    spe = runPCA(spe)
    reducedDim(spe, "PCA_harmony") <- harmony::RunHarmony(reducedDim(spe, "PCA"),
                                                          colData(spe), 'sample')
  }

  # generating h5ad file for each sample
  samplePaths_vec = c()
  sampleNames_vec = c()
  samples = unique(spe$sample)
  for (i in samples) {
    speX = spe[, spe$sample == i]
    sample = paste0("sample_", i)
    spCoords = as.data.frame(spatialCoords(speX))
    rownames(spCoords) = NULL
    colnames(spCoords) = NULL
    count_val = t(counts(speX))
    log_val = t(logcounts(speX))
    metaData = data.frame(domainAnnotations = as.character(colData(speX)[, annots]),
                          sampleID = rep(sample, times = nrow(spCoords)),
                          row.names = rownames(colData(speX)))
    metaData$domainAnnotations <- metaData$domainAnnotations %>% replace_na("Unknown")
    # in anndata object
    # X is a matrix
    # obs is a data frame
    # obsm and layers are lists of matrices
    ad = AnnData(
      X = t(counts(speX)),
      obs = metaData,
      obsm = list(
        spatial = as.matrix(spCoords),
        X_pca = reducedDim(speX, "PCA"),
        X_pca_harmony = reducedDim(speX, "PCA_harmony")
      ),
      layers = list(
        logcounts = t(logcounts(speX))))
    write_h5ad(ad, paste0(fileSavepath, sample, ".h5ad"))
    samplePaths_vec = append(samplePaths_vec, paste0(fileSavepath, sample, ".h5ad"))
    sampleNames_vec = append(sampleNames_vec, sample)
  }
  # return a list of vectors containing
  # anndata filepaths and sample names
  return(list("samplePaths" = samplePaths_vec,
              "sampleNames" = sampleNames_vec))
}



# Function for running SLAT
runSLAT <- function(python_path, pyscript_path,
                    samplePaths_vec, sampleNames_vec,
                    sampleIDs = 'sampleID',
                    domains = 'domainAnnotations',
                    cos = 0.3) {

  # IMPORTANT: requires python v3.8

  # create a virtual/conda environment and install
  # python v3.8 to run slat.
  # more information here: https://github.com/gao-lab/SLAT
  # and https://slat.readthedocs.io/en/latest/install.html
  #
  #### SLAT options explanation
  # python_path = path to python library in virtual environment;
  #   usually /venv/bin/python
  # pyscript_path = path to script with Slat python function
  # samplePaths_vec = vector of anndata filepaths
  # sampleNames_vec = vector of sample names in anndata file
  # The file and sample names in the samplePaths_vec and
  # sampleNames_vec vectors should be ordered from first to
  # last slice. This is important for slat multislice runs,
  # because it performs sequential alignment.
  # sampleIDs = name of anndata obs column containing sample names,
  #   the values in the column should be same as in sampleNames_vec
  # domains = name of anndata obs column containing domain
  #   annotations
  # Cos = cut-off for cosine similarity when aligning slices.
  #   By default, slat uses high cutoff of 0.85, which leads
  #   to very few alignments and likely errors during label
  #   borrowing.

  use_python(python_path)
  print(paste("Python environment set up.", Sys.time()))
  source_python(paste0(pyscript_path, "benchmark_slatFunc.py"))
  print(paste("Python script sourced.", Sys.time()))
  # Performing slat multislice runs, where every sample is
  # used as reference once.
  mega_List = list()
  for (filenum in 1:length(samplePaths_vec)){
    if (filenum == 1){
      # for first run use the default samplePaths_vec
      # as input for slat
      ad_list = samplePaths_vec
      ad_names = paste0("Run", filenum, "-", sampleNames_vec)
    } else {
      # for all other runs, bring reference to first place and
      # move the previous non-reference samples to the end of
      # sample list
      ad_list = samplePaths_vec[c(filenum:length(samplePaths_vec), 1:(filenum-1))]
      ad_names = paste0("Run", filenum, "-",
                        sampleNames_vec[c(filenum:length(samplePaths_vec), 1:(filenum-1))])
    }
    ad_list_out = slat_run(ad_list, cos)
    names(ad_list_out) = ad_names
    print(names(ad_list_out))
    mega_List = append(mega_List, ad_list_out)
    print(paste("SLAT run", filenum, "completed.",  Sys.time()))
  }

  #### clusters, ari, nmi
  clusters_df = data.frame(matrix(nrow = 0, ncol = 0))
  ari_df = data.frame(matrix(nrow = 0, ncol = length(sampleNames_vec)))
  nmi_df = data.frame(matrix(nrow = 0, ncol = length(sampleNames_vec)))
  for (sl in sampleNames_vec) {
    ari_vec = c()
    nmi_vec = c()
    run1_names = paste0("Run1", "-", sl)
    tmp_df = data.frame(X = mega_List[[run1_names]]$obsm$spatial[,1],
                        Y = mega_List[[run1_names]]$obsm$spatial[,2],
                        Cells_barcode = mega_List[[run1_names]]$obs_names,
                        Sample = mega_List[[run1_names]]$obs[[sampleIDs]],
                        Annotation = mega_List[[run1_names]]$obs[[domains]])
    for (i in 1:length(sampleNames_vec)) {
      slice_name = paste0("Run", i, "-", sl)
      tmp_df = cbind(tmp_df, cluster = mega_List[[slice_name]]$obs$clust_ref)
      ari_vec = c(ari_vec, as.numeric(ARI(mega_List[[slice_name]]$obs[[domains]],
                                          mega_List[[slice_name]]$obs$clust_ref)))
      nmi_vec = c(nmi_vec, as.numeric(NMI(mega_List[[slice_name]]$obs[[domains]],
                                          mega_List[[slice_name]]$obs$clust_ref)))
    }
    clusters_df = rbind(clusters_df, tmp_df)
    ari_df = rbind(ari_df, ari_vec)
    nmi_df = rbind(nmi_df, nmi_vec)
    print(paste0("Cluster and metrics information gathered for ", sl, ". ", Sys.time()))
  }
  colnames(clusters_df) <- c("X", "Y", "Cells_barcode", "Sample",
                                  "Annotation", paste0("SLATcluster-ref",
                                                       sampleNames_vec))
  colnames(ari_df) <- c(paste0("SLAT_ref", sampleNames_vec))
  colnames(nmi_df) <- c(paste0("SLAT_ref", sampleNames_vec))
  row.names(ari_df) <- sampleNames_vec
  row.names(nmi_df) <- sampleNames_vec
  # return a list containing data frames of clusters,
  # ari, and nmi from different slat runs, and a
  # list anndata objects from all runs
  return(list("clusters" = clusters_df,
              "ari" = ari_df,
              "nmi" = nmi_df,
              "slat_output" = mega_List))
}



# Function for running MENDER
runMENDER <- function(python_path, pyscript_path,
                      samplePaths_vec, sampleNames_vec,
                      sampleIDs = 'sampleID',
                      domains = 'domainAnnotations',
                      scale = 6, mode = 'ring', radius = 6,
                      seed = 101, batch = 'True', msm_res = 8){

  # IMPORTANT: requires python v3.10
  #
  # create a virtual/conda environment and install
  # python v3.10 to run mender.
  # more information here: https://github.com/yuanzhiyuan/MENDER
  #
  #### MENDER options explanation
  # python_path = path to python library in virtual environment;
  #   usually /venv/bin/python3.10
  # pyscript_path = path to script with Mender python function
  # samplePaths_vec = vector of anndata filepaths
  # sampleNames_vec = vector of sample names in anndata file
  # sampleIDs = name of anndata obs column containing sample names,
  #   the values in the column should be same as in sampleNames_vec
  # domains = name of anndata obs column containing domain
  #   annotations
  # Default value of scale = 6.
  # For single cell,
  #     mode = 'radius', radius = 15.
  # For spot,
  #     mode = 'ring', radius = 6 (visium) or 4 (ST).
  # nSeed is value for random_seed for reproducibility.
  # batch = True, when batch effect needs to be corrected.
  # msm_res is resolution for final leiden clustering step.
  #     positive values for expected number of domains,
  #     negative values for clustering resolution,
  #     e.g. when msm_res = 8, mender will search for appropriate
  #     resolution value that will result in 8 leiden clusters; and
  #     when msm = -0.2, mender will perform leiden with resolution 0.2.

  use_python(python_path)
  print(paste("Python environment set up.", Sys.time()))
  source_python(paste0(pyscript_path, "benchmark_menderFunc.py"))
  print(paste("Python script sourced.", Sys.time()))
  # The anndata objects should contain PCA values ('X_pca' for
  # batch = False option) and batch-corrected PCA values
  # ('X_pca_harmony' to use batch = True option).
  mender_out = mender_run(samplePaths_vec, scale, mode, radius,
                          seed, batch, msm_res)
  print(paste("MENDER run completed.",  Sys.time()))
  #### clusters, ari, nmi
  cells = paste0(str_split_i(as.vector(mender_out$obs_names), pattern = "-", i = 1),
                 "-",
                 str_split_i(as.vector(mender_out$obs_names), pattern = "-", i = 2))
  clusters_df = data.frame(X = mender_out$obsm$spatial[, 1],
                                  Y = mender_out$obsm$spatial[, 2],
                                  Cells_barcode = cells,
                                  Sample = mender_out$obs[[sampleIDs]],
                                  Annotation = mender_out$obs[[domains]],
                                  Initial_leiden = mender_out$obs$ct,
                                  MENDER_leiden = mender_out$obs$MENDER)
  metrics_df = clusters_df %>%
    group_by(Sample = as.character(clusters_df$Sample)) %>%
    summarise(ARI = ARI(Annotation, MENDER_leiden),
              NMI = NMI(Annotation, MENDER_leiden)) %>%
    as.data.frame()
  print(paste("Cluster and metrics information gathered.", Sys.time()))
  # return a list containing data frames of clusters
  # and metrics like ari and nmi, and an anndata
  # object containing mender output
  return(list("clusters" = clusters_df,
              "metrics" = metrics_df,
              "mender_output" = mender_out))
}

# Function for running BANKSY
runBANKSY <- function(spe, batch = TRUE, sample_info, annots_label = NULL,
                      sample_label, k_geom, lambda, res, npcs = 20,
                      SEED = 1000, use_agf = TRUE, compute_agf = TRUE){
  #### Banksy parameters
  # compute_agf = TRUE computes both weighted neighborhood mean (H_0) and
  #   the azimuthal Gabor filter (H_1).
  # lambda is mixing parameter, ranges from 0-1.
  #   Smaller lambda for cell-typing mode (recommended value is 0.2),
  #   Higher lambda for domain-finding mode (recommended value is 0.8).
  #   Recommended lambda = 0.2 for visium data.
  # k_geom is neighbourhood size.
  #   k_geom = 6 corresponds to first-order neighbors in visium data,
  #   whereas k_geom = 18 corresponds to first and second-order neighbors.
  #   Recommended k_geom = 18 for Visium spots, and
  #   and k_geom = c(15, 30) for other data.
  # res is Leiden clustering resolution.
  #   Higher value gives more clusters.
  #   Used res = 0.55 for visium data.
  # SEED to set.seed() for reproducibility
  # npcs is the number of PCA dimensions to calculate
  #
  #### Other parameters
  # batch = TRUE indicates batch correction is required.
  # sample_info is a data frame where first column contains sample names
  #   and second column contains patient/subject/group names.
  # annots_label is name of column in colData(spe) containing cell/spot
  #   annotations.
  # sample_label is name of column in colData(spe) conatining sample names.

  #### Preprocessing
  print(paste("Preprocessing started.", Sys.time()))
  if (is.null(annots_label) == FALSE) {
    # replacing 'NA' annotated labels with 'Unknown'
    annots = as.data.frame(as.character(spe[[annots_label]]))
    annots[, 1] = replace_na(annots[, 1], "Unknown")
    spe[[annots_label]] = factor(as.character(annots[, 1]))
    print(paste("Annotated data checked - NA replaced with Unknown.", Sys.time()))
  }
  # Removing NA spots
  # na_id <- which(is.na(spe$layer_guess_reordered_short))
  # spe <- spe[, -na_id]
  # Trimming dataset
  imgData(spe) <- NULL
  assay(spe, "logcounts") <- NULL
  reducedDims(spe) <- NULL
  rowData(spe) <- NULL
  print(paste("Trimmed spe object.", Sys.time()))
  # Grouping samples by source
  spe$subject = factor(as.character(lapply(spe[[sample_label]], function(x) {
    sample_info[sample_info[, 1] == x, 2]})))
  print(paste("Added sample group information.", Sys.time()))

  if (batch == TRUE) {
    print(paste("Multisample run with batch correction.", Sys.time()))
    colnames(spe) <- paste0(colnames(spe), "_", spe[[sample_label]])
    # Staggering spatial coordinates
    locs = spatialCoords(spe)
    locs = cbind(locs, sample = factor(spe[[sample_label]]))
    locs_dt = data.table(locs)
    colnames(locs_dt) <- c("sdimx", "sdimy", "group")
    locs_dt[, sdimx := sdimx - min(sdimx), by = group]
    global_max = max(locs_dt$sdimx) * 1.5
    locs_dt[, sdimx := sdimx + group * global_max]
    locs = as.matrix(locs_dt[, 1:2])
    rownames(locs) <- colnames(spe)
    spatialCoords(spe) <- locs
    print(paste("Spatial coordinates of samples staggered.", Sys.time()))
    # Seurat
    # Identifying HVGs
    seu = as.Seurat(spe, data = NULL)
    seu = FindVariableFeatures(seu, nfeatures = 2000)
    # Normalizing data
    scale_factor = median(colSums(assay(spe, "counts")))
    seu = NormalizeData(seu, scale.factor = scale_factor,
                        normalization.method = "RC")
    # Adding data to spe object and subsetting to HVGs
    assay(spe, "normcounts") <- GetAssayData(seu)
    spe = spe[VariableFeatures(seu), ]
    print(paste("Seurat feature selection and normalisation complete.", Sys.time()))

    #### Running BANKSY
    print(paste("BANKSY run started.", Sys.time()))
    spe <- computeBanksy(spe, assay_name = "normcounts",
                         compute_agf = compute_agf,
                         k_geom = k_geom)
    spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda,
                        npcs = npcs, seed = SEED)
    # Harmony batch correction
    PCA_label = paste0("PCA_M", as.numeric(use_agf), "_lam", lambda)
    set.seed(SEED)
    harmony_embedding = RunHarmony(data_mat = reducedDim(spe, PCA_label),
                                   meta_data = colData(spe),
                                   vars_use = c(sample_label, "subject"),
                                   max.iter = 20,
                                   verbose = FALSE)
    reducedDim(spe, "PCA_harmony") <- harmony_embedding
    print(paste("Batch correction completed.", Sys.time()))
    # Banksy clustering
    spe = clusterBanksy(spe, dimred = "PCA_harmony", use_agf = use_agf,
                         lambda = lambda, resolution = res, seed = SEED)
    print(paste("BANKSY clustering completed.", Sys.time()))
  }
  else {
    # separating samples into individual spe objects
    sample_names = unique(spe[[sample_label]])
    spe_list = lapply(sample_names, function(x) spe[, spe[[sample_label]] == x])
    # Seurat
    # Normalizing data
    seu_list = lapply(spe_list, function(x) {
      x = as.Seurat(x, data = NULL)
      NormalizeData(x, scale.factor = 5000, normalization.method = 'RC')})
    # Identifying HVGs
    hvgs = lapply(seu_list, function(x) {
      VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))})
    hvgs = Reduce(union, hvgs)
    # Adding data to spe object and subsetting to HVGs
    spe_list <- Map(function(spe, seu) {
      assay(spe, "normcounts") <- GetAssayData(seu)
      spe[hvgs,]},
      spe_list, seu_list)
    print(paste("Seurat feature selection and normalisation complete.", Sys.time()))

    #### Running BANKSY
    print(paste("BANKSY run started.", Sys.time()))
    spe_list = lapply(spe_list, computeBanksy, assay_name = "normcounts",
                      compute_agf = compute_agf, k_geom = k_geom)
    # merging samples for downstream steps
    spe <- do.call(cbind, spe_list)
    spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda,
                            group = sample_label, seed = SEED)
    spe <- clusterBanksy(spe, use_agf = use_agf, lambda = lambda,
                             resolution = res, seed = SEED)
    print(paste("BANKSY clustering completed.", Sys.time()))
  }
  clust_label = names(colData(spe))[startsWith(names(colData(spe)), "clust")]
  metrics_df = as.data.frame(colData(spe)) %>%
    group_by(!!sym(sample_label)) %>%
    summarise(ARI = ARI(layer_guess_reordered_short, !!sym(clust_label)),
              NMI = NMI(layer_guess_reordered_short, !!sym(clust_label)))
  print(paste("ARI and NMI metrics calculated.", Sys.time()))
  clusters_df = data.frame(X = spatialCoords(spe)[, 1],
                           Y = spatialCoords(spe)[, 2],
                           Sample = spe[[sample_label]],
                           Subject = spe$subject,
                           Annotation = spe[[annots_label]],
                           Banksy_cluster = spe[[clust_label]])
  return(list("clusters" = clusters_df,
              "metrics" = metrics_df,
              "banksy_output" = spe))
}
