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


# Function for creating anndata objects for
# each sample in spe object and saving as h5ad files.
create_annFiles <- function(spe, annots, filepath) {

  # ####  Options explanation
  # spe = spe object containing gene expressions of
  # cells from all samples
  # annots = name of colData column with domain annotations
  # filepath = path where anndata files should be saved to

  spe = logNormCounts(spe)
  spe = runPCA(spe)
  reducedDim(spe, "PCA_harmony") <- harmony::RunHarmony(reducedDim(spe, "PCA"),
                                                        colData(spe), 'sample')
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
    write_h5ad(ad, paste0(filepath, sample, ".h5ad"))
    samplePaths_vec = append(samplePaths_vec, paste0(filepath, sample, ".h5ad"))
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
