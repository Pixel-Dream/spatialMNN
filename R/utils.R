#' Preprocess
#'
#' Do column scale
#'
#' @param data_mat Expression Matrix
#' @param scale_coef scale factor
#'
#' @return Void
#'
#' @examples
#' # No Example

myPreprocess <- function(data_mat, scale_coef = 10000){
  data_mat <- data_mat/matrix(rep(colSums(data_mat), nrow(data_mat)),
                              nrow = nrow(data_mat),
                              byrow = F) * scale_coef
  data_mat
}

#' Swap
#'
#' Swap elements
#'
#' @param vec Vector with two elements
#'
#' @return Void
#'
#' @examples
#' # No Example

swap <- function(vec){
  c(vec[2],vec[1])
}


#' Convert spe to seurat list
#'
#' Convert spatialExperiment object to a list of Seurat objects based on sample id
#'
#' @param spe Spatial Experiment or Singel Cell Experiment Object
#' @param sample_id Column name contain sample IDs
#' @param sel_assay Select assay to
#' @param sel_col Columns to be included in the `meta.data`
#' @param col_name Names of the selected columns
#'
#' @return A list of Seurat Objects
#'
#' @import Seurat
#' @import SpatialExperiment
#'
#' @export
#'
#' @examples
#' No Example

spe2SeuList <- function(spe,
                       sel_assay="logcounts",
                       sample_id = "sample_id",
                       sel_col = c("layer_guess_reordered_short","spatialLIBD"),
                       col_name = c("layer","spatialLIBD")){
  seu_ls <- list()
  #hvg_ls <- c()
  #sel_assay <- "logcounts" # counts logcounts
  stopifnot(is.null(sel_col) | is.null(col_name) | length(sel_col) == length(col_name))

  if(!is.null(sel_col) & is.null(col_name))col_name = sel_col

  for(i in unique(spe[[sample_id]])){
    idx <- spe[[sample_id]] == i
    meta_df <- data.frame(barcode = spe@colData@rownames[idx],
                          coord_x = spe@int_colData@listData[["spatialCoords"]][idx,1],
                          coord_y = spe@int_colData@listData[["spatialCoords"]][idx,2],
                          row = spe$array_row[idx],
                          col = spe$array_col[idx],
                          row.names = colnames(spe@assays@data@listData[[sel_assay]])[idx])
    if(!is.null(sel_col)){
      for(j in seq_along(sel_col)){
        meta_df[[col_name[j]]] <- spe@colData[idx,sel_col[j]]
      }
    }

    seu_ls[[i]] <- SeuratObject::CreateSeuratObject(counts = spe@assays@data@listData[[sel_assay]][,idx],
                                      project = paste0("spe_",i),
                                      meta.data = meta_df)
    seu_ls[[i]] <- Seurat::FindVariableFeatures(seu_ls[[i]], verbose = F)
    #hvg_ls <- unique(c(hvg_ls,VariableFeatures(seurat_ls[[i]])))
  }
  seu_ls
}


#' Louvain with Correlation Matrix
#'
#'
#'
#' @param cor_mat_ Correlation Matrix
#' @param nn_ ?
#' @param res_ ?
#'
#' @return Clustered Data
#'
#' @import igraph
#'
#'
#' #' @examples
#' # No Example
louvain_w_cor <- function(cor_mat_, nn_=10, res_ = 1){
  for(j in seq_len(nrow(cor_mat_))){
    not_nn_vec <- sort(cor_mat_[,j],decreasing = T)[(nn_+2):ncol(cor_mat_)] %>% names()
    cor_mat_[j,not_nn_vec] <- 0
    cor_mat_[j,j] <- 0
  }
  g <- graph.adjacency(cor_mat_, mode = "directed",
                       weighted = TRUE, diag = TRUE)
  g <- as.undirected(g,mode = "mutual")

  louvain_res <- cluster_louvain(g, resolution = res_)
  louvain_res
}


#' Rasterize Single cell datasets
#'
#' Convert spatialExperiment object to a list of Seurat objects based on sample id
#'
#' @param spe Spatial Experiment or Singel Cell Experiment Object
#' @param sample_id Column name contain sample IDs
#' @param sel_assay Select assay to
#' @param sel_col Columns to be included in the `meta.data`
#' @param col_name Names of the selected columns
#'
#' @return A list of Seurat Objects
#'
#' @import Seurat
#' @import dplyr
#' @import Matrix
#' @import DescTools
#'
#'
#' @examples
#' No Example

rasterSeuratList = function(seu_ls,
                            preserveAnno = NULL,
                            tile_width = 3,
                            verbose = F) {

  pss_list = list()

  for (i in names(seu_ls)){
    if(verbose) message(paste("Rasterizing sample",i))
    xr = range(seu_ls[[i]]$coord_x)
    yr = range(seu_ls[[i]]$coord_y)
    by = tile_width

    # Add pseudo-spot info
    x_location = seq(from = xr[1] - by/2, to = xr[2] + by, by = by)
    y_location = seq(from = yr[1] - by/2, to = yr[2] + by, by = by)


    seu_ls[[i]]$psd_cell <- apply(seu_ls[[i]]@meta.data,1,
                                  function(x){
                                    sum(as.numeric(x[["coord_x"]]) > x_location) +
                                      (length(x_location)-1)*(sum(as.numeric(x[["coord_y"]]) > y_location)-1)
                                  })
    # Merge expression
    psd_cell_idx <- unique(seu_ls[[i]]$psd_cell)
    expr_mat <- seu_ls[[i]]@assays$RNA$counts
    idx_mat <- matrix(0,ncol = length(psd_cell_idx), nrow = ncol(expr_mat)) %>% as("dgCMatrix")
    for(c in 1:length(psd_cell_idx)){
      idx_mat[seu_ls[[i]]$psd_cell == psd_cell_idx[c],c] <- 1
    }
    expr_mat <- expr_mat %*% idx_mat
    colnames(expr_mat) <- paste0("s_",psd_cell_idx)
    meta_df <- data.frame(spot_id = psd_cell_idx, row.names = paste0("s_",psd_cell_idx)) %>%
      mutate(coord_x = x_location[1+((spot_id-1) %% (length(x_location)-1))] + by/2,
             coord_y = y_location[ceiling(spot_id/(length(x_location)-1))] + by/2)
    if(!is.null(preserveAnno)){
      for(anno in preserveAnno){
        if(seu_ls[[i]][[anno]] %>% unlist() %>% class() == "factor"){
          meta_df[[anno]] <- sapply(psd_cell_idx,
                                    function(idx){
                                      anno_vec <- seu_ls[[i]]@meta.data[[anno]][seu_ls[[i]]$psd_cell == as.numeric(idx)]
                                      mode_ = DescTools::Mode(anno_vec)[1] %>% as.character()
                                      if(is.na(mode_)) mode_ = as.character(anno_vec[1])
                                      mode_
                                    }) %>% as.factor()
        }else warning(paste("Discard non-factor annotation:",anno))

      }
    }
    pss_list[[i]] <- CreateSeuratObject(counts = expr_mat,
                                        meta.data = meta_df)

  }
  gc()
  pss_list
}

