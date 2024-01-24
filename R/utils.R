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
#' #' @examples
#' # No Example

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
#' @export
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


#' ZZZ
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
#' @export
#'
#' #' @examples
#' # No Example
mat_cor <- function(mat_x, mat_y){
  cor_mat <- matrix(0,nrow = nrow(mat_x),ncol = nrow(mat_y))
  for(i in seq_len(nrow(mat_x))){
    cor_mat[i,] <- sapply(seq_len(nrow(mat_y)),
                          function(x){
                            cor(x=mat_x[i,],y=mat_y[x,],method = "pearson")
                          })
  }
  row.names(cor_mat) <- row.names(mat_x)
  colnames(cor_mat) <- row.names(mat_y)
  cor_mat
}


