#' Package Initialization
#'
#' Load Required package
#' @param verbose Show initialization information
#' @return Void
#' @examples 
#' init()
#' @export
init <- function(verbose=F){
  #suppressMessages(library(spatialLIBD))
  suppressMessages(library(Seurat))
  suppressMessages(library(ggpubr))
  suppressMessages(library(tidyverse))
  #suppressMessages(library(PRECAST))
  suppressMessages(library(igraph))
  suppressMessages(library(psych))
  suppressMessages(library(ggnewscale))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(circlize))
  #suppressMessages(library(SpotClean))
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(require(circlize))
  suppressMessages(library(reshape2))
  suppressMessages(library(graphlayouts))
  suppressMessages(require(scater))
  suppressMessages(require(scry))
  if(verbose) message("Initialize done!")
}

#' Preprocess
#'
#' Do column scale
#' @param data_mat Expression Matrix
#' @param scale_coef scale factor
#' @return Void
#' @examples 
#' # No Example
myPreprocess <- function(data_mat, scale_coef = 10000){
  data_mat <- data_mat/matrix(rep(colSums(data_mat), nrow(data_mat)), 
                              nrow = nrow(data_mat),
                              byrow = F) * scale_coef
  data_mat
}


swap <- function(vec){
  c(vec[2],vec[1])
}

#' Convert spe to seurat list
#'
#' Convert spatialExperiment object to a list of Seurat objects based on sample id
#' @param data_mat Expression Matrix
#' @param scale_coef scale factor
#' @return Void
#' @examples 
#' # No Example

new_seu_ls <- function(sel_assay="logcounts",
                       sel_col = c("layer","spatialLIBD"),
                       col_name = c("layer_guess_reordered_short","spatialLIBD")){
  seu_ls <- list()
  #hvg_ls <- c()
  #sel_assay <- "logcounts" # counts logcounts
  stopifnot(is.null(sel_col) | length(sel_col) == length(col_name))
  
  for(i in unique(spe$sample_id)){
    idx <- spe$sample_id == i
    meta_df <- data.frame(barcode = spe@colData@rownames[idx],
                          coord_x = spe@int_colData@listData[["spatialCoords"]][idx,1],
                          coord_y = spe@int_colData@listData[["spatialCoords"]][idx,2],
                          row = spe$array_row[idx],
                          col = spe$array_col[idx],
                          row.names = colnames(spe@assays@data@listData[[sel_assay]])[idx])
    if(!is.null(sel_col)){
      for(j in seq_along(sel_col)){
        meta_df[[col_name[j]]] <- spe@colData[idx,sel_col[i]]
      }
    }
    
    seu_ls[[i]] <- CreateSeuratObject(counts = spe@assays@data@listData[[sel_assay]][,idx],
                                      project = paste0("spe_",i),
                                      meta.data = meta_df)
    seu_ls[[i]] <- FindVariableFeatures(seu_ls[[i]], verbose = F)
    #hvg_ls <- unique(c(hvg_ls,VariableFeatures(seurat_ls[[i]])))
  }
  seu_ls
}


