
#' Core function for Niche Identification
#'
#' @param seu_obj A  Seurat Objects, columns `coord_x,coord_y` must be included in the `meta.data`
#' @param cor_threshold Threshold for edge pruning
#' @param nn Number of NN for graph building
#' @param nn_2 For SNN correlation only
#' @param cl_resolution Louvain resolution for clustering
#' @param top_pcs Number of selected top PCs
#' @param cl_min Minimal nodes in identified niches
#' @param find_HVG Whether find HVG for downstream analysis or use full gene list
#' @param hvg HVG numbers, default: `2000`
#' @param cor_met Correlation calculation method, available methods:`("PC","HVG")`
#' @param edge_smoothing Perform smoothed edge detection
#' @param use_glmpca Use GLMPCA or regular PCA
#' @param verbose Output clustering information
#'
#' @return A Seurat Objects
#'
#' @import Seurat
#' @import scater
#' @import scry
#' @import parallel
#' @import ggnewscale
#' @import igraph
#' @import DescTools
#' @import HiClimR
#' @importFrom dbscan kNN
#'
#'
#' @examples
#' TBD

stg1_func <- function(seu_obj, cor_threshold = 0.2, nn = 12, nn_2=20, cl_resolution = 10,
                       top_pcs = 30, cl_min=5, find_HVG = T, hvg = 2000, cor_met = "PC",
                       edge_smoothing = T, use_glmpca = T, verbose = F){
  obj_name <- seu_obj@meta.data[["orig.ident"]][1]
  if(find_HVG){
    seu_obj <- NormalizeData(seu_obj)
    stopifnot("HVG# Exceeded" = hvg <= nrow(seu_obj))
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = hvg, verbose = F)
  }else{
    VariableFeatures(object = seu_obj) <- row.names(seu_obj)
    hvg <- nrow(seu_obj)
  }

  if(use_glmpca == T){
    mat <- as.matrix(seu_obj@assays$RNA$counts[VariableFeatures(object = seu_obj),]) # CH: changed @ to $ before "counts"
    mat <- nullResiduals(mat, type="deviance")
    PCA_res <- suppressWarnings(calculatePCA(mat,ncomponents = top_pcs, scale = TRUE, BSPARAM = BiocSingular::RandomParam()))
    seu_obj@misc[["glmpca"]] <- PCA_res %>% as.matrix()
    #seu_obj@reductions[["pca"]]@feature.loadings <- res1$loadings %>% as.matrix()
    rm(mat)
    gc()
  }else{
    seu_obj <- ScaleData(seu_obj, features = rownames(seu_obj),verbose = F)
    seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj),verbose = F)
    PCA_res <- seu_obj@reductions[["pca"]]@cell.embeddings
  }

  dist_knn <- dbscan::kNN(seu_obj@meta.data %>% select(coord_x,coord_y),k=nn)
  if(cor_met == "SNN") expr_knn <- dbscan::kNN(PCA_res[,1:top_pcs],k=nn_2)
  tic<-Sys.time()
  nn_df <- data.frame(i = pmin(rep(1:nrow(dist_knn[["id"]]), nn),as.numeric(dist_knn[["id"]])),
                      j = pmax(rep(1:nrow(dist_knn[["id"]]), nn),as.numeric(dist_knn[["id"]])),
                      x = rep(0,nrow(dist_knn[["id"]]) * nn)) %>% distinct()

  if(!edge_smoothing){
    nn_df$x <- apply(nn_df, 1,
                     FUN = function(x){
                       i=as.numeric(x[["i"]]);j=as.numeric(x[["j"]]);
                       seq_vec <- seq_len(ncol(seu_obj))
                       if(cor_met == "SNN"){
                         length(intersect(unlist(expr_knn$id[j,]),
                                          unlist(expr_knn$id[k,])))
                       }
                       else{
                         if(cor_met == "PC"){
                           cor(x=PCA_res[i,1:top_pcs],x=PCA_res[j,1:top_pcs],method = "pearson")
                         }
                         else if(cor_met == "HVG"){
                           cor_mat <- cor(x=seu_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_obj),i],
                                          y=seu_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_obj),j],method = "pearson")
                         }
                       }
                     })
  }else{# Enable edge_smoothing
    nn_df$x <- apply(nn_df, 1,
                     FUN = function(x){
                       i=as.numeric(x[["i"]]);j=as.numeric(x[["j"]]);
                       nn_vec1 <- unlist(dist_knn$id[i,])
                       nn_vec2 <- unlist(dist_knn$id[j,])

                       nn_common <- c(i,j,intersect(nn_vec1, nn_vec2))

                       nn_vec1 <- c(i,nn_vec1[!nn_vec1 %in% nn_common])
                       nn_vec2 <- c(j,nn_vec2[!nn_vec2 %in% nn_common])

                       if(cor_met == "PC"){
                         cor_val <- cor(x = colMeans(matrix(PCA_res[nn_vec1,1:top_pcs], ncol = top_pcs)),
                                        y = colMeans(matrix(PCA_res[nn_vec2,1:top_pcs], ncol = top_pcs)),
                                        method = "pearson")
                       }else if(cor_met == "HVG") {
                         cor_mat <- cor(seu_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_obj),],method = "pearson")
                         cor_val <- cor(x = rowMeans(matrix(seu_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_obj),nn_vec1], nrow = hvg)),
                                        y = rowMeans(matrix(seu_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_obj),nn_vec2], nrow = hvg)),
                                        method = "pearson")
                       }else if(cor_met == "SNN"){
                         cor_val <- map2(.x = rep(nn_vec1,length(nn_vec2)),
                                         .y = rep(nn_vec2,each = length(nn_vec1)),
                                         .f = function(x,y){
                                           length(intersect(unlist(expr_knn$id[x,]),
                                                            unlist(expr_knn$id[y,])))
                                         }) %>% unlist() %>% median()
                       }
                     })
  }

  nn_df <- nn_df[nn_df$x > cor_threshold,]

  cor_mat <- Matrix::sparseMatrix(i = nn_df$i,
                                  j = nn_df$j,
                                  x = nn_df$x,
                                  dims = c(ncol(seu_obj),ncol(seu_obj)),
                                  dimnames = list(colnames(seu_obj), colnames(seu_obj)))

  g <- graph.adjacency(cor_mat, mode = "directed",
                       weighted = TRUE, diag = TRUE)
  g <- as.undirected(g,mode = "collapse")

  seu_obj@misc[["edges"]] <- igraph::as_data_frame(g,"edges")

  cl_res <- cluster_louvain(g, resolution = cl_resolution)
  seu_obj@meta.data[["cluster"]] <- as.factor(cl_res$membership)

  if(verbose) message(paste(obj_name,"cluster number:",length(levels(as.factor(cl_res$membership)))))
  seu_obj@misc[["graph_plot_cluster"]] <-
    draw_slide_graph(seu_obj@meta.data,seu_obj@misc[["edges"]],
                     cor_threshold,"cluster") + theme(legend.position = "none")

  # Cluster merging
  seu_obj@meta.data[["merged_cluster"]] <- seu_obj@meta.data[["cluster"]]
  while(sum(table(seu_obj@meta.data[["merged_cluster"]]) < cl_min) > 0){
    sml_cl_idx <- names(table(seu_obj@meta.data[["merged_cluster"]]))[table(seu_obj@meta.data[["merged_cluster"]]) < cl_min]

    for(j in sml_cl_idx){
      node_ls <- which(seu_obj@meta.data[["merged_cluster"]]==j)
      nn_ls <- lapply(node_ls,function(x){unlist(dist_knn$id[x,])}) %>% unlist %>% unique
      nn_ls <- nn_ls[!nn_ls %in% node_ls]
      nn_cl <- seu_obj@meta.data[["merged_cluster"]][nn_ls]
      seu_obj@meta.data[["merged_cluster"]][node_ls] <- names(sort(table(nn_cl),decreasing=TRUE)[1])
    }
    seu_obj@meta.data[["merged_cluster"]] <- droplevels(seu_obj@meta.data[["merged_cluster"]])
    if(verbose) message(paste(obj_name,"merged cluster number:",length(levels(seu_obj@meta.data[["merged_cluster"]]))))
  }


  seu_obj@misc[["graph_plot_cluster_merged"]] <-
    draw_slide_graph(seu_obj@meta.data,seu_obj@misc[["edges"]],
                     cor_threshold,"merged_cluster") + theme(legend.position = "none")
  seu_obj

}




#' Niche Identification
#' Main function ??
#'
#' @param seu_ls A list of Seurat Objects, columns `coord_x,coord_y` must be included in the `meta.data`
#' @param cor_threshold Threshold for edge pruning
#' @param nn Number of NN for graph building
#' @param nn_2 For SNN correlation only
#' @param cl_resolution Louvain resolution for clustering
#' @param top_pcs Number of selected top PCs
#' @param cl_min Minimal nodes in identified niches
#' @param find_HVG Whether find HVG for downstream analysis or use full gene list
#' @param hvg HVG numbers, default: `2000`
#' @param cor_met Correlation calculation method, available methods:`("PC","HVG")`
#' @param edge_smoothing Perform smoothed edge detection
#' @param use_glmpca Use GLMPCA or regular PCA
#' @param verbose Output clustering information
#'
#' @return A list of Seurat Objects
#'
#' @import Seurat
#' @import scater
#' @import scry
#' @import parallel
#' @import ggnewscale
#' @import igraph
#' @importFrom dbscan kNN
#'
#' @export
#'
#' @examples
#' TBD
stage_1 <- function(seu_ls, cor_threshold = 0.2, nn = 12, nn_2=20, cl_resolution = 10,
                    top_pcs = 30, cl_min=5, find_HVG = T, hvg = 2000, cor_met = "PC",
                    edge_smoothing = T, use_glmpca = T, num_core = 1, verbose = F){
  stopifnot("Not supported Cor Calc Method" = cor_met %in% c("PC","HVG","SNN"))



  tic <- Sys.time()

  if(num_core>1){
    cl <- makeCluster(num_core)
    pkg_ls <- clusterEvalQ(cl, {
      library(atlasClustering)
      library(ggplot2)
      library(Seurat)
      library(scater)
      library(scry)
      library(dbscan)
      library(dplyr)
      library(igraph)
      library(ggnewscale)
    })
    clusterExport(cl, c("cor_threshold","nn","nn_2","cl_resolution","top_pcs",
                        "cl_min","find_HVG","hvg","cor_met","edge_smoothing",
                        "use_glmpca", "verbose","stg1_func"), envir=environment())

    num_batch = ceiling(length(seu_ls)/num_core)
    for(i in 1:num_batch){
      exec_seq <- (1+(i-1)*num_core):(min(i*num_core,length(seu_ls)))
      for(j in seq_along(exec_seq)){
        seu_obj <- seu_ls[[exec_seq[j]]]
        clusterExport(cl[j],"seu_obj", envir=environment())
      }
      seu_ls[exec_seq] <- clusterEvalQ(cl[1:length(exec_seq)],{
        stg1_func(seu_obj, cor_threshold, nn, nn_2, cl_resolution,
                  top_pcs, cl_min, find_HVG, hvg, cor_met,
                  edge_smoothing, use_glmpca, verbose)
      })
    }

    stopCluster(cl)
  }else{
    for(i in names(seu_ls)){
      seu_ls[[i]] <- stg1_func(seu_ls[[i]], cor_threshold, nn, nn_2, cl_resolution,
                top_pcs, cl_min, find_HVG, hvg, cor_met,
                edge_smoothing, use_glmpca, verbose)
    }
  }

  toc <- Sys.time()
  if(verbose){
    message("Elasped Time:")
    message(format(difftime(toc, tic, units='secs'),digits = 3))
  }

  seu_ls
}



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

#' Integrative Clustering
#'
#' Main function 2
#' @param seu_ls A list of Seurat Objects from `stage_1` output
#' @param nn_2 MNN parameter
#' @param top_pcs Number of selected top PCs
#' @param cl_key Column contains primary cluster labels
#' @param rtn_seurat Return a Seurat object with averaged expression profile and integrative clustering results
#' @param method Secondary Clustering algorithm, available methods:`("Louvain","MNN")`
#' @param hly_cor Threshold for highly correlated clusters
#' @param find_HVG Whether find HVG for downstream analysis or use full gene list
#' @param hvg HVG number for secondary clustering
#' @param cor_met Correlation calculation method, available methods:`("PC","HVG")`
#' @param use_glmpca Use GLMPCA or regular PCA
#' @param resolution Secondary clustering resolution
#' @param rare_ct Rare cell type detection method, available methods:`("a","m","none")`
#' @param verbose Output clustering information
#'
#' @return A list of Seurat Objects
#'
#' @import scater
#' @import scry
#' @import HiClimR
#' @import igraph
#' @export
#'
#' @examples
#' TBD
stage_2 <- function(seu_ls, top_pcs = 30, nn_2=10, cl_key = "merged_cluster",
                    rtn_seurat = F, method="louvain", hly_cor = 0.9, find_HVG = T, hvg = 2000, cor_met = "PC",
                    use_glmpca = F, resolution = 1,
                    rare_ct="none", verbose = F){
  tic <- Sys.time()
  hvg_union <- c()
  gene_intersect <- row.names(seu_ls[[names(seu_ls)[1]]])
  cl_num <- c()

  stopifnot("Not supported Cor Calc Method" = cor_met %in% c("PC","HVG"))

  for(i in names(seu_ls)){
    hvg_union <- union(hvg_union,VariableFeatures(object = seu_ls[[i]]))
    gene_intersect <- intersect(gene_intersect, row.names(seu_ls[[i]]))
    cl_num <- append(cl_num,
                     c(length(levels(droplevels(as.factor(seu_ls[[i]]@meta.data[[cl_key]]))))))
  }

  hvg_union <- intersect(hvg_union, gene_intersect)

  cl_expr <- matrix(0, nrow = sum(cl_num),ncol = length(hvg_union))
  cl_df <- data.frame(sample=rep(names(seu_ls),cl_num),cluster=0)

  for(i in names(seu_ls)){
    cl_df[["cluster"]][cl_df$sample==i] <- levels(droplevels(as.factor(seu_ls[[i]]@meta.data[[cl_key]])))
  }
  row.names(cl_expr) <- paste0("X",cl_df$sample,"_",cl_df$cluster)
  colnames(cl_expr) <- hvg_union

  for(i in seq_len(nrow(cl_df))){
    #if(i%%1000==0)message(i)
    idx <- which(seu_ls[[cl_df$sample[i]]]@meta.data[[cl_key]] == cl_df$cluster[i])
    if(length(idx)>1) cl_expr[i,] <-  seu_ls[[cl_df$sample[i]]]@assays[["RNA"]]$counts[hvg_union,idx] %>% # CH: change @data to counts
        as.matrix(.) %>% rowSums(.)/length(idx)
    else cl_expr[i,] <-  seu_ls[[cl_df$sample[i]]]@assays[["RNA"]]$counts[hvg_union,idx] # CH: change @data to counts
  }

  cl_expr_obj <- CreateSeuratObject(t(cl_expr), meta.data = data.frame(row.names = row.names(cl_expr), cl_df))

  if(find_HVG){
    cl_expr_obj <- NormalizeData(cl_expr_obj)
    stopifnot("HVG# Exceeded" = hvg <= nrow(cl_expr_obj))
    cl_expr_obj <- FindVariableFeatures(cl_expr_obj, selection.method = "vst", nfeatures = hvg, verbose = F)
  }else{
    VariableFeatures(object = cl_expr_obj) <- row.names(cl_expr_obj)
  }

  if(use_glmpca){
    mat <- as.matrix(cl_expr_obj@assays$RNA$counts[VariableFeatures(object = cl_expr_obj),])
    mat <- scry::nullResiduals(mat, type="deviance")
    PCA_res <- scater::calculatePCA(mat,ncomponents = top_pcs, scale = TRUE, BSPARAM = BiocSingular::RandomParam())
    #cl_expr_obj@reductions[["pca"]]@cell.embeddings <- res1 %>% as.matrix()
    cl_expr_obj @misc[["glmpca"]] <- PCA_res %>% as.matrix()
    #cl_expr_obj@reductions[["pca"]]@feature.loadings <- gp_res$loadings %>% as.matrix()
    rm(mat)
    gc()
  }else{
    cl_expr_obj <- ScaleData(cl_expr_obj, features = rownames(cl_expr_obj),verbose = F)
    cl_expr_obj <- RunPCA(cl_expr_obj, features = VariableFeatures(object = cl_expr_obj),verbose = F)
    cl_expr_obj <- RunUMAP(cl_expr_obj, dims = 1:top_pcs,verbose = F)
    PCA_res <- cl_expr_obj@reductions[["pca"]]@cell.embeddings
  }

  if(method == "louvain"){
    if(cor_met == "PC"){
      cor_mat <- HiClimR::fastCor(t(as.matrix(PCA_res[,1:top_pcs])), upperTri = F)
    }else if(cor_met == "HVG"){
      cor_mat <- HiClimR::fastCor(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj),], upperTri = F)
    }

    louvain_res <- louvain_w_cor(cor_mat)
    #table(louvain_res$membership)
    cl_df[["louvain"]] <- as.character(louvain_res$membership)
    sml_cl_idx <- names(table(cl_df[["louvain"]]))[table(cl_df[["louvain"]]) < 2]
    cl_df[["louvain"]][cl_df[["louvain"]] %in% sml_cl_idx] <- NA
  }else if(method == "MNN"){
    g <- igraph::make_empty_graph(directed = F)
    node_df <- data.frame(ID=seq_along(colnames(cl_expr_obj)),
                          cl = colnames(cl_expr_obj),
                          sample = cl_expr_obj$orig.ident,
                          row.names = colnames(cl_expr_obj))
    g <- igraph::add_vertices(g, nrow(node_df))
    igraph::V(g)$cl <- node_df$cl
    igraph::V(g)$sample <- node_df$sample

    all_samples <- unique(cl_expr_obj@meta.data[["sample"]])
    # Enable MNN within sample
    if(rare_ct == "m"){
      mnn_seq <- seq_along(all_samples)
    }else{
      mnn_seq <- seq_along(all_samples)[-1]
    }

    for(i in mnn_seq){
      if(rare_ct == "m"){
        mnn_seq2 <- 1:i
      }else{
        mnn_seq2 <- 1:(i-1)
      }
      for(j in mnn_seq2){
        #message(paste(i,j))
        idx_i <- cl_expr_obj@meta.data[["sample"]] == all_samples[i]

        idx_j <- cl_expr_obj@meta.data[["sample"]] == all_samples[j]

        if(cor_met == "PC"){
          mat_i <- as.matrix(PCA_res[idx_i,1:top_pcs])
          mat_j <- as.matrix(PCA_res[idx_j,1:top_pcs])
        }else if(cor_met == "HVG"){
          mat_i <- as.matrix(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj), idx_i]) %>% t()
          mat_j <- as.matrix(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj), idx_j]) %>% t()
        }

        cor_mat <- HiClimR::fastCor(t(rbind(mat_i,mat_j)), upperTri = T)
        cor_mat <- t(cor_mat[(nrow(mat_i)+1):nrow(cor_mat),1:nrow(mat_i)])

        hly_cor_mat <- cor_mat > hly_cor
        i_nn_mat <- matrix(0,nrow = nrow(mat_i),ncol = nrow(mat_j))
        j_nn_mat <- matrix(0,nrow = nrow(mat_i),ncol = nrow(mat_j))
        for(k in seq_len(nrow(mat_i))){
          i_nn_mat[k,cor_mat[k,] %>% order(decreasing = T) %>% .[1:nn_2]] <- T
        }
        for(k in seq_len(nrow(mat_j))){
          j_nn_mat[cor_mat[,k] %>% order(decreasing = T) %>% .[1:nn_2],k] <- T
        }
        graph_mat <- hly_cor_mat | (i_nn_mat & j_nn_mat)
        # Add edges
        row_id <- node_df[row.names(mat_i),1]
        col_id <- node_df[row.names(mat_j),1]
        edge_vec <- c()
        for(k in seq_len(nrow(graph_mat))){
          tmp_vec <- col_id[which(graph_mat[k,])]
          if(length(tmp_vec)>0){
            edge_vec <- c(edge_vec,
                          paste(row_id[k],tmp_vec) %>% str_split(pattern = " ") %>% unlist() %>% as.numeric())
          }
        }
        g <- igraph::add.edges(g,edges = edge_vec)

      }
    }

    louvain_res <- igraph::cluster_louvain(g, resolution = resolution)

    cl_df[["louvain"]] <- as.character(louvain_res$membership)
    cl_df[["raw_louvain"]] <- cl_df[["louvain"]]
    sml_cl_idx <- names(table(cl_df[["louvain"]]))[table(cl_df[["louvain"]]) < 2]
    cl_df[["louvain"]][cl_df[["louvain"]] %in% sml_cl_idx] <- NA
    if(rare_ct == "a"){
      if(cor_met == "PC"){
        cor_mat <- cor(t(as.matrix(PCA_res[is.na(cl_df[["louvain"]]),1:top_pcs])),method = "pearson")
      }else if(cor_met == "HVG") cor_mat <- cor(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj),is.na(cl_df[["louvain"]])],method = "pearson")

      louvain_res <- louvain_w_cor(cor_mat, nn_=20,res_ = 1)
      #table(louvain_res$membership)
      cl_df[["louvain_2"]] <- NA
      cl_df[["louvain_2"]][is.na(cl_df[["louvain"]])] <- as.character(louvain_res$membership)
      cl_df[["louvain"]] <- apply(cl_df,1,
                                  FUN = function(x){
                                    ifelse(is.na(x[["louvain"]]),
                                           paste0("R2_",x[["louvain_2"]]),
                                           paste0("R1_",x[["louvain"]]))
                                  })

    }


  }

  #cl_df[["umap_1"]] <- cl_expr_obj@reductions[["umap"]]@cell.embeddings[,1]
  #cl_df[["umap_2"]] <- cl_expr_obj@reductions[["umap"]]@cell.embeddings[,2]
  #cl_df[["layer"]] <- apply(cl_df,1,
  #                          FUN = function(x){
  #                            idx <- seu_ls[[x[["sample"]] ]]@meta.data[[cl_key]] == x[["cluster"]]
  #                            layer_vec <- seu_ls[[x[["sample"]] ]]@meta.data[["layer"]][idx]
  #                            table(layer_vec) %>% sort(decreasing = T) %>% names() %>% .[1]
  #
  #                          })
  toc <- Sys.time()
  if(verbose){
    message("Elasped Time:")
    message(format(difftime(toc, tic, units='secs'),digits = 3))
  }
  if(rtn_seurat){
    if(method == "MNN") list(seurat_obj = cl_expr_obj, cl_df = cl_df, g=g)
    else list(seurat_obj = cl_expr_obj, cl_df = cl_df)
  }else{
    cl_df
  }
}


#' Assign labels
#'
#' Assign labels to each spot based on integrated clustering results
#' @param seu_ls A list of Seurat Objects from `stage_1` results
#' @param cl_df A Dataframe from `stage_2`, containing the clustering results
#' @param anno Annotation string
#' @param cor_threshold Threshold for edge pruning in stage 1
#' @param cl_key Column for clustering labels
#'
#' @return A list of Seurat Objects
#'
#' @import scater
#' @import scry
#'
#' @export
#'
#' @examples
#' TBD
assign_label <- function(seu_ls, cl_df, anno, cor_threshold,
                         cl_key = "merged_cluster"){
  # Assign secondary clustering label
  for(i in names(seu_ls)){
    #message(i)
    seu_ls[[i]]@meta.data[[paste0("sec_cluster_",anno)]] <- apply(seu_ls[[i]]@meta.data,1,
                                                                FUN = function(x){
                                                                  cl_df$louvain[cl_df$sample==i & cl_df$cluster==x[[cl_key]]]
                                                                })
    seu_ls[[i]]@misc[[paste0("graph_plot_cluster_sec_",anno)]] <-
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["edges"]],
                       cor_threshold,paste0("sec_cluster_",anno)) +
      labs(title = paste("sample:", i))
  }
  seu_ls
}
