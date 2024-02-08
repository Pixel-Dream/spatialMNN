##### Helper functions #####
library(SpatialExperiment)
library(STexampleData)
library(scater)             # log-transformation
library(here)
# * Smooth & Circular -------------------------------------------------------
smth_circ_fun <- function(x, y){ (x-0.5)^2 + (y-0.5)^2 }

# * Smooth & Linear -------------------------------------------------------
smth_lnr_fun <- function(x, y) { 0.5-x + 0.5-y }


# * Layered & Circular -------------------------------------------------------
# Calcualte center of the graph
lyr_circ_fun <- function(x,y){
  # ret <- NA_real_
  dplyr::case_when(
    (0.5-x)^2 + (0.5-y)^2 < (0.125)^2 ~ 40,
    (0.5-x)^2 + (0.5-y)^2 < (0.25)^2 ~ 30,
    (0.5-x)^2 + (0.5-y)^2 < (0.375)^2 ~ 20,
    (0.5-x)^2 + (0.5-y)^2 < (1)^2 ~ 10,
    TRUE ~ NA_real_
  )
  # return(ret)
}

# * Layered & Linear -------------------------------------------------------
lyr_lnr_fun <- function(x,y){
  # ret <- NA_real_
  dplyr::case_when(
    1.5 - x - y <= 0 ~ 40,
    1 - x - y < 0 ~ 30,
    0.5 - x - y < 0 ~ 20,
    0 - x - y < 0 ~ 10,
    TRUE ~ NA_real_
  )
}

# Corresponding Functions to Create Covariates ----------------------------
# * Layered Pattern -------------------------------------------------------
# Caliberate

lnr_group_fun <- function(x,y){
  dplyr::case_when(
    1.5 - x - y <= 0 ~ "Group 1",
    1 - x - y < 0 ~ "Group 2",
    0.5 - x - y < 0 ~ "Group 3",
    0 - x - y < 0 ~  "Group 4",
    TRUE ~ NA_character_
  )
}

# # * Circular Pattern -------------------------------------------------------
circ_group_fun <- function(x,y){
  dplyr::case_when(
    (0.5-x)^2 + (0.5-y)^2 < (0.125)^2 ~ "Group 1",
    (0.5-x)^2 + (0.5-y)^2 < (0.25)^2 ~ "Group 2",
    (0.5-x)^2 + (0.5-y)^2 < (0.375)^2 ~ "Group 3",
    (0.5-x)^2 + (0.5-y)^2 < (1)^2 ~ "Group 4",
    TRUE ~ NA_character_
  )
}

semicirc_group_fun <- function(x,y){
  dplyr::case_when(
    (0.5-x)^2 + y^2 < (0.3)^2 ~ "Group 1",
    (0.5-x)^2 + y^2 < (0.5)^2 ~ "Group 2",
    (0.5-x)^2 + y^2 < (0.7)^2 ~ "Group 3",
    (0.5-x)^2 + y^2 < 1.25  ~ "Group 4",
    #(0.5-x)^2 + y^2 < 1.25 ~ "Group 5",
    TRUE ~ NA_character_
  )
}

qurtcirc_group_fun <- function(x,y){
  dplyr::case_when(
    x^2 + y^2 < (0.4)^2 ~ "Group 1",
    x^2 + y^2 < (0.8)^2 ~ "Group 2",
    x^2 + y^2 < (1)^2 ~ "Group 3",
    x^2 + y^2 < 2 ~ "Group 4",
    #x^2 + y^2 < 2 ~ "Group 6",
    TRUE ~ NA_character_
  )
}

norm <- function(x){
  (x - min(x))/ (max(x) - min(x))
}

assign_celltype <- function(probs, cell_max = 1){
  stopifnot(sum(probs) == 1)
  sample.int(length(probs), sample(seq_len(cell_max),1),
             prob = probs, replace = T)
}

##### Main Simulation #####
spe_template <- STexampleData::Visium_humanDLPFC()

my_sim <- function(it = 1, tmp_spe=spe_template, n_gene = 200, noise = 0.2, ig_ratio = 1, top_pcs = 30, 
                   n_group = 4, n_celltype = 4, map_mat = NULL, cell_max = 1, n_sample = 4, 
                   segmentation = T, integration = T){
  n_spots <- ncol(tmp_spe)
  stopifnot(ig_ratio<=1)
  stopifnot(ncol(map_mat)==n_celltype)
  seu_ls <- list()
  set.seed(it)
  
  if(is.null(map_mat)){
    map_mat <- matrix(0,nrow = n_group, ncol = n_celltype, byrow = T)
    for(i in 1:n_group) map_mat[i,i] <- 1
  }
  
  row.names(map_mat) <- paste("Group",1:n_group)
  colnames(map_mat) <- paste0("Celltype_",1:n_celltype)
  
  total_gene_num <- round(n_gene*n_celltype/ig_ratio)
  # Ground truth ranking
  rowData_df <- data.frame(
    gene_idx  = 1:total_gene_num, 
    mu_shift = runif(n = total_gene_num, min = 2, max = 4),# Different mean expression level
    var_scale = runif(n = total_gene_num, min = 1, max = 4) # Different effect size
  ) |> 
    mutate(gene_name = paste0("gene_", gene_idx),
           marker=ifelse(ceiling(gene_idx/n_gene) <= n_celltype,
                         paste0("Celltype_",ceiling(gene_idx/n_gene)),
                         "noise")) |> 
    column_to_rownames("gene_name")
  
  if(n_sample!=4){
    pattern_vec <- sample(c("lnr","circ","semicirc","qurtcirc"),n_sample,replace = T)
  }else pattern_vec <- c("lnr","circ","semicirc","qurtcirc")
  
  # Choose one of the spatial pattern
  for(i in seq_along(pattern_vec)) {
    pattern <- pattern_vec[i]
    message(paste("Simulating Pattern #",i,":",pattern))
    str_func <- switch(pattern,
                       "lnr"=lnr_group_fun,
                       "circ"=circ_group_fun,
                       "semicirc"=semicirc_group_fun,
                       "qurtcirc"=qurtcirc_group_fun) 
    set.seed(it)
    spa_str_df <- spatialCoords(tmp_spe) |> 
      data.frame() |> 
      mutate(
        x = norm(pxl_col_in_fullres),
        y = norm(pxl_row_in_fullres),
        z = str_func(
          x, y
        ),
        cell_type = sapply(z,
                           FUN = function(x){
                             paste0("Celltype_",
                                    assign_celltype(map_mat[x,], cell_max)) %>%
                               paste(collapse = ",")
                           }),
        cell_num = sapply(cell_type,
                          FUN = function(x){
                            str_split(x,pattern=',') %>% unlist() %>% length()
                          }),
        
        #std_z = scale(z) # standardized
        coord_x = tmp_spe@colData@listData[["array_col"]],
        coord_y = tmp_spe@colData@listData[["array_row"]]
      )
    # Check cell type ratio
    
    for(j in 1:n_celltype){
      spa_str_df[[paste0("Celltype_",j)]] <- sapply(spa_str_df$cell_type,
                                                    FUN = function(x){
                                                      str_count(x,pattern = paste0("Celltype_",j))
                                                    })
    }
    
    #stopifnot(all(!is.na(spa_str_df$std_z)))
    # Check region
    #ggplot(spa_str_df,aes(x=pxl_col_in_fullres,y=pxl_row_in_fullres, color=as.factor(z))) + 
    #  geom_point() + 
    #  theme_classic()
    # Simulated raw counts following Poisson noise
    # gene by spot
    gene_count_mat <- 
      map2(.x = rowData_df$mu_shift,
           #.y = rowData_df$var_scale,
           .y = rowData_df$gene_idx,
           .f = function(shift, gene_idx){
             #eta_vec <- spa_str_df$std_z*scale + shift
             markers <- rowData_df$marker[rowData_df$gene_idx == gene_idx]
             if(markers=="noise"){
               eta_vec <- rep(shift,n_spots)
               ret_vec <- rnorm(
                 n = n_spots,
                 mean = eta_vec,
                 sd = noise
               )
               
             }else{
               cell_type_vec <- str_split(spa_str_df$cell_type,pattern = ",") %>% unlist()
               eta_vec <- (cell_type_vec == markers)*shift
               tmp_vec <- rnorm(
                 n = length(eta_vec),
                 mean = eta_vec,
                 sd = noise
               )
               l_vec <- r_vec <- cumsum(spa_str_df$cell_num)
               l_vec[2:length(r_vec)] <- r_vec[2:length(r_vec) - 1]+1
               l_vec[1] = 1
               ret_vec <- sapply(1:nrow(spa_str_df),
                                 FUN = function(x){
                                   #pmax(round(mean(tmp_vec[l_vec[x]:r_vec[x]])),0)
                                   mean(tmp_vec[l_vec[x]:r_vec[x]])
                                 })
             }
             
             stopifnot(length(ret_vec) == n_spots)
             return(data.frame(ret_vec) %>% `colnames<-`(paste0("gene_",gene_idx)))
           }) |> 
      list_cbind() |>  
      t()
    
    # Cap
    gene_count_mat[gene_count_mat<0] <- 0
    
    rownames(gene_count_mat) <- rownames(rowData_df)
    colnames(gene_count_mat) <- rownames(spa_str_df)
    
    seu_obj <- CreateSeuratObject(counts = gene_count_mat,
                                  project = "tmp_seurat",
                                  meta.data = spa_str_df)
    seu_obj <- ScaleData(seu_obj, features = row.names(seu_obj), do.scale = F, do.center = F)
    #DoHeatmap(seu_obj, features = row.names(seu_obj)[c(1:10,201:210,401:410,601:610,801:810)],group.by = "z") + NoLegend()
    seu_ls[[paste0(pattern,"_",i)]] <- seu_obj
  }
  # Create Seurat
  #ggplot(seu_obj@meta.data,aes(x = as.factor(cell_num))) + 
  #  geom_histogram(stat="count") +
  #  theme_classic()
  for(pattern in names(seu_ls)){
    seu_ls[[pattern]]@meta.data[["barcode"]] = row.names(seu_ls[[pattern]]@meta.data)
    seu_ls[[pattern]]@meta.data[["batch"]] = pattern
    seu_ls[[pattern]]@meta.data[["layer"]] = seu_ls[[pattern]]@meta.data[["z"]]
    seu_ls[[pattern]]@meta.data[["row"]] = seu_ls[[pattern]]@meta.data[["coord_x"]]
    seu_ls[[pattern]]@meta.data[["col"]] = seu_ls[[pattern]]@meta.data[["coord_y"]]
    seu_ls[[pattern]]@misc[["group"]] <- draw_slide_graph(seu_ls[[pattern]]@meta.data,col_sel = "z")
  }
  
  if(segmentation) seu_ls <- stage_1(seu_ls,preprocess = F, top_pcs = top_pcs)
  if(segmentation & integration){
    rtn_ls <- stage_2(seu_ls,cl_key = "merged_cluster",rtn_seurat = T,nn_2 = 1,method = "MNN", top_pcs = top_pcs)
    seu_ls <- assign_label(seu_ls,rtn_ls[["cl_df"]], "MNN", cl_key = "merged_cluster")
  }
  seu_ls
}
