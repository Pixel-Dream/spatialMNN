#' Plot Confusion Matrix
#'
#' Draw Confusion Matrix using ComplexHeatmap
#' @param x,y vector with discrete values
#' @param col_title Heatmap title
#'
#' @return A Heatmap-class object.
#'
#' @import ComplexHeatmap
#' @import circlize
#'
#' @export
#'
#' @examples
#' x_ <- sample(1:4,100,replace = T)
#' y_ <- sample(1:5,100,replace = T)
#' plotConfusionMatrix(x_,y_,col_title = "Confusion Matrix")
#'
plotConfusionMatrix <- function(x,y,col_title = ""){
  na_index <- (is.na(x) | is.na(y))
  x <- x[!na_index]
  y <- y[!na_index]
  u_x <- unique(x)
  u_y <- unique(y)
  hm_mat <- matrix(0, nrow = length(u_x), ncol = length(u_y))
  row.names(hm_mat) <- u_x
  colnames(hm_mat) <- u_y
  for(i in seq_len(length(u_x))){
    for(j in seq_len(length(u_y))){
      hm_mat[i,j] = sum(x==u_x[i] & y ==u_y[j])
    }
  }
  hm_mat <- hm_mat/matrix(rep(rowSums(hm_mat),ncol(hm_mat)),
                          byrow = F,nrow = nrow(hm_mat), ncol = ncol(hm_mat))
  Heatmap(hm_mat,cluster_columns = F, cluster_rows = F,
          col = colorRamp2(seq(0, 1,length.out=5), viridis::viridis(5)),
          rect_gp = gpar(col = "white", lwd = 1),column_title = col_title)
}



#' Plot SRT Slide
#'
#' Draw SRT spot array with customized labels
#' @param meta_data A Data Frame contains at least coordinates
#' @param edge_df A Data Frame contains the edges information, at least 5 columns: `x, y, xend, yend, weight`
#' @param threshold Threshold of edge weight for displaying edges
#' @param col_sel column name in `meta_data`
#' @param pal Palette
#' @param flip Whether swap x-y coordinates
#'
#' @return A Heatmap-class object.
#'
#' @export
#'
#' @examples
#' #TBD

draw_slide_graph <- function(meta_data,
                             edge_df=NULL,
                             threshold=NULL,
                             col_sel,
                             pal = NULL,
                             flip = T){
  g <- ggplot()
  if(!is.null(edge_df)){
    edge_df[["x"]] <- meta_data[edge_df$from,'coord_x']
    edge_df[["y"]] <- meta_data[edge_df$from,'coord_y']
    edge_df[["xend"]] <- meta_data[edge_df$to,'coord_x']
    edge_df[["yend"]] <- meta_data[edge_df$to,'coord_y']
    if(is.null(threshold)){
      g <- g + geom_segment(mapping = aes(x=x,y=y,xend=xend,yend=yend,
                                          color = weight),
                            size=abs(edge_df$weight),
                            data = edge_df) +
        scale_colour_gradientn(colours = rev(viridis::rocket(n=10)))
    }else{
      edge_df[["filtered"]] <- as.numeric(edge_df$weight < threshold)
      g <- g + geom_segment(mapping = aes(x=x,y=y,xend=xend,yend=yend,
                                          color = weight<threshold),
                            size=abs(edge_df$filtered+0.1),
                            data = edge_df) +
        scale_color_manual(values = c("grey70","red"))
    }
    g <- g + new_scale_colour()
  }
  g <- g + geom_point(mapping = aes_string(x = "coord_x", y = "coord_y", color=col_sel),
                      data = meta_data)
  if(is.numeric(meta_data[[col_sel]])){
    require(viridis)
    g <- g + scale_colour_gradientn(colours = viridis::viridis(n=10))
  }else{
    if(is.null(pal)) g <- g + scale_color_discrete()
    else g <- g + scale_colour_manual(values = pal)
  }
  g <- g + theme_classic()

  if(flip) g + coord_flip()
  else g
}
