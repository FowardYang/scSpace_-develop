library(Seurat)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(progress)
library(ggsignif)
library(pheatmap)
library(viridisLite)

# data processing
data_process <- function(sc_data,
                         sc_meta){
  rownames(sc_meta) <- colnames(sc_data)
  sc_seu <- CreateSeuratObject(sc_data, meta.data = sc_meta)
  sc_seu <- NormalizeData(sc_seu)
  sc_seu <- FindVariableFeatures(sc_seu)
  sc_seu <- ScaleData(sc_seu)
  sc_seu <- RunPCA(sc_seu)
  sc_seu <- RunTSNE(sc_seu)
  
  return(sc_seu)
}


# clustering
get_louvain_clu <- function(sc_seu,
                            res){
  sc_seu <- FindNeighbors(sc_seu)
  sc_seu <- FindClusters(sc_seu, resolution = res)
  sc_seu$louvain <- sc_seu$seurat_clusters
  
  return(sc_seu)
}


get_kmeans_clu <- function(sc_seu,
                           target_num){
  sc_pca_ori <- sc_seu@reductions$pca@cell.embeddings
  kmeans_label <- kmeans(sc_pca_ori, target_num)$cluster
  sc_seu$kmeans <- kmeans_label
  
  return(sc_seu)
}


# calculate dist
calculate_dist <- function(scspace_meta,
                           group_by,
                           selected_type,
                           ignore_select_type = FALSE){
  dist <- as.matrix(dist(scspace_meta[,c('Pseudo_space1', 'Pseudo_space2')]))
  
  ct_type <- sort(unique(scspace_meta[[group_by]]))
  
  if(ignore_select_type){
    ct_type <- setdiff(ct_type, selected_type)
  }
  
  select_ct_index <- which(scspace_meta[[group_by]] == selected_type)
  
  pseudo_dist_table <- data.frame()
  
  message('Beginning normalized distance calculating...')
  pb <- progress::progress_bar$new(format = ' Calculating [:bar] :percent eta: :eta', 
                                   total = length(ct_type), clear = FALSE, width = 60, 
                                   complete = "+", incomplete = "-")
  
  # dist
  for (i in ct_type) {
    ct_index <- which(scspace_meta[[group_by]] == i)
    dist_list <- c()
    if(i == selected_type){
      for (j in ct_index) {
        col_index <- j
        dist_list <- c(dist_list, dist[ct_index, col_index])
        ct_index <- setdiff(ct_index, col_index)
      }
    }
    else{
      for (j in select_ct_index) {
        col_index <- j
        dist_list <- c(dist_list, dist[ct_index, col_index])
      }
    }
    
    dist_table <- data.frame(dist = dist_list,
                             group = rep(paste0(selected_type, '_', i), length(dist_list)))
    
    pseudo_dist_table <- rbind(pseudo_dist_table, dist_table)
    
    pb$tick()
  }
  print('Normalized distance calculating done.')
  
  pseudo_dist_table$dist <- pseudo_dist_table$dist / max(pseudo_dist_table$dist)
  return(pseudo_dist_table)
  
}


# plot
# pseudo space
plot_pseudo_space <- function(scspace_meta,
                              group_by,
                              show_ct = NULL){
  if(is.null(show_ct)){
    ggplot(scspace_meta) +
      geom_point(aes(Pseudo_space1, Pseudo_space2, color = .data[[group_by]])) +
      theme_base()
  }else{
    ggplot() +
      geom_point(scspace_meta[!(scspace_meta[[group_by]] %in% show_ct), ], 
                 mapping = aes(Pseudo_space1, Pseudo_space2,color = .data[[group_by]]), color = 'grey90') +
      geom_point(scspace_meta[scspace_meta[[group_by]] %in% show_ct, ], 
                 mapping = aes(Pseudo_space1, Pseudo_space2,color = .data[[group_by]])) +
      theme_base()
  }
}


# tsne
plot_tsne <- function(scspace_meta,
                      group_by){
  scspace_meta[[group_by]] <- as.factor(scspace_meta[[group_by]])
  ggplot(scspace_meta) +
    geom_point(aes(tsne1, tsne2, color = .data[[group_by]])) +
    theme_base()
}


plot_dist <- function(dist_obj){
  plot_obj <- data.frame(aggregate(normalized_dist$dist, by=list(normalized_dist$group), median))
  colnames(plot_obj) <- c('group', 'median')
  plot_obj$sd <-  aggregate(normalized_dist$dist, by=list(normalized_dist$group), sd)$x
  
  ggplot(plot_obj, aes(x = group, y = median, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = (median - sd), ymax = (median + sd)), width = 0.4, size = 1) +
    scale_color_manual(values = c(rep('black', nrow(plot_obj)))) +
    theme_base()
}


plot_density <- function(scspace_meta, 
                         group_by,
                         selected_type,
                         color_high){
  ct_type <- sort(unique(scspace_meta[[group_by]]))
  p <- ggplot(scspace_meta ,mapping = aes(x = Pseudo_space1, y = Pseudo_space2,color = .data[[group_by]])) +
    geom_point(size = 3) +
    stat_density2d(scspace_meta[scspace_meta[[group_by]] == selected_type, ], mapping = aes(fill = ..density..), geom = "raster", contour = F, n = 200) +
    scale_fill_gradient(low = NA, high = color_high) +
    scale_color_manual(values = rep('grey90', length(ct_type))) +
    theme_base() +
    ggtitle(selected_type) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}




