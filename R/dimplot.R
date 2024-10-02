#' Beautify Seurat's DimPlot Function
#'
#' This function enhances the DimPlot function from the Seurat package by adding
#' additional customization options for plot aesthetics and layout.
#'
#' @param object A Seurat object.
#' @param groupBy Grouping variable.
#' @param style Style of the plot. Defaults to 1.
#' @param plot.title Title of the plot. Defaults to NA.
#' @param legend.point.size Size of the legend points. Defaults to 4.
#' @param reduction Reduction method to use. Defaults to "umap".
#' @param legend.position Position of the legend. Defaults to "bottom".
#' @param label Logical; whether to label the plot. Defaults to FALSE.
#' @param label.size Size of the labels. Defaults to 4.
#' @param point.size Size of the points. Defaults to 0.5.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' dimPlotMisc(seurat_object, groupBy = "celltype", style = 1, label = TRUE)
#' }
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom ggplot2 ggplot geom_point theme_bw theme element_text element_blank element_rect labs scale_fill_manual scale_color_manual guides guide_legend alpha
#' @importFrom paletteer paletteer_d
#' @importFrom ggrepel geom_label_repel
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate inner_join select rename
#' @importFrom tidydr theme_dr
#' @importFrom grid arrow unit
#' @importFrom Seurat DimPlot
#'
#' @export
dimPlotMisc <- function(object,
                           groupBy,
                           style=1,
                           plot.title = NA,
                           legend.point.size = 4,
                           reduction = "umap",
                           legend.position = "bottom",
                           label = FALSE,
                           label.size = 4,point.size = 0.5){
  suppressPackageStartupMessages({
    library(tibble)
    library(ggplot2)
    library(paletteer)
    library(ggrepel)
    library(ggthemes)
    library(stringr)
  })
  # groupBy = 'celltype'
  # object = sce

  if(!is.null(plot.title)){
    if(is.na(plot.title)){
      plot.title = groupBy
    }
  }
  # (1) 获取非线性降维坐标，可以选择tsne或者umap，前提是存在哦：
  plot_data = object@reductions[[reduction]]@cell.embeddings %>%
    as.data.frame() %>%
    dplyr::mutate(Barcode = rownames(.)) %>%
    dplyr::inner_join(object@meta.data %>% tibble::rownames_to_column('Barcode'), by = 'Barcode') %>%
    tibble::column_to_rownames('Barcode') %>%
    dplyr::select(1,2, {{groupBy}}) %>%
    dplyr::rename(Group = as.name(groupBy))
  colnames(plot_data) = gsub("_"," ",colnames(plot_data))

  if(grepl(x = reduction,pattern = '*tsne*' )) {
    colnames(plot_data)[1] = 'tSNE 1'
    colnames(plot_data)[2] = 'tSNE 2'
  } else {
    colnames(plot_data)[1] = 'UMAP 1'
    colnames(plot_data)[2] = 'UMAP 2'
  }

  # (2) 生成聚类中心坐标
  centroids = aggregate(as.matrix(plot_data[,c(1,2)]) ~ Group,
                        data = plot_data,
                        FUN = mean)

  # (3) 两种风格的绘图
  if(grepl(x = reduction,pattern = '*tsne*' )) {
    x = 'tSNE 1'; y = 'tSNE 2'
    segment.df=data.frame(x=c(-55,-55),
                          xend=c(-30,-55),
                          y=c(-44,-44),
                          yend=c(-44,-19))
  } else {
    x = 'UMAP 1'; y = 'UMAP 2'
    segment.df=data.frame(x=c(-15,-15),
                          xend=c(-5,-15),
                          y=c(-15,-15),
                          yend=c(-15,-5))
  }

  if(style == 1){
    p1 =  DimPlot(object, reduction = reduction,group.by = groupBy,
                  label = label,label.box = T,label.size = label.size,repel = T) + theme_bw() +
      labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
      theme(panel.grid=element_blank(), # 去网格线
            plot.title = element_text(size = 12,color="black",hjust = 0.5),
            axis.text.x = element_text(size = 10, color = 'black'),
            axis.text.y = element_text(size = 10, color = 'black'),
            axis.title.x = element_text(size = 10, color = 'black'),
            axis.title.y = element_text(size = 10, color = 'black'),
            axis.ticks = element_line(color = 'black', lineend = 'round'),
            legend.position = legend.position,
            legend.text = element_text(size = 10, color = 'black'),
            legend.title = element_text(size = 10, color = 'black'),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      guides(color=guide_legend(override.aes = list(size=legend.point.size)))

    if(length(unique(plot_data$Group))<21){
      p1 = p1 +
        scale_fill_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 1)) +
        scale_color_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 1))
    }
    return(p1)
  }

  if(style > 1){
    if(style == 2){
      p2 = ggplot(data = plot_data, mapping = aes(x = !!as.name(x), y = !!as.name(y))) +
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = legend.point.size) +
        geom_point(mapping = aes(fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = 3,show.legend = F) +
        theme_bw() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
        theme(
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12,color="black",hjust = 0.5),
          axis.text.x = element_text(size = 10, color = 'black'),
          axis.text.y = element_text(size = 10, color = 'black'),
          axis.title.x = element_text(size = 10, color = 'black'),
          axis.title.y = element_text(size = 10, color = 'black'),
          axis.ticks = element_line(color = 'black', lineend = 'round'),
          legend.position = legend.position,
          legend.text = element_text(size = 10, color = 'black'),
          legend.title = element_text(size = 13, color = 'black'),
          panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
    }

    if(style == 3){
      p2 = ggplot(data = plot_data,mapping = aes(x = !!as.name(x), y = !!as.name(y),color=Group))+
        geom_point(mapping = aes(color = Group),size = point.size,show.legend = F) +
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 color = Group),alpha = 1,size = 1,show.legend = T)+
        theme_classic() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
        tidydr::theme_dr(arrow = grid::arrow(length = unit(0.3, "cm"), type = "open")) +
        # geom_segment(data = segment.df,mapping = aes(x=x,xend=xend,y=y,yend=yend),
        #              arrow = arrow(length=unit(0.3, "cm")),color = "black",size=0.6,show.legend = F) +
        #zfm:这个主题更简洁一点
        coord_cartesian(clip = "off")+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = legend.position,
              plot.title = element_text(size = 12,color="black",hjust = 0.5),
              # axis.title.x.bottom = element_text(hjust = 0.12,size = 12,margin=margin(-5,0,0,0)),
              # axis.title.y.left = element_text(hjust = 0.12,size = 12,margin=margin(0,-8,0,0)),
              axis.title.x.bottom = element_text(hjust = 0.12,size = 12),
              axis.title.y.left = element_text(hjust = 0.12,size = 12),
              axis.title.x = element_text(size = 12, color = 'black'),
              axis.title.y = element_text(size = 12, color = 'black'),
              # legend.text = element_text(size = 10, color = 'black'),
              legend.title = element_text(size = 10, color = 'black'),
              # plot.margin = margin(50,50,50,50),
              legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(size=legend.point.size)))
    }
    if(style == 4){
      p2 = ggplot(data = plot_data,mapping = aes(x = !!as.name(x), y = !!as.name(y),color=Group))+
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = legend.point.size) +
        geom_point(mapping = aes(fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = 3,show.legend = F) +
        theme_classic() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
        tidydr::theme_dr(arrow = grid::arrow(length = unit(0.3, "cm"), type = "open")) +
        # geom_segment(data = segment.df,mapping = aes(x=x,xend=xend,y=y,yend=yend),
        #              arrow = arrow(length=unit(0.3, "cm")),color = "black",size=0.6,show.legend = F) +
        #zfm:这个主题更简洁一点
        coord_cartesian(clip = "off")+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12,color="black",hjust = 0.5),
              legend.position = legend.position,
              # axis.title.x.bottom = element_text(hjust = 0.12,size = 12,margin=margin(-5,0,0,0)),
              # axis.title.y.left = element_text(hjust = 0.12,size = 12,margin=margin(0,-8,0,0)),
              axis.title.x.bottom = element_text(hjust = 0.12,size = 12),
              axis.title.y.left = element_text(hjust = 0.12,size = 12),
              axis.title.x = element_text(size = 12, color = 'black'),
              axis.title.y = element_text(size = 12, color = 'black'),
              # legend.text = element_text(size = 10, color = 'black'),
              legend.title = element_text(size = 10, color = 'black'),
              # plot.margin = margin(50,50,50,50),
              legend.background = element_blank())
    }
    if(label == T){
      if(length(unique(plot_data$Group))<21){
        p3 = p2 + ggrepel::geom_label_repel(data = centroids,
                                            mapping = aes(fill = Group,
                                                          label= Group),
                                            size = label.size,
                                            fontface = 'plain', #bold
                                            color="black",
                                            # family = 'Arial',
                                            show.legend = FALSE) +
          scale_fill_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 1)) +
          scale_color_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 1))
      }else{
        warning("Too many label box (> 20). Show label is off!")
        p3 = p2
        #   p3 = p2 + ggrepel::geom_label_repel(data = centroids,
        #                                       mapping = aes(fill = Group,
        #                                                     label= Group),
        #                                       size = label.size,
        #                                       fontface = 'plain', #bold
        #                                       color="black",
        #                                       family = 'Arial',
        #                                       show.legend = FALSE)
        #   # scale_fill_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 0.65)) +
        #   # scale_color_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 0.65))
      }
      return(p3)
    }
    if(label == F){
      if(length(unique(plot_data$Group))<21){
        p3 = p2 +
          scale_color_manual(values = alpha(paletteer::paletteer_d('ggthemes::Classic_20'), 1))

      }else{
        p3 = p2
      }
      return(p3)
    }
  }
}
