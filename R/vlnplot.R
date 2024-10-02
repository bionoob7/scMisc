#' Stacked Violin Plot for Seurat Object
#'
#' This function generates a stacked violin plot for a Seurat object, allowing for
#' customization of plot aesthetics and layout.
#'
#' @param object A Seurat object.
#' @param features A vector of features (genes) to plot.
#' @param group Grouping variable.
#' @param text.size Size of the text. Defaults to 8.
#' @param text.angle Angle of the text. Defaults to 45.
#' @param text.hjust Horizontal justification of the text. Defaults to 1.
#' @param legend.position Position of the legend. Defaults to "right".
#' @param switch Switch position for facet labels. Defaults to "left".
#' @param fill.cols Colors for violin plot fill. Defaults to NULL.
#' @param cols Colors for cell type annotation. Defaults to NULL.
#' @param widths Widths for plot layout. Defaults to c(3, 0.08).
#' @param heights Heights for plot layout. Defaults to c(3, 0.08).
#' @param legend.title Title for the legend. Defaults to "Ave.exp".
#' @param x.lab Label for the x-axis. Defaults to NULL.
#' @param y.lab Label for the y-axis. Defaults to NULL.
#' @param title Title of the plot. Defaults to NULL.
#' @param ... Additional parameters passed to FetchData.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' vlnPlotStacked(object = seurat.data,
#'                features = marker_gene,
#'                group = "celltype",
#'                legend.key.size = 0.3,
#'                combine = TRUE, stack = TRUE, flip = TRUE,
#'                cols = paletteer::paletteer_d('ggthemes::Classic_20'))
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr select mutate group_by rename
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_violin scale_fill_gradientn facet_grid theme element_text element_blank labs scale_x_discrete scale_y_discrete geom_tile scale_fill_manual unit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork plot_layout
#' @importFrom Seurat FetchData
#' @importFrom paletteer paletteer_d
#'
#' @export
vlnPlotStacked <- function(object,features,
                     group,
                     text.size = 8,
                     text.angle = 45,
                     text.hjust = 1,
                     legend.position = "right",
                     switch = "left",
                     fill.cols = NULL,
                     cols=NULL,
                     widths  = c(3,0.08),
                     heights = c(3,0.08),
                     legend.title = "Ave.exp",
                     x.lab=NULL,y.lab=NULL,title =NULL,...){
  suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(patchwork)
    library(grDevices)
    library(dplyr)
  })

  if(T){
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"),
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                     #axis.line = element_line(color = "black"),
                     #axis.ticks = element_line(color = "black"),
                     #panel.grid.minor.y = element_blank(),
                     #panel.grid.minor.x = element_blank(),
                     panel.grid=element_blank(), # 去网格线
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size),
                     # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
                     strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07",
    )
  }
  #从Seurat对象中提取细胞注释以及基因表达量
  features = features[features%in%row.names(object)]
  vln.dat=FetchData(object,c(features,group))
  colnames(vln.dat)[length(features)+1] = "celltype.sub"
  # 定义因子顺序，防止画图时对细胞注释进行重排
  # vln.dat$celltype.sub
  vln.dat=vln.dat[order(vln.dat$celltype.sub),]
  vln.dat.melt=vln.dat %>%
    reshape2::melt(,features) %>%
    dplyr::rename("Gene"="variable") %>%
    dplyr::group_by(celltype.sub,Gene) %>%
    dplyr::mutate(fillcolor=mean(value))

  if(switch == "right"){
    switch = "x"
  }else{
    switch = "y"
  }
  if (is.null(fill.cols)) {
    # 小提琴图的填充颜色
    pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)
    # 堆积小提琴图
    p1 = ggplot(vln.dat.melt,aes(x=celltype.sub,y=value,fill=fillcolor))+
      # 把小提琴图的外缘轮廓去除
      geom_violin(linetype="blank",scale = "width")+
      scale_fill_gradientn(colors=pal,name=legend.title)+
      facet_grid(Gene~.,switch = switch)+mytheme+
      theme(panel.grid = element_blank(),
            strip.text.y.left = element_text(angle = 0,hjust=1),
            strip.text.y = element_text(angle = 0,hjust=1),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            legend.position = legend.position

      )
    # p1
  }else{
    p1 = ggplot(vln.dat.melt,aes(x=celltype.sub,y=value,fill=celltype.sub))+
      # 把小提琴图的外缘轮廓去除
      geom_violin(linetype="blank",scale = "width",aes(fill=celltype.sub))+
      scale_fill_manual(values = fill.cols ,name=legend.title)+
      facet_grid(Gene~.,switch = switch)+mytheme+
      theme(panel.grid = element_blank(),
            strip.text.y.left = element_text(angle = 0,hjust=1),
            strip.text.y = element_text(angle = 0,hjust=1),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            legend.position = legend.position

      )
  }

  if(switch == "x"){
    p1 = p1 + theme(strip.text.y = element_text(angle = 0,hjust=0))
  }

  p1 = p1+labs(x=x.lab,y=y.lab,title =title)


  # 我们用geom_tile()完成细胞注释，当然也可以向上面一样用geom_segment或者geom_bar来完成
  p3=ggplot(vln.dat%>%
              dplyr::select(celltype.sub)%>%
              unique()%>%
              dplyr::mutate(value='A'),
            aes(x=celltype.sub,y=value,fill=celltype.sub))+
    geom_tile()+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    # facet_grid(.~celltype.sub,scales = "free",space = "free",switch = 'x')+
    mytheme+
    theme(panel.background = element_blank(),
          strip.text = element_text(angle = 0,hjust = 0.5,vjust = 1),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          legend.position = "none")
  p3
  if(!is.null(cols)){
    p3 = p3+scale_fill_manual(values = cols)
  }
  (p1+ theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
    (p3 + theme(plot.margin = unit(c(0,30,0,0), "pt")))+
    plot_layout(ncol = 1, widths  = widths,heights = heights)
}



#' Create a beautified violin plot for Seurat objects
#'
#' This function enhances Seurat's VlnPlot by adding mean points, error bars, and custom theming.
#'
#' @param object A Seurat object
#' @param features Features to plot (e.g. genes)
#' @param cols Colors for the plot (optional)
#' @param idents Identity classes to include in the plot
#' @param pt.size Size of the points in the plot (default: 0)
#' @param sort Whether to sort the identity classes (default: FALSE)
#' @param assay Assay to use
#' @param group.by Variable to group cells by
#' @param split.by Variable to split plots by
#' @param adjust Adjust parameter for density estimation (default: 1)
#' @param y.max Maximum y-axis value
#' @param same.y.lims Use the same y-axis limits for all plots (default: FALSE)
#' @param log Use log scale for y-axis (default: FALSE)
#' @param ncol Number of columns in the plot grid
#' @param split.plot Whether to split the plot (default: FALSE)
#' @param stack Whether to stack the plots (default: FALSE)
#' @param combine Whether to combine the plots (default: TRUE)
#' @param fill.by How to fill the violin plots (default: "feature")
#' @param flip Flip the plot coordinates (default: FALSE)
#' @param add.noise Add jitter to the plot (default: TRUE)
#' @param raster Convert points to raster format
#' @param x.lab Label for x-axis (optional)
#' @param y.lab Label for y-axis (default: "Expression level")
#' @param legend.position Position of the legend (default: "none")
#' @param ... Additional arguments passed to Seurat's VlnPlot function
#'
#' @return A ggplot object representing the beautified violin plot
#'
#' @examples
#' \dontrun{
#' VlnPlotMisc(object = seurat.data, features = c("CD79A", "MS4A1"), group.by = "celltype",pt.size = 0)
#' }
#'
#' @import ggplot2
#' @import ggpubr
#'
#' @export
vlnPlotMisc <- function(object,features, cols = NULL,idents = NULL, pt.size = 0,
                      sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL,
                      adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE,
                      ncol = NULL,  split.plot = FALSE, stack = FALSE,
                      combine = TRUE, fill.by = "feature", flip = FALSE, add.noise = TRUE,
                      raster = NULL,x.lab=NULL,y.lab="Expression level",
                      legend.position = "none",...
){
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
  })
  VlnPlot(object = object, features = features, cols = cols,
          pt.size = pt.size,idents = idents,
          sort = sort, assay = assay, group.by = group.by, split.by = split.by,
          adjust = adjust, y.max = y.max, same.y.lims = same.y.lims, log = log,
          ncol = ncol, split.plot = split.plot, stack = stack,
          combine = combine, fill.by = fill.by, flip = flip, add.noise = add.noise, raster = raster,...
  )&
    labs(x=x.lab,y=y.lab)&
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.9)) &
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,position = position_dodge(0.9))&
    mytheme&theme(legend.position = legend.position)
}
