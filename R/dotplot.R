#' Beautify Seurat's DotPlot Function
#'
#' This function enhances the DotPlot function from the Seurat package by adding
#' additional customization options for plot aesthetics and layout.
#'
#' @param object A Seurat object.
#' @param features A vector of features to plot.
#' @param assay The assay to use; defaults to NULL.
#' @param scale Logical; whether to scale the data. Defaults to TRUE.
#' @param dot.range.min Minimum value for the dot size scale. Defaults to 0.
#' @param dot.range.max Maximum value for the dot size scale. Defaults to 3.5.
#' @param Combine Logical; whether to combine plots. Defaults to FALSE.
#' @param legend.position Position of the legend. Defaults to "right".
#' @param label.size Size of the labels. Defaults to 4.
#' @param label_widths Width of the label panel. Defaults to 0.1.
#' @param x.lab Label for the x-axis. Defaults to NULL.
#' @param y.lab Label for the y-axis. Defaults to NULL.
#' @param title Title of the plot. Defaults to NULL.
#' @param text.size Size of the text. Defaults to 8.
#' @param text.angle Angle of the text. Defaults to 90.
#' @param text.vjust Vertical justification of the text. Defaults to 0.5.
#' @param text.hjust Horizontal justification of the text. Defaults to 1.
#' @param group.by Grouping variable. Defaults to NULL.
#' @param color.use Colors to use for the plot. Defaults to NULL.
#' @param cols Colors for the dot plot. Defaults to c("lightgrey", "blue").
#' @param legend.key.size Size of the legend keys. Defaults to 0.5.
#' @param col.min Minimum value for color scale. Defaults to -2.5.
#' @param col.max Maximum value for color scale. Defaults to 2.5.
#' @param dot.min Minimum value for dot size. Defaults to 0.
#' @param dot.scale Scaling factor for dot size. Defaults to 6.
#' @param idents Identities to include. Defaults to NULL.
#' @param split.by Variable to split by. Defaults to NULL.
#' @param cluster.idents Logical; whether to cluster identities. Defaults to FALSE.
#' @param scale.by Scaling method. Defaults to "radius".
#' @param scale.min Minimum scale value. Defaults to NA.
#' @param scale.max Maximum scale value. Defaults to NA.
#' @param ... Additional parameters passed to DotPlot.
#'
#' @return A ggplot object or a combined plot if Combine is TRUE.
#'
#' @examples
#' \dontrun{
#' dotPlotMisc(seurat_object, features = c("Gene1", "Gene2"), Combine = TRUE)
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot geom_point scale_color_manual theme_classic scale_x_continuous element_blank element_text theme labs scale_size guide_legend guide_colorbar unit alpha
#' @importFrom aplot insert_left
#' @importFrom Seurat DotPlot
#' @importFrom paletteer paletteer_d
#'
#' @export
dotPlotMisc <- function(object,
                      features,
                      assay = NULL,
                      scale = T,
                      dot.range.min = 0,
                      dot.range.max = 3.5,
                      Combine=F,
                      legend.position = "right",
                      label.size = 4,
                      label_widths = 0.1,
                      x.lab = NULL,
                      y.lab = NULL,
                      title = NULL,
                      text.size = 8,
                      text.angle = 90,
                      text.vjust = 0.5,
                      text.hjust = 1,
                      group.by = NULL,
                      color.use = NULL,
                      cols = c("lightgrey",
                               "blue"),
                      legend.key.size = 0.5,
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                      idents = NULL, split.by = NULL, cluster.idents = FALSE,
                      scale.by = "radius", scale.min = NA, scale.max = NA,
                      ...
){
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(aplot)
  })
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"),
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust,
                                              vjust = text.vjust), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )

  if(!is.null(group.by)){
    Idents(object) = group.by
  }

  if(Combine){
    if(is.null(color.use)){
      color.use <- alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)
    }
    df <- data.frame(x = 0, y = levels(object), stringsAsFactors = F )
    df$y <- factor(df$y, levels = df$y )
    p1 <- ggplot(df, aes(x, y, color = factor(y))) +
      geom_point(size = label.size, show.legend = F) +
      scale_color_manual(values = color.use) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) + mytheme +
      theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 0,color="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  }

  p2 = DotPlot(object = object,cols = cols,
               assay = assay,
               col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale,
               idents = idents,split.by = split.by, cluster.idents = cluster.idents,
               scale.by = scale.by, scale.min = scale.min, scale.max = scale.max,
               features = features,scale = scale,...)+theme_bw()+
    mytheme + labs(x = x.lab,y = y.lab,title = title)  +
    scale_size(range = c(dot.range.min,dot.range.max))+
    theme(legend.key.size = unit(legend.key.size, "cm"))+
    guides(size = guide_legend(title = "Per.Exp"),
           color = guide_colorbar(title = paste("Aver.Exp.","Scaled",sep = "\n")))
  if(Combine){
    p3 = p2 %>% insert_left(p1, width=label_widths)
    return(p3)
  }else{
    return(p2)
  }
}
