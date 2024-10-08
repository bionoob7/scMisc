% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vlnplot.R
\name{vlnPlotStacked}
\alias{vlnPlotStacked}
\title{Stacked Violin Plot for Seurat Object}
\usage{
vlnPlotStacked(
  object,
  features,
  group,
  text.size = 8,
  text.angle = 45,
  text.hjust = 1,
  legend.position = "right",
  switch = "left",
  fill.cols = NULL,
  cols = NULL,
  widths = c(3, 0.08),
  heights = c(3, 0.08),
  legend.title = "Ave.exp",
  x.lab = NULL,
  y.lab = NULL,
  title = NULL,
  ...
)
}
\arguments{
\item{object}{A Seurat object.}

\item{features}{A vector of features (genes) to plot.}

\item{group}{Grouping variable.}

\item{text.size}{Size of the text. Defaults to 8.}

\item{text.angle}{Angle of the text. Defaults to 45.}

\item{text.hjust}{Horizontal justification of the text. Defaults to 1.}

\item{legend.position}{Position of the legend. Defaults to "right".}

\item{switch}{Switch position for facet labels. Defaults to "left".}

\item{fill.cols}{Colors for violin plot fill. Defaults to NULL.}

\item{cols}{Colors for cell type annotation. Defaults to NULL.}

\item{widths}{Widths for plot layout. Defaults to c(3, 0.08).}

\item{heights}{Heights for plot layout. Defaults to c(3, 0.08).}

\item{legend.title}{Title for the legend. Defaults to "Ave.exp".}

\item{x.lab}{Label for the x-axis. Defaults to NULL.}

\item{y.lab}{Label for the y-axis. Defaults to NULL.}

\item{title}{Title of the plot. Defaults to NULL.}

\item{...}{Additional parameters passed to FetchData.}
}
\value{
A ggplot object.
}
\description{
This function generates a stacked violin plot for a Seurat object, allowing for
customization of plot aesthetics and layout.
}
\examples{
\dontrun{
vlnPlotStacked(object = seurat.data,
               features = marker_gene,
               group = "celltype",
               legend.key.size = 0.3,
               combine = TRUE, stack = TRUE, flip = TRUE,
               cols = paletteer::paletteer_d('ggthemes::Classic_20'))
}

}
