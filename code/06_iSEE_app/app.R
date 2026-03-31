library("SingleCellExperiment")
library("iSEE")
library("shiny")

## Load data
sce <- readRDS("sce_ERC_iSEE.rds")
sn_colors <- readRDS("sn_colors.rds")

sce <- registerAppOptions(sce, color.maxlevels = 38)

## iSEE configuration
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, YAxis = 2L, 
                                          FacetRowByColData = "sample_id", FacetColumnByColData = "sample_id", 
                                          ColorByColumnData = "cell_type_anno", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_id", 
                                          SizeByColumnData = "sum", TooltipColumnData = character(0), 
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "DDX11L2", 
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                          ColorBySampleName = "1", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                          ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                          ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE, 
                                          VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                          FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                          PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                          CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1, 
                                          LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                          LabelCenters = FALSE, LabelCentersBy = "sample_id", LabelCentersColor = "#000000", 
                                          VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version", 
                                                                                                              "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                          SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "GFAP", Search = "", SearchColumns = c("", 
                                                                                                    "", "", "", "", "", "", "", "", "", "", "", "", ""), HiddenColumns = character(0), 
                                  VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version", 
                                                                                                      "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L, 
                                  PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "cell_type_anno", XAxisFeatureName = "DDX11L2", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "APOE", YAxisFeatureSource = "RowDataTable1", 
                                      YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "sample_id", 
                                      FacetColumnByColData = "sample_id", ColorByColumnData = "cell_type_anno", 
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                      ShapeByColumnData = "sample_id", SizeByColumnData = "sum", 
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                      ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                      ColorByFeatureName = "DDX11L2", ColorByFeatureSource = "---", 
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "1", 
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE, 
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                      FixAspectRatio = FALSE, ViolinAdd = TRUE, PointSize = 1, 
                                      PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "1", FontSize = 1, 
                                      LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                      LabelCenters = FALSE, LabelCentersBy = "sample_id", LabelCentersColor = "#000000", 
                                      VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version", 
                                                                                                          "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L), 
                                      PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                      RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                      SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "AQP4\nCD68\nPDGFRA\nMOBP\nCLDN5\nSLC17A7\nGAD1", 
                                        CapRowSelection = 200L, ClusterRows = FALSE, ClusterRowsDistance = "spearman", 
                                        ClusterRowsMethod = "ward.D2", DataBoxOpen = FALSE, VisualChoices = "Annotations", 
                                        ColumnData = "cell_type_anno", RowData = character(0), CustomBounds = FALSE, 
                                        LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE, 
                                        AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow", 
                                        ShowDimNames = "Rows", LegendPosition = "Bottom", LegendDirection = "Horizontal", 
                                        VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10, 
                                        ShowColumnSelection = TRUE, OrderColumnSelection = TRUE, 
                                        VersionInfo = list(iSEE = structure(list(c(2L, 20L, 0L)), class = c("package_version", 
                                                                                                            "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 6L, 
                                        SelectionBoxOpen = FALSE, RowSelectionSource = "RowDataTable1", 
                                        ColumnSelectionSource = "ReducedDimensionPlot1", RowSelectionDynamicSource = FALSE, 
                                        ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                        ColumnSelectionRestrict = FALSE, SelectionHistory = list())
## Build the iSEE app
# Warning in .coerce_type(v, f, max_levels = max_levels, df = df) :
#     covariate has too many unique values, coercing to numeric

iSEE(sce,
        appTitle = "LFF_ERC_snRNA-seq", 
        initial = initial,
        colormap = ExperimentColorMap(
            colData = list(
                 cell_type_anno= function(n) {
                    return(sn_colors$anno)
                } ,
                cell_type_broad = function(n) {
                    return(sn_colors$broad)
                }
            )
        )
    )
