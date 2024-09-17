library("SingleCellExperiment")
library("iSEE")
library("shiny")

## Load data
sce <- readRDS("sce_ERC_iSEE.rds")
# sn_colors <- readRDS("sn_colors.rds")

## iSEE configuration
initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "TSNE", XAxis = 1L, YAxis = 2L, 
                                          FacetRowByColData = "sample_id", FacetColumnByColData = "sample_id", 
                                          ColorByColumnData = "snn_k15", ColorByFeatureNameAssay = "logcounts", 
                                          ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "sample_id", 
                                          SizeByColumnData = "sum", TooltipColumnData = character(0), 
                                          FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
                                          ColorByDefaultColor = "#000000", ColorByFeatureName = "ENSG00000290825", 
                                          ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
                                          ColorBySampleName = "cell1", ColorBySampleSource = "---", 
                                          ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
                                          SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
                                          VisualBoxOpen = TRUE, VisualChoices = "Color", ContourAdd = FALSE, 
                                          ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
                                          Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
                                          CustomLabelsText = "cell1", FontSize = 1, LegendPointSize = 1, 
                                          LegendPosition = "Bottom", HoverInfo = TRUE, LabelCenters = FALSE, 
                                          LabelCentersBy = "sample_id", LabelCentersColor = "#000000", 
                                          VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version", 
                                                                                                              "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), 
                                          PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                          RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                          DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                          RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                          SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
                                        CustomRowsText = "POU4F1\nGPR151\nCHRNB4\nHTR2C\nLYPD6B\nADARB2\nRORB\nSYT1\nSLC17A6\nGAD1\nMOBP\nPDGFRA\nAQP4\nITIH5\nCSF1R\n\n# ONECUT2\n# CRH\n# MCOLN3\n# TLE2\n# SEMA3D\n# ESRP1\n# CCK\n# CHAT\n# EBF3", 
                                        ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
                                        DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "snn_k15", 
                                        RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_, 
                                        UpperBound = NA_real_, AssayCenterRows = TRUE, AssayScaleRows = FALSE, 
                                        DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows", 
                                        LegendPosition = "Right", LegendDirection = "Vertical", VisualBoxOpen = FALSE, 
                                        NamesRowFontSize = 10, NamesColumnFontSize = 10, ShowColumnSelection = TRUE, 
                                        OrderColumnSelection = TRUE, VersionInfo = list(iSEE = structure(list(
                                            c(2L, 16L, 0L)), class = c("package_version", "numeric_version"
                                            ))), PanelId = c(ComplexHeatmapPlot = 1L), PanelHeight = 500L, 
                                        PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "RowDataTable1", 
                                        ColumnSelectionSource = "FeatureAssayPlot1", RowSelectionDynamicSource = FALSE, 
                                        ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                        ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "ENSG00000197971", Search = "MBP", 
                                  SearchColumns = c("", "", "", ""), HiddenColumns = "Type", 
                                  VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version", 
                                                                                                      "numeric_version"))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L, 
                                  PanelWidth = 6L, SelectionBoxOpen = FALSE, RowSelectionSource = "---", 
                                  ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, 
                                  ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE, 
                                  ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
                                      XAxisColumnData = "snn_k15", XAxisFeatureName = "ENSG00000290825", 
                                      XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
                                      YAxisFeatureName = "ENSG00000197971", YAxisFeatureSource = "RowDataTable1", 
                                      YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "sample_id", 
                                      FacetColumnByColData = "sample_id", ColorByColumnData = "snn_k15", 
                                      ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000", 
                                      ShapeByColumnData = "sample_id", SizeByColumnData = "sum", 
                                      TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None", 
                                      ColorBy = "Column data", ColorByDefaultColor = "#000000", 
                                      ColorByFeatureName = "ENSG00000290825", ColorByFeatureSource = "---", 
                                      ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "cell1", 
                                      ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE, 
                                      ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1, 
                                      ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = TRUE, 
                                      VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF", 
                                      PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
                                      CustomLabels = FALSE, CustomLabelsText = "cell1", FontSize = 1, 
                                      LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
                                      LabelCenters = FALSE, LabelCentersBy = "sample_id", LabelCentersColor = "#000000", 
                                      VersionInfo = list(iSEE = structure(list(c(2L, 16L, 0L)), class = c("package_version", 
                                                                                                          "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L), 
                                      PanelHeight = 500L, PanelWidth = 6L, SelectionBoxOpen = FALSE, 
                                      RowSelectionSource = "---", ColumnSelectionSource = "---", 
                                      DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
                                      RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
                                      SelectionHistory = list())
## Build the iSEE app
iSEE(sce,
        appTitle = "LFF_ERC - snRNA-seq", initial = initial,
        # colormap = ExperimentColorMap(
        #     colData = list(
        #         final_Annotations = function(n) {
        #             return(sn_colors)
        #         },
        #         broad_Annotations = function(n) {
        #             return(bulk_colors)
        #         }
        #     )
        # )
    )
