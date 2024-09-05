library("rsconnect")
library("here")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appDir = here("code", "16_iSEE", "02_snRNAseq"),
    appFiles = c("app.R", "sce.rds", "sn_colors.rds", "bulk_colors.rds"),
    appName = "habenulaPilot_snRNAseq",
    account = "libd",
    server = "shinyapps.io"
)
