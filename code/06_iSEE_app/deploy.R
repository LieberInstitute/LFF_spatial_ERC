library("rsconnect")
library("here")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appDir = here("code", "06_iSEE_app"),
    appFiles = c("app.R", "sce_ERC_iSEE.rds"),
    appName = "habenulaPilot_snRNAseq",
    account = "libd",
    server = "shinyapps.io"
)
