library("rsconnect")
library("here")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appDir = here("code", "06_iSEE_app"),
    appFiles = c("app.R", "sce_ERC_iSEE.rds","sn_colors.rds"),
    appName = "LFF_ERC_snRNA-seq",
    account = "libd",
    server = "shinyapps.io"
)
