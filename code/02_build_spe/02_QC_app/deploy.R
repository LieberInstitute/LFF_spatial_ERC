library("rsconnect")
library("here")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "03_spatialLIBD_app", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())
# getOption("repos") # 'getOption("repos")' replaces Bioconductor standard repositories

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "02_build_spe","02_QC_app"),
    appFiles = c(
        "app.R",
        "spe_raw.rds",
        "spe_qc_anno_clean.csv",
        # withr::with_dir(here("code", "03_spatialLIBD_app"), dir("clusters_BayesSpace", full.names = TRUE)),
        withr::with_dir(here("code", "02_build_spe", "02_QC_app"), dir("www", full.names = TRUE))
    ),
    appName = "LFF_ERC_Visium",
    account = "libd",
    server = "shinyapps.io"
)