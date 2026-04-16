## Louise Huuki-Myers, April 2026
## prep xeniumranger run

library(here)

list.files(here("raw-data", "xenium"))


slurmjobs::job_single('07_xenium_ranger', create_shell = TRUE, memory = '25G', command = "xeniumranger")
