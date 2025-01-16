library(getopt)
library(sessioninfo)

# Import command-line parameters
spec <- matrix(
    c(
        c("ssn"),
        c("s"),
        rep("1", 1),
        rep("character", 1),
        rep("Add variable description here", 1)
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

session_info()

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
