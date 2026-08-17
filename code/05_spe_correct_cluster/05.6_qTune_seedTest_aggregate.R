## Aggregate the 05.5_qTune_seedTest.R array job outputs and check whether the
## q=10 vs q=11 gap (and the neighboring plateau) is stable across seeds or
## within run-to-run MCMC noise.
##
## Run this after all array tasks from 05.5_qTune_seedTest.sh have finished.

library("here")
library("ggplot2")
library("dplyr")

data_dir <- here("processed-data", "05_spe_correct_cluster", "05.5_qTune", "seedTest")

files <- list.files(data_dir, pattern = "^qTune_logliks_SVGm_seed.*\\.csv$", full.names = TRUE)
message("Found ", length(files), " seed-test result files")
stopifnot(length(files) > 0)

read_one <- function(f) {
    df <- read.csv(f, row.names = 1)
    ## pull seed/task back out of the filename for bookkeeping
    m <- regmatches(basename(f), regexec("seed([0-9]+)_task([0-9]+)", basename(f)))[[1]]
    df$seed <- as.integer(m[2])
    df$task <- as.integer(m[3])
    df
}

all_runs <- do.call(rbind, lapply(files, read_one))

## per-q summary across seeds
summ <- all_runs |>
    group_by(q) |>
    summarize(
        n_seeds = n(),
        mean_loglik = mean(loglik),
        sd_loglik = sd(loglik),
        min_loglik = min(loglik),
        max_loglik = max(loglik),
        .groups = "drop"
    )

print(summ, n = 30)

## Is the q=10 vs q=11 difference bigger than the run-to-run SD at those q's?
q10 <- summ$mean_loglik[summ$q == 10]
q11 <- summ$mean_loglik[summ$q == 11]
sd10 <- summ$sd_loglik[summ$q == 10]
sd11 <- summ$sd_loglik[summ$q == 11]
message(sprintf(
    "q10 mean=%.0f (sd=%.0f), q11 mean=%.0f (sd=%.0f), delta=%.0f, delta/pooled_sd=%.2f",
    q10, sd10, q11, sd11, q11 - q10, (q11 - q10) / sqrt((sd10^2 + sd11^2) / 2)
))
## Rule of thumb: |delta| / pooled_sd >> 1 (e.g. > ~2-3) suggests a real, seed-stable
## effect; a ratio near/below 1 suggests the q10-q11 gap seen in the original single
## run is within the noise band and shouldn't be leaned on for choosing q.

## save combined table
write.csv(all_runs, here(data_dir, "qTune_logliks_SVGm_seedTest_combined.csv"), row.names = FALSE)
write.csv(summ, here(data_dir, "qTune_logliks_SVGm_seedTest_summary.csv"), row.names = FALSE)

## plot: individual seed curves (thin) + mean +/- SD ribbon
p <- ggplot(all_runs, aes(x = q, y = loglik)) +
    geom_line(aes(group = seed), alpha = 0.25, color = "steelblue") +
    geom_ribbon(
        data = summ,
        aes(x = q, ymin = mean_loglik - sd_loglik, ymax = mean_loglik + sd_loglik, y = mean_loglik),
        inherit.aes = FALSE, alpha = 0.2
    ) +
    geom_line(data = summ, aes(x = q, y = mean_loglik), inherit.aes = FALSE, size = 1, color = "black") +
    geom_vline(xintercept = c(10, 11), linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
        title = "qTune log-likelihood across seeds (thin blue = individual seeds)",
        subtitle = "black = mean, ribbon = +/-1 SD, dashed red = q10/q11",
        x = "q", y = "log-likelihood"
    ) +
    theme_bw()

ggsave(here(data_dir, "qTune_seedTest_spread.pdf"), p, width = 8, height = 5)
message("Saved plot to ", here(data_dir, "qTune_seedTest_spread.pdf"))
