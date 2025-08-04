#   Create boxplots by sample of each type of nuclear count

library(here)
library(tidyverse)
library(sessioninfo)

sample_1_path = here(
    'code', '01_spaceranger', 'spaceranger_parameters_1v-16v.txt'
)
sample_2_path = here(
    'code', '01_spaceranger', 'spaceranger_parameters_17v-31v_trimmed.txt'
)
count_path = here('processed-data', '16_nuclear_counts', '%s.csv.gz')
plot_dir = here('plots', '16_nuclear_counts')

sample_boxplot = function(count_df, var_name, filename) {
    #   Get a sensible upper bound for the y-axis
    upper_bound = count_df |>
        group_by(sample_id) |>
        summarize(
            q1 = quantile(!!sym(var_name), 0.25),
            q3 = quantile(!!sym(var_name), 0.75),
            whisker_max = 2.5 * q3 - 1.5 * q1
        ) |>
        pull(whisker_max) |>
        max()
    
    p = ggplot(count_df, aes(x = sample_id, y = !!sym(var_name))) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        coord_cartesian(ylim = c(0, upper_bound)) +
        labs(x = 'Sample ID')
    
    pdf(file.path(plot_dir, filename), width = 20, height = 5)
    print(p)
    dev.off()
}

sample_ids = rbind(
    read_csv(sample_1_path, col_names = FALSE, show_col_types = FALSE),
    read_csv(sample_2_path, col_names = FALSE, show_col_types = FALSE)
)[[1]]

#   Read in nuclear counts for all samples as one tibble
count_df_list = list()
for (sample_id in sample_ids) {
    count_df_list[[sample_id]] = sprintf(count_path, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)
}
count_df = do.call(rbind, count_df_list)

#   Boxplots for both types of nuclear counts
sample_boxplot(count_df, 'num_nuclei_within', 'nuclei_within_boxplot.pdf')
sample_boxplot(count_df, 'num_nuclei_intersect', 'nuclei_intersect_boxplot.pdf')

session_info()
