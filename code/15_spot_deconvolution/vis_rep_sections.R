library(patchwork)

rep_sections_tb <- read.csv(here("processed-data", "05_spe_correct_cluster", "22_SpD_clean_plots", "rep_section.csv")) |>
    filter(rep_section)

vis_rep_sections <- function(spe, 
                             geneid, 
                             assayname = "logcounts",
                             point_size = 1){

    ct_plots <- map(c("AA", "EA"), function(anc) {
        samples <- rep_sections_tb |> filter(Ancestry == anc) |> arrange(APOE)
        
        cluster_row_plots <- map(samples$sample_id, function(s) {
            vis_clus_plot <- vis_gene(
                spe = spe,
                geneid = geneid,
                assayname = assayname,
                point_size = point_size,
                sampleid = s,
                cont_colors = viridisLite::rocket(10, direction = -1)
            ) 
            return(vis_clus_plot)
        })
        cluster_row <- Reduce("+", cluster_row_plots) + plot_layout(nrow = 1)
        return(cluster_row)
    })
    
    ct_grid <- Reduce("/", ct_plots)
    return(ct_grid)
}