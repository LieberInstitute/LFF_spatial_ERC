## Louise Huuki-Myers, Feb 2026
## Calculate local ancestry of ERC donors
## Original code from Bernie Mulvey 2025-02-21 
## https://github.com/LieberInstitute/LFF_spatial_LC/blob/21567b9445291adf9f6d35ffc22900c8c331970d/code/02_build_spe/02-calc_APOE_haplolocus_ancestry.Rmd#L19

#### Set up ####
library("data.table")
library("tidyverse")
library("here")
library("sessioninfo")

plot_dir <- here("plots", "00_project_prep", "10_local_ancestry")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

data_dir <- here("processed-data", "00_project_prep", "10_local_ancestry")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### Load local ancestry data ####
## load chr19 local ancestry estimates from /dcs04/lieber/shared/statsgen/local_ancestry/flare_array_data/merged-R.9_local_ancestry_chr19.anc.vcf.gz

sample_info <- read.csv(here('processed-data', '00_project_prep', '05_pathology', 'sample_taupathy.csv')) |>
    dplyr::filter(BrNum != "Br1289")

brs <- sample_info$BrNum

anc <- fread("/dcs04/lieber/shared/statsgen/local_ancestry/flare_array_data/merged-R.9_local_ancestry_chr19.anc.vcf.gz")
keepn <- c("POS","REF","ALT",brs)
anc <- anc[,..keepn]
gc(full=T)

## coords in the vcf are in hg38. based on Aguet 2023 guide to molQTLs in Nat Reviews Methods Primers, ± 1Mb is a standard window for local ancestry estimates. since APOE is defined by two SNPs very close to each other, get the 1MB window around the pair (in hg38 coords): 44908684-1mil to 44908822+1mil
anc <- anc[POS>=(44908684-1000000)&POS<=(44908822+1000000)]

## extract haplotype 1 and haplotype 2 into separate columns for each donor
posn <- c("POS",brs)
anc1 <- copy(anc[,..posn])
anc1 <- melt.data.table(anc1,id.vars="POS")
anc1[,value:=as.numeric(gsub(value,pattern="^.\\|.:(.):.$",replacement="\\1"))]
## drop the occassional CHB genotypes (treat as 0)
anc1[value==2,value:=0]
hap1 <- anc1[,mean(value),by="variable"]
setnames(hap1,"V1","apoe_hap1_estPropYRI")

anc2 <- copy(anc[,..posn])
anc2 <- melt.data.table(anc2,id.vars="POS")
anc2[,value:=as.numeric(gsub(value,pattern="^.\\|.:.:(.)$",replacement="\\1"))]
## drop the occassional CHB genotypes (treat as 0)
anc2[value==2,value:=0]
hap2 <- anc2[,mean(value),by="variable"]
setnames(hap2,"V1","apoe_hap2_estPropYRI")

## get haplotype 1 APOE isoform and hap2 APOE isoform
# 44908684:rs429358
# 44908822:rs7412

# rs42 rs74 apoe
# T 	T 	ε2
# T 	C 	ε3 
# C 	C	ε4

hapg1 <- copy(anc)
hapg1.42 <- hapg1[POS==44908684]
hapg1.42 <- melt.data.table(hapg1.42,id.vars=c("POS","REF","ALT"))
hapg1.42[,hap1gt42:=gsub(value,pattern="^(.)\\|.*$",replacement="\\1")]
hapg1.42[hap1gt42=="0",gt42:=REF]
hapg1.42[hap1gt42=="1",gt42:=ALT]

hapg1.74 <- hapg1[POS==44908822]
hapg1.74 <- melt.data.table(hapg1.74,id.vars=c("POS","REF","ALT"))
hapg1.74[,hap1gt74:=gsub(value,pattern="^(.)\\|.*$",replacement="\\1")]
hapg1.74[hap1gt74=="0",gt74:=REF]
hapg1.74[hap1gt74=="1",gt74:=ALT]

hapg1 <- merge(hapg1.42[,.(variable,gt42)],hapg1.74[,.(variable,gt74)],by="variable")
hapg1[,hap1gt:=paste0(gt42,gt74)]

## haplotype 2 APOE isoform
hapg2 <- copy(anc)
hapg2.42 <- hapg2[POS==44908684]
hapg2.42 <- melt.data.table(hapg2.42,id.vars=c("POS","REF","ALT"))
hapg2.42[,hap2gt42:=gsub(value,pattern="^.\\|(.):.*$",replacement="\\1")]
hapg2.42[hap2gt42=="0",gt42:=REF]
hapg2.42[hap2gt42=="1",gt42:=ALT]

hapg2.74 <- hapg2[POS==44908822]
hapg2.74 <- melt.data.table(hapg2.74,id.vars=c("POS","REF","ALT"))
hapg2.74[,hap2gt74:=gsub(value,pattern="^.\\|(.):.*$",replacement="\\1")]
hapg2.74[hap2gt74=="0",gt74:=REF]
hapg2.74[hap2gt74=="1",gt74:=ALT]

hapg2 <- merge(hapg2.42[,.(variable,gt42)],hapg2.74[,.(variable,gt74)],by="variable")
hapg2[,hap2gt:=paste0(gt42,gt74)]

## merge haplotype 1 and haplotype 2 APOE isoform
hapg <- merge(hapg1[,.(variable,hap1gt)],hapg2[,.(variable,hap2gt)],by="variable")
## append local ancestry estimates per haplotype
hapg <- merge(hapg,hap1[,.(variable,apoe_hap1_estPropYRI)],by="variable")
hapg <- merge(hapg,hap2[,.(variable,apoe_hap2_estPropYRI)],by="variable")

hapg[hap1gt=="TT",hap1gt:="E2"]
hapg[hap1gt=="TC",hap1gt:="E3"]
hapg[hap1gt=="CC",hap1gt:="E4"]

hapg[hap2gt=="TT",hap2gt:="E2"]
hapg[hap2gt=="TC",hap2gt:="E3"]
hapg[hap2gt=="CC",hap2gt:="E4"]

## make sure the genotypes we're getting here agree with the source metadata
check <- merge.data.table(unique(sample_info[,c("BrNum","APOE")]),hapg,by.x="BrNum",by.y="variable")
## looks good.

hapg <- sample_info |> 
    select(BrNum, APOE, Ancestry, Anc_Afr) |> 
    left_join(
        as.data.frame(hapg) |>
            dplyr::rename(BrNum = variable) )


write.csv(hapg, file = here(data_dir, "ERC_apoe_haplotypes_w_localYRI_est.txt"))

hapg |> arrange(Anc_Afr)

#### plot local ancestries ####
load(here("processed-data", "project_colors.Rdata"), verbose = TRUE)

hapg_anc_long <- hapg |>
    select(-starts_with("hap")) |>
    pivot_longer(!c("BrNum" ,"APOE" ,"Ancestry" ,"Anc_Afr"),
                 names_to = "hap", values_to = "estPropYRI") |>
    mutate(hap = gsub("apoe_|_estPropYRI", "", hap)) |> 
    left_join(hapg |>
                  select(-starts_with("apoe_hap")) |>
                  pivot_longer(!c("BrNum" ,"APOE" ,"Ancestry" ,"Anc_Afr"),
                               names_to = "hap", values_to = "APOE_gt") |>
                  mutate(hap = gsub("gt", "", hap))
    ) |>
    mutate(BrNum = fct_reorder(BrNum, -Anc_Afr))


local_ancestry_tile <- hapg_anc_long |> 
    group_by(APOE, Ancestry) |>
    mutate(APOE_anno = paste(APOE, Ancestry, "n=", n()/2)) |>
    ggplot(aes(x = BrNum, y = hap, fill = estPropYRI)) +
    geom_tile(color = "white") +
    geom_text(aes(label = APOE_gt)) +
    facet_wrap(~APOE_anno, scales = "free_x", ncol = 2)   +
    # scale_fill_manual(values = ancestry_colors2) +
    scale_fill_gradient(low = ancestry_colors[["EA"]], high = ancestry_colors[["AA"]]) +
    theme_bw() 

ggsave(local_ancestry_tile, filename = here(plot_dir, "ERC_local_ancestry_tile.png"), width = 8)

