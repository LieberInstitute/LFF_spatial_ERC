## Louise Huuki-Myers, June 2025
## models for dreamlet DGE

dreamlet_models_sn <- list(
    carrier_n0 = ~ APOE_carrier,
    carrier = ~ APOE_carrier + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round)  + subsets_Mito_percent,
    apoe = ~ APOE + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    e4e4 = ~ APOE_E4E4 + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    contrast = ~0 + APOE_syn  + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    # interaction syntax x + (x | g)
    carrier_i = ~ APOE_carrier*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    apoe_i = ~ APOE*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent,
    e4e4_i = ~ APOE_E4E4*Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round) + subsets_Mito_percent
)

# dreamlet_models_sn <- list(
#     carrier = ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
#     apoe = ~ (1 | APOE) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
#     e4e4 = ~ (1 | APOE_E4E4) + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round) ,
#     # interaction syntax x + (x | g)
#     carrier_i = ~ Anc_Afr + (Anc_Afr | APOE_carrier) + (1 | Sex) + Age + Rin + (1 | exp_round),
#     apoe_i = ~ Anc_Afr + (Anc_Afr | APOE) + (1 | Sex) + Age + Rin + (1 | exp_round),
#     e4e4_i = ~ Anc_Afr + (Anc_Afr | APOE_E4E4) + (1 | Sex) + Age + Rin + (1 | exp_round)
# )