## Louise Huuki-Myers, June 2025
## models for dreamlet DGE

dreamlet_models_sn <- list(
    carrier = ~ (1 | APOE_carrier) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    apoe = ~ (1 | APOE) + Anc_Afr + (1 | Sex) + Age + Rin + (1 | exp_round),
    e4e4 = ~ (1 | APOE_E4E4) + Anc_Afr  + (1 | Sex) + Age + Rin + (1 | exp_round) ,
    # interaction syntax x + (x | g)
    carrier_i = ~ Anc_Afr + (Anc_Afr | APOE_carrier) + (1 | Sex) + Age + Rin + (1 | exp_round),
    apoe_i = ~ Anc_Afr + (Anc_Afr | APOE) + (1 | Sex) + Age + Rin + (1 | exp_round),
    e4e4_i = ~ Anc_Afr + (Anc_Afr | APOE_E4E4) + (1 | Sex) + Age + Rin + (1 | exp_round)
)