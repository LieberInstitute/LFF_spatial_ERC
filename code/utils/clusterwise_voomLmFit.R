
clusterwise_voomLmFit <- function(sce_pb, clus, batch = "exp_round"){

        dge <- sce_pb[,sce_pb$registration_variable == clus]
        
        des <- model.matrix(~0 + APOE_syn + Sex + Age + Anc_Afr + pseudo_expr_chrM_ratio, data = colData(dge))
        des <- as.data.frame(des)
        
        # filter low expression genes
        dge <- edgeR::calcNormFactors(dge)
        keep <- edgeR::filterByExpr.DGEList(dge,design=des)
        dge <- dge[keep,,keep.lib.sizes=FALSE]
        dge <- edgeR::calcNormFactors(dge)
        
        message(Sys.time(), sprintf(" - voomLmFit - cluster: %s, block= '%s', ncol: %s, ngene: %i", clus, batch, ncol(dge), nrow(dge$genes)))
        
        # make these more readable
        colnames(des) <- gsub(colnames(des),pattern="_syn",replacement="_")
        
        ## run voomLmFit for the pseudobulked data, referring donor to duplicateCorrelation; 
        ## using an adaptive span (number of genes, based on the number of genes in the dge) for smoothing the mean-variance trend
        v.swt <- voomLmFit(dge,design = des,block = as.factor(dge$samples[[batch]]),adaptive.span = T,sample.weights = T)
        
        cont <- makeContrasts(
            ## main
            carrier = "-0.5*(APOE_E2.E2 + APOE_E2.E3) + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
            # E4E4 = "-APOE_E4.E4 + (APOE_E2.E2 + APOE_E2.E3 + APOE_E3.E4)/3",
            # ## apoe pairwise
            # apoe_E2E2_E4E4 = "-APOE_E2.E2 + APOE_E4.E4",
            # apoe_E3E4_E4E4 = "-APOE_E3.E4 + APOE_E4.E4",
            # apoe_E2E3_E4E4 = "-APOE_E2.E3 + APOE_E4.E4",
            # apoe_E2E2_E3E4 = "-APOE_E2.E2 + APOE_E3.E4",
            # apoe_E2E2_E2E3 = "-APOE_E2.E2 + APOE_E2.E3",
            # apoe_E2E3_E3E4 = "-APOE_E2.E3 + APOE_E3.E4",
            # # heterozygous vs. homozygous
            # anyE2_E4E4 = "- 0.5*(APOE_E2.E3 + APOE_E2.E2) + APOE_E4.E4",
            # E2E2_anyE4 = "-APOE_E2.E2 + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
            # E2E3_anyE4 = "-APOE_E2.E3 + 0.5*(APOE_E3.E4 + APOE_E4.E4)",
            # ## other
            # Sex="SexM",
            # Anc="Anc_Afr",
            levels=des
        )
        
        v.swt.fit <- contrasts.fit(v.swt,contrasts=cont)
        v.swt.fit.e <- eBayes(v.swt.fit)
        
        ## run top table over contrasts
        v.swt.e.tt <- topTable(v.swt.fit.e, coef = 'carrier', number=Inf, adjust.method = "BH") |>
                                     mutate(cluster = clus,
                                            contrast = "carrier", 
                                            .before = 1) |>
                                     arrange(adj.P.Val) 
        
        message("Done - Save data")
        saveRDS(v.swt.e.tt, file = here(data_dir, sprintf("voomLmFit_%s.rds", clus)))
        return(v.swt.e.tt$)

}