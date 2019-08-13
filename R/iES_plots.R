########## The main demonstration plots #############
water_fall = function(iES_mat, gs_str, title = T){
    group_colors = c(tumor = "Brown", normal = "#56B4E9")
    gs_ind = which(rownames(iES_mat[[1]]) == gs_str)
    normal = iES_mat[[1]][gs_ind,]; tumor = iES_mat[[2]][gs_ind, ]

    n_gap = round((length(tumor) + length(normal)) * 0.01)
    tmp_iES_mat = iES_mat.frame(value = c(tumor[order(tumor)], rep(0, n_gap), normal[order(normal)]),
                          type = c(rep("tumor", length(tumor)), rep("tumor", n_gap), rep("normal", length(normal))),
                          fill = c(rep("fill", length(tumor)), rep(NA, n_gap), rep("fill", length(normal))))
    if (title ==T){
        p = ggplot(tmp_iES_mat, aes(x = 1:nrow(tmp_iES_mat), y = value)) +
            geom_area(aes(fill = type)) +
            theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            labs(x="Samples", y="Enrichment Score") +
            scale_fill_manual(values=group_colors) +
            ggtitle(gs_str)
    }else{
        p = ggplot(tmp_iES_mat, aes(x = 1:nrow(tmp_iES_mat), y = value)) +
            geom_area(aes(fill = type)) +
            theme(legend.position="top", legend.direction="horizontal", panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            labs(x="Samples", y="Enrichment Score") +
            scale_fill_manual(values=group_colors)
    }
    return(p)
}


density_fall = function(iES_mat, gs_str, title = T){
    gs_ind = which(rownames(iES_mat[[1]]) == gs_str)
    normal = iES_mat[[1]][gs_ind,]; tumor = iES_mat[[2]][gs_ind, ]
    tmp_m = Mclust(normal, parameter = T, modelNames = "V")
    id = which.max(tmp_m$parameters$pro)
    tmp_mean = tmp_m$parameters$mean[id]
    tmp_iES_mat = iES_mat.frame(value = c(normal, tumor),
                          type = c(rep("normal", length(normal)), rep("tumor", length(tumor))))
    if (title ==T){
        p = ggdensity(tmp_iES_mat, x = "value",rug = TRUE,
                      color = "type", fill = "type",
                      palette = c("#56B4E9", "Brown"),
                      main = gs_str, legend.title = "") +
            geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
    }else{
        p = ggdensity(tmp_iES_mat, x = "value", rug = TRUE,
                      color = "type", fill = "type",
                      palette = c("#56B4E9", "Brown"),
                      legend.title = "")+
            geom_vline(xintercept= c(tmp_mean, mean(tumor)), linetype="dashed", color = c("#56B4E9", "Brown"))
    }
    return(p)

}


## One survival plot
surv_one = function(GSDB, iES_mat, cli, gs_str, samp_thre=9, title = T){
    data1 = iES_mat$normal_patients; data2 = iES_mat$tumor_patients
    rem_ids = as.vector(which(unlist(apply(data1, 1, sd))!=0))
    tum_samps = Reduce(intersect, list(colnames(iES_mat$tumor_patients), cli$bcr_patient_barcode))
    cli = cli[cli$bcr_patient_barcode %in% tum_samps, ]
    data2 = data2[, tum_samps]
    normal = data1[gs_str,]
    tumor = data2[gs_str,]
    tmp_m = Mclust(normal, parameter = T, modelNames = "V")
    id = which.max(tmp_m$parameters$pro)
    tmp_mean = tmp_m$parameters$mean[id]
    tmp_sd = sqrt(tmp_m$parameters$variance$sigmasq[id])
    # UP Regulate
    if (tmp_mean < mean(tumor)){
        act = names(tumor)[which(tumor > tmp_mean + 2*tmp_sd)]
    }else{
        act = names(tumor)[which(tumor < tmp_mean - 2*tmp_sd)]
    }
    non = setdiff(names(tumor), act)
    if (length(act) > samp_thre & length(non) > samp_thre){
        new_cli = rbind(cli[cli$bcr_patient_barcode %in% act,],
                             cli[cli$bcr_patient_barcode %in% non,])
        new_cli$activition = c(rep("act", length(act)), rep("non", length(non)))
        tmp_cox = coxph(Surv(times, patient.vital_status) ~ activition, data=new_cli, method = "breslow")
        nacts = length(act)
        gs_len = length(GSDB$genesets[[which(GSDB$geneset.names==gs_str)]])
        c = round(summary(tmp_cox)$concordance[1], 4)
        coef = round(tmp_cox$coefficients, 4)
        pval = summary(tmp_cox)$sctest["pvalue"]
        gs_gs= paste0(GSDB$genesets[[which(GSDB$geneset.names==gs_str)]], collapse = "|")
        if (pval < 0.05){
            sfit <- survfit(Surv(times, patient.vital_status)~activition, data = new_cli, conf.type="log-log")
            if (title ==T){
                p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=F,
                               legend.labs=c("perturbed", "normal-like"), legend.title="",
                               palette=c("darkred", "darksalmon"),
                               linetype = "strata",
                               title= paste0("Kaplan-Meier Curve for ", gs_str), data = new_cli)
            }else{
                p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=F,
                               legend.labs=c("perturbed", "normal-like"), legend.title="",
                               palette=c("darkred", "darksalmon"),
                               linetype = "strata", data = new_cli)
            }
        }
        else{
            p = NULL
        }
    }
    sur_rslt = data.frame(gs_str = gs_str, nacts = nacts,
                          gs_len = gs_len, ci=c,
                          coef = coef,
                          pval = as.numeric(as.character(pval)),
                          genes = gs_gs)
    return(list(surv = sur_rslt, one_plot = p))
}

