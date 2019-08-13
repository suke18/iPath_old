########## survival analysis for two groups ##########

#' iES calculation Function
#'
#' This function allows you to express your love of cats.
#' @param GSDB,iES_mat is the gene set database and iES_mat with tumor and normal
#' @keywords iPath survival analysis for two groups of patients: perturbed and normal-like
#' @export
#' @examples
#' iES_surv()

iES_surv = function(GSDB, iES_mat, cli, gs_thre = 500, samp_thre=9, qval = T){
    gs_lens = sapply(1:length(GSDB$geneset.names), function(i) length(GSDB$genesets[[i]]))
    gs_ids = which(gs_lens < gs_thre)
    data1 = iES_mat$normal_patients; data2 = iES_mat$tumor_patients
    rem_ids = as.vector(which(unlist(apply(data1, 1, sd))!=0))
    gss = intersect(rownames(data1)[rem_ids], GSDB$geneset.names[gs_ids])
    tum_samps = Reduce(intersect, list(colnames(iES_mat$tumor_patients), cli$bcr_patient_barcode))
    cli = cli[cli$bcr_patient_barcode %in% tum_samps, ]
    data2 = data2[, tum_samps]
    rslt_mat = NULL; p_plots = NULL
    for(gs in gss){
        gs_str = gs
        normal = data1[gs,]
        tumor = data2[gs,]
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
            rslt_mat = rbind(rslt_mat, c(gs_str, nacts, gs_len, c, coef, pval, gs_gs))
            if (pval < 0.01){
                sfit <- survfit(Surv(times, patient.vital_status)~activition, data = new_cli, conf.type="log-log")
                p = ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=F,
                               legend.labs=c("Perturbed", "Normal-like"), legend.title="",
                               palette=c("red", "green1"),
                               linetype = "strata",
                               title=paste0("Kaplan-Meier Curve for ", gs_str), data = new_cli)
                p_plots[[gs_str]] = p
            }
        }
    }
    ord_inds = order(as.numeric(as.character(rslt_mat[,6])))
    rslt_mat = rslt_mat[ord_inds, ]
    if (qval == T){
        qval = qvalue(as.numeric(rslt_mat[,6]))$qvalues
        sur_rslt = data.frame(gs_str = rslt_mat[,1], nacts = rslt_mat[,2],
                              gs_len = rslt_mat[,3], ci=rslt_mat[,4],
                              coef = rslt_mat[,5],
                              pval = as.numeric(as.character(rslt_mat[,6])),
                              qval = as.numeric(as.character(qval)),
                              genes = rslt_mat[,7])
    }else{
        sur_rslt = data.frame(gs_str = rslt_mat[,1], nacts = rslt_mat[,2],
                              gs_len = rslt_mat[,3], ci=rslt_mat[,4],
                              coef = rslt_mat[,5],
                              pval = as.numeric(as.character(rslt_mat[,6])),
                              genes = rslt_mat[,7])
    }

    return(list(surv = sur_rslt, plots = p_plots))
}




