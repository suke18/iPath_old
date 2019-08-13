########## iES score matrix calculation ##########

rem_data = function(x){
    rem_ids = which(apply(x, 1, sd) == 0)
    if (is_empty(rem_ids)){
        return(x)
    }else{
        return(x[-rem_ids, ])
    }
}
f_colnames = function(x) gsub("\\.","-",substr(x,1,12))

GSEA = function(gene_list, gene_set, stats_vector){
    # Input genelist must be ordered.
    tag_indicator = sign(match(gene_list, gene_set, nomatch = 0))
    no_tag_indicator = 1- tag_indicator
    N = length(gene_list); Nh = length(gene_set); Nm =  N - Nh

    sum_rank_tag = sum(stats_vector[tag_indicator ==1])
    if (sum_rank_tag == 0){
        ES = -1
    }else{
        norm_tag = 1.0/ sum_rank_tag; norm_no_tag = 1/Nm
        RES = cumsum(tag_indicator*stats_vector*norm_tag - no_tag_indicator*norm_no_tag)
        max.ES = max(RES)
        min.ES = min(RES)
        if (max.ES > - min.ES) {
            ES = signif(max.ES, digits=5)
        } else {
            ES = signif(min.ES, digits=5)
        }
    }
    return(ES)
}

#' iES calculation Function
#'
#' This function allows you to express your love of cats.
#' @param Y is the expression matrix, GSDB is the gene set database
#' @keywords iES statistics calculatioin
#' @export
#' @examples
#' iES_cal()

iES_cal = function(Y, GSDB){
    if (is(Y, "ExpressionSet")){
        pheno = phenoData(Y)
        n_index = which(pheno$phenotypes == "normal")
        t_index = which(pheno$phenotypes == "tumor")
        Y = exprs(Y)
    }else{Y = as.matrix(Y)
        n_index = which(substr(colnames(Y),14,14) == "1")
        t_index = which(substr(colnames(Y),14,14) == "0")
        colnames(Y) = f_colnames(colnames(Y))
    }
    Y = rem_data(Y)
    row_mean = rowMeans(Y)
    row_sd = rowSds(Y)
    tmp_norm = apply(Y, 2, function(x) abs((x - row_mean)/(row_sd)))

    n1 = length(GSDB$genesets); n2 = ncol(tmp_norm)
    es_Y = matrix(rep(0, n1*n2), ncol = n2)

    order_array = apply(tmp_norm, 2, function(x) rev(order(x)))
    for (j in 1:n1){
        tmp_gs = GSDB$genesets[[j]]
        gs_score = sapply(1:n2, function(i){
            tmp_glist = tmp_norm[,i][order_array[,i]]
            return(GSEA(names(tmp_glist), tmp_gs, as.numeric(tmp_glist)))
        })
        es_Y[j,] = gs_score
    }
    dimnames(es_Y) = list(GSDB$geneset.names, colnames(Y))
    nor_patients = es_Y[, n_index]; tum_patients = es_Y[, t_index]
    return(list(normal_patients = nor_patients, tumor_patients = tum_patients))
}





