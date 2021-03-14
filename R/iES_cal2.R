
#' iES calculation Function
#'
#' This function allows you to express your love of cats.
#' @param Y is the expression matrix, GSDB is the gene set database
#' @keywords iES statistics calculatioin
#' @export
#' @examples

iES_cal_2 = function(Y, GSDB){
    Y = as.matrix(Y)
    colnames(Y) = f_colnames(colnames(Y))
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
    return(es_Y)
}





