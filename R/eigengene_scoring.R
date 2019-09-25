#' @title scoreEigengenes
#'
#' @param object Either a matrix of expression values, a DESeqDataSet object, or a Seurat object
#' @param module_list A named lists of lists where each sublist are the genes that make up the
#' module and each sublist is named.
#' @param md Optional meta data to add to the scores (for plotting purposes). Default: NULL
#' @param assay Seurat assay object from which to pull data
#' @param slot Assay slot to use
#' @param return_self return a form of `object` with the scores added to it. Default: TRUE
#' @param score_func Factorization method to use in scoring. Currently only 'rsvd', 'svd', and 'nmf' are allowed. Default: 'rsvd'
#' @param ... Not used
#'
#' @return
#' @export
#'
scoreEigengenes <- function(object,...){
  UseMethod("scoreEigengenes")
}


#' @rdname scoreEigengenes
#' @method scoreEigengenes default
#'
#' @importFrom dplyr inner_join intersect
#' @importFrom tibble rownames_to_column column_to_rownames as_tibble
#' @importFrom magrittr %<>% %>%
#'
#' @return
#' @export
scoreEigengenes.default <- function(object,
                                    module_list,
                                    md = NULL,
                                    score_func = "rsvd",
                                    ...){
  scores <- score_matrix(object = object,
                         module_list = module_list,
                         score_func = score_func)
  scores %<>% as.matrix() %>% as_tibble()
  names(scores) <- names(module_list)
  scores[["sample"]] <- colnames(object)
  if (!is.null(md)){
    md %<>%
      as_tibble(rownames = "sample") %>%
      select(
        unique(
          c("sample",
            colnames(md)[which(!colnames(md) %in% colnames(scores))]))) %>%
      inner_join(scores)
    return(md)
  } else {
    return(scores)
  }
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes DESeqDataSet
#'
#' @importFrom DESeq2 vst
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom S4Vectors DataFrame
#'
#' @return
#' @export
scoreEigengenes.DESeqDataSet <- function(object,
                                         module_list,
                                         return_self = TRUE,
                                         score_func = "rsvd",
                                         ...){
  exprs <- vst(object) %>% assay() %>% as.matrix()
  md <- colData(object)
  scores <- scoreEigengenes.default(object = exprs,
                                    module_list = module_list,
                                    md = md,
                                    score_func = score_func)

  if(isTRUE(return_self)){
    scores %<>% DataFrame(row.names = .[["sample"]])
    scores[["sample"]] <- NULL
    colData(object) <- scores
    return(object)
  } else {
    return(scores)
  }
}

#' @title score_matrix
#' @rdname score_matrix
#'
#' @param object matrix of gene expression to assign an eigenvalue
#' @param module_list list of gene modules
#' @param score_func function to use to assign an eigenvalue. Functions allowed include 'svd', 'rsvd', and 'nmf'
#'
#' @importFrom rsvd rsvd
#' @importFrom NMF nmf
#' @importFrom furrr future_map_dfc
#' @importFrom purrr map_dfc
#' @importFrom dplyr intersect
#'
#' @return
#' @export
#'
#' @examples
score_matrix <- function(object,
                         module_list,
                         score_func = "rsvd"){
  switch(
    score_func,
    "svd" = future_map_dfc(names(module_list), function(j) {
      modgenes <- intersect(module_list[[j]], rownames(object))
      if(length(modgenes) > 1){
        exprDat <- object[modgenes,]
        expr <- svd(x = exprDat,
                     nv = 1,
                     nu = 1)
        expr <- expr$v
        return(expr)
      } else {
        return(matrix(rep(0,ncol(object))))
      }
    }),
    "rsvd" = future_map_dfc(names(module_list), function(j) {
      modgenes <- intersect(module_list[[j]], rownames(object))
      if(length(modgenes) > 1){
        exprDat <- object[modgenes,]
        expr <- rsvd(A = exprDat,
                     k = 1)
        expr <- expr$v
        return(expr)
      } else {
        return(matrix(rep(0,ncol(object))))
      }
    }),
    "nmf" = map_dfc(names(module_list), function(j) {
      modgenes <- intersect(module_list[[j]], rownames(object))
      if(length(modgenes) > 1){
        exprDat <- object[modgenes,]
        expr <- nmf(exprDat, rank=1) %>%
          slot("fit") %>%
          slot("H") %>%
          t()
        as.matrix(expr)
      } else {
        matrix(rep(0,ncol(object)))
      }
    })
    )
  }
