#' @title scoreEigengenes
#'
#' @param object Either a matrix of expression values, a DESeqDataSet object, or a Seurat object
#' @param module_list A named lists of lists where each sublist are the genes that make up the
#' module and each sublist is named.
#' @param md Optional meta data to add to the scores (for plotting purposes). Default: NULL
#' @param assay Seurat assay object from which to pull data
#' @param slot Assay slot to use
#' @param return_self return a form of `object` with the scores added to it. Default: TRUE
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
#' @importFrom rsvd rsvd
#' @importFrom furrr future_map_dfc
#' @importFrom dplyr inner_join intersect
#' @importFrom tibble rownames_to_column column_to_rownames as_tibble
#' @importFrom magrittr %<>% %>%
#'
#' @return
#' @export
scoreEigengenes.default <- function(object,
                                    module_list,
                                    md = NULL,
                                    ...){
  scores <- future_map_dfc(names(module_list), function(j) {
    modgenes <- intersect(module_list[[j]], rownames(object))
    exprDat <- object[modgenes,]
    expr <- rsvd(exprDat, k = 1)
    expr <- expr$v
    expr
  })
  scores %<>% as.matrix() %>% as_tibble()
  names(scores) <- names(module_list)
  scores[["sample"]] <- colnames(object)
  if (!is.null(md)){
    md %<>% as_tibble(rownames = "sample") %>% inner_join(scores)
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
                                         ...){
  exprs <- vst(object) %>% assay() %>% as.matrix()
  md <- colData(object)
  scores <- scoreEigengenes.default(object = exprs, module_list = module_list, md = md)

  if(isTRUE(return_self)){
    scores %<>% DataFrame(row.names = .[["sample"]])
    scores[["sample"]] <- NULL
    colData(object) <- scores
    return(object)
  } else {
    return(scores)
  }
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes Seurat
#'
#' @importFrom Seurat FetchData
#' @importFrom glue glue
#' @importFrom methods slot
#' @importFrom furrr future_map_dfc
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr inner_join
#'
#' @return
#' @export
scoreEigengenes.Seurat <- function(object,
                                   module_list,
                                   return_self = TRUE,
                                   ...){
  intersecting_ml <- map(names(module_list), function(x){
    intersect(module_list[[x]], rownames(object))
  })
  names(intersecting_ml) <- names(module_list)
  module_list <- intersecting_ml[lapply(intersecting_ml, length)>0]
  scores <- future_map_dfc(.x = names(module_list),
                           .progress = TRUE,
                           .f = function(j) {
    exprDat <- FetchData(object = object,
                         vars = module_list[[j]]) %>%
      t()
      eigen <- rsvd(exprDat, k = 1)
      eigen[["v"]]
    }) %>%
    as.data.frame()

  colnames(scores) <- names(module_list)
  scores[["cells"]] <- colnames(object)


  if(isTRUE(return_self)){
    object@meta.data %<>%
      rownames_to_column("cells") %>%
      inner_join(scores) %>%
      column_to_rownames("cells")
    return(object)
  } else {
    return(scores)
  }
}

