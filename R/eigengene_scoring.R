#' @title scoreEigengenes
#'
#' @param object Either a matrix of expression values, a DESeqDataSet object, or a Seurat object
#' @param module_list A named lists of lists where each sublist are the genes that make up the
#' module and each sublist is named.
#' @param md Optional meta data to add to the scores (for plotting purposes). Default: NULL
#' @param assay Seurat assay object from which to pull data
#' @param slot Assay slot to use
#'
#' @return
#' @export
#'
#' @examples
scoreEigengenes <- function(object,...){
  UseMethod("scoreEigengenes")
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes default
#'
#' @importFrom rsvd rsvd
#' @importFrom furrr future_map_dfc
#' @importFrom dplyr inner_join intersect
#' @importFrom tibble rownames_to_column
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
  names(scores) <- names(module_list)
  scores %<>% as.data.frame()
  scores$sample <- colnames(object)
  if (!is.null(md)){
    md %<>% as.data.frame() %>% rownames_to_column("sample")
    scores %<>% inner_join(md)
  }
  return(scores)
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes DESeqDataSet
#'
#' @importFrom DESeq2 vst
#' @importFrom SummarizedExperiment assay colData
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
    tmp <- colData(object)
    tmp %>%
      rownames_to_column("sample") %>%
      inner_join(scores) %>%
      column_to_rownames("sample")
    colData(object) <- tmp
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

