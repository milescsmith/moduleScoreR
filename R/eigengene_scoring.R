#' @title score_module_eigengenes
#'
#' @param object Either a matrix of expression values of a DESeqDataSet
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
score_module_eigengenes <- function(object,...){
  UseMethod("score_module_eigengenes")
}

#' @rdname score_module_eigengenes
#' @method score_module_eigengenes default
#'
#' @importFrom rsvd rsvd
#' @importFrom furrr future_map_dfc
#' @importFrom dplyr inner_join intersect
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %<>% %>%
#'
#' @return
#' @export
score_module_eigengenes.default <- function(object,
                                            module_list,
                                            md = NULL){
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
    md %<>% as.data.frame() %>% rownames_to_column('sample')
    scores %<>% inner_join(md)
  }
  return(scores)
}

#' @rdname score_module_eigengenes
#' @method score_module_eigengenes DESeqDataSet
#'
#' @importFrom DESeq2 vst
#' @importFrom SummarizedExperiment assay colData
#'
#' @return
#' @export
score_module_eigengenes.DESeqDataSet <- function(object, module_list){
  exprs <- vst(object) %>% assay() %>% as.matrix()
  md <- colData(object)
  scores <- score_module_eigengenes.default(object = exprs, module_list = module_list, md = md)
  return(scores)
}

#' @rdname score_module_eigengenes
#' @method score_module_eigengenes Seurat
#'
#' @importFrom Seurat FetchData
#' @importFrom glue glue
#'
#' @return
#' @export
score_module_eigengenes.Seurat <- function(object,
                                           module_list,
                                           assay = "RNA",
                                           slot = "data"){
  scores <- future_map_dfc(names(module_list), function(j) {
    i <- glue("{tolower(assay)}_{module_list[[j]]}") %>% as.character()
    exprDat <- FetchData(object = object,
                         vars = i,
                         slot = slot) %>%
      t()
    eigen <- rsvd(exprDat, k = 1)
    eigen[["v"]]
  }) %>%
    as.matrix()
  rownames(scores) <- colnames(object)
  colnames(scores) <- names(module_list)
  #object[["eigengenes"]] <- CreateDimReducObject(embeddings = scores,
  #                                            key = "eigen_")
  return(scores)
}
