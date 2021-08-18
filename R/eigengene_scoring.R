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
#' @importFrom dplyr inner_join intersect mutate
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom tidyselect any_of
#' @importFrom rlang set_names
#'
#' @return
#' @export
scoreEigengenes.default <-
  function(
    object,
    module_list,
    md = NULL,
    score_func = "rsvd",
    ...){

  scores <-
    score_matrix(
      object = object,
      module_list = module_list,
      score_func = score_func
      )

  scores <-
    as.matrix(scores) |>
    tibble::as_tibble() |>
    rlang::set_names(nm=names(module_list)) |>
    dplyr::mutate(sample = colnames(object))

  if (!is.null(md)){
    out <-
      tibble::as_tibble(
        x = md,
        rownames = "sample"
        ) |>
      dplyr::select(
        -tidyselect::any_of(colnames(scores)),
        sample
        ) |>
      dplyr::inner_join(scores)
  } else {
    out <- scores
  }

  out
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes DESeqDataSet
#'
#' @importFrom DESeq2 vst rlog
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom tibble column_to_rownames
#'
#' @return DESeqDataSet with scores added to the colData
#' @export
scoreEigengenes.DESeqDataSet <-
  function(
    object,
    module_list,
    return_self = TRUE,
    score_func = "rsvd",
    normalize_func = c("vst", "rlog"),
    ...
    ){

    normalize_func <- match.arg(normalize_func)

    exprs <-
      switch(
        normalize_fun,
        vst  = DESeq2::vst(object),
        rlog = DESeq2::rlog(object)
      ) |>
      SummarizedExperiment::assay() |>
      as.matrix()

    md <- SummarizedExperiment::colData(object)

    scores <-
      scoreEigengenes.default(
        object = exprs,
        module_list = module_list,
        md = md,
        score_func = score_func
        )

    if(isTRUE(return_self)){
      scores <-
        tibble::column_to_rownames(
          scores,
          var = "sample"
        ) |>
        DataFrame()

      colData(object) <- scores
      out <- object
    } else {
      out <- scores
    }

    out
}

#' @rdname scoreEigengenes
#' @method scoreEigengenes DGEList
#'
#' @importFrom edgeR rpkm cpm
#' @importFrom tibble column_to_rownames
#'
#' @return DGEList with module scores added to "samples"
#' list component
#' @export
scoreEigengenes.DGEList <-
  function(
    object,
    module_list,
    return_self = TRUE,
    score_func = "rsvd",
    normalize_func = c("cpm", "rpkm"),
    ...
  ){

    normalize_func = match.arg(normalize_func)

    exprs <- switch(
      normalize_func,
      rpkm = edgeR::rpkm(object, log = TRUE, normalized.lib.sizes = TRUE),
      cpm = edgeR::cpm(object, log = TRUE, normalized.lib.sizes = TRUE)
    )

    md <- object[["samples"]]

    scores <-
      scoreEigengenes.default(
        object      = exprs,
        module_list = module_list,
        md          = md,
        score_func  = score_func
      )

    if(isTRUE(return_self)){
      object[["samples"]] <-
        tibble::column_to_rownames(
          .data = scores,
          var   = "sample"
        )

      out <- object
    } else {
      out <- scores
    }
    out
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
#' @return Matrix of module scores
#' @export
#'
#' @examples
score_matrix <- function(
  object,
  module_list,
  score_func = c("rsvd", "svd", "nmf")
  ){

  score_func <- match.arg(score_func)

  out <- switch(
    score_func,
    svd =
      furrr::future_map_dfc(
        .x = names(module_list),
        .f =  function(j) {
          modgenes <- intersect(module_list[[j]], rownames(object))
          if(length(modgenes) > 1){
            svd(
              x = object[modgenes,],
              nv = 1,
              nu = 1
              ) |>
              purrr::chuck("v")
          } else {
            return(matrix(rep(0,ncol(object))))
          }
        }
      ),
    rsvd =
      furrr::future_map_dfc(
        .x = names(module_list),
        .f = function(j) {
          modgenes <- intersect(module_list[[j]], rownames(object))
          if(length(modgenes) > 1){
            rsvd::rsvd(
              A = object[modgenes,],
              k = 1
              ) |>
              purrr::chuck("v")
          } else {
            return(matrix(rep(0,ncol(object))))
          }
        }
      ),
    nmf =
      purrr::map_dfc(
        .x = names(module_list),
        .f =  function(j) {
          modgenes <- intersect(module_list[[j]], rownames(object))
          if(length(modgenes) > 1){
            NMF::nmf(
              x = object[modgenes,], rank=1
              ) |>
              slot("fit") |>
              slot("H") |>
              t() |>
              as.matrix()
          } else {
            return(matrix(rep(0,ncol(object))))
          }
        }
      )
    )

  out
}
