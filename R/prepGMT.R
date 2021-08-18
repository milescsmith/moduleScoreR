#' @title prepGMT
#'
#' @description Prepare a GMT file for use in score_module_eigengenes
#'
#' @param pathway_list A data.frame consisting of two columns (`term` and `gene`) or the
#' path to a gmt file
#' @param ... Not used
#'
#' @return
#' @export
#'
prepGMT <- function(pathway_list, ...){
  UseMethod("prepGMT")
}

#' @rdname prepGMT
#' @method prepGMT data.frame
#'
#' @importFrom purrr map
#' @importFrom dplyr pull
#'
#' @return Named list of modules and the genes constituting
#' that module
#' @export
#'
prepGMT.data.frame <- function(pathway_list, ...){

  split(
    x = pathway_list,
    f = pathway_list[["term"]]
    ) |>
  purrr::map(\(x) dplyr::pull(.data = x, var = "gene"))
}

#' @rdname prepGMT
#' @method prepGMT character
#'
#' @importFrom clusterProfiler read.gmt
#' @importFrom glue glue
#'
#' @return Named list of modules and the genes constituting
#' that module
#' @export
#'
prepGMT.character <- function(pathway_list, ...){
  if(!file.exists(pathway_list)){
    stop(glue::glue("Can't find {pathway_list}"))
  }

  clusterProfiler::read.gmt(pathway_list) |>
    prepGMT.data.frame()

}
