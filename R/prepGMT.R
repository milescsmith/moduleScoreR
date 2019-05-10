#' @title prepGMT
#'
#' @description Prepare a GMT file for use in score_module_eigengenes
#'
#' @param pathway_list A data.frame consisting of two columns (`ont`` and `gene`) or the
#' path to a gmt file
#' @param ... Not used
#'
#' @return
#' @export
#'
#' @examples
#'
prepGMT <- function(pathway_list,...){
  UseMethod("prepGMT")
}

#' @rdname prepGMT
#' @method prepGMT data.frame
#' @importFrom purrr map
#' @return
#' @export
#'
prepGMT.data.frame <- function(pathway_list){
  z <- split(x = pathway_list, f = pathway_list[["ont"]])
  z %<>% map(function(x){
    y <- x %>% pull(gene)
  })
  names(z) <- unique(pathway_list$ont)
  return(z)
}

#' @rdname prepGMT
#' @method prepGMT character
#' @importFrom clusterProfiler read.gmt
#' @return
#' @export
#'
prepGMT.character <- function(pathway_list){
  if(!file.exists(pathway_list)){
    stop(glue("Can't find {pathway_list}"))
  }
  z <- read.gmt(pathway_list)
  z <- prepGMT.data.frame(z)
  return(z)
}
