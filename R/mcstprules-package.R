#' @keywords internal
#' @description
#' A comprehensive set of tools to compute, analyze, and visualize cost allocation rules
#' for Minimum Cost Spanning Tree Problems (MCSTP).
#'
#' @details
#' This package implements two major approaches for cost allocation found in the literature:
#' algorithmic rules and rules defined through cooperative games. Additionally, it provides
#' a robust suite of analytical tools to evaluate core stability, geometric properties,
#' and cost sensitivity.
#'
#' For a detailed overview, classification, and references of all
#' implemented rules and analysis tools, please see the central documentation page:
#' \code{\link{mcstRules}}.
#'
#' @docType package
#' @name mcstprules-package
"_PACKAGE"

## usethis namespace: start
#' @importFrom graphics
#' abline barplot contour legend matplot points polygon
#' segments text title mtext par
#'
#' @importFrom grDevices
#' adjustcolor chull hcl.colors rainbow rgb devAskNewPage
#'
#' @importFrom stats
#' setNames
#'
#' @importFrom utils
#' read.csv read.csv2 combn
## usethis namespace: end
NULL

