#' Gene Expression Values for PDAC Cancer Cell Lines exposed to Hypoxia
#'
#' @docType data
#' @usage data(expected_score_output)
#'
#' @source Derived Data
#' @format A data frame with columns:
#' \describe{
#'  \item{sample_id}{String. The name of the sample. Samples with "hyp"
#'  or "norm" in the sample id are cell lines that were exposed to
#'  hypoxic or normoxic conditions respectively. Samples with "ctrl" or
#'  "noHIF" were samples that were able to produce a HIF-mediated hypoxic
#'  response or not, respectively.}
#'  \item{pathway_score}{Float. The estimated hypoxia score for this sample.}
#' }
#' @examples
#' \dontrun{
#'  expected_score_output
#' }
"expected_score_output"
