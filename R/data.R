
#' Example of proper dataset
#'
#' A dataset as an example containing randomly generated numbers and fitting the requirements of `est_func()`.
#'
#' @format A data frame with 50 rows and 14 variables:
#' \describe{
#'   \item{A}{the treatment arm}
#'   \item{Z}{the negative control exposure}
#'   \item{W}{the negative control outcome}
#'   \item{Y}{the outcome}
#'   \item{V2-V10}{known covariates}
#' }
#'
"bin_data"

#' Example of proper arguments
#'
#' A dataset as an example containing arguments and fitting the requirements of `est_func()`.
#'
#' @format A list with 4 variables, of which `MW.wrong`, `MR.wrong`, `MY.wrong`, `MZ.wrong` are all `FALSE`.
"bin_args"
