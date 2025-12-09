

# función egm -------------------------------------------------------------

#' Extinct Generation Method (EGM) from Lexis-square deaths
#'
#' Compute cohort population counts (the extinct generation method) by
#' accumulating deaths along cohort lines in a Lexis diagram, locating estimates at
#' at the beginning of each year. Inputs are
#' calendar year \code{y}, age at the Lexis square \code{x}, and deaths
#' \code{d} in each square. Cohort (year of birth, \code{yb}) can be supplied
#' directly, or derived by splitting each Lexis square into upper/lower
#' triangles using \code{alpha}.
#'
#' @details
#' The extinct generation method reconstructs cohort population size at age
#' \eqn{x} as the (reverse) cumulative sum of cohort deaths from age \eqn{x}
#' onward.  Implementation note: this uses \code{dplyr::mutate(..., .by = yb)}, which
#' requires \strong{dplyr >= 1.1.0}.
#'
#' @param y Integer or numeric vector. Calendar year (Lexis square).
#' @param x Integer or numeric vector. Age (Lexis square).
#' @param d Numeric vector. Deaths in the Lexis square \eqn{(y, x)}.
#' @param yb Integer or numeric vector (optional). Year of birth (cohort).
#'   If supplied, it must be the same length as \code{y}/\code{x}/\code{d} and
#'   will be used directly (no triangle split). If \code{NULL}, \code{yb} is
#'   derived from \code{y} and \code{x} using \code{alpha}.
#' @param alpha Numeric scalar in \eqn{[0, 1]} (optional). Upper triangle share
#'   used to split Lexis squares when \code{yb} is \code{NULL}. Default is
#'   \code{0.5}. Ignored when \code{yb} is provided.
#'
#' @return
#' A data frame with one row per cohort–age–year of population stock at the beggining of the year:
#' \describe{
#'   \item{yb}{Year of birth (cohort).}
#'   \item{x}{Age (exact age on vertical segment).}
#'   \item{y}{Calendar year corresponding to age \eqn{x}.}
#'   \item{AgeDeath105}{Inform if max age with data is 105 or more. A flag for extinct interpretation.}
#'   \item{pop}{Cohort population at exact age \eqn{x}, reconstructed by reverse
#'   cumulative deaths.}
#' }
#'
#' @section Assumptions:
#' \itemize{
#' \item Death registration is complete for extinct (or nearly extinct) cohorts.
#' \item When \code{yb} is not provided, the within-square timing of deaths is
#'       approximated by a constant split \code{alpha} between triangles.
#' }
#' @export

egm <- function(y, x, d, yb = NULL, alpha = NULL){

  # initial check
  stopifnot(length(y)==length(x))
  stopifnot(length(y)==length(d))

  # set cohort deaths
  if(!is.null(yb)){
    deaths <- data.frame(y, x, d, yb)
  }else{
    if(is.null(alpha)){
      alpha <- .5
    }
    deaths <- data.frame(y, x, d, alpha)
    deaths_lower <- deaths %>%
      dplyr::mutate(d = (1-alpha) * d, triangle = "l", yb = y - x)
    deaths_upper <- deaths %>%
      dplyr::mutate(d = alpha * d, triangle = "u", yb = y - x - 1)
    deaths <- dplyr::bind_rows(
      deaths_lower, deaths_upper) %>%
      dplyr::select(yb, x, y, triangle, d) %>%
      dplyr::arrange(yb, y, desc(triangle))
  }

  # compute population at the vertical segment (aged x) and horizontal (exact age)
  pop_egm_Lexis <- deaths %>%
    dplyr::mutate(pop = rev(cumsum(rev(d))),.by = yb) %>%
    dplyr::mutate(lexis_segment = ifelse(triangle == "u", "vertical", "horizontal"))

  # return aged x for stock counts. mark as probably not extincted
  pop_egm <- pop_egm_Lexis %>%
    left_join(pop_egm_Lexis %>%
                dplyr::summarise(max_x = max(x), .by = yb)) %>%
    dplyr::mutate(AgeDeath105 = max_x >= 105) %>%
    dplyr::filter(lexis_segment == "vertical") %>%
    dplyr::select(yb, x, y, AgeDeath105, pop)

  return(pop_egm)
}


# ejemplo -----------------------------------------------------------------

# con planilla "CHL_GenExtintas.xlsx"
library(tidyverse)
matrix_period_deaths <- readxl::read_xlsx("CHL_GenExtintas.xlsx", range = "A2:AI122")
table_period_deaths <- matrix_period_deaths %>%
  pivot_longer(`1990`:`2023`, names_to = "year", values_to = "deaths") %>%
  rename(age = Edad)
pop_egm <- egm(y = table_period_deaths$year %>% as.integer(),
               x = table_period_deaths$age %>% as.integer(),
               d = table_period_deaths$deaths %>% as.numeric())
pop_egm %>%
  filter(y == 1992, AgeDeath105) %>%
  ggplot(aes(x, pop)) +
  geom_line()
