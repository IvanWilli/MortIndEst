
#' Extinct Generation Method (EGM)
#'
#' Compute cohort population counts by accumulating deaths along cohort lines
#' in a Lexis diagram. Population stocks are located at the beginning of each
#' calendar year. See vignette for an example.
#'
#' @param y Numeric vector. Calendar year.
#' @param x Numeric vector. Age in the Lexis square.
#' @param d Numeric vector. Deaths in the Lexis square.
#' @param yb Numeric vector (optional). Year of birth. If NULL, cohort deaths
#'   are obtained by splitting each square using \code{alpha}.
#' @param alpha Numeric between 0 and 1. Upper-triangle share. Default 0.5.
#'
#' @return A data frame with \code{yb}, \code{x}, \code{y},
#'   and reconstructed population \code{pop}.
#'
#' @export
egm <- function(y, x, d, yb = NULL, alpha = NULL){

  deaths <- .lexis_cohort_deaths(y, x, d, yb, alpha)

  # Reverse cumulative cohort deaths along the Lexis trajectory
  pop_lexis <- deaths %>%
    dplyr::mutate(pop = rev(cumsum(rev(d))), .by = yb)

  # Maximum observed age at death by cohort
  max_age <- deaths %>%
    dplyr::summarise(max_x = max(x), .by = yb)

  # Vertical segment = population aged x at beginning of year
  pop_lexis %>%
    dplyr::filter(triangle == "u") %>%
    dplyr::left_join(max_age, by = "yb") %>%
    # dplyr::mutate(AgeDeath105 = max_x >= 105) %>%
    dplyr::select(yb, x, y, pop)
}

#' Nearly-Extinct Cohort Methods
#'
#' Estimate population stocks of nearly-extinct cohorts using the Survivor
#' Ratio (SR), Das Gupta (DG), or Survivor Ratio Advanced (SA) method. See vignette for an example.
#'
#' Deaths are first assigned to cohorts from Lexis triangles. The method
#' estimates the population remaining at the end of the observed death series;
#' earlier population stocks are then reconstructed backwards by adding
#' observed cohort deaths.
#'
#' @param y Numeric vector. Calendar year.
#' @param x Numeric vector. Age in the Lexis square.
#' @param d Numeric vector. Deaths in the Lexis square.
#' @param yb Numeric vector (optional). Year of birth. If NULL, cohort deaths
#'   are obtained by splitting each square using \code{alpha}.
#' @param alpha Numeric scalar in between 0 and 1. Upper-triangle share. Default 0.5.
#' @param method Character. \code{"SR"}, \code{"DG"}, or \code{"SA"}.
#' @param k Number of years over which survivor ratios are calculated.
#'   Used by SR and SA. Default 5.
#' @param m Number of older cohorts used to estimate survivor/death ratios.
#'   Default 5.
#' @param n Number of older cohorts used to estimate survival improvement
#'   in SA. Default 10.
#' @param min_age Youngest age to estimate. Default 85.
#' @param omega Extinction age. Cohorts aged \code{omega} or older at the
#'   calculation date are assumed extinct. Default is max observed age + 1.
#'
#' @return A data frame with \code{yb}, \code{x}, \code{y}
#'   and estimated population \code{pop}.
#'
#' @references
#' Terblanche, W. and Wilson, T. (2015). An Evaluation of Nearly-Extinct
#' Cohort Methods for Estimating the Very Elderly Populations of Australia
#' and New Zealand. PLoS ONE 10(4): e0123692.
#'
#' @export
nec <- function(y, x, d, yb = NULL, alpha = NULL,
                method = c("SR", "DG", "SA"),
                k = 5, m = 5, n = 10,
                min_age = 85, omega = NULL){
  method <- match.arg(method)
  deaths <- .lexis_cohort_deaths(y, x, d, yb, alpha)
  Y <- max(y)
  if(is.null(omega)) {
    omega <- max(x) + 1
    message("omega==NULL. Set to max(x) + 1")
  }
  # Observed future deaths from each vertical Lexis segment
  obs <- deaths %>%
    dplyr::mutate(O = rev(cumsum(rev(d))), .by = yb) %>%
    dplyr::filter(triangle == "u") %>%
    dplyr::select(yb, x, O)
  # O(c,a): observed deaths after cohort c reaches age a
  O <- function(c, a)
    sum(obs$O[obs$yb == c & obs$x == a], na.rm = TRUE)
  # Cohorts currently aged min_age,..., omega-1; oldest first
  cohorts <- Y - ((omega-1):min_age)
  pcur <- stats::setNames(rep(NA_real_, length(cohorts)), cohorts)
  # Population at an earlier vertical segment:
  # current survivors + subsequent observed deaths
  P <- function(c, a){
    pc <- pcur[as.character(c)]
    if(length(pc) == 0 || is.na(pc)) pc <- 0
    pc + O(c, a)
  }
  # Cohort deaths between age a and a+1
  D <- function(c, a) O(c, a) - O(c, a+1)
  for(c in cohorts){
    a <- Y-c  # age at beginning of year Y+1
    # -------------------------------------------------------
    # Survivor Ratio / Survivor Ratio Advanced
    # -------------------------------------------------------
    if(method %in% c("SR", "SA")){
      old <- c - seq_len(m)
      # Pooled survivor ratio among m older cohorts
      R <- sum(sapply(old, \(z) P(z, a))) /
           sum(sapply(old, \(z) P(z, a-k)))
      if(method == "SA"){
        # Cohort-specific survivor ratios, oldest -> newest
        old_n <- sort(c - seq_len(n))
        Rn <- sapply(old_n, \(z) P(z, a) / P(z, a-k))
        Rn <- Rn[is.finite(Rn) & Rn > 0]
        # Mean improvement in survivor ratios
        r <- if(length(Rn) > 1)
          mean(Rn[-1] / Rn[-length(Rn)] - 1) else 0
        R <- R * (1+r)^((m+1)/2)
        # Avoid impossible survivor ratios

        R <- min(R, .999999)
      }
      # Deaths observed during the previous k years
      Dk <- O(c, a-k) - O(c, a)
      pcur[as.character(c)] <-
        if(is.finite(R) && R > 0 && R < 1)
          R/(1-R) * Dk else NA_real_
    }
    # -------------------------------------------------------
    # Das Gupta
    # -------------------------------------------------------
    if(method == "DG"){
      # Deaths during the last observed cohort-year
      Dh <- D(c, a-1)
      pc <- 0
      # Project successive future cohort deaths
      for(A in a:(omega-1)){
        h <- A-a
        # Cohorts old enough to have already experienced age A
        old <- c - h - 1 - 0:(m-1)
        den <- sum(sapply(old, \(z) D(z, A-1)))
        num <- sum(sapply(old, \(z) D(z, A)))
        dr <- if(den > 0) num/den else 0
        Dh <- Dh * dr
        pc <- pc + Dh
      }
      pcur[as.character(c)] <- pc
    }
  }
  # Reconstruct population stocks backwards from current estimate
  out <- dplyr::bind_rows(lapply(cohorts, function(c){
    amax <- Y-c
    aa <- min_age:amax
    data.frame(
      yb  = c,
      x   = aa,
      y   = c + aa + 1,
      pop = sapply(aa, \(a) P(c, a))
    )
  }))
  # Same diagnostic as egm()
  max_age <- deaths %>%
    dplyr::summarise(max_x = max(x), .by = yb)
  # return
  out %>%
    dplyr::left_join(max_age, by = "yb") %>%
    # dplyr::mutate(AgeDeath105 = max_x >= 105) %>%
    dplyr::select(yb, x, y, pop)
}

# Internal: convert Lexis-square deaths to cohort-triangle deaths
.lexis_cohort_deaths <- function(y, x, d, yb = NULL, alpha = NULL){
  stopifnot(length(y) == length(x), length(y) == length(d))
  if(is.null(yb)){
    if(is.null(alpha)) alpha <- .5
    stopifnot(length(alpha) == 1, alpha >= 0, alpha <= 1)
    deaths <- dplyr::bind_rows(
      data.frame(y, x, d = (1-alpha)*d, triangle = "l", yb = y-x),
      data.frame(y, x, d = alpha*d,     triangle = "u", yb = y-x-1)
    )
  }else{
    stopifnot(length(yb) == length(y))
    # Infer triangle from cohort membership
    triangle <- ifelse(yb == y-x-1, "u",
                ifelse(yb == y-x, "l", NA))
    if(anyNA(triangle))
      stop("yb must correspond to either y-x or y-x-1.")
    deaths <- data.frame(y, x, d, triangle, yb)
  }
  deaths %>%
    dplyr::select(yb, x, y, triangle, d) %>%
    dplyr::arrange(yb, y, dplyr::desc(triangle))
}
