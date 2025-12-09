#' Intercensal survival estimates of 15q60 based on Li-Gerland (2012).
#' @description Adapted from matlab code provided from PG.
#' @param c1 numeric vector. Population at observation 1.
#' @param c2 numeric vector. Population at observation 2.
#' @param age integer vector. Lower bound of age groups from first census. Last age is assumed as lower age open age group.
#' @param date1 Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for the source (typically census or survey).
#' @param date2 Same than date1 but for observation 2.
#' @param shortcut logical. If reply matlab code from author. `FALSE` is the default, following the paper.
#' @export
#' @examples
#' \dontrun{
#' # Albania. Replicate examples from spreadsheet "OldMortPG.xls"
#' alb_males <- census_method(c(1326536, 566784, 687820),
#'                            c(1441428,848609,762845),
#'                            c(60, 65, 70), "1981/3/1", "1998/3/2")
#' alb_females <- census_method(c(929401,435774,487321),
#'                              c(1241433,698000,624473),
#'                              c(60, 65, 70), "1981/3/1", "1998/3/2")
#' c(alb_males$q60_15, alb_females$q60_15)
#' }
census_method<-function(c1, c2, age, date1, date2, shortcut = F){

  # check input
  ages_scope <- c(60, 65, 70)
  if(any(!age %in% ages_scope)) stop("needs pop in ages 60, 65 and 70")
  ages <- length(age)
  c1 <- c1[age %in% ages_scope]
  c2 <- c2[age %in% ages_scope]
  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2),2)
  interc_t <- (date2-date1)

  # step 1: r-variable stable pop
  r <- log(c2/c1)/interc_t
  cum_r <- c(r[1]*2.5, r[1]*5 + r[2]*2.5, (r[1]+r[2])*5 + r[3]*2.5)
  N <- (c1*c2)^.5
  L <- N * exp(cum_r)

  # step 2: age heaping correction: heaping ratio at age 60 equals to that at age 70
  S <- c(L[2]/L[1], L[3]/L[2])
  a <- -.28; b <- 1.27 # (estimated from paper)
  S65_line <- a + b * S[1]
  L60 <- L[1]; L65 <- L[2]; L70 <- L[3]

  # adjust age heaping
  if(L60*L70>L65^2){ # from matlab code (basically that S[1]<S[2]). The condition in paper is S[2]>S65_line
    alpha <- L60/L70
    A <- b - a * alpha - alpha
    B <- a*(L60-alpha*L65)+2*b*L65+L60+alpha*L70 #?
    C <- L65*(a*L60+b*L65)-L60*L70
    delta <- (-B + (B^2 - 4*A*C)^.5)/2/A
    L_hat <- c(L60-L60/L70*delta, L65+delta, L70-delta)
  }else{
  # minimal adjustment
    S_hat <- (-a * b + L65/L60 + b * L70/L65) / (1+b^2)
    S_hat[2] <- a + b * S_hat[1]
    if(shortcut){
      ### in matlab code (more simplified, with w=1)
      L_hat <- c()
      L_hat[1] <- (L60+S_hat[1]*L65+S_hat[1]*S_hat[2]*L70) / (1+ S_hat[1]^2+S_hat[1]^2*S_hat[2]^2)
      L_hat[2] <- L_hat[1] * S_hat[1]
      L_hat[3] <- L_hat[2] * S_hat[2]
    }else{
      ### in paper
      w <- .5
      L_hat[1] <- w * (L60+S_hat[1]*L65+S_hat[1]*S_hat[2]*L70) / (1+ S_hat[1]^2+S_hat[1]^2*S_hat[2]^2) + (1-w) * L_hat[1]
      L_hat[2] <- w * S_hat[1] * (L60+S_hat[1]*L65+S_hat[1]*S_hat[2]*L70) / (1+ S_hat[1]^2+S_hat[1]^2*S_hat[2]^2) + (1-w) * L_hat[2]
      L_hat[3] <- w * S_hat[1] * S_hat[] * (L60+S_hat[1]*L65+S_hat[1]*S_hat[2]*L70) / (1+ S_hat[1]^2+S_hat[1]^2*S_hat[2]^2) + (1-w) * L_hat[3]
    }
  }
  # step 3: gompertz extension
  initial_values <- c(L_hat[1]/5, .02, .1) # from matlab code
  y_star <- pracma::fminsearch(min_L_gomp, initial_values, L = L_hat)$xmin
  mx <- c(y_star[2] * exp(2.5 * y_star[3]), y_star[2] * exp(7.5 * y_star[3]), y_star[2] * exp(12.5 * y_star[3]))
  lx <- get_Lx_gomp(y_star)$lx
  q60_15 <- 1 - lx[4]/lx[1]

  # out
  return(list(q60_15 = q60_15,
              par = y_star,
              lt = data.frame(lx = lx, Lx = c(L_hat, NA), mx = c(mx, NA)))
         )
}

# compute Lx and lx (x = 60, 65, 70) given parameters l(60), mu(60) and b in a Gompertz function. Kind of numerical integration
get_Lx_gomp <- function(x){
  lx <- rep(0, 1500) # 15 years interval
  for (i in 1:1500) {
    u <- i/100
    lx[i] <- x[1] * exp((-x[2]/x[3]) * (exp(x[3]*u) - 1)) # gompertz fun
  }
  Lm <- rep(0, 3)
  for (i in 1:500) {
    Lm[1] <- Lm[1] + lx[i]/100
    Lm[2] <- Lm[2] + lx[i+500]/100
    Lm[3] <- Lm[3] + lx[i+1000]/100
  }
  return(list(Lx = Lm, lx = lx[c(1, 500, 1000, 1500)]))
}

# minimize function for parameters l(60), mu(60) and b.
min_L_gomp <- function(x, L){
  Lm <- get_Lx_gomp(x)$L
  y <- 0
  for (i in 1:3) {
    y <- y + (L[i] - Lm[i])^2
  }
  return(y)
}
