#' Extrapolate mortality by age, constrained to a life expectancy value.
#' @description Instead of fitting parameters from some mortality law in previous ages to an observed OAG, find the parameters for the mortality function that replicates e(OAG) (Ediev, 2016).
#' @param nMx numeric. Vector of mortality rates in abridged or single age classes.
#' @param Age integer. Vector with lower bound for each age class (could be integer or abridged). Last age is assumed as lower age from open age group.
#' @param Sex character. Either male \code{"m"}, or female \code{"f"}.
#' @param method character. This indicates how to calculate `e_{OAG}` with methods: `classical` (default), `H-C` or `Mitra` (Ediev, 2014).
#' @param extrapLaw character. Available options: `Gompertz` or `Kannisto`.
#' @param eOAG numeric. An estimate of life expectancy in the input OAG age. If this has a value, then replace `method`.
#' @param OAnew integer. After extrapolating, pick a new OAG.
#' @param alpha numeric. Parameter for `H-C` method. Default values are automatically selected for input OAG.
#' @param beta numeric. Parameter for `H-C` method. Default values are automatically selected for input OAG.
#' @param r numeric. Parameter for `H-C` and `Mitra` methods. Annual growth rate of the population in the open age interval.
#' @param x_hat numeric. Parameter for `Mitra` method. Mean age of the population in the open age interval.
#' @export
#' @examples
#' \dontrun{
#' # Pakistan - 1968-1971 - males - OAG=65 - HLD (2020)
#' nMx <- c(0.13328, 0.01539, 0.0031, 0.00155, 0.00169, 0.00185, 0.00201,
#'          0.0024, 0.00289, 0.00367, 0.00498, 0.00736, 0.01214, 0.02301, 0.0879)
#' Age <- c(0L, 1L, 5L, 10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L, 65L)
#' fit_extrap_constr <- lt_extrap_constrained(nMx, Age, OAnew = 100, extrapLaw = "Gompertz")
#' fit_extrap_constr$ex[fit_extrap_constr$Age == 65] - (1/nMx[Age == 65])
#' # if some other value is required on e(65), like 10. The function forces extrapolation to fit that value.
#' nMx <- c(0.13328, 0.01539, 0.0031, 0.00155, 0.00169, 0.00185, 0.00201,
#'          0.0024, 0.00289, 0.00367, 0.00498, 0.00736, 0.01214, 0.02301, 0.0879)
#' Age <- c(0L, 1L, 5L, 10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L, 65L)
#' fit_extrap_constr <- lt_extrap_constrained(nMx, Age, OAnew = 100, extrapLaw = "Gompertz", eOAG = 10)$lt
#' fit_extrap_constr$ex[fit_extrap_constr$Age == 65] - 10
#' }

lt_extrap_constrained <- function(nMx, Age,
                                  Sex = "m",
                                  method = "classical",
                                  extrapLaw = "Gompertz",
                                  eOAG = NULL,
                                  OAnew = 100,
                                  alpha = NULL, # H-C recommended for OAG 65
                                  beta = NULL, # H-C recommended for OAG 65
                                  r = NULL, # for non classical e(OAG) computation
                                  x_hat = NULL # for non classical e(OAG) computation
){

  # initial settings
  if(is_abridged(Age)) complete = FALSE else complete = TRUE
  extrapFrom <- max(Age)
  Sex <- match.arg(Sex, c("m", "f"))
  method <- match.arg(method, c("classical", "H-C", "Mitra"))
  OAG <- max(Age)
  this_Sex <- Sex

  # coefficients for non-classical variants
  coefficients <- data.frame(
    c(rep("f",6),rep("m",6),rep("b",6)),
    matrix(c(
      40, 1.0, 0.283, 0.321, 50.045, 0.241, -4.918,
      55, 1.1, 0.207, 0.241, 61.025, 0.303, -4.503,
      65, 1.4, 0.095, 0.100, 69.200, 0.335, -3.670,
      75, 1.4, 0.095, 0.109, 77.701, 0.380, -2.676,
      85, 1.4, 0.095, 0.104, 86.460, 0.470, -1.883,
      95, 1.4, 0.095, 0.062, 95.591, 0.626, -0.867,
      40, 1.0, 0.283, 0.330, 50.924, 0.196, -3.919,
      55, 1.1, 0.207, 0.236, 61.406, 0.269, -3.722,
      65, 1.4, 0.095, 0.102, 69.229, 0.318, -3.180,
      75, 1.4, 0.095, 0.108, 77.563, 0.379, -2.398,
      85, 1.4, 0.095, 0.102, 86.355, 0.482, -1.863,
      95, 1.4, 0.095, 0.058, 95.633, 0.609, -0.914,
      40, 1.0, 0.283, 0.308, 50.839, 0.206, -3.849,
      55, 1.1, 0.207, 0.234, 61.115, 0.293, -4.030,
      65, 1.4, 0.095, 0.099, 69.117, 0.335, -3.324,
      75, 1.4, 0.095, 0.108, 77.583, 0.387, -2.481,
      85, 1.4, 0.095, 0.102, 86.405, 0.477, -1.803,
      95, 1.4, 0.095, 0.061, 95.518, 0.658, -0.929), ncol = 7, byrow = T)) %>%
    setNames(c("Sex","a","alpha","beta","beta.hmd","C","k1","k2")) %>%
    mutate(diff = abs(a-OAG)) %>%
    filter(Sex == this_Sex) %>%
    arrange(diff) %>%
    slice(1)

  # open age group rate
  Ma <- last(nMx)

  # life expectancy computation depending variant
  if(method == "classical"){
    ex_obj <- 1/Ma
  }
  if(method == "H-C"){
    if(is.null(r)) {r <- .01; message("was assumed r=.01")}
    if(is.null(alpha))  alpha <- coefficients$alpha
    if(is.null(beta))  beta <- coefficients$beta
    ex_obj <- 1/Ma*exp(-beta*r*Ma^(-alpha))
  }
  if(method == "Mitra"){
    if(is.null(r)) {r <- .01; message("was assumed r=.01")}
    if(is.null(x_hat)) {
      C <- coefficients$C
      k1 <- coefficients$k1
      k2 <- coefficients$k2
      message("mean age was approximated using regression formula")
      x_hat <- C + 1/Ma*k1 + 1/Ma*k2*r
    }
    ex_obj <- 1/Ma*exp(-r*(1/Ma-(1+r*1/Ma)*(x_hat-extrapFrom)))
  }
  if(!is.null(eOAG)){
    ex_obj <- eOAG
  }

  # choose b for fitting ex_obj
  if(complete){
    Age_extrap <- 0:OAnew
  }else{
    Age_extrap <- c(0,1,seq(5,OAnew,5))
    }
  nMx_extrap <- c(nMx,rep(NA,length(Age_extrap)-length(Age)))
  b_optim <- optimise(fo_extrap, interval=c(.0001,1),
                      Age=Age_extrap, nMx = nMx_extrap, Sex = Sex,
                      extrapFrom = extrapFrom, ex_obj = ex_obj,
                      extrapLaw = extrapLaw)$minimum

  # extrapolate rates
  nMx_prev <- nMx_extrap[Age_extrap<extrapFrom]
  nMx_extrap <- law_extrap(nMx = nMx_extrap, Age = Age_extrap,
                           extrapFrom = extrapFrom,
                           b = b_optim, extrapLaw = extrapLaw)
  nMx_hat <- c(nMx_prev, nMx_extrap)

  # compute final lt
  lt_out <- lt_ambiguous(nMx_or_nqx_or_lx = nMx_hat, type = "m", Sex = Sex,
                         Age = Age_extrap, OAnew=OAnew, Single = complete, OAG = FALSE)

  # converged?
  if(round(lt_out$ex[lt_out$Age==OAG],1)!=round(ex_obj,1)){
    mess <- paste0("Not converged. Objective e(a) is ", round(ex_obj,1),
                   " and got ",round(lt_out$ex[lt_out$Age==OAG],1),".")
  }else{
    mess <- "Converged"
  }
  out <- list(conv = mess, lt = lt_out)
  return(out)
}

# function to minimize
fo_extrap <- function(b, Age, Sex = Sex, extrapFrom, ex_obj, nMx, extrapLaw){
  nMx_extrap <- law_extrap(nMx = nMx, Age = Age, extrapFrom = extrapFrom, b = b, extrapLaw = extrapLaw)
  nMx_prev <- nMx[Age<extrapFrom]
  nMx_hat <- c(nMx_prev, nMx_extrap)
  if(is_abridged(Age)) complete = FALSE else complete = TRUE
  ex_hat <- lt_ambiguous(nMx_or_nqx_or_lx = nMx_hat, type = "m", Age = Age, Sex = Sex, OAnew=extrapFrom, OAG = FALSE, Single = complete) %>%
    filter(Age==extrapFrom) %>%
    pull(ex)
  # quadratic relative diff
  quad_dif <- (ex_hat/ex_obj-1)^2
  return(quad_dif)
}

# given parameter "b", apply extrapolation law from some age "extrapFrom".
law_extrap <- function(nMx, Age, extrapFrom, b, extrapLaw = "Gompertz", k=1){
  # where is the starting age for extrapolate
  id <- which(Age==extrapFrom)
  # index of observations to average
  id_mean <- seq(id - 1,length.out = k, by=-1)
  # abr or complete
  if(is_abridged(Age)) AgeSpan = 5 else AgeSpan = 1

  # apply some law
  if(extrapLaw == "Gompertz"){
    C_mean <- exp(mean(log(nMx[id_mean])))
    nMx_hat <- C_mean * exp(b*(Age-extrapFrom+AgeSpan))
  }
  if(extrapLaw == "Kannisto"){
    C_mean <- exp(mean(log(nMx[id_mean])))/(1-exp(mean(log(nMx[id_mean]))))
    nMx_hat <- C_mean * exp(b*(Age-extrapFrom+AgeSpan))/(1 + C_mean * exp(b*(Age-extrapFrom+AgeSpan)))
  }
  # return only ages that matters
  return(nMx_hat[Age>=extrapFrom])
}
