#' coherent rate extension in older ages
#' @description Based on cokannisto function from Mostcast package (Ševčíková, 2023), extend to gompertz law and with time.
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

# co-gompertz (simple code before turning into a function)
coherent_mx <- function(mx,
                        law = "kannisto",
                        fit.ages = NULL, estim.ages = NULL,
                        method = "coherent",
                        age_conv,
                        fit.years = NULL, estim.years = NULL){

  # handle cathegories
  sexs <- data.frame(sex = sort(unique(mx$sex)),
                     sex_int = as.integer(as.factor(sort(unique(mx$sex)))) - 1)
  mx <- merge(mx, sexs, by = "sex", all.x = TRUE)
  ages <- sort(unique(mx$age))

  # handle errors
  if(is.null(fit.ages)){
    fit.ages <- ages
  }
  if(is.null(estim.ages)){
    estim.ages <- ages
  }

  # law
  if(law == "kannisto"){
    mx$y <- log(mx$mx) - log(1 - mx$mx)
  }else if(law == "gompertz"){
    mx$y <- log(mx$mx)
  }

  # type model
  if(length(unique(mx$year)) > 1){
    mx$year <- as.numeric(mx$year)
    years <- sort(unique(mx$year))
    if(is.null(fit.years)){
      fit.years <- years
    }
    if(is.null(estim.years)){
      estim.years <- years
    }
    # browser()
    fit.model <- lm(y ~ age + sex_int + year,
                    data = mx[mx$age %in% fit.ages & mx$year %in% fit.years, ])
    new_data <- expand.grid(age = estim.ages, sex_int = sexs$sex_int, year  = estim.years)
    pred.model <- predict(fit.model, newdata = new_data)
    mx.estim <- cbind(new_data, y = pred.model)
  }else{
    mx$year <- NULL
    if(method == "convergent"){
      fit.model.0 <- lm(y ~ age,
                        data = mx[mx$age %in% fit.ages & mx$sex_int == 0, ])
      new_data.0 <- data.frame(age = estim.ages)
      pred.model.0 <- predict(fit.model.0, newdata = new_data.0)
      init_age <- fit.ages[which.min(estim.ages[1] - fit.ages)]
      fit.model.1 <- lm(y ~ age,
                        data = data.frame(age = c(init_age, age_conv),
                                          y = c(mx$y[mx$age %in% init_age & mx$sex_int == 1],
                                                pred.model.0[estim.ages == age_conv])))
      pred.model.1 <- predict(fit.model.1, newdata = new_data.0)
      pred.model.1[estim.ages>age_conv] <- pred.model.0[estim.ages>age_conv]
      mx.estim <- data.frame(y = c(pred.model.0, pred.model.1),
                             age = rep(estim.ages, 2),
                             sex_int = c(rep(0, length(estim.ages)), rep(1, length(estim.ages)))
      )
    }else{
      fit.model <- lm(y ~ age + sex_int,
                      data = mx[mx$age %in% fit.ages, ])
      new_data <- expand.grid(age = estim.ages, sex_int = sexs$sex_int)
      pred.model <- predict(fit.model, newdata = new_data)
      mx.estim <- cbind(new_data, y = pred.model)
    }
  }

  # return
  if(law == "gompertz"){
    mx.estim$mx <- exp(mx.estim$y)
  }else if(law == "kannisto"){
    mx.estim$mx <- exp(mx.estim$y)/(1+exp(mx.estim$y))
  }
  mx.estim <- merge(mx.estim, sexs, by = "sex_int", all.x = TRUE)
  mx.estim$y <- NULL; mx.estim$sex_int <- NULL
  mx$y <- NULL; mx$sex_int <- NULL
  # browser()
  mx.out <- rbind(mx[mx$age < min(estim.ages),], mx.estim)
  return(mx.out)
}
