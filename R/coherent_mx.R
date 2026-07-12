#' Extrapolate mortality rates coherently at older ages
#'
#' Fit a simple old-age mortality model to observed rates and predict rates for
#' selected ages, sexes, and optionally years. Rates are transformed with either
#' a Gompertz log transform or a Kannisto logit transform before fitting a
#' linear model.
#'
#' @description
#' Based on the coherent Kannisto idea used in the \pkg{MortCast} package,
#' extended here to allow a Gompertz law and an optional period trend.
#'
#' @param mx Data frame with mortality rates. It must contain \code{age},
#'   \code{nMx}, and \code{sex}. It may also contain \code{year}. The
#'   \code{age} column gives age lower bounds, \code{nMx} gives mortality
#'   rates, \code{sex} identifies sex groups, and \code{year} identifies the
#'   period when a time trend is fitted.
#' @param law Character. Mortality law used for the transformed rates. Options
#'   are \code{"kannisto"} and \code{"gompertz"}.
#' @param fit.ages Integer or numeric vector. Ages used to fit the regression.
#'   If \code{NULL}, all ages in \code{mx} are used.
#' @param estim.ages Integer or numeric vector. Ages to return after
#'   extrapolation. If \code{NULL}, all ages in \code{mx} are returned.
#' @param method Character. Extrapolation method when only one year is supplied.
#'   Use \code{"coherent"} for a common age slope by sex, or
#'   \code{"convergent"} to force the second sex group to converge to the first
#'   at \code{age_conv}.
#' @param age_conv Numeric. Age where the second sex group converges to the
#'   first when \code{method = "convergent"}.
#' @param fit.years Integer or numeric vector. Years used to fit the regression
#'   when \code{mx} contains more than one year. If \code{NULL}, all years are
#'   used.
#' @param estim.years Integer or numeric vector. Years to return when \code{mx}
#'   contains more than one year. If \code{NULL}, all years are returned.
#'
#' @return
#' A data frame with observed rates below the first estimated age and fitted or
#' extrapolated rates for \code{estim.ages}. Columns include \code{age},
#' \code{nMx}, \code{sex}, and \code{year} when a year dimension is supplied.
#'
#' @export
#'
#' @examples
#' mx <- data.frame(
#'   age = rep(seq(60, 80, 5), 2),
#'   sex = rep(c("f", "m"), each = 5),
#'   nMx = c(0.020, 0.035, 0.060, 0.100, 0.160,
#'           0.030, 0.050, 0.080, 0.130, 0.200)
#' )
#'
#' coherent_mx(
#'   mx,
#'   law = "gompertz",
#'   fit.ages = seq(60, 80, 5),
#'   estim.ages = seq(60, 100, 5)
#' )
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
    mx$y <- log(mx$nMx) - log(1 - mx$nMx)
  }else if(law == "gompertz"){
    mx$y <- log(mx$nMx)
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
      # init_age <- fit.ages[which.min(estim.ages[1] - fit.ages)]
      init_age <- estim.ages[1]
      estim.ages.all <- init_age:120
      pred.model.0.all <- predict(fit.model.0, newdata = data.frame(age = estim.ages.all))
      pred.model.0 <- predict(fit.model.0, newdata = new_data.0)
      fit.model.1 <- lm(y ~ age,
                        data = data.frame(age = c(init_age, age_conv),
                                          y = c(mx$y[mx$age %in% init_age & mx$sex_int == 1],
                                                pred.model.0.all[estim.ages.all == age_conv])))
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
    mx.estim$nMx <- exp(mx.estim$y)
  }else if(law == "kannisto"){
    mx.estim$nMx <- exp(mx.estim$y)/(1+exp(mx.estim$y))
  }
  mx.estim <- merge(mx.estim, sexs, by = "sex_int", all.x = TRUE)
  mx.estim$y <- NULL; mx.estim$sex_int <- NULL
  mx$y <- NULL; mx$sex_int <- NULL
  # browser()
  mx.out <- rbind(mx[mx$age < min(estim.ages),], mx.estim)
  return(mx.out)
}
