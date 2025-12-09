# library(DemoTools)
# library(tidyverse)

# compute error -----------------------------------------------------------
get_error_sq <- function(param, x, pop1, pop2, n,
                         m_age_fit = NULL, age_fit = 80,
                         law, OA){

  # browser()
  b <- param[1]
  ages <- 1:(length(age_fit:120)-1)
  if(law == "gompertz"){
    a <- m_age_fit
    m_hat <- a * exp(b*ages)
  }
  if(law == "kannisto"){
    a <- m_age_fit/(1-m_age_fit)
    m_hat <- (a * exp(b * ages))/(1 + a * exp(b * ages))
  }
  lt <- lt_single_mx(Age = age_fit:120, nMx = c(m_age_fit, m_hat), OAnew = OA)
  pop_proj_hat <- surv_pop(x[x >= age_fit],
                           pop0 = pop1[x >= age_fit],
                           S = c(lt$Sx[-1], lt$Sx[lt$Age == OA]),
                           t = n)
  # browser()
  pop2_cum <- rev(cumsum(rev(pop2)))[x >= (age_fit+n)]
  pop2_hat <- pop_proj_hat[,ncol(pop_proj_hat)][(n+1):nrow(pop_proj_hat)]
  pop2_cum_hat <- rev(cumsum(rev(pop2_hat)))
  stopifnot(length(pop2_cum_hat) == length(pop2_cum))
  error <- sum(((pop2_cum - pop2_cum_hat)/pop2_cum)^2)
  return(error)
}


# optim -------------------------------------------------------------------
# optimize parameters to fit pop2 from pop 1, given a specific law and age
law_fit_old_pop <- function(x, pop1, pop2, n = 10,
                            age_fit = 80,
                            m_age_fit = NULL,
                            law = "gompertz",
                            lower = .2, upper = .6){

  # browser()
  # set OA with data
  OA <- max(x)

  # optimize
  get_optim <- optimize(f = get_error_sq,
                        lower = lower, upper = upper,
                        law = law,
                        x = x,
                        pop1 = pop1,
                        pop2 = pop2,
                        n = n,
                        age_fit = age_fit,
                        m_age_fit = m_age_fit,
                        OA = OA)

  # what pop arrives
  b_hat <- get_optim$minimum
  ages <- 1:(length(age_fit:120)-1)
  if(law == "gompertz"){
    a <- m_age_fit
    m_hat <- a * exp(b_hat*ages)
  }
  if(law == "kannisto"){
    a <- m_age_fit/(1-m_age_fit)
    m_hat <- (a * exp(b_hat * ages))/(1 + a * exp(b_hat * ages))
  }
  lt <- lt_single_mx(Age = age_fit:120, nMx = c(m_age_fit, m_hat), OAnew = OA)
  pop_proj_hat <- surv_pop(x[x >= age_fit],
                           pop0 = pop1[x >= age_fit],
                           S = c(lt$Sx[-1], lt$Sx[lt$Age == OA]),
                           t = n)
  pop2_cum <- rev(cumsum(rev(pop2)))[x %in% (age_fit+n)]
  pop2_hat <- pop_proj_hat[,ncol(pop_proj_hat)]
  pop2_cum_hat <- rev(cumsum(rev(pop2_hat)))[n+1]
  pop2_result <- data.frame(age = age_fit,
                            pop2 = pop2_cum,
                            pop2_fit = pop2_cum_hat,
                            ex_fit = lt$ex[lt$Age==age_fit])

  # output
  return(list(params = c(b_hat),
              pop2_fit_cum = pop2_result,
              pop2_fit = data.frame(x = x[x >= age_fit],
                                    pop2 = pop2[x >= age_fit],
                                    pop2_hat = pop2_hat),
              lt = lt,
              optim_info = get_optim)
  )
}

# función de sobreviviencia 60+ -------------------------------------------

surv_pop <- function(x, pop0, S, t = NULL){
  if(!is.matrix(S)){S <- matrix(S, nrow = length(S), ncol = t)}
  xs <- length(x)
  n <- ncol(S)
  if(nrow(S)!=length(pop0) | length(pop0)!=xs) stop("dimension issue")
  pop <- matrix(0, xs, n+1)
  pop[,1] <- pop0
  for(t in 1:n){
    # t = 1; print(t)
    U <- matrix(0, xs, xs)
    U[row(U)-1 == col(U)] <- S[-xs, t]
    U[xs, xs] <- S[xs, t]
    pop[,t+1] <- U %*% pop[,t]
  }
  return(cbind(x, pop))
}

# example -----------------------------------------------------------------

# x <- 80:100
# age_fit <- 80
# m_age_fit_2001 <- 0.04070319
# pop2001 <- c(11450.5, 10300, 9296, 8732, 8228.5, 7716.5, 7257.5, 6487.5,
#              5614.5, 4785, 3866.5, 3110.5, 2439.5, 1845.5, 1341.5, 936, 681.5,
#              481, 324, 204, 282)
# pop2010 <- c(12338, 10812, 10561, 9928, 9101, 8482, 7370, 6636, 5823, 4828,
#              4352, 2975, 2463, 1948, 1581, 1270, 959, 708, 477, 313, 556)
# pop2022 <- c(11384, 10436, 10143, 9208, 8288, 7648, 7228, 6632, 5826, 5715,
#              5136, 4418, 3826, 3202, 2563, 1834, 1496, 1182, 747, 580, 785)
#
#
# ############# 2001-2010
# pop1 <- pop2001
# pop2 <- pop2010
# m_age_fit <- m_age_fit_2001
# n <- 2010-2001
# law <- "kannisto"
# OA <- 100
# get_optim <- optimize(
#   f = get_error_sq,
#   lower = .001, upper = 1,
#   law = "kannisto",
#   x = x,
#   pop1 = pop1,
#   pop2 = pop2,
#   n = n,
#   age_fit = 80,
#   m_age_fit = m_age_fit_2001,
#   OA = 100)
# b_hat <- get_optim$minimum
# ages <- 1:(length(age_fit:120)-1)
# if(law == "gompertz"){
#   a <- m_age_fit
#   m_hat <- a * exp(b_hat*ages)
# }
# if(law == "kannisto"){
#   a <- m_age_fit/(1-m_age_fit)
#   m_hat <- (a * exp(b_hat * ages))/(1 + a * exp(b_hat * ages))
# }
# lt <- lt_single_mx(Age = age_fit:120, nMx = c(m_age_fit, m_hat), OAnew = OA)
# pop_proj_hat <- surv_pop(x[x >= age_fit],
#                          pop0 = pop1[x >= age_fit],
#                          S = c(lt$Sx[-1], lt$Sx[lt$Age == OA]),
#                          t = n)
# pop2_cum <- rev(cumsum(rev(pop2)))[x %in% (age_fit+n)]
# pop2_hat <- pop_proj_hat[,ncol(pop_proj_hat)]
# pop2_cum_hat <- rev(cumsum(rev(pop2_hat)))[n+1]
# plot(pop2); lines(pop2_hat)
#
# ############# 2010-2022
# pop1 <- pop2010
# pop2 <- pop2022
# n <- 2010-2001
# law <- "kannisto"
# get_optim <- optimize(
#   f = get_error_sq,
#   lower = .001, upper = 1,
#   law = "kannisto",
#   x = x,
#   pop1 = pop1,
#   pop2 = pop2,
#   n = n,
#   age_fit = 80,
#   m_age_fit = m_age_fit_2001,
#   OA = 100)
# b_hat <- get_optim$minimum
# ages <- 1:(length(age_fit:120)-1)
# if(law == "gompertz"){
#   a <- m_age_fit
#   m_hat <- a * exp(b_hat*ages)
# }
# if(law == "kannisto"){
#   a <- m_age_fit/(1-m_age_fit)
#   m_hat <- (a * exp(b_hat * ages))/(1 + a * exp(b_hat * ages))
# }
# lt <- lt_single_mx(Age = age_fit:120, nMx = c(m_age_fit, m_hat), OAnew = OA)
# pop_proj_hat <- surv_pop(x[x >= age_fit],
#                          pop0 = pop1[x >= age_fit],
#                          S = c(lt$Sx[-1], lt$Sx[lt$Age == OA]),
#                          t = n)
# pop2_cum <- rev(cumsum(rev(pop2)))[x %in% (age_fit+n)]
# pop2_hat <- pop_proj_hat[,ncol(pop_proj_hat)]
# pop2_cum_hat <- rev(cumsum(rev(pop2_hat)))[n+1]
# plot(pop2); lines(pop2_hat)
