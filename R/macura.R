### macura´s implementation
# main differences:
  # - instead using mlt levels, I use level from logit transformation
  # - I use a first alpha parameter for right previous level of the survey/census, and then drift parameters that constrain increasing mortality level with time within bounds
  # - instead of looking by 10 year interval first, I do it by five directly
  # - can be done with a log-quad first iteration of alpha parameters
  # - can be done using period approach (Macura´s intention) or cohort one (each optimized alpha for each parent´s age)

# warnings:
  # - Lot of local optimums
  
# functions -------------------------------------------------------------------------

# main function
orphanhood_macura_proxy <- function(NO, age_NO, date, pi, 
                                    mlt_family = "CD_West",
                                    mlt_e0 = 60, 
                                    mlt_5q0 = NULL,
                                    period = TRUE,
                                    sex = "f"){
  age <- seq(0, max(age_NO), 5)
  ages <- length(age)
  if(mlt_family == "log-quad"){
    load("orphanhood/log_quad_coeff_hmd1719.rdata")
    mlt_e0_orig <- mlt_e0
    mlt_e0 <- c(mlt_e0_orig, rep(min(mlt_e0_orig), 10))[1:ages]
    mlt_5q0_orig <- mlt_5q0
    mlt_5q0 <- c(mlt_5q0_orig, rep(min(mlt_5q0_orig), 10))[1:ages]
    if(sex == "f"){
      lx <- sapply(1:ages, function(y) MortalityEstimate::wilmothLT(Wf, q0_5 = mlt_5q0[y], e0 = mlt_e0[y])$lt$lx/1e5)
    }else{
      lx <- sapply(1:ages, function(y) MortalityEstimate::wilmothLT(Wm, q0_5 = mlt_5q0[y], e0 = mlt_e0[y])$lt$lx/1e5)
    }
    alphas_init <- apply(lx, 2, function(y) logit(1-y[14]/y[5]) - logit(1-lx[14,1]/lx[5,1]))
    lx <- lx[,1]
  }else{
    mlt_e0 <- round(mlt_e0/2.5,0) * 2.5
    mlt_data <- MortCast::MLTlookup %>% 
      filter(type %in% mlt_family, sex == ifelse(sex == "f", 2, 1), e0 == mlt_e0) %>% 
      mutate(lx = lx/1e5)
    lx <- mlt_data$lx[mlt_data$age %in% c(0, 1, seq(5, 110,5))]
    alphas_init <-rep(0,ages)
  }
  l2.5 <- splinefun(x = c(0, 1, seq(5, 110, 5)), lx, method = "monoH.FC")(seq(2.5, 107.5, 5))
  Sx <- pmax(0.01, c(l2.5[-1]/l2.5[-length(l2.5)]))
  
  # optimize alphas
  optim_star_NO <- optim(alphas_init, 
                         obj_fun, 
                         lower = c(-2, rep(0, ages-1)), 
                         upper = c(2,  rep(.5, ages-1)), 
                         Sx = Sx, pi = pi, NO = NO, age_NO = age_NO, 
                         period = period, method = "L-BFGS-B")
  alphas_optim <- cumsum(optim_star_NO$par) 
  lt_out <- map_df(1:ages, function(i){
    lx_hat <- 1 - logit_inv(alphas_optim[i] + logit(1-lx))
    this_age <- age[i]
    DemoTools::lt_abridged(lx = lx_hat[1:22], Age = c(0,1,seq(5,100,5)), Sex = sex) %>% 
      mutate(age = this_age, time_location = date - (i-1)*5 - 2.5)
  }) 
  out <- list(alphas = alphas_optim,
              pars = optim_star_NO$par,
              NO_fit = get_NO(optim_star_NO$par, Sx, pi, age, period),
              lt = lt_out,
              adult = lt_out %>% group_by(age, time_location) %>% summarise(q15_45 = 1-lx[Age==60]/lx[Age==15]))
  return(out)
}

# get expected NO based on survival levels and pi
get_NO <- function(alphas, Sx, pi, age, period = TRUE){
  ages <- length(age)
  alpha_this <- cumsum(alphas)
  Sx_hat <- lapply(1:ages,  function(x) 1 - logit_inv(alpha_this[x] + logit(1-Sx)))
  U <- lapply(1:ages, function(x) do_U(Sx_hat[[x]]))
  U[[1]] <- U[[1]]^.5 # first mothers survive half period 
  if(!period){
    U_cum <- U %>% map2(1:length(U), function(x, y) expm::`%^%`(x, y))
  }else{
    U_cum <- U %>% accumulate(function(a,b) a %*% b) 
  }
  if(length(pi)<length(U_cum)){
    for(i in (length(pi)+1):length(U_cum)){pi[[i]] <- pi[[length(pi)]]}
  }
  NO_hat <- lapply(1:ages, function(x) sum(U_cum[[x]] %*% pi[[x]])) %>% unlist()
  data.frame(age, NO_hat)
}

# do U matrix for a set of probabilities, no mater length
do_U <- function(lx){
  ages <- length(lx)
  U <- matrix(0, ages, ages)
  U[row(U)-1 == col(U)] <- lx[-ages]
  U[ages, ages] <- lx[ages]
  U
}

# cost to minimize
obj_fun <- function(alphas, Sx, pi, NO, age_NO, period){
  age <- seq(0, max(age_NO), 5)
  NO_hat <- get_NO(alphas, Sx, pi, age, period)
  NO_hat <- NO_hat$NO_hat[NO_hat$age %in% age_NO]
  # sum((NO/NO_hat-1)^2)
  sum(abs(NO/NO_hat-1))
  # sum(abs(NO-NO_hat))
}