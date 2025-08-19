#' Estimate adult mortality with one survey/census from orphanhood data.
#' @description Follows IUSSP (2012) implementation for estimate adult mortality from orphanhood data, based in Timæus (1992).
#' The function also include the possibility to match the closest model life table (not using Brass logit given a pattern as default method), for some specific family (if no family is set, then return results for all CD and UN).
#' Two options for considering HIV populations: one is following Timæus and Nunn (1997) (implemented in IUSSP (2012)), and other is using the Spectrum model (Stover and others, 2012). The last needs additional data on `HIV_prop_partner`, `HIV_art`.
#' @param prop_not_orphan numeric vector. Population proportion of respondents with mother/father alive, by age groups.
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for the source (typically census or survey).
#' @param age integer vector. Lower bound of age groups from first census. Last age is assumed as lower age from the open age group.
#' @param maternal logical. Either `TRUE` for maternal orphanhood (default) or `FALSE` for paternal orphanhood.
#' @param mac numeric. Mean age at childbearing for mothers/fathers.
#' @param mlt_data_input data.frame. If is `NULL` then model life tables is used, available from \code{Morcast} package. But specific pattern can be included (`l_x` must be included, see examples).
#' @param mlt_family character. Options: "CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian". If `NULL` returns results for all.
#' @param brass_logit logical. Doing a level smoothing with 1-parameter logit Brass.
#' @param brass_logit_e0 numeric. Life expectancy at birth when `brass_logit` is `TRUE`.
#' @param brass_logit_5q0 numeric. Find best level when `brass_logit_e0` is `NULL`. If this is not `NULL`, then implied level replaces `brass_logit_e0`.
#' @param HIV_prev_antenatal numeric. Estimates of the prevalence by age of HIV infection among women attending antenatal clinics. If some value is assigned, then will be assumed AIDS variation of the method.
#' @param HIV_vert_transm numeric. Estimate of proportion of infants who acquire HIV from their mothers.
#' @param HIV_fert_ratio numeric. Estimate of relative level of fertility among HIV-positive women compared with HIV-negative women.
#' @param HIV_prev_males numeric. Estimates of population-based estimate of HIV prevalence among men by age. If some value is assigned, then will be assumed AIDS variation of the method.
#' @param HIV_prop_partner numeric. Estimate of the proportion of men with infected female partner.
#' @param e0_accept integer vector. Range acceptable when calculating the median of implied level by age (avoid non-possible extrapolations). By deafult between 20 and 100.
#' @export
#' @examples
#' \dontrun{
#' # Examples from IUUSP Tools for Demographic Estimation
#' # Orhpanhood 1 census
#' # Maternal orphanhood - Male respondents - Irak ("AM_Orphanhood_OneCensus_Basic - Iraq_v2.xlsx")
#' standard_iusssp <- data.frame(age = seq(15, 85, 5),
#'                               lx = c(0.8868,0.8763,0.8621,0.8467,0.8296,0.8094,0.7845, 0.7527,0.7096,0.6515,0.5724,0.4700,0.3443,0.2101,0.0962))
#' irak_1997 <- orphanhood_one(prop_not_orphan = c(0.9918,0.9819,0.9665,0.9415,0.9005,0.8374,0.7678,0.6380,0.4941),
#'                age = seq(5, 45, 5),
#'                mac = 28.2807464,
#'                date = 1997.79,
#'                mlt_data_input = standard_iusssp,
#'                brass_logit = TRUE)
#' publ_rr <- data.frame(
#'   age = seq(5, 45, 5),
#'   q15_45 = c(0.072, 0.085, 0.096, 0.110, 0.126, 0.143, 0.137, 0.159, 0.167),
#'   Date     = c(1994.2, 1992.1, 1990.1, 1988.5, 1987.0, 1985.9, 1985.5, NA, NA))
#' round(publ_rr$q15_45-irak_1997$adult_mort_index$q15_45, 2)
#' # Orphanhood HIV
#' # maternal Kenya HIV 1999 ("AM_Orphanhood_OneCensus_AIDS Kenya_v2.xlsx")
#' standard_iusssp_aids <- data.frame(age = seq(15, 85, 5),
#'                                    lx = c(0.8847,0.8603, 0.8173,0.7677, 0.7247,0.6913, 0.6645,0.6281, 0.5783,0.5110, 0.4225,0.3149, 0.2012,0.1021, 0.0377))
#' kenya_1999_maternal_aids <- orphanhood_one(
#'   prop_not_orphan = c(0.9732, 0.9560, 0.9336, 0.9080, 0.8771, 0.8244, 0.7691, 0.6685, 0.5653),
#'   age = seq(5, 45, 5),
#'   mac = 26.75048387,
#'   maternal = TRUE,
#'   date = 1999.64,
#'   mlt_data_input = standard_iusssp_aids,
#'   brass_logit = TRUE,
#'   HIV_prev_antenatal = c(.07, .01, rep(0,7)),
#'   HIV_vert_transm = 1/3,
#'   HIV_fert_ratio = 3/4)$adult_mort_index
#' rr_45q15 <- c(0.368, 0.179, 0.183, 0.189, 0.183, 0.185, 0.165, 0.161, 0.139)
#' round((kenya_1999_maternal_aids$q15_45/rr_45q15)-1,2)
#' }

orphanhood_one <- function(prop_not_orphan,
                           age,
                           maternal = TRUE,
                           mac = NULL,
                           date = NULL,
                           mlt_data_input = NULL,
                           mlt_family = "CD_West",
                           brass_logit = FALSE,
                           brass_logit_e0 = 60,
                           brass_logit_5q0 = NULL,
                           HIV_prev_antenatal = NULL,
                           HIV_vert_transm = 1/3,
                           HIV_fert_ratio = 3/4,
                           HIV_parents_surv = c(.5, .75, rep(1, length(age)-2)),
                           HIV_prop_partner = 1/2,
                           HIV_prev = NULL,
                           HIV_art = NULL,
                           e0_accept = c(20, 100),
                           verbose = TRUE){

  # initial argument checks
  if(any(is.na(prop_not_orphan) | prop_not_orphan<0 | prop_not_orphan>1)) stop("Not possible values in prop_not_orphan")
  if((!maternal & !is.null(HIV_prev_antenatal))) stop("No consistency between sex for AIDS adjustment and proportion alive by age")
  if(!is.null(mlt_data_input) & !brass_logit) stop("With custom mlt is not possible to interpolate over different UN/CD levels. Use brass_logit=TRUE")

  # get coefficients if maternal or parental
  if(maternal){
    sex <- "f"
    y <- 25
    timaes_coeff <- matrix(c(
      10,	-0.2894, 	0.00125, 	1.2559,
      15,	-0.1718,	0.00222, 	1.1123,
      20,	-0.1513,	0.00372,	1.0525,
      25,	-0.1808,	0.00586,	1.0267,
      30,	-0.2511,	0.00885,	1.0219,
      35,	-0.3644,	0.01287,	1.0380,
      40,	-0.5181,	0.01795,	1.0753,
      45,	-0.6880,	0.02342,	1.1276,
      50,	-0.8054,	0.02721,	1.1678),
      ncol=4, byrow = T,
      dimnames = list(c(), c("age", "a", "b", "c"))) %>%
      as.data.frame()

    # HIV settings
    is.HIV <- FALSE
    timaes_coeff_aids <- matrix(c(
      10, -0.3611, 	0.00125, 	1.2974,
      15, -0.4030, 	0.00222, 	1.3732,
      20, -0.2120, 	0.00372, 	1.1342,
      25, -0.2389, 	0.00586, 	1.1131,
      30, -0.2513, 	0.00885, 	1.0223,
      35,	-0.3644, 	0.01287, 	1.0380,
      40,	-0.5181, 	0.01795, 	1.0753,
      45,	-0.6880, 	0.02342, 	1.1276,
      50,	-0.8054, 	0.02721, 	1.1678),
      ncol=4, byrow = T,
      dimnames = list(c(), c("age","a", "b", "c")))

    # use pre natal clinic prev
    if(!is.null(HIV_prev_antenatal)){
      if(length(HIV_prev_antenatal) != length(prop_not_orphan)) stop("diff length between proportion alive and prevalence data")
      timaes_coeff[HIV_prev_antenatal>.05,] <- timaes_coeff_aids[HIV_prev_antenatal>.05,]
      HIV_factor_adj <- (1-HIV_vert_transm*HIV_prev_antenatal)/(1+(1-HIV_fert_ratio)/HIV_fert_ratio*HIV_prev_antenatal)
      HIV_factor_adj <- 1-(1-HIV_factor_adj) * HIV_parents_surv
      prop_not_orphan <- HIV_factor_adj * prop_not_orphan
      is.HIV <- TRUE
      # use general prev
    }else if(!is.null(HIV_prev)){
      if(length(HIV_prev) != length(prop_not_orphan)) stop("diff length between proportion alive and prevalence data")
      HIV_factor_adj <- 1-(1-(1-HIV_vert_transm)*HIV_fert_ratio)*HIV_prev*HIV_parents_surv
      prop_not_orphan <- HIV_factor_adj * prop_not_orphan
      is.HIV <- TRUE
    }
  }else{
    sex <- "m"
    y <- 35
    timaes_coeff <- matrix(c(
      10,	-0.5578, 	0.00040, 	1.4708, 	0.0698,
      15,	-0.4013, 	0.00576, 	1.5602, 	-0.3522,
      20,	-0.3329, 	0.01031, 	0.6656, 	0.3419,
      25,	-0.4726, 	0.01559, 	0.2161, 	0.7896,
      30,	-0.7056, 	0.02076, 	0.1997, 	0.9066,
      35,	-0.9153, 	0.02493, 	0.3484, 	0.8631,
      40,	-0.9950, 	0.02635, 	0.4269, 	0.8263),
      ncol=5, byrow = T,
      dimnames = list(c(), c("age", "a", "b", "c", "d"))) %>%
      as.data.frame()
    # HIV settings
    is.HIV <- FALSE
    # use general prev
    if(!is.null(HIV_prev)){
      if(length(HIV_prev) != length(prop_not_orphan)) stop("diff length between proportion alive and prevalence data")
      HIV_factor_adj <- 1-(1-(1-HIV_vert_transm)*HIV_fert_ratio)*HIV_prop_partner*HIV_prev*HIV_parents_surv
      prop_not_orphan <- HIV_factor_adj * prop_not_orphan
      is.HIV <- TRUE
    }
  }

  # check family and method (not HIV)
  mlt_families <- c("CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian")
  if(!is.null(mlt_family)){
    mlt_family <- match.arg(mlt_family, mlt_families, several.ok = TRUE)
  }else{
    mlt_family <- mlt_families
  }

  # age details
  age_int <- max(diff(age))
  OAG <- max(age)
  ages <- length(age)

  # MLT to use. If no user data as input, then use mortcast for maternal/paternal. If HIV then use Spectrum model
  if(is.null(mlt_data_input)){
    # No HIV
    if(!is.HIV){
      mlt_data <- MortCast::MLTlookup %>%
        filter(type %in% mlt_family, sex == ifelse(sex == "f", 2, 1)) %>%
        mutate(lx = lx/1e5, age_resp = "all")
      # e0 should be rounded to proximate available level
      brass_logit_e0 <- mlt_data %>% filter(e0 == round(brass_logit_e0/2.5,0)*2.5)
    # HIV
    }else{
      if(is.null(brass_logit_5q0) | is.null(HIV_art)) stop("You need 5q0 and/or HIV ART")
      if(length(brass_logit_5q0)!=length(age) | length(HIV_art)!=length(age) | length(HIV_prev)!=length(age)) stop("You need same length 5q0, prev, and art")
      #  create a range of patterns, given prev, art and 5q0.
      this_sex <- ifelse(maternal, "female", "male")
      comb <- expand.grid(age = age, q15_45 = seq(.1,.9,.1))
      mlt_data <- map2_df(comb$age, comb$q15_45, function(x, y){
        hiv_svd_comp_x <- predictNQX(this_sex,
                                     cm = brass_logit_5q0[age == x],
                                     am = y,
                                     hiv = HIV_prev[age == x],
                                     art = HIV_art[age == x],
                                     adult = "q45") %>% pull()
        lx_hiv_svd_comp_x <- lt_id_q_l(expit(hiv_svd_comp_x))
        lt_abridged(lx = lx_hiv_svd_comp_x[c(0,1,seq(5,100,5))+1], Age = c(0,1,seq(5,100,5))) %>%
          mutate(type = "HIVSpectrum", age_resp = x) %>%
          group_by(type, age_resp) %>%
          mutate(e0 = ex[Age == 0])}) %>%
          select(age_resp, type, e0, age = Age, lx) %>%
          mutate(lx = lx/1e5)
        # if apply Brass, then look closer level
        if(brass_logit){
          brass_logit_e0 <- mlt_data %>%
            group_by(type, age_resp) %>%
            filter(abs(e0-brass_logit_e0) == min(abs(e0-brass_logit_e0))) %>%
            distinct(type, age_resp, e0)
      }
    }
  # if using custom mlt
  }else{
    if(!all(colnames(mlt_data_input) %in% colnames(MortCast::MLTlookup))) stop("Be sure to have same col names and cathegories than MortCast::MLTlookup")
    mlt_data <- mlt_data_input %>% mutate(type = "user", e0 = "user", age_resp = "all")
    brass_logit_e0 <- data.frame(age = "all", age_resp = "all", type = "user", e0 = "user")
  }

  # compute probabilities to be matched against observed data
  this_mlt_family <- mlt_data %>%
    group_by(age_resp, type, e0) %>%
    mutate(ly = lx[age==y], py_n = pmin(lx/ly, 1)) %>%
    select(age_resp, type, age, e0, lx_mlt = lx, ly_mlt = ly, py_n_mlt = py_n) %>%
    ungroup() %>% arrange(age_resp, type, e0, age)

  # estimate conditional probabilities of survival using timaes coefficients
  orp_data <- data.frame(age = age, prop_not_orphan = prop_not_orphan, mac =  mac, n = age + 5) %>%
    left_join(timaes_coeff, by = c("n" = "age" )) %>%
    mutate(age_n = age + n,
           age_mean = age + age_int/2,
           maternal = maternal,
           age_yplus_n = ifelse(!maternal, pmax(50, n + y), n + y), # page 226
           prop_not_orphan_next = lead(prop_not_orphan),
           py_n = ifelse(maternal,
                         a + b * mac + c * prop_not_orphan,
                         a + b * mac + c * prop_not_orphan + d * prop_not_orphan_next)) %>%
    filter(between(py_n,0,1)) %>%
    mutate(# seems two different versions in IUSSP (2012), excel template or pdf
           time = ifelse(maternal,
                         age_mean/2*(1-log(prop_not_orphan/((1-(mac+age_mean)/80)/(1-mac/80)))/3),
                         # age_mean/2*(1-log(prop_not_orphan)/3)+log((80-M-age_mean)/(80-M))/3,
                         (n+.75)/2*(1-log(sqrt(prop_not_orphan_next*prop_not_orphan)/((1-(mac+n)/80)/(1-(mac-.75)/80)))/3)
                         # (age_mean+.75)/2*(1-log(sqrt(lead(prop_not_orphan)*prop_not_orphan))/3)+log((80-M-age_mean)/(80-M+.75))/3),
                         ),
           time_location = date - time)

  # Check probability btwn 0 and 1
  if(nrow(orp_data)==0) stop("py_n out of [0,1] range")

  # find adult mortality depending if brass smoothing
  if(brass_logit){

    # find alpha
    # HIV and no custom input mlt
    if(is.HIV & is.null(mlt_data_input)){
      mlt_closest <- orp_data %>%
        left_join(this_mlt_family %>%
                    inner_join(brass_logit_e0, by = c("age_resp", "type", "e0")),
                  by = c("age" = "age_resp", "age_yplus_n" = "age")) %>%
        ungroup %>% arrange(age, type, e0) %>%
        mutate(alpha = -.5*log(1+(py_n/lx_mlt-1/ly_mlt)/(1-py_n)))

      # lx estimate
      lx_out <- mlt_closest %>% select(-lx_mlt, -ly_mlt, -py_n_mlt) %>%
        left_join(this_mlt_family %>%
                    inner_join(brass_logit_e0, by = c("age_resp", "type", "e0")) %>%
                    select(age_resp, type, e0, age_mlt = age, lx_mlt),
                  by = c("age" = "age_resp", "type", "e0")) %>%
        mutate(lx_interp = 1 - logit_inv(alpha + logit(1-lx_mlt)))
    # no HIV or custom input mlt
    }else{
      mlt_closest <- orp_data %>%
        left_join(this_mlt_family %>%
                    filter(e0 == unique(brass_logit_e0$e0)), by = c("age_yplus_n" = "age")) %>%
        ungroup %>% arrange(type, e0, age) %>%
        mutate(alpha = -.5*log(1+(py_n/lx_mlt-1/ly_mlt)/(1-py_n)))

      # lx estimate
      lx_out <- mlt_closest %>%
        select(-lx_mlt, -ly_mlt, -py_n_mlt) %>%
        left_join(this_mlt_family %>%
                    filter(e0 == unique(brass_logit_e0$e0)) %>%
                    select(type, e0, age_mlt = age, lx_mlt), by = c("type", "e0")) %>%
        mutate(lx_interp = 1 - logit_inv(alpha + logit(1-lx_mlt)))
    }

  # no brass, interpolate with MLT
  }else{
    # HIV
    if(is.HIV){
      # find level for each npy
      obs_data <- orp_data %>% select(age_resp = age, age = age_yplus_n, py_n)
      mlt_closest <- map_dfr(obs_data$age_resp, function(x){
        obs_data_x <- obs_data %>% filter(age_resp == x) %>% select(-age_resp)
        mlt_data_x <- this_mlt_family %>% filter(age_resp == x) %>% select(type,  age, e0, py_n = py_n_mlt)
        interp_level_mlt(obs_data_x, mlt_data_x, "e0", "py_n") %>% mutate(age_resp = x)
      })

      # interpolated lx for each respondent age
      lx_out <- orp_data %>% filter(!is.na(py_n)) %>%
        left_join(mlt_closest %>% select(type, age_resp, e0_interp), by = c("age" = "age_resp")) %>%
        split(list(.$age, .$type), drop = TRUE) %>%
        map_df(function(X){
          obs_data <- data.frame(age = c(0,1,seq(5,100,5)), e0 = X$e0_interp)
          mlt_data <- this_mlt_family %>%
            filter(age_resp == unique(X$age)) %>%
            select(type,  age, e0, lx = lx_mlt)
          mlt_closest <- interp_level_mlt(obs_data, mlt_data, "lx", "e0")
          return(data.frame(X, mlt_closest %>% select(age_mlt = age, lx_interp), row.names = NULL))
        })

    }else{
      # find level for each npy
      obs_data <- orp_data %>% select(age_resp = age, age = age_yplus_n, py_n)
      mlt_data <- this_mlt_family %>% select(type,  age, e0, py_n = py_n_mlt)
      mlt_closest <- map_dfr(obs_data$age_resp, function(x){
        obs_data_x <- obs_data %>% filter(age_resp == x) %>% select(-age_resp)
        interp_level_mlt(obs_data_x, mlt_data, "e0", "py_n") %>% mutate(age_resp = x)
      })

      # interpolated lx for each respondent age
      lx_out <- orp_data %>% filter(!is.na(py_n)) %>%
        left_join(mlt_closest %>% select(type, age_resp, e0_interp), by = c("age" = "age_resp")) %>%
        split(list(.$age, .$type), drop = TRUE) %>%
        map_df(function(X){
          obs_data <- data.frame(age = c(0,1,seq(5,100,5)), e0 = X$e0_interp)
          mlt_data <- this_mlt_family %>% select(type,  age, e0, lx = lx_mlt) %>% filter(type == unique(X$type))
          mlt_closest <- interp_level_mlt(obs_data, mlt_data, "lx", "e0")
          return(data.frame(X, mlt_closest %>% select(age_mlt = age, lx_interp), row.names = NULL))
        })
    }

  }

  # return adult index
  adult_mort <- lx_out %>% select(type, age, time_location, age_mlt, lx_interp) %>%
    group_by_at(vars(-c(age_mlt, lx_interp))) %>%
    summarise(q15_45 = 1-lx_interp[age_mlt==60]/lx_interp[age_mlt==15],
              q15_35 = 1-lx_interp[age_mlt==50]/lx_interp[age_mlt==15], .groups = "keep") %>%
    ungroup()

  # list of results to return
  return(list(adult_mort_index = adult_mort,
              lx_estimates = lx_out))
}

# orphanhood two census/surveys -------------------------------------------

#' Estimate adult mortality with two survey/census observations from orphanhood data.
#' @description Follows IUSSP(2012) template implementation of Timæus (1992) for estimate adult mortality using two observations.
#' The other option is using r-variable method (Preston, 2001, p. 239).
#' @param prop1_not_orphan numeric vector. Population proportion of respondents with mother/father alive, by age groups, from source 1.
#' @param prop2_not_orphan numeric vector. Population proportion of respondents with mother/father alive, by age groups, from source 2.
#' @param date1 Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for the source 1 (typically census or survey).
#' @param date2 Same than `date1` but for the second source.
#' @param age1 integer vector. Lower bound of age groups from first census. Last age is assumed as lower age open age group.
#' @param age2 Same than `age1` but for the second source.
#' @param maternal logical. Either `TRUE` for maternal orphanhood (default) or `FALSE` for paternal orphanhood.
#' @param mac1 numeric. Mean age at childbearing for mothers/fathers from source 1.
#' @param mac2 Same than `mac1` but for the second source.
#' @param method character. Choose between `timæus` (default) (Timæus, 1991) or `preston_chen` using r-variable method (Preston, 2001, p. 239).
#' @param mlt_data_input data.frame. If is `NULL` then model life tables is used, available from \code{Morcast} package. But specific pattern can be included (`l_x` must be included see examples).
#' @param mlt_family character. Options: "CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian". If `NULL` returns for all.
#' @param brass_logit logical. Doing a level smoothing with 1-parameter logit Brass.
#' @param brass_logit_e0 numeric. Life expectancy at birth when `brass_logit` is `TRUE`.
#' @param e0_accept integer vector. Range acceptable when calculating the median of implied level by age (avoid non-possible extrapolations). By deafult between 20 and 100.
#' @export
#' @examples
#' \dontrun{
#' # Example from IUUSP Tools for Demographic Estimation
#' # Orhpanhood 2 census
#' # Paternal orphanhood - Salomons 1987 - 1999 ("AM_Orphanhood_TwoCensus - Solomons_2.xlsx")
#' salomons_1986_1999_paternal <- orphanhood_two(
#'   prop1_not_orphan = c(0.9644 ,0.9340 ,0.8862 ,0.8088 ,0.7060 ,0.5812 ,0.4620 ,0.3317),
#'   prop2_not_orphan = c(0.9710 ,0.9454 ,0.9071 ,0.8493 ,0.7717 ,0.6568 ,0.5375 ,0.3891),
#'   age1 = seq(5, 40, 5), age2 = seq(5, 40, 5),
#'   date1 = 1986.9, date2 = 1999.9,
#'   maternal = FALSE,
#'   mac1 = 34.08, mac2 = 34.08,
#'   mlt_data_input = standard_iusssp,
#'   brass_logit = TRUE)$adult_mort
#' rr_45q15 <- c(0.175, 0.185,0.187,0.187)
#' round((salomons_1986_1999_paternal$q15_45/rr_45q15)-1,2)
#' # Preston_chen method
#' # table 11.3 Preston (2001) r-var
#' panama_1977_1980 <- orphanhood_two(
#'   prop1_not_orphan = c(.996,.9844,.9753,.9556,.9195,.8908,.8157,.7618,.6166),
#'   prop2_not_orphan = c(.9949,.9903,.9794,.9608,.9316,.8967,.8450,.7746,.6781),
#'   age1 = seq(0, 40, 5),
#'   age2 = seq(0, 40, 5),
#'   date1 = "1977-07-01",
#'   date2 = "1980-05-11",
#'   maternal = TRUE,
#'   mac1 = 27,
#'   mac2 = 27,
#'   method = "preston_chen",
#'   mlt_data_input = NULL,
#'   mlt_family = "CD_West",
#'   brass_logit = FALSE,
#'   brass_logit_e0 = 60)
#' rr_levels <- c(22.4, 23, 22.7, 22.6, 22.9, 23, 23.5, 23.9)
#' }

orphanhood_two <- function(prop1_not_orphan, prop2_not_orphan,
                           age1, age2,
                           date1, date2,
                           maternal = TRUE,
                           method = "timæus",
                           mac1 = NULL, mac2 = NULL,
                           mlt_data_input = NULL,
                           mlt_family = NULL,
                           brass_logit = FALSE,
                           brass_logit_e0 = 60,
                           HIV_prev = NULL,
                           HIV_art = NULL,
                           HIV_5q0 = NULL,
                           e0_accept = c(20, 100),
                           verbose = TRUE){

  # initial checks
  if(any(is.na(prop1_not_orphan) | prop1_not_orphan<0 | prop1_not_orphan>1)) stop("Not possible values in prop1_not_orphan")
  if(any(is.na(prop2_not_orphan) | prop2_not_orphan<0 | prop2_not_orphan>1)) stop("Not possible values in prop2_not_orphan")
  if(!is.null(mlt_data_input) & !brass_logit) stop("With custom mlt is not possible to interpolate. Use brass_logit=TRUE")

  # find a close date for date1, so difference is a multiple of age_int
  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2),2)
  interc_t <- date2-date1

  # maternal or parental
  if(maternal){
    sex <- "f"
    y <- 25
    z <- 45
    timaes_coeff <- matrix(c(
      25,	-0.8623, 	0.00292, 	1.7861,
      30,	-0.3822, 	0.00679, 	1.2062,
      35,	-0.4355, 	0.01197, 	1.1310,
      40,	-0.5995, 	0.01847, 	1.1419,
      45,	-0.7984, 	0.02547, 	1.1866,
      50,	-0.9360, 	0.03039, 	1.2226),
      ncol=4, byrow = T,
      dimnames = list(c(), c("age","a", "b", "c"))) %>%
      as.data.frame()
  }else{
    sex <- "m"
    y <- 35
    z <- 55
    timaes_coeff <- matrix(c(
      25, 	-0.0554, 	0.00757, 	0.0239, 	0.8080,
      30, 	-0.7539, 	0.01558, 	0.6452, 	0.6498,
      35, 	-1.0809, 	0.02273, 	0.9289, 	0.4807,
      40, 	-1.1726, 	0.02647, 	0.9381, 	0.4372),
      ncol=5, byrow = T,
      dimnames = list(c(), c("age","a", "b", "c", "d"))) %>%
      as.data.frame()
  }

  # check family and method
  mlt_families <- c("CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean",
                    "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian")
  # which families consider
  mlt_family_input <- mlt_family
  if(!is.null(mlt_family)){
    mlt_family <- match.arg(mlt_family, mlt_families, several.ok = TRUE)
  }else{
    mlt_family <- mlt_families
  }

  # mlt data to use. If no input data, then use mortcast for maternal/paternal
  if((is.null(HIV_prev) + is.null(HIV_art)) == 0) is.HIV = TRUE
  if((is.null(HIV_prev) + is.null(HIV_art)) == 2) is.HIV = FALSE
  if((is.null(HIV_prev) + is.null(HIV_art)) == 1) stop("for hiv you need prevalence and art")
  if(is.null(mlt_data_input)){
    if(!is.HIV){
      mlt_data <- MortCast::MLTlookup %>%
        filter(type %in% mlt_family, sex == ifelse(sex == "f", 2, 1)) %>%
        mutate(lx = lx/1e5)
      # e0 should be rounded to proximate available level
      brass_logit_e0 <- round(brass_logit_e0/2.5,0)*2.5
    }else{
      if(is.null(HIV_5q0)) stop("You need 5q0 for spectrum in HIV setting")
      this_sex <- ifelse(maternal, "female", "male")
      mlt_data <- lapply(seq(.1,.9,.1), function(x){
        hiv_svd_comp_x <- predictNQX(this_sex,
                                     cm = HIV_5q0,
                                     am = x,
                                     hiv = HIV_prev,
                                     art = HIV_art,
                                     adult = "q45") %>% pull()
        lx_hiv_svd_comp_x <- lt_id_q_l(expit(hiv_svd_comp_x))
        lt_abridged(lx = lx_hiv_svd_comp_x[c(0,1,seq(5,100,5))+1], Age = c(0,1,seq(5,100,5))) %>%
          mutate(type = "HIVSpectrum") %>%
          group_by(type) %>%
          mutate(e0 = ex[Age == 0])}) %>%
        bind_rows() %>%
        select(type, e0, age = Age, lx) %>%
        mutate(lx = lx/1e5)
      # If brass
      if(brass_logit){
        actual_levels <- unique(mlt_data$e0)
        closer_level <- which(abs(actual_levels-brass_logit_e0)==min(abs(actual_levels-brass_logit_e0)))
        brass_logit_e0 <- actual_levels[closer_level]
      }
    }
  }else{
    mlt_data <- mlt_data_input %>% mutate(type = "user", e0 = "user")
    brass_logit_e0 <- "user"
  }

  # select family and standardize for same age_int, OAG and risk interval
  this_mlt_family <- mlt_data %>%
    group_by(type, e0) %>%
    mutate(lz = lx[age==z], pz_yplusn = pmin(lx/lz, 1)) %>%
    select(type, age, e0, lx_mlt = lx, lz_mlt = lz, pz_yplusn_mlt = pz_yplusn) %>%
    ungroup() %>% arrange(type, e0, age)

  # manage ages: fit to more wider group and minor OAG.
  # Always assume last age is an OAG
  if(length(age1)!=length(prop1_not_orphan) | length(age2)!=length(prop2_not_orphan)) stop("Not same interval between pop and age")
  age_int <- 5
  first_age <- max(min(age1), min(age2))
  last_age <- min(max(age1), max(age2))
  age <- seq(first_age, last_age, age_int)
  prop1_not_orphan <- prop1_not_orphan[age1 %in% age]
  prop2_not_orphan <- prop2_not_orphan[age2 %in% age]
  ages <- length(age)
  M <- (mac1 + mac2)/2

  # method selection
  if(method == "timæus"){

    # need ages
    if(!any((age+5) %in% timaes_coeff$age)) stop("For applying this method you need adult ages for respondants (IUSSP, 2012)")

    # org data
    orp_data <- data.frame(age = age,
                           prop1_not_orphan = prop1_not_orphan,
                           prop2_not_orphan = prop2_not_orphan,
                           mac1 =  mac1, mac2 =  mac2, n = age + 5) %>%
      left_join(timaes_coeff, by = c("n" = "age" )) %>%
      mutate(maternal = maternal,
             age_n = age + n,
             age_yplus_n = y + n,
             prop_avg_not_orphan = sqrt(prop1_not_orphan * prop2_not_orphan),
             mac_avg = (mac1 + mac2)/2,
             r = log(prop2_not_orphan/prop1_not_orphan)/(date2-date1),
             S20_adj = sqrt(prop_avg_not_orphan[age==15] * prop_avg_not_orphan[age==20])) %>%
      filter(age >= 20) %>%
      mutate(r_cum = cumsum(r),
             r_cum_5 = 5 * lag(r_cum, default = 0) + 2.5 * r,
             Sn_adj = prop_avg_not_orphan/S20_adj*exp(r_cum_5),
             pz_yplusn = ifelse(maternal,
                                a + b * mac_avg + c * Sn_adj,
                                a + b * mac_avg + c * Sn_adj + d * lead(Sn_adj)))

    # find adult mortality depending if brass smoothing
    if(brass_logit | !is.null(mlt_data_input)){

      mlt_closest <- orp_data %>%
        left_join(this_mlt_family %>% filter(e0 == brass_logit_e0), by = c("age_yplus_n" = "age")) %>%
        ungroup %>% arrange(type, e0, age)

      # only level
      lx_out <- mlt_closest  %>%
        mutate(alpha = -.5*log(1+(pz_yplusn/lx_mlt-1/lz_mlt)/(1-pz_yplusn))) %>%
        select(-lx_mlt) %>%
        left_join(this_mlt_family%>% filter(e0 == brass_logit_e0) %>%
                    select(type, e0, age_mlt = age, lx_mlt), by = c("type", "e0")) %>%
        mutate(logit_mlt_lx = logit(1-lx_mlt),
               lx_interp =  1 - logit_inv( alpha + logit(1-lx_mlt))) %>%
        select(age, prop1_not_orphan, prop2_not_orphan, mac1, mac2, type, alpha, age_mlt, lx_interp)

    }else{
      # find level for each np25
      obs_data <- orp_data %>% select(age = age_yplus_n, pz_yplusn)
      mlt_data <- this_mlt_family %>% select(type, age, e0, pz_yplusn = pz_yplusn_mlt)
      mlt_closest <- interp_level_mlt(obs_data, mlt_data, "e0", "pz_yplusn")

      # not possible level values
      ages_filter_out <- mlt_closest %>% filter(!between(e0_interp, e0_accept[1], e0_accept[2])) %>% pull(age)
      if(length(ages_filter_out)>0 & verbose) message(paste0("Age ", ages_filter_out, " filtered out because not acceptable mortality level."))

      # interpolated lx for each respondent age
      lx_out <- orp_data %>%
        inner_join(mlt_closest %>%
                     filter(!age %in% ages_filter_out) %>%
                     select(type, age_yplus_n = age, e0_interp))%>%
        split(list(.$type, .$age_n)) %>%
        map_df(function(X){
          obs_data <- data.frame(age = c(0,1,seq(5,100,5)), e0 = X$e0_interp)
          mlt_data <- this_mlt_family %>% select(type,  age, e0, lx = lx_mlt) %>% filter(type == unique(X$type))
          mlt_closest <- interp_level_mlt(obs_data, mlt_data, "lx", "e0")
          return(data.frame(X, mlt_closest %>% select(age_mlt = age, lx_interp), row.names = NULL))
        })
    }
  }

  # variable r method (Preston, 2001, p. 239)
  if(method == "preston_chen"){

    # preston et al (2001), page 239
    # adj prop: the proportion of non-orphaned at age x in a stationary population
    # based on the force of mortality of mothers during the intersurvey period.
    # the proportion of non-orphanedat age x in a stationary population
    # based on the force of mortality of mothers during the intersurvey period.
    orp_data <- data.frame(age = age,
                           prop1_not_orphan = prop1_not_orphan,
                           prop2_not_orphan = prop2_not_orphan,
                           mac1 =  mac1, mac2 =  mac2, n = age + 5) %>%
      mutate(maternal = maternal,
             prop_avg_not_orphan = sqrt(prop1_not_orphan * prop2_not_orphan),
             mac_avg = (mac1 + mac2)/2,
             r = log(prop2_not_orphan/prop1_not_orphan)/(date2-date1),
             r_cum = cumsum(r),
             r_cum_5 = 5 * lag(r_cum, default = 0) + 2.5 * r,
             prop_avg_not_orphan_adj = prop_avg_not_orphan*exp(r_cum_5)) %>%
      filter(prop_avg_not_orphan_adj>=0)

    # add hiv period data if is HIV
    if(is.HIV){
      HIV_prev = rep(HIV_prev, ages)
      HIV_art = rep(HIV_art, ages)
      HIV_5q0 = rep(HIV_5q0, ages)
      brass_logit_e0 = rep(brass_logit_e0, ages)
    }
    lx_out <- orphanhood_one(prop_not_orphan = orp_data$prop_avg_not_orphan_adj,
                             age = orp_data$age,
                             maternal = maternal,
                             mac = M, date = (date1+date2)/2,
                             mlt_data_input = mlt_data_input,
                             mlt_family = mlt_family_input,
                             brass_logit = brass_logit,
                             brass_logit_e0 = brass_logit_e0,
                             brass_logit_5q0 = HIV_5q0,
                             HIV_prev = HIV_prev,
                             HIV_art = HIV_art)$lx_estimates %>%
      mutate(time_location = (date1 + date2)/2) %>%
      select(age, prop_not_orphan, time_location, mac, type, age_mlt, lx_interp)
  }

  # return adult index
  lx_out <- lx_out[!is.na(lx_out$lx_interp),]
  adult_mort <- lx_out %>% ungroup %>% select(type, age, age_mlt, lx_interp) %>%
    group_by_at(vars(-c(age_mlt, lx_interp))) %>%
    summarise(q15_45 = 1-lx_interp[age_mlt==60]/lx_interp[age_mlt==15],
              q15_35 = 1-lx_interp[age_mlt==50]/lx_interp[age_mlt==15],
              .groups = "keep") %>%
    ungroup()

  # list of results to return
  return(list(adult_mort_index = adult_mort,
              lx_estimates = lx_out))
}

# hiv svd.comp model. Coded by Sarah Hertog based on model from Sam Clark (https://github.com/sinafala/svd-comp)

## sex="female"; cm=hiv.countries.f$q5; hiv=hiv.countries.f$hiv; art=hiv.countries.f$art; adult="q35"; am=NULL
predictNQX <- function(sex, cm, am=NULL, hiv, art, adult, im=NULL) {
  # sex: "female" or "male"
  # cm: vector of 5q0 values
  # am: vector of 45q15 or 35q15 values
  # hiv: vector of HIV values (%)
  # art: vector of ART values (%)
  # adult: "q45" or "q35"
  # im: vector of 1q0 values

  # Logit transform
  cml <- logit(cm)

  # predict 1q0
  if(missing(im)) {
    cmls <- cml^2
    preds.q0 <- data.frame(
      cml = as.numeric(cml),
      cmls = as.numeric(cmls)
    )
    iml <- predict(mods.r[[sex]]$q0, newdata=preds.q0)
  } else{
    iml <- logit(im)
  }

  # predict adult mx
  if(missing(am)) {
    preds.aml <- data.frame(
      cm = as.numeric(cm),
      cml = as.numeric(cml),
      hiv = as.numeric(hiv),
      art = as.numeric(art)
    )
    aml <- predict(mods.r[[sex]][[adult]]$aml, newdata=preds.aml)
  } else {
    aml <- logit(am)
  }

  preds.vs <- data.frame(
    cml = as.numeric(cml),
    aml = as.numeric(aml),
    hiv = as.numeric(hiv),
    art = as.numeric(art)
  )

  v1 <- predict(mods.r[[sex]][[adult]]$v1, newdata=preds.vs)
  v2 <- predict(mods.r[[sex]][[adult]]$v2, newdata=preds.vs)
  v3 <- predict(mods.r[[sex]][[adult]]$v3, newdata=preds.vs)
  v4 <- predict(mods.r[[sex]][[adult]]$v4, newdata=preds.vs)

  # adjust weights
  input.list <- list()
  for (i in 1:length(v1)) {
    input.list[[i]] <- list(ws = c(v1[i], v2[i], v3[i], v4[i]),
                            q5.ref = cm[i], sex = sex, qadult.ref = am[i],
                            adult = adult, q0.ref = iml[i])
    if(missing(am)) input.list[[i]]$qadult.ref <- -9999
    #            if(missing(am)) input.list[[i]]$qadult.ref <- expit(aml[i])
  }

  opt.out <- lapply(input.list, function (x) { optim(x$ws, fn = error.svd,
                                                     q5.ref = x$q5.ref, sex = x$sex, qadult.ref = x$qadult.ref,
                                                     adult = x$adult, q0.ref = x$q0.ref,
                                                     lower=c((x$ws[1] - abs(0.0001*x$ws[1])), (x$ws[2] - abs(1*x$ws[2])), (x$ws[3] - abs(.25*x$ws[3])), (x$ws[4] - abs(.25*x$ws[4]))),
                                                     upper=c((x$ws[1] + abs(0.0001*x$ws[1])), (x$ws[2] + abs(1*x$ws[2])), (x$ws[3] + abs(.25*x$ws[3])), (x$ws[4] + abs(.25*x$ws[4]))),
                                                     method="L-BFGS-B")$par})
  names(opt.out) <- 1:length(opt.out)

  # Predict qx values
  v <- matrix(unlist(opt.out), ncol = length(cml))
  r.p <- t(mods.r[[sex]]$components) %*% v + mods.r[[sex]][[adult]]$offset
  r.p <- data.frame(r.p)
  # splice in q0, predicted or original
  r.p[1,] <- iml

  return(r.p)

}

error.svd <- function(weights, sex, q5.ref, qadult.ref, adult, q0.ref) {

  adult_q <- "none"
  if (qadult.ref!=-9999 & adult == "q45") adult_q <- "q45"
  if (qadult.ref!=-9999 & adult == "q35") adult_q <- "q35"

  # Predict qx values
  r.p <- matrix(data = 0, ncol = 1, nrow = dim(mods.r[[sex]]$components)[2])
  for (z in 1:4) {
    r.p <- r.p + mods.r[[sex]]$components[z,] * weights[z]
    # r.p <- r.p + mods.r[[sex]]$components[z,] %*% t(weights[,z])
  }
  r.p <- r.p + mods.r[[sex]][[adult]]$offset

  r.p <- data.frame(r.p)

  # Splice in q0
  r.p[1,] <- q0.ref

  r.p.exp <- expit(r.p)

  # Predict q5
  q5.pred <- 1-sapply(1-r.p.exp[1:5,,drop=FALSE], prod)

  # Predict qAdult
  if (adult_q == "q45") {
    qadult.pred <- 1-sapply(1-r.p.exp[16:60,,drop=FALSE], prod)
  } else if (adult_q == "q35" ) {
    qadult.pred <- 1-sapply(1-r.p.exp[16:50,,drop=FALSE], prod)
  } else {
    qadult.pred <- qadult.ref <- 0
  }
  # Calculate and return the errors
  return((abs(q5.ref-q5.pred)*.1)+(abs(qadult.ref-qadult.pred)*.01))
}
