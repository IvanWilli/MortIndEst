#' Estimate adult mortality with one survey/census from sibling survival data.
#' @description Follows IUSSP(2012) template implementation for estimate adult mortality from siblings data, based on Timæus, Zaba and Ali (2001).
#' The function also include the possibility to match the closest model life table (not using Brass logit given a pattern as default method), for some specific family (if no family is set, then return results for all CD and UN).
#' An option for considering HIV populations is using the Spectrum model (Stover and others, 2012).
#' #' @param prop_not_dead numeric vector. Population proportion of respondents with mother/father alive, by age groups.
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for the source (typically census or survey).
#' @param age integer vector. Lower bound of age groups from first census. Last age is assumed as lower age open age group.
#' @param sex_siblings logical. Either `TRUE` for maternal orphanhood (default) or `FALSE` for paternal orphanhood.
#' @param mlt_data_input data.frame. If is `NULL` then model life tables is used, available from \code{Morcast} package. But specific pattern can be included (`l_x` must be included see examples).
#' @param mlt_family character. Options: "CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian". If `NULL` returns for all.
#' @param brass_logit logical. Doing a level smoothing with 1-parameter logit Brass.
#' @param brass_logit_e0 numeric. Life expectancy at birth when `brass_logit` is `TRUE`.
#' @param brass_logit_5q0 numeric. Find best level when `brass_logit_e0` is `NULL`. If is not `NULL`, then implied level replaces `brass_logit_e0`.
#' @param HIV_prev numeric. Estimates of population-based estimate of HIV prevalence among men by age. If some value is assigned, then will be assumed AIDS variation of the method.
#' @param HIV_art numeric. Estimate of the proportion of men with infected female partner.
#' @param e0_accept integer vector. Range acceptable when calculating the median of implied level by age (avoid non-possible extrapolations). By deafult between 20 and 100.
#' @export
#' @examples
#' \dontrun{
#' # Examples from IUUSP (2012) Tools for Demographic Estimation
#' # Bangladesh 2003, all brothers. "AM_Siblings_Indirect - Bangladesh_2_5.xlsx"
#' standard_iusssp <- c(-0.8899, -0.8600,-0.8202,-0.7801, -0.7371, -0.6902, -0.6345, -0.5682, -0.4824, -0.3737)
#' standard_iusssp <- data.frame(age = seq(15, 60, 5), lx = 1/(1+exp(2*standard_iusssp)))
#' iussp_all_brothers <- siblings_one(prop_not_dead = c(0.946425,0.943374,0.941222,0.935335,0.928475,0.882786,0.832533),
#'                                    age = seq(15, 45, 5),
#'                                    sex_siblings = "m",
#'                                    date = 2003.4,
#'                                    mlt_data_input = standard_iusssp,
#'                                    brass_logit = TRUE)
#' iussp_results <- c(0.417, 0.344, 0.288, 0.242, 0.285, 0.306)
#' round(iussp_all_brothers$adult_mort_index$q15_45 - iussp_results,3)
#'
#' # Bangladesh 2003, all sisters. "AM_Siblings_Indirect - Bangladesh_2_5.xlsx"
#' iussp_all_sisters <- siblings_one(prop_not_dead = c(0.9708,0.9447,0.9313,0.9180,0.9010,0.8748,0.8189),
#'                                   age = seq(15, 45, 5),
#'                                   sex_siblings = "f",
#'                                   date = 2003.4,
#'                                   mlt_data_input = standard_iusssp,
#'                                   brass_logit = TRUE)
#' iussp_results <- c(0.411, 0.384, 0.345, 0.315, 0.301, 0.327)
#' round(iussp_all_sisters$adult_mort_index$q15_45 - iussp_results,3)
#' }

siblings_one <- function(prop_not_dead,
                           age,
                           sex_siblings = "f",
                           date = NULL,
                           mlt_data_input = NULL,
                           mlt_family = "CD_West",
                           brass_logit = FALSE,
                           brass_logit_e0 = 60,
                           brass_logit_5q0 = NULL,
                           HIV_prev = NULL,
                           HIV_art = NULL,
                           e0_accept = c(20, 100),
                           verbose = TRUE){

  # coefficents - Timaeus et al., 2001
    timaeus_coeff <- matrix(c(
      25,	-0.0003, 1.0011,	3.23,  1.12,
      30,	-0.1546, 1.1560,	5.46,  1.95,
      35,	-0.1645, 1.1660,	7.52,  2.78,
      40,	-0.1388, 1.1406,	9.38,  3.62,
      45,	-0.1140, 1.1168,	11,    4.45,
      50,	-0.1018, 1.1066,	12.32, 5.28),
      ncol=5, byrow = T,
      dimnames = list(c(), c("age","aS", "bS", "aT", "bT"))) %>%
      as.data.frame()

  # is hiv
  is.HIV <- FALSE
  if(!is.null(HIV_prev)){
    if(length(HIV_prev) != length(prop_not_dead)) stop("diff length between proportion alive and prevalence data")
    is.HIV <- TRUE
  }

  # initial argument checks
  if(any(is.na(prop_not_dead) | prop_not_dead<0 | prop_not_dead>1)) stop("Not possible values in prop_not_dead")
  if(!is.null(mlt_data_input) & !brass_logit) stop("With custom mlt is not possible to interpolate over different UN/CD levels. Use brass_logit=TRUE")

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

  # mlt data to use. If no input data, then use mortcast. If HIV then use Spectrum model
  if(is.null(mlt_data_input)){
    # No HIV
    if(!is.HIV){
      mlt_data <- MortCast::MLTlookup %>%
        filter(type %in% mlt_family, sex == ifelse(sex_siblings == "f", 2, 1)) %>%
        mutate(lx = lx/1e5, age_resp = "all")
      # e0 should be rounded to proximate available level
      brass_logit_e0 <- mlt_data %>% filter(e0 == round(brass_logit_e0/2.5,0)*2.5)
    # HIV
    }else{
      if(is.null(brass_logit_5q0) | is.null(HIV_art)) stop("needs 5q0 and/or HIV ART")
      if(length(brass_logit_5q0)!=length(age) | length(HIV_art)!=length(age) | length(HIV_prev)!=length(age)) stop("needs same length 5q0, prev, and art")
      #  create a range of patterns, given prev, art and 5q0.
      this_sex <- ifelse(sex_siblings == "f", "female", "male")
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
  # custom mlt
  }else{
    if(!all(colnames(mlt_data_input) %in% colnames(MortCast::MLTlookup))) stop("Be sure to have same col names and cathegories than MortCast::MLTlookup")
    mlt_data <- mlt_data_input %>% mutate(type = "user", e0 = "user", age_resp = "all")
    brass_logit_e0 <- data.frame(age = "all", age_resp = "all", type = "user", e0 = "user")
  }

  # compute probabilities to be matched against observed data
  y <- 15
  this_mlt_family <- mlt_data %>%
    group_by(age_resp, type, e0) %>%
    mutate(ly = lx[age==y], py_n = pmin(lx/ly, 1)) %>%
    select(age_resp, type, age, e0, lx_mlt = lx, ly_mlt = ly, py_n_mlt = py_n) %>%
    ungroup() %>% arrange(age_resp, type, e0, age)

  # estimate conditional probabilities of survival using coefficients
  sibl_data <- data.frame(age = age, prop_not_dead = prop_not_dead, n = age + 5) %>%
    left_join(timaeus_coeff, by = c("n" = "age" )) %>%
    mutate(age_mean = age + age_int/2,
           sex_siblings = sex_siblings,
           py_n = aS + bS * prop_not_dead) %>%
    filter(between(py_n, 0, 1)) %>%
    mutate(time = aT - bT * log(prop_not_dead),
           time_location = date - time)

  # no guarantee probability btwn 0 and 1
  if(nrow(sibl_data)==0) stop("py_n out of [0,1] range for all ages")

  # get mlt level for each age resp
  if(brass_logit){
    # find alpha
    # HIV and no custom input mlt
    if(is.HIV & is.null(mlt_data_input)){
      mlt_closest <- sibl_data %>%
        left_join(this_mlt_family %>% inner_join(brass_logit_e0, by = c("age_resp", "type", "e0")),
                  by = c("age" = "age_resp", "n" = "age")) %>%
        ungroup %>% arrange(age, type, e0) %>%
        mutate(alpha = -.5*log(1+(py_n/lx_mlt-1/ly_mlt)/(1-py_n)))

      # lx estimate
      lx_out <- mlt_closest %>% select(-lx_mlt, -ly_mlt, -py_n_mlt) %>%
        left_join(this_mlt_family %>% inner_join(brass_logit_e0, by = c("age_resp", "type", "e0")) %>%
                    select(age_resp, type, e0, age_mlt = age, lx_mlt),
                  by = c("age" = "age_resp", "type", "e0")) %>%
        mutate(lx_interp = 1 - logit_inv(alpha + logit(1-lx_mlt)))

    # no HIV or custom input mlt
    }else{
    mlt_closest <- sibl_data %>%
      left_join(this_mlt_family %>%
                  filter(e0 == unique(brass_logit_e0$e0)), by = c("n" = "age")) %>%
      ungroup %>% arrange(type, e0, age) %>%
      mutate(alpha = -.5*log(1+(py_n/lx_mlt-1/ly_mlt)/(1-py_n)))

    # lx estimate
    lx_out <- mlt_closest %>% select(-lx_mlt, -ly_mlt, -py_n_mlt) %>%
      left_join(this_mlt_family %>%
                  filter(e0 == unique(brass_logit_e0$e0)) %>%
                  select(type, e0, age_mlt = age, lx_mlt), by = c("type", "e0")) %>%
      mutate(lx_interp = 1 - logit_inv(alpha + logit(1-lx_mlt)))
    }
  }else{
    # find level for each npy
    obs_data <- sibl_data %>% select(age_resp = age, age = n, py_n)
    mlt_data <- this_mlt_family %>% select(type,  age, e0, py_n = py_n_mlt)
    mlt_closest <- map_dfr(obs_data$age_resp, function(x){
      obs_data_x <- obs_data %>% filter(age_resp == x) %>% select(-age_resp)
      interp_level_mlt(obs_data_x, mlt_data, "e0", "py_n") %>% mutate(age_resp = x)
    })

    # interpolated lx for each respondent age
    lx_out <- sibl_data %>% filter(!is.na(py_n)) %>%
      left_join(mlt_closest %>% select(type, age_resp, e0_interp), by = c("age" = "age_resp")) %>%
      split(list(.$age, .$type), drop = TRUE) %>%
      map_df(function(X){
        obs_data <- data.frame(age = c(0,1,seq(5,100,5)), e0 = X$e0_interp)
        mlt_data <- this_mlt_family %>% select(type,  age, e0, lx = lx_mlt) %>% filter(type == unique(X$type))
        mlt_closest <- interp_level_mlt(obs_data, mlt_data, "lx", "e0")
        return(data.frame(X, mlt_closest %>% select(age_mlt = age, lx_interp), row.names = NULL))
      })
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

