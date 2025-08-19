#' Estimate adult mortality with one survey/census from widowhood data.
#' @description Follows Hill and Trusell (1977) revised regression coefficients, included in Manual X from UN.
#' The function also include the possibility of using Brass logit given a pattern as default method.
#' An option for considering HIV populations is using the Spectrum model (Stover and others, 2012).
#' @param prop_not_widowed numeric vector. Population proportion of ever-married (first spouse assumption), by age groups.
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for the source (typically census or survey).
#' @param age integer vector. Lower bound of age groups. Last age is assumed as lower age open age group.
#' @param sex_spouse character. "f" for male respondent (estimate female mortality), and viceversa.
#' @param smam_f numeric. Singulate mean age at marriage for females.
#' @param smam_m numeric. Singulate mean age at marriage for males.
#' @param mlt_data_input data.frame. If is `NULL` then model life tables is used, available from \code{Morcast} package. But specific pattern can be included (`l_x` must be included see examples).
#' @param mlt_family character. Options: "CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian". If `NULL` returns for all.
#' @param brass_logit logical. Doing a level smoothing with 1-parameter logit Brass.
#' @param brass_logit_e0 numeric. Life expectancy at birth when `brass_logit` is `TRUE`.
#' @param brass_logit_5q0 numeric. Find best level when `brass_logit_e0` is `NULL`. If is not `NULL`, then implied level replaces `brass_logit_e0`.
#' @param HIV_prev numeric. Estimates of population-based estimate of HIV prevalence among men by age. If some value is assigned, then will be assumed AIDS variation of the method.
#' @param HIV_art numeric. Estimate of the proportion of men with infected female partner.
#' @export
#' @examples
#' \dontrun{

#' # Bolivia 1975. Example from UN Manual X.
#' manualX_bolivia_male <- widowhood_one(prop_not_widowed = c(0.9798, 0.9729, 0.9514, 0.9170, 0.8735, 0.8195, 0.7054, 0.6520),
#'                                       age = seq(20, 55, 5),
#'                                       sex_spouse = "m",
#'                                       smam_f = 23.2 ,
#'                                       smam_m = 25.3 ,
#'                                       date = 1975.6,
#'                                       mlt_family = "CD_West",
#'                                       brass_logit = FALSE)
#' un_results_py_n <- c(.956, .931, .894, .85, .798, .687, .638)
#' un_results_level <- c(16.3, 16.6, 15.7, 15.4, 15.3, 13.3, 14.7)
#' round(unique(manualX_bolivia_male$lx_estimates$py_n) - un_results_py_n,3)
#' un_results_e0 <- approx(13:17, c(47.082, 49.546, 51.816, 54.122, 56.45), un_results_level)$y
#' round(unique(manualX_bolivia_male$lx_estimates$e0_interp) - un_results_e0,1) # ok with MORTPACK TOO
#'}
widowhood_one <- function(prop_not_widowed,
                           age,
                           sex_spouse = "f",
                           smam_f = NULL,
                           smam_m = NULL,
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

  # coefficents
  if(sex_spouse == "m"){
    # table 97
    hill_trussel_coeff <- matrix(c(
      25,	 0.1082,  -0.00209,	0.00072, 0.9136,
      30,	-0.0284,	-0.00465,	0.00157, 1.0822,
      35,	-0.0159,	-0.00638,	0.00253, 1.0831,
      40,	 0.0041,	-0.00784,	0.00395, 1.0596,
      45,	 0.0152,	-0.00953,	0.00611, 1.0324,
      50,	 0.0087,	-0.01189,	0.00925, 1.0144,
      55, -0.0169,  -0.01515, 0.01353, 1.0111,
      60, -0.0590,  -0.01940, 0.01880, 1.0291),
      ncol=5, byrow = T,
      dimnames = list(c(), c("age","a", "b", "c", "d"))) %>%
      as.data.frame()
  }else{
    # table 98
    hill_trussel_coeff <- matrix(c(
      25,	-0.0208,  0.00052,	-0.00137, 1.0451,
      30,	-0.2135,	0.00104,	-0.00329, 1.2791,
      35,	-0.1896,	0.00162,	-0.00492, 1.2884,
      40,	-0.1290,	0.00236,	-0.00624, 1.2483,
      45,	-0.0713,	0.00340,	-0.00624, 1.2005,
      50,	-0.0327,	0.00502,  -0.00860, 1.1590,
      55, -0.0139,  0.00749,  -0.01019, 1.1297),
      ncol=5, byrow = T,
      dimnames = list(c(), c("age","a", "b", "c", "d"))) %>%
      as.data.frame()
  }

  # is hiv
  is.HIV <- FALSE
  if(!is.null(HIV_prev)){
    if(length(HIV_prev) != length(prop_not_widowed)) stop("diff length between proportion alive and prevalence data")
    is.HIV <- TRUE
  }

  # initial argument checks
  if(any(is.na(prop_not_widowed) | prop_not_widowed<0 | prop_not_widowed>1)) stop("Not possible values in prop_not_widowed")
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
        filter(type %in% mlt_family, sex == ifelse(sex_spouse == "f", 2, 1)) %>%
        mutate(lx = lx/1e5, age_resp = "all")
      # e0 should be rounded to proximate available level
      brass_logit_e0 <- mlt_data %>% filter(e0 == round(brass_logit_e0/2.5,0)*2.5)
      # HIV
    }else{
      if(is.null(brass_logit_5q0) | is.null(HIV_art)) stop("needs 5q0 and/or HIV ART")
      if(length(brass_logit_5q0)!=length(age) | length(HIV_art)!=length(age) | length(HIV_prev)!=length(age)) stop("needs same length 5q0, prev, and art")
      #  create a range of patterns, given prev, art and 5q0.
      this_sex <- ifelse(sex_spouse == "f", "female", "male")
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
  y <- 20
  this_mlt_family <- mlt_data %>%
    group_by(age_resp, type, e0) %>%
    mutate(ly = lx[age==y], py_n = pmin(lx/ly, 1)) %>%
    select(age_resp, type, age, e0, lx_mlt = lx, ly_mlt = ly, py_n_mlt = py_n) %>%
    ungroup() %>% arrange(age_resp, type, e0, age)

  # estimate conditional probabilities of survival using coefficients
  wid_data <- data.frame(age = age, prop_not_widowed = prop_not_widowed,
                         smam_f =  smam_f, smam_m =  smam_m, n = age + 5) %>%
    left_join(hill_trussel_coeff, by = c("n" = "age" )) %>%
    mutate(age_mean = age + age_int/2,
           sex_spouse = sex_spouse,
           prop_not_widowed_next = lead(prop_not_widowed),
           py_n = ifelse(sex_spouse == "f",
                         a + b * smam_f + c * smam_m + d * prop_not_widowed_next,
                         a + b * smam_f + c * smam_m + d * prop_not_widowed)) %>%
    filter(between(py_n, 0, 1)) %>%
    mutate(u = ifelse(sex_spouse == "f",
                      1/3*log(prop_not_widowed_next) + interp_Z(smam_f+n+2.5-smam_m) + .0037*(27-smam_f),
                      1/3*log(prop_not_widowed)      + interp_Z(smam_m+n-2.5-smam_f) + .0037*(27-smam_m)),
           time = ifelse(sex_spouse == "f",
                         (n+2.5-smam_m)*(1-u)/2,
                         (n-2.5-smam_f)*(1-u)/2),
           time_location = date - time) %>%
    filter(time > 0) # avoid not possible time locations from extrapolate from Z

  # no guarantee probability btwn 0 and 1
  if(nrow(wid_data)==0) stop("py_n out of [0,1] range for all ages")

  # get mlt level for each age resp
  if(brass_logit){
    # find alpha
    # HIV and no custom input mlt
    if(is.HIV & is.null(mlt_data_input)){
      mlt_closest <- wid_data %>%
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
      mlt_closest <- wid_data %>%
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
    obs_data <- wid_data %>% select(age_resp = age, age = n, py_n)
    mlt_data <- this_mlt_family %>% select(type,  age, e0, py_n = py_n_mlt)
    mlt_closest <- map_dfr(obs_data$age_resp, function(x){
      obs_data_x <- obs_data %>% filter(age_resp == x) %>% select(-age_resp)
      interp_level_mlt(obs_data_x, mlt_data, "e0", "py_n") %>% mutate(age_resp = x)
    })

    # interpolated lx for each respondent age
    lx_out <- wid_data %>% filter(!is.na(py_n)) %>%
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


# get_widowhood_coeff <- function(sex = "f"){
#   if(sex == "m"){
#     # table 97
#     hill_trussel_coeff <- matrix(c(
#       25,	 0.1082,  -0.00209,	0.00072, 0.9136,
#       30,	-0.0284,	-0.00465,	0.00157, 1.0822,
#       35,	-0.0159,	-0.00638,	0.00253, 1.0831,
#       40,	 0.0041,	-0.00784,	0.00395, 1.0596,
#       45,	 0.0152,	-0.00953,	0.00611, 1.0324,
#       50,	 0.0087,	-0.01189,	0.00925, 1.0144,
#       55, -0.0169,  -0.01515, 0.01353, 1.0111,
#       60, -0.0590,  -0.01940, 0.01880, 1.0291),
#       ncol=5, byrow = T,
#       dimnames = list(c(), c("age","a", "b", "c", "d"))) %>%
#       as.data.frame()
#   }else{
#     # table 98
#     hill_trussel_coeff <- matrix(c(
#       25,	-0.0208,  0.00052,	-0.00137, 1.0451,
#       30,	-0.2135,	0.00104,	-0.00329, 1.2791,
#       35,	-0.1896,	0.00162,	-0.00492, 1.2884,
#       40,	-0.1290,	0.00236,	-0.00624, 1.2483,
#       45,	-0.0713,	0.00340,	-0.00624, 1.2005,
#       50,	-0.0327,	0.00502,  -0.00860, 1.1590,
#       55, -0.0139,  0.00749,  -0.01019, 1.1297),
#       ncol=5, byrow = T,
#       dimnames = list(c(), c("age","a", "b", "c", "d"))) %>%
#       as.data.frame()
#   }
#   return(hill_trussel_coeff)
# }
