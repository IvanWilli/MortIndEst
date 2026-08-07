# orphanhood_two_b: two-observation maternal orphanhood with HIV/ART correction
#
# This is intentionally kept separate from orphanhood_two(). It implements the
# Masquelier and Timaeus HIV/AIDS adaptation without changing the existing
# package function or NAMESPACE.

orphanhood_two_b <- function(prop1_not_orphan, prop2_not_orphan,
                             age1, age2,
                             date1, date2,
                             maternal = TRUE,
                             mac1 = NULL, mac2 = NULL,
                             mlt_data_input = NULL,
                             mlt_family = "CD_West",
                             brass_logit = FALSE,
                             brass_logit_e0 = 60,
                             HIV_prev = NULL,
                             HIV_art = NULL,
                             HIV_pmtct = NULL,
                             HIV_5q0 = NULL,
                             HIV_period_covariates = c("midpoint", "average"),
                             HIV_age_above30 = c("identity_if_low_exposure", "drop", "identity"),
                             HIV_birth_exposure_tol = 1e-3,
                             e0_accept = c(20, 100),
                             verbose = TRUE){

  # initial checks -----------------------------------------------------------
  if(any(is.na(prop1_not_orphan) | prop1_not_orphan < 0 | prop1_not_orphan > 1)) stop("Not possible values in prop1_not_orphan")
  if(any(is.na(prop2_not_orphan) | prop2_not_orphan < 0 | prop2_not_orphan > 1)) stop("Not possible values in prop2_not_orphan")

  if(!maternal) stop("Masquelier-Timaeus HIV adaptation is only for maternal orphanhood")
  if(is.null(mac1) | is.null(mac2)) stop("mac1 and mac2 are required")
  if(is.null(HIV_prev) | is.null(HIV_art) | is.null(HIV_pmtct)) stop("HIV_prev, HIV_art, and HIV_pmtct are required")
  if(!is.null(mlt_data_input) & !brass_logit) stop("With custom mlt is not possible to interpolate. Use brass_logit=TRUE")

  HIV_period_covariates <- match.arg(HIV_period_covariates)
  HIV_age_above30 <- match.arg(HIV_age_above30)

  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1), 2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2), 2)
  interc_t <- date2 - date1
  if(interc_t <= 0) stop("date2 must be greater than date1")
  if((interc_t < 3 | interc_t >= 11) & verbose) {
    message("Masquelier and Timaeus combine inquiries separated by at least 3 and less than 11 years.")
  }

  # Small local reader for scalar, vector, named-vector, or year/value data.
  value_at <- function(x, years, nm){
    if(is.data.frame(x)){
      if(!all(c("year", "value") %in% names(x))) stop(nm, " data.frame must have columns year and value")
      x <- x[order(x$year), ]
      return(stats::approx(x$year, x$value, xout = years, rule = 2)$y)
    }
    if(!is.null(names(x))){
      x_year <- suppressWarnings(as.numeric(names(x)))
      if(!any(is.na(x_year))){
        ord <- order(x_year)
        return(stats::approx(x_year[ord], as.numeric(x)[ord], xout = years, rule = 2)$y)
      }
    }
    if(length(x) == 1) return(rep(as.numeric(x), length(years)))
    if(length(x) == length(years)) return(as.numeric(x))
    stop(nm, " must be scalar, same length as requested years, named by year, or data.frame(year, value)")
  }

  check_prop <- function(x, nm){
    if(any(is.na(x) | x < 0 | x > 1)) stop(nm, " must contain proportions in [0, 1], not percentages")
  }

  # Coefficients from Masquelier and Timaeus, Tables 2 and 3.
  hiv_adj_coeff <- data.frame(
    age = c(0, 5, 10, 15, 20, 25, 30),
    b0 = c(1.0000, 1.0000, 1.0000, 1.0000, 0.9999, 0.9999, 1.0000),
    b1 = c(-0.0813, -0.2296, -0.3594, -0.4672, -0.5200, -0.5500, -0.5550),
    b2 = c(-0.0013, 0.0032, 0.0164, 0.0264, 0.0154, 0.0054, 0.0007)
  )

  hiv_surv_coeff <- data.frame(
    n = c(10, 15, 20, 25, 30, 35, 40, 45),
    b0 = c(-0.274189617817127, -0.148848703703217, -0.112309508695488, -0.133634622472498,
           -0.205065486280316, -0.307392285323096, -0.457986086528175, -0.624366557287962),
    b1 = c(0.000993598298186803, 0.00155894922231065, 0.00269536818964186, 0.00466106703035042,
           0.00753906204322899, 0.0112726481404247, 0.0162149450809221, 0.0214859022382269),
    b2 = c(1.24553666114657, 1.10531464967079, 1.03873329415231, 1.00822655497049,
           1.00780192924595, 1.01820155115521, 1.05299881802666, 1.10311550541285),
    b3 = c(-0.203921717036762, -0.261120474591671, -0.274723448434756, -0.284092135558322,
           -0.270108805561241, -0.221259863398156, -0.132314934149978, -0.0342885657390863),
    b4 = c(0.260879782530258, 0.353923879422056, 0.371063447646154, 0.33857804831828,
           0.29949901737767, 0.264892594865085, 0.176954263183001, 0.0670584923547602),
    b5 = c(-0.353182881308647, -0.239318327740746, -0.169271163226226, -0.155896002995966,
           -0.191526099927397, -0.270792703796941, -0.380139154781085, -0.511864034685874),
    b6 = c(0.00207138273062763, 0.00258063107334759, 0.00311851288042298, 0.00449210831105194,
           0.00686496899967131, 0.00996070146629463, 0.0135930556543234, 0.0175501628270404),
    b7 = c(1.30296729042681, 1.17315156049862, 1.0853955190426, 1.0441897859521,
           1.02228671726323, 1.02410122967073, 1.05360033966495, 1.1132092238738),
    b8 = c(-0.133622307716433, -0.105560680415059, -0.107118206258426, -0.18663730365779,
           -0.246632431328696, -0.221979415795482, -0.170765341251885, -0.13182557837591),
    b9 = c(-0.00173920887716492, 0.0183429690620448, -0.0329664262975753, -0.114712444607727,
           -0.0750762170019715, -0.0228384661714024, -0.0690907462489522, -0.238188567817781)
  )

  # check family
  mlt_families <- c("CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean",
                    "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian")
  if(!is.null(mlt_family)){
    mlt_family <- match.arg(mlt_family, mlt_families, several.ok = TRUE)
  }else{
    mlt_family <- mlt_families
  }

  # manage ages; keep age 0 for the synthetic-cohort chaining, but Table 3
  # only converts estimates for respondent ages 5-40 (n = 10-45).
  if(length(age1) != length(prop1_not_orphan) | length(age2) != length(prop2_not_orphan)) stop("Not same interval between pop and age")
  age <- seq(max(min(age1), min(age2)), min(max(age1), max(age2)), 5)
  if(!all(age %in% age1) | !all(age %in% age2)) stop("age1 and age2 must contain common 5-year age groups")

  prop1_not_orphan <- prop1_not_orphan[age1 %in% age]
  prop2_not_orphan <- prop2_not_orphan[age2 %in% age]

  if(is.data.frame(mac1)){
    mac1 <- value_at(mac1, date1 - (age + 2.5), "mac1")
  }else if(length(mac1) == 1){
    mac1 <- rep(mac1, length(age))
  }else if(length(mac1) == length(age1)){
    mac1 <- mac1[age1 %in% age]
  }else if(length(mac1) != length(age)){
    stop("mac1 must be scalar, same length as age1/common age, or data.frame(year, value)")
  }

  if(is.data.frame(mac2)){
    mac2 <- value_at(mac2, date2 - (age + 2.5), "mac2")
  }else if(length(mac2) == 1){
    mac2 <- rep(mac2, length(age))
  }else if(length(mac2) == length(age2)){
    mac2 <- mac2[age2 %in% age]
  }else if(length(mac2) != length(age)){
    stop("mac2 must be scalar, same length as age2/common age, or data.frame(year, value)")
  }

  # HIV-related correction to proportions of mothers surviving (Table 2).
  adj_long <- rbind(
    data.frame(source = "one",
               age = age,
               prop_not_orphan = prop1_not_orphan,
               date = date1,
               birth_year = date1 - (age + 2.5),
               mac = mac1),
    data.frame(source = "two",
               age = age,
               prop_not_orphan = prop2_not_orphan,
               date = date2,
               birth_year = date2 - (age + 2.5),
               mac = mac2)
  ) %>%
    dplyr::mutate(HIV_birth = value_at(HIV_prev, birth_year, "HIV_prev"),
                  PMTCT_birth = value_at(HIV_pmtct, birth_year, "HIV_pmtct"),
                  ART_survey = value_at(HIV_art, date, "HIV_art")) %>%
    dplyr::left_join(hiv_adj_coeff, by = "age")

  check_prop(adj_long$HIV_birth, "HIV_prev")
  check_prop(adj_long$PMTCT_birth, "HIV_pmtct")
  check_prop(adj_long$ART_survey, "HIV_art")

  no_adj_coeff <- is.na(adj_long$b0)
  if(any(no_adj_coeff)){
    if(HIV_age_above30 == "drop"){
      if(verbose) message("Dropping ages without Table 2 HIV adjustment coefficients: ",
                          paste(sort(unique(adj_long$age[no_adj_coeff])), collapse = ", "))
      adj_long <- adj_long[!no_adj_coeff, ]
    }else if(HIV_age_above30 == "identity"){
      adj_long$b0[no_adj_coeff] <- 1
      adj_long$b1[no_adj_coeff] <- 0
      adj_long$b2[no_adj_coeff] <- 0
    }else{
      low_exposure <- adj_long$HIV_birth * (1 - adj_long$PMTCT_birth) <= HIV_birth_exposure_tol
      drop_rows <- no_adj_coeff & !low_exposure
      if(any(drop_rows)){
        if(verbose) message("Dropping ages without Table 2 coefficients and non-negligible birth HIV exposure: ",
                            paste(sort(unique(adj_long$age[drop_rows])), collapse = ", "))
        adj_long <- adj_long[!drop_rows, ]
        no_adj_coeff <- is.na(adj_long$b0)
      }
      adj_long$b0[no_adj_coeff] <- 1
      adj_long$b1[no_adj_coeff] <- 0
      adj_long$b2[no_adj_coeff] <- 0
    }
  }

  adj_long <- adj_long %>%
    dplyr::mutate(HIV_birth_exposure = HIV_birth * (1 - PMTCT_birth),
                  HIV_adj_ratio = pmin(b0 + b1 * HIV_birth_exposure + b2 * ART_survey, 1),
                  prop_not_orphan_hiv_adj = pmin(pmax(prop_not_orphan * HIV_adj_ratio, 0), 1))

  orp_data <- adj_long %>%
    dplyr::select(source, age, prop_not_orphan, prop_not_orphan_hiv_adj, HIV_adj_ratio,
                  HIV_birth, PMTCT_birth, ART_survey, mac) %>%
    tidyr::pivot_wider(names_from = source,
                       values_from = c(prop_not_orphan, prop_not_orphan_hiv_adj,
                                       HIV_adj_ratio, HIV_birth, PMTCT_birth,
                                       ART_survey, mac),
                       names_sep = "") %>%
    dplyr::arrange(age) %>%
    dplyr::mutate(n = age + 5,
                  mac_avg = (macone + mactwo) / 2,
                  prop_avg_not_orphan_hiv_adj = sqrt(prop_not_orphan_hiv_adjone * prop_not_orphan_hiv_adjtwo),
                  r = log(prop_not_orphan_hiv_adjtwo / prop_not_orphan_hiv_adjone) / interc_t,
                  r_cum = cumsum(r),
                  r_cum_5 = 5 * dplyr::lag(r_cum, default = 0) + 2.5 * r,
                  S_synthetic = pmin(pmax(prop_avg_not_orphan_hiv_adj * exp(r_cum_5), 0), 1))

  # Table 3 conversion of synthetic proportions into n_p_25.
  hiv1 <- value_at(HIV_prev, date1, "HIV_prev")
  hiv2 <- value_at(HIV_prev, date2, "HIV_prev")
  art1 <- value_at(HIV_art, date1, "HIV_art")
  art2 <- value_at(HIV_art, date2, "HIV_art")
  hiv_period <- if(HIV_period_covariates == "midpoint"){
    value_at(HIV_prev, (date1 + date2) / 2, "HIV_prev")
  }else{
    mean(c(hiv1, hiv2))
  }
  art_period <- if(HIV_period_covariates == "midpoint"){
    value_at(HIV_art, (date1 + date2) / 2, "HIV_art")
  }else{
    mean(c(art1, art2))
  }

  check_prop(c(hiv1, hiv2, hiv_period), "HIV_prev")
  check_prop(c(art1, art2, art_period), "HIV_art")

  d_hiv <- abs(hiv2 - hiv1)
  hiv_untreated_period <- hiv_period * (1 - art_period)
  # The supplementary workbook uses the signed change for this term.
  d_hiv_untreated <- hiv2 * (1 - art2) - hiv1 * (1 - art1)

  p25_data <- orp_data %>%
    dplyr::filter(age >= 5, age <= 40) %>%
    dplyr::left_join(hiv_surv_coeff, by = "n")

  if(any(is.na(p25_data$b0))) stop("No Table 3 coefficients for at least one respondent age")

  p25_data <- p25_data %>%
    dplyr::mutate(HIV_period = hiv_period,
                  ART_period = art_period,
                  dHIV = d_hiv,
                  HIV_untreated_period = hiv_untreated_period,
                  dHIV_untreated = d_hiv_untreated,
                  p25_n = dplyr::if_else(
                    ART_period <= 0,
                    b0 + b1 * mac_avg + b2 * S_synthetic + b3 * HIV_period + b4 * dHIV,
                    b5 + b6 * mac_avg + b7 * S_synthetic +
                      b8 * HIV_untreated_period + b9 * dHIV_untreated
                  ),
                  age_25n = 25 + n,
                  time_location = (date1 + date2) / 2) %>%
    dplyr::filter(dplyr::between(p25_n, 0, 1))

  if(nrow(p25_data) == 0) stop("p25_n out of [0,1] range")

  # MLT selection and conversion to full lx, following orphanhood_two style.
  if(is.null(mlt_data_input)){
    if(is.null(HIV_5q0)){
      mlt_data <- MortCast::MLTlookup %>%
        dplyr::filter(type %in% mlt_family, sex == 2) %>%
        dplyr::mutate(lx = lx / 1e5)
      brass_logit_e0 <- round(brass_logit_e0 / 2.5, 0) * 2.5
    }else{
      mlt_data <- lapply(seq(.1, .9, .1), function(x){
        hiv_svd_comp_x <- predictNQX("female",
                                     cm = HIV_5q0,
                                     am = x,
                                     hiv = hiv_period,
                                     art = art_period,
                                     adult = "q45") %>% dplyr::pull()
        lx_hiv_svd_comp_x <- lt_id_q_l(expit(hiv_svd_comp_x))
        DemoTools::lt_abridged(lx = lx_hiv_svd_comp_x[c(0, 1, seq(5, 100, 5)) + 1],
                               Age = c(0, 1, seq(5, 100, 5))) %>%
          dplyr::mutate(type = "HIVSpectrum") %>%
          dplyr::group_by(type) %>%
          dplyr::mutate(e0 = ex[Age == 0])
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::select(type, e0, age = Age, lx) %>%
        dplyr::mutate(lx = lx / 1e5)
      if(brass_logit){
        actual_levels <- unique(mlt_data$e0)
        brass_logit_e0 <- actual_levels[which(abs(actual_levels - brass_logit_e0) == min(abs(actual_levels - brass_logit_e0)))[1]]
      }
    }
  }else{
    mlt_data <- mlt_data_input %>% dplyr::mutate(type = "user", e0 = "user")
    brass_logit_e0 <- "user"
  }

  this_mlt_family <- mlt_data %>%
    dplyr::group_by(type, e0) %>%
    dplyr::mutate(l25 = lx[age == 25],
                  p25_n_mlt = pmin(lx / l25, 1)) %>%
    dplyr::select(type, age, e0, lx_mlt = lx, l25_mlt = l25, p25_n_mlt) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(type, e0, age)

  if(brass_logit | !is.null(mlt_data_input)){
    mlt_closest <- p25_data %>%
      dplyr::left_join(this_mlt_family %>% dplyr::filter(e0 == brass_logit_e0),
                       by = c("age_25n" = "age")) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(type, e0, age)

    lx_out <- mlt_closest %>%
      dplyr::mutate(alpha = -.5 * log(1 + (p25_n / lx_mlt - 1 / l25_mlt) / (1 - p25_n))) %>%
      dplyr::select(-lx_mlt) %>%
      dplyr::left_join(this_mlt_family %>%
                         dplyr::filter(e0 == brass_logit_e0) %>%
                         dplyr::select(type, e0, age_mlt = age, lx_mlt),
                       by = c("type", "e0")) %>%
      dplyr::mutate(logit_mlt_lx = logit(1 - lx_mlt),
                    lx_interp = 1 - logit_inv(alpha + logit(1 - lx_mlt)))
  }else{
    obs_data <- p25_data %>% dplyr::select(age = age_25n, p25_n)
    mlt_level_data <- this_mlt_family %>% dplyr::select(type, age, e0, p25_n = p25_n_mlt)
    mlt_closest <- interp_level_mlt(obs_data, mlt_level_data, "e0", "p25_n")

    ages_filter_out <- mlt_closest %>%
      dplyr::filter(!dplyr::between(e0_interp, e0_accept[1], e0_accept[2])) %>%
      dplyr::pull(age)
    if(length(ages_filter_out) > 0 & verbose) {
      message(paste0("Age ", ages_filter_out, " filtered out because not acceptable mortality level."))
    }

    lx_out <- p25_data %>%
      dplyr::inner_join(mlt_closest %>%
                          dplyr::filter(!age %in% ages_filter_out) %>%
                          dplyr::select(type, age_25n = age, e0_interp),
                        by = "age_25n") %>%
      split(list(.$type, .$age), drop = TRUE) %>%
      purrr::map_df(function(X){
        obs_lx <- data.frame(age = c(0, 1, seq(5, 100, 5)), e0 = X$e0_interp)
        mlt_lx <- this_mlt_family %>%
          dplyr::select(type, age, e0, lx = lx_mlt) %>%
          dplyr::filter(type == unique(X$type))
        mlt_closest_lx <- interp_level_mlt(obs_lx, mlt_lx, "lx", "e0")
        return(data.frame(X, mlt_closest_lx %>% dplyr::select(age_mlt = age, lx_interp), row.names = NULL))
      })
  }

  lx_out <- lx_out[!is.na(lx_out$lx_interp), ]
  adult_mort <- lx_out %>%
    dplyr::ungroup() %>%
    dplyr::select(type, age, age_mlt, lx_interp) %>%
    dplyr::group_by_at(dplyr::vars(-c(age_mlt, lx_interp))) %>%
    dplyr::summarise(q15_45 = 1 - lx_interp[age_mlt == 60] / lx_interp[age_mlt == 15],
                     q15_35 = 1 - lx_interp[age_mlt == 50] / lx_interp[age_mlt == 15],
                     .groups = "keep") %>%
    dplyr::ungroup()

  return(list(adult_mort_index = adult_mort,
              lx_estimates = lx_out,
              p25_estimates = p25_data,
              synthetic_maternal_survival = orp_data,
              hiv_adjusted_proportions = adj_long))
}

# Toy/example using the South Africa workbook distributed with the
# Masquelier-Timaeus HIV/AIDS orphanhood materials. Proportions are from the
# 2007 and 2016 inquiries in that workbook; HIV/PMTCT/ART inputs are the
# workbook's already-interpolated values at the exact dates used below.
#
sa_age <- seq(0, 40, 5)
sa_2007 <- c(0.976752231834, 0.938522218426, 0.904016522616,
             0.879178529507, 0.862578576842, 0.833214037245,
             0.782326632973, 0.694392148056, 0.591888426780)
sa_2016 <- c(0.990523726582, 0.961273816529, 0.915876052936,
             0.865710780852, 0.827517688102, 0.802972857395,
             0.779306963601, 0.744475915327, 0.667449129474)
sa_mac <- c(NA, 27.093958904110, 27.326510958904, 27.331552638671,
            27.419768979714, 27.457958881653, 27.471149644434,
            27.701852062280, 27.838196781196)
sa_hiv <- data.frame(
  year = c(1964.620548, 1969.620548, 1973.773224, 1974.620548,
           1978.773224, 1979.620548, 1983.773224, 1984.620548,
           1988.773224, 1989.620548, 1993.773224, 1994.620548,
           1998.773224, 1999.620548, 2003.773224, 2004.620548,
           2007.120548, 2008.773224, 2011.696886, 2013.773224,
           2016.273224),
  value = c(0, 0, 0.0001, 0.000112055, 0.000427322, 0.000524110,
            0.001154645, 0.001312055, 0.002655191, 0.003937534,
            0.033433880, 0.044920274, 0.110798361, 0.123339452,
            0.167385246, 0.173147397, 0.188409452, 0.198893989,
            0.216464741, 0.224919672, 0.231201093)
)
sa_pmtct <- data.frame(
  year = c(1964.620548, 1969.620548, 1973.773224, 1974.620548,
           1978.773224, 1979.620548, 1983.773224, 1984.620548,
           1988.773224, 1989.620548, 1993.773224, 1994.620548,
           1998.773224, 1999.620548, 2003.773224, 2004.620548,
           2008.773224, 2013.773224),
  value = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.000373699, 0.101726776, 0.165495342,
            0.630328415, 0.995784699)
)
sa_art <- data.frame(
  year = c(2007.120548, 2011.696886, 2016.273224),
  value = c(0.139831644, 0.449831643, 0.576669126)
)
sa_out <- orphanhood_two_b(
  prop1_not_orphan = sa_2007,
  prop2_not_orphan = sa_2016,
  age1 = sa_age,
  age2 = sa_age,
  date1 = 2007.120548,
  date2 = 2016.273224,
  mac1 = sa_mac,
  mac2 = sa_mac,
  HIV_prev = sa_hiv,
  HIV_art = sa_art,
  HIV_pmtct = sa_pmtct,
  mlt_family = "CD_West"
)
round(sa_out$p25_estimates$p25_n, 6)
# Expected from the workbook, approximately:
# 0.927425 0.884010 0.843594 0.814588
# 0.773862 0.735641 0.711444 0.687329
