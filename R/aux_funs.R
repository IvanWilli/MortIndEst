# aux funs
logit <- function(q){
  log(q/(1-q))/2
}
logit_inv <- function(logit_q){
  exp(logit_q * 2) / (1 + exp(logit_q * 2))
}
expit <- function(x) {
  return(exp(x)/(1 + exp(x)))
}


# level interpolation for any lt function ------------------------------------

interp_level_mlt <- function(obs_data, mlt_data, var_interp = "e0", var_ref = "nSx"){

  # match var names
  if(!all(colnames(obs_data) %in% colnames(mlt_data))) stop("names in obs_data are not in mlt_data")
  match.arg(var_interp, colnames(mlt_data))

  # rename
  var_reference_index_obs <- which(colnames(obs_data)==var_ref)
  colnames(obs_data)[var_reference_index_obs] <- "ref"
  var_reference_index_mlt <- which(colnames(mlt_data)==var_ref)
  colnames(mlt_data)[var_reference_index_mlt] <- "var_ref"
  var_interp_index_mlt <- which(colnames(mlt_data)==var_interp)
  colnames(mlt_data)[var_interp_index_mlt] <- "var_int"

  # join and find closest
  data <- dplyr::inner_join(obs_data, mlt_data,
                            by = if(!"type" %in% colnames(obs_data)) "age" else {c("type", "age")}) %>%
    dplyr::mutate(diff = abs(var_ref/ref-1)) %>%
    dplyr::arrange(type, age, var_int) %>%
    dplyr::group_by(type, age) %>%
    dplyr::arrange(diff) %>%
    dplyr::slice(1:2) %>%
    dplyr::mutate(levels = c("left", "right")) %>%
    tidyr::pivot_wider(id_cols = -diff, names_from = levels, values_from = c(var_int, var_ref)) %>%
    dplyr::mutate(diff = min(abs(ref - var_ref_left), abs(ref - var_ref_right)),
                  interp = var_int_left + (ref - var_ref_left)/(var_ref_right - var_ref_left) * (var_int_right - var_int_left)) %>%
    dplyr::select(type, age, ref, var_ref_left, var_ref_right, var_int_left, var_int_right, interp) %>%
    ungroup()

  # rename
  colnames(data) <- c("type", "age", paste0(var_ref,"_obs"),
                      paste0(var_ref,"_left"), paste0(var_ref,"_right"),
                      paste0(var_interp,"_left"), paste0(var_interp,"_right"),
                      paste0(var_interp,"_interp"))
  return(data)
}

# interpolate Z
interp_Z <- function(x){
  # table 88
  brass_time_coeff <- data.frame(age = 26:75,
                                 Z = c(rep(0.09, 9), 0.091, 0.092, 0.093, 0.095, 0.099, 0.104, .109,
                                       .115, .122, .13, .139, .149, .160, .171, .182, .193, .205, .218, .231,
                                       .245, .259, .274, .289, .305, .321, .338, .356, .374, .392, .411, .431,
                                       .452, .473, .495, .518, .542, .568, .595, .622, .650, .678))
  stats::approx(brass_time_coeff$age,
                brass_time_coeff$Z,
                x)$y
}

# get mlt e0
get_cd_females_level_e0 <- function(level = NULL, e0 = NULL){
  out <- if(!is.null(e0)) approx(seq(20,100,2.5), 1:25, e0)
  out <- if(!is.null(level)) 20 + (level - 1) * 2.5
  out
}
