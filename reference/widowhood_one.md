# Estimate adult mortality with one survey/census from widowhood data.

Follows Hill and Trusell (1977) revised regression coefficients,
included in Manual X from UN. The function also include the possibility
of using Brass logit given a pattern as default method. An option for
considering HIV populations is using the Spectrum model (Stover and
others, 2012).

## Usage

``` r
widowhood_one(
  prop_not_widowed,
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
  verbose = TRUE
)
```

## Arguments

- prop_not_widowed:

  numeric vector. Population proportion of ever-married (first spouse
  assumption), by age groups.

- age:

  integer vector. Lower bound of age groups. Last age is assumed as
  lower age open age group.

- sex_spouse:

  character. "f" for male respondent (estimate female mortality), and
  viceversa.

- smam_f:

  numeric. Singulate mean age at marriage for females.

- smam_m:

  numeric. Singulate mean age at marriage for males.

- date:

  Either a `Date` class object or an unambiguous character string in the
  format `"YYYY-MM-DD"`. Reference date for the source (typically census
  or survey).

- mlt_data_input:

  data.frame. If is `NULL` then model life tables is used, available
  from `Morcast` package. But specific pattern can be included (`l_x`
  must be included see examples).

- mlt_family:

  character. Options: "CD_East", "CD_North", "CD_South", "CD_West",
  "UN_Chilean", "UN_Far_Eastern", "UN_General", "UN_Latin_American",
  "UN_South_Asian". If `NULL` returns for all.

- brass_logit:

  logical. Doing a level smoothing with 1-parameter logit Brass.

- brass_logit_e0:

  numeric. Life expectancy at birth when `brass_logit` is `TRUE`.

- brass_logit_5q0:

  numeric. Find best level when `brass_logit_e0` is `NULL`. If is not
  `NULL`, then implied level replaces `brass_logit_e0`.

- HIV_prev:

  numeric. Estimates of population-based estimate of HIV prevalence
  among men by age. If some value is assigned, then will be assumed AIDS
  variation of the method.

- HIV_art:

  numeric. Estimate of the proportion of men with infected female
  partner.

## Examples

``` r
if (FALSE) { # \dontrun{
# Bolivia 1975. Example from UN Manual X.
manualX_bolivia_male <- widowhood_one(prop_not_widowed = c(0.9798, 0.9729, 0.9514, 0.9170, 0.8735, 0.8195, 0.7054, 0.6520),
                                      age = seq(20, 55, 5),
                                      sex_spouse = "m",
                                      smam_f = 23.2 ,
                                      smam_m = 25.3 ,
                                      date = 1975.6,
                                      mlt_family = "CD_West",
                                      brass_logit = FALSE)
un_results_py_n <- c(.956, .931, .894, .85, .798, .687, .638)
un_results_level <- c(16.3, 16.6, 15.7, 15.4, 15.3, 13.3, 14.7)
round(unique(manualX_bolivia_male$lx_estimates$py_n) - un_results_py_n,3)
un_results_e0 <- approx(13:17, c(47.082, 49.546, 51.816, 54.122, 56.45), un_results_level)$y
round(unique(manualX_bolivia_male$lx_estimates$e0_interp) - un_results_e0,1) # ok with MORTPACK TOO
} # }
```
