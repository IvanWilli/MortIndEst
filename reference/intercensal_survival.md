# Estimate adult mortality using cohort relationship between two censuses.

This function implements a variety of methods. See `method` argument for
more details. An option for considering HIV populations is using the
Spectrum model (Stover and others, 2012).

## Usage

``` r
intercensal_survival(
  c1,
  c2,
  date1,
  date2,
  age1 = seq(0, length(c1) * 5, 5),
  age2 = seq(0, length(c2) * 5, 5),
  sex = "f",
  mlt_family = "CD_West",
  method = "match",
  ages_fit = seq(10, 70, 5),
  e0_accept = c(20, 100),
  q01_q05 = NULL,
  mlt_e0_logit_feeney = NULL,
  span_pre_smooth = NULL,
  mlt_input_data = NULL,
  HIV_prev = NULL,
  HIV_art = NULL,
  extrapLaw = "makeham",
  verbose = TRUE
)
```

## Arguments

- c1:

  numeric vector. Population counts in five-age groups, from first
  census with an exact reference date.

- c2:

  numeric vector. Population counts in five-age groups, from second
  census with an exact reference date.

- date1:

  Either a `Date` class object or an unambiguous character string in the
  format `"YYYY-MM-DD"`. Reference date for first census.

- date2:

  Either a `Date` class object or an unambiguous character string in the
  format `"YYYY-MM-DD"`. Reference date for second census.

- sex:

  character string. Either `"m"` for males, `"f"` for females (default).

- mlt_family:

  character string. Model life table family: `"Chilean"`, `"East"`,
  `"Far_East_Asian"`, `"General"`, `"Latin"`, `"North"`, `"South"`,
  `"South_Asian"` or `"West"` (default).

- method:

  character. Options include:

  - `"match"` use survival ratios for finding best mlt level (default).

  - `"bproj"` back-projection from second census (UN, 1983).

  - `"fproj"` forward-projection from first census (UN, 1983).

  - `"var-r"` variable-r method (Preston & Benet, 1986).

  - `"feeney"` smoothing survival ratios into survival function
    implementation by Feeney (UN, 2002).

  - `logit` based on Preston (2001, section 11.5.2)

- ages_fit:

  integer vector. Ages to be considered when calculating the median of
  implied level by age. By default 10 to 70, but depends data and
  method.

- e0_accept:

  integer vector. Range acceptable when calculating the median of
  implied level by age (avoid non-possible extrapolations). By deafult
  between 20 and 100.

- q01_q05:

  numeric. Both values in case final life table be computed considering
  child mortality input values.

- mlt_e0_logit_feeney:

  numeric. Level as reference in `logit` and `feeney` methods.

- span_pre_smooth:

  numeric. smooth is applied to survival ratios by age in case is not
  null. Value is `span` parameter from `loess` function (between 0 and
  1). Default `FALSE`.

- mlt_input_data:

  data.frame. By default is used data from `Morcast` package. But
  specific model can be included (be sure same columns are included with
  same names), ant that `mlt_family` be included as category in `type`
  column.

- HIV_prev:

  numeric. Estimates of population-based estimate of HIV prevalence by
  age. If some value is assigned, then will be assumed AIDS variation of
  the method based on Spectrum patterns.

- HIV_art:

  numeric. Estimate of proportion using antiretroviral therapy (ART).

- extrapLaw:

  character. Which mortality law is chosen for extension to 100 age in
  lt output. See `Demotools::lt_abridged` for mode details

- age:

  integer vector. Lower bound of age groups from second census. Last age
  is assumed OAG.

## Examples

``` r
if (FALSE) { # \dontrun{
# match: PANAMA 1960-1970 from UNX (Table 174)
panama_match <- intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "match")

# UNX removes outliers from ages 10, 20, 35, 55 and 65 (implying the use of the initial OAG)
# and get the mean 16.1 which is between 57.5 and 60 of e0.
# by default the function only removes age 10 because nSx and give a result a higher one.
# But is possible to reconstruct the result using the output where is found the best level for each age:
panama_match$rank_match %>% mutate(color = ifelse(age %in% c(10, 20, 35, 55, 65), "outlier", "used")) %>%
panama_match$rank_match %>% filter(!age %in% c(10, 20, 35, 55, 65)) %>% summarise(mean(e0_interp))
# to replicate this during the process, parameters ages_fit or e0_range_accept should be used
panama_match_modif_ages <- intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "match",
  ages_fit = seq(0,75,5)[!seq(0,75,5) %in% c(10, 20, 35, 55, 65)])$selected_level

# projection. Indonesia Males 1980-1990, SA pattern (Preston, 2001; Box 11.4) - OK
india_fproj <- intercensal_survival(
  c1 = c(10815974, 10832383, 9131871, 7512541, 5978576, 5612684, 4022625, 4190944, 3644053, 3012756, 2717883, 1720501, 1559230, 811113, 689074, 688422),
  c2 = c(10760859,11928095,11044127,9520440,7583305,7457150,6584325,5788441,4010254,3723922,3289190,2321621,2219069,1329162,945876,867636),
  date1 = "1980/07/01",
  date2 = "1990/07/01",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "m",
  mlt_family = "UN_South_Asian",
  method = "fproj",
  e0_accept = c(51,57),
  ages_fit = seq(0,40,5))

# Feeney (2002, ESA/P/WP.175): japan females 1960-1970. OK - level at age 5, 68.7 (page 27)
# level 72.5 has ~ e80=6.06. It´s not exact but very accurate
japan_feeney <- intercensal_survival(
  c1 <- c(3831870,4502304,5397061,4630775,4193184,4114704,3770907,3274822,
          2744786,2559755,2160716,1839025,1494043,1133409,870238,577972,313781, 131547),
  c2 <- c(4292503,3988292,3852101,4492096,5347327,4571868,4190340,4085338,
          3674127,3198934,2648360,2382691,1970485,1584699,1172155,736258,408191, 53116),
  date1 = "1960/08/18",
  date2 = "1970/08/18",
  age1 = seq(0,85,5),
  age2 = seq(0,85,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "feeney",
  ages_fit = 5,
  mlt_e0_logit_feeney = 72.5)

# Preston (1983) logit - r. Indian Females 1961-1971 - OK e5 = 51
india_preston <- intercensal_surv_Preston_logit_var_r(
        c1 = c(32889570, 31569832, 23003038, 17262494, 19113102, 18031562,
               14837941, 11846139, 10758321, 8310282, 7966324, 4540060, 5521764,
               2374069, 4434064),
       c2 = c(39050792, 40165344, 32193430, 22241264, 21524600, 20477162,
             17892874, 15632259, 13264020, 10381364, 9489160, 5875138, 6893570,
             3268044, 5664190),
       date1 = 1961.164,
       date2 = 1971.249,
       age1 = seq(0,70,5),
       age2 = seq(0,70,5),
       sex = "f",
       mlt_family = "CD_South",
       mlt_e0 = 50,
       ages_fit = seq(5,65,5),
       q0_5 = 1 - .776)

# variable-r (Preston-Bennet): PANAMA 1960-1970 from UNX (Table 186). - OK
# It´s picked ages_fit 10:30 as Manual choose, because growing too much below that
# promedio e(10) 17.1 ~ 60.3, en comparación con 60.1
panama_var_r <- intercensal_survival(
  c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
  c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
  date1 = "1960/12/11",
  date2 = "1970/05/10",
  age1 = seq(0,75,5),
  age2 = seq(0,75,5),
  sex = "f",
  mlt_family = "CD_West",
  method = "var-r",
  ages_fit = 10:30)
} # }
```
