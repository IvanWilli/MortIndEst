# Estimate life expectancy at age 5 (`e(5)`) using cohort relationship between two censuses, based on r-variable methods.

Method developed by Preston (2001, section 11.5.2)

## Usage

``` r
intercensal_surv_Preston_logit_var_r(
  c1,
  c2,
  date1,
  date2,
  age1 = 1:length(c1) - 1,
  age2 = 1:length(c2) - 1,
  sex = "f",
  mlt_family = "CD_West",
  mlt_e0 = 60,
  ages_fit = NULL,
  q0_5 = NULL,
  verbose = FALSE
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

- age1:

  integer vector. Lower bound of age groups from first census. Last age
  is assumed OAG.

- age2:

  integer vector. Lower bound of age groups from second census. Last age
  is assumed OAG.

- sex:

  character string. Either `"m"` for males, `"f"` for females (default).

- mlt_family:

  character string. Model life table family: `"Chilean"`, `"East"`,
  `"Far_East_Asian"`, `"General"`, `"Latin"`, `"North"`, `"South"`,
  `"South_Asian"` or `"West"` (default).

- ages_fit:

  integer vector. Ages to be considered when calculating the median of
  implied level by age. By default 10 to 70, but depends data and
  method.

- q05:

  numeric. Probability of dying from birth to age 5.
