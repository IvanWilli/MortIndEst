# Extrapolate mortality rates coherently at older ages

Based on the coherent Kannisto idea used in the MortCast package,
extended here to allow a Gompertz law and an optional period trend.

## Usage

``` r
coherent_mx(
  mx,
  law = "kannisto",
  fit.ages = NULL,
  estim.ages = NULL,
  method = "coherent",
  age_conv,
  fit.years = NULL,
  estim.years = NULL
)
```

## Arguments

- mx:

  Data frame with mortality rates. It must contain `age`, `nMx`, and
  `sex`. It may also contain `year`. The `age` column gives age lower
  bounds, `nMx` gives mortality rates, `sex` identifies sex groups, and
  `year` identifies the period when a time trend is fitted.

- law:

  Character. Mortality law used for the transformed rates. Options are
  `"kannisto"` and `"gompertz"`.

- fit.ages:

  Integer or numeric vector. Ages used to fit the regression. If `NULL`,
  all ages in `mx` are used.

- estim.ages:

  Integer or numeric vector. Ages to return after extrapolation. If
  `NULL`, all ages in `mx` are returned.

- method:

  Character. Extrapolation method when only one year is supplied. Use
  `"coherent"` for a common age slope by sex, or `"convergent"` to force
  the second sex group to converge to the first at `age_conv`.

- age_conv:

  Numeric. Age where the second sex group converges to the first when
  `method = "convergent"`.

- fit.years:

  Integer or numeric vector. Years used to fit the regression when `mx`
  contains more than one year. If `NULL`, all years are used.

- estim.years:

  Integer or numeric vector. Years to return when `mx` contains more
  than one year. If `NULL`, all years are returned.

## Value

A data frame with observed rates below the first estimated age and
fitted or extrapolated rates for `estim.ages`. Columns include `age`,
`nMx`, `sex`, and `year` when a year dimension is supplied.

## Details

Fit a simple old-age mortality model to observed rates and predict rates
for selected ages, sexes, and optionally years. Rates are transformed
with either a Gompertz log transform or a Kannisto logit transform
before fitting a linear model.

## Examples

``` r
mx <- data.frame(
  age = rep(seq(60, 80, 5), 2),
  sex = rep(c("f", "m"), each = 5),
  nMx = c(0.020, 0.035, 0.060, 0.100, 0.160,
          0.030, 0.050, 0.080, 0.130, 0.200)
)

coherent_mx(
  mx,
  law = "gompertz",
  fit.ages = seq(60, 80, 5),
  estim.ages = seq(60, 100, 5)
)
#>    age        nMx sex
#> 1   60 0.02152700   f
#> 2   65 0.03541838   f
#> 3   70 0.05827387   f
#> 4   75 0.09587801   f
#> 5   80 0.15774810   f
#> 6   85 0.25954297   f
#> 7   90 0.42702609   f
#> 8   95 0.70258609   f
#> 9  100 1.15596499   f
#> 10  60 0.02926447   m
#> 11  65 0.04814883   m
#> 12  70 0.07921929   m
#> 13  75 0.13033950   m
#> 14  80 0.21444760   m
#> 15  85 0.35283067   m
#> 16  90 0.58051236   m
#> 17  95 0.95511707   m
#> 18 100 1.57145427   m
```
