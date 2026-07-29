# Estimate child mortality with the Brass indirect method

Estimate probabilities of dying in childhood from children ever born and
children dead reported by women by age group. The function applies
either the Trussell or Palloni-Helligman multipliers, locates each
estimate in time, and interpolates the matching model life table level.

## Usage

``` r
im_brass(
  CEB,
  CD,
  W,
  age = seq(15, 45, 5),
  variant = "trussell",
  mlt_model = NULL,
  Sex = NULL,
  date_obs = 0,
  mam = NULL
)
```

## Arguments

- CEB:

  Numeric vector. Children ever born by women's age group.

- CD:

  Numeric vector. Children dead by women's age group.

- W:

  Numeric vector. Women by age group.

- age:

  Numeric vector. Lower bounds of women's five-year age groups. Default
  is `seq(15, 45, 5)`.

- variant:

  Character. Multiplier variant to use. Options are `"trussell"` and
  `"palloni_helligman"`.

- mlt_model:

  Character. Model life table family. For `"trussell"`, use one of
  `"North"`, `"South"`, `"East"`, or `"West"`. For
  `"palloni_helligman"`, use one of `"Latin"`, `"Chilean"`,
  `"South_Asian"`, `"Far_East_Asian"`, or `"General"`. If `NULL`,
  `"West"` is assumed.

- Sex:

  Character. Sex of the child mortality estimates in the model life
  table, usually `"male"` or `"female"`. If `NULL`, `"male"` is assumed.

- date_obs:

  Numeric. Reference date of the census or survey. Used to compute
  `date_t`, the time location of each estimate.

- mam:

  Numeric. Mean age at maternity. Used only with
  `variant = "palloni_helligman"`. If `NULL`, 27 is assumed.

## Value

A list with two data frames:

- mi_brass:

  Brass estimates by women's age group, including the input age,
  estimated survivorship `lx`, childhood mortality age `x`, time lag
  `t`, and reference date `date_t`.

- mlts_interp:

  Interpolated model life table rows matched to each Brass estimate,
  including `mx_mlt`, `qx_mlt`, `lx_mlt`, and `ex_mlt`.

## References

United Nations. 1983. *Manual X: Indirect Techniques for Demographic
Estimation*. United Nations, New York.

## Examples

``` r
if (FALSE) { # \dontrun{
# Panama example from UN Manual X (1983)
panama <- data.frame(
  age = seq(15, 45, 5),
  W = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
  CEB = c(278, 1380, 2395, 3097, 3444, 3274, 2682),
  CD = c(24, 77, 172, 236, 348, 394, 354)
)
im_brass(
  CEB = panama$CEB,
  CD = panama$CD,
  W = panama$W,
  age = panama$age,
  mlt_model = "West",
  Sex = "male",
  date_obs = 1970
)
} # }
```
