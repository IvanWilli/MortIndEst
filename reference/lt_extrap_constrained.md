# Extrapolate mortality by age, constrained to a life expectancy value.

Instead of fitting parameters from some mortality law in previous ages
to an observed OAG, find the parameters for the mortality function that
replicates e(OAG) (Ediev, 2016).

## Usage

``` r
lt_extrap_constrained(
  nMx,
  Age,
  Sex = "m",
  method = "classical",
  extrapLaw = "Gompertz",
  eOAG = NULL,
  OAnew = 100,
  alpha = NULL,
  beta = NULL,
  r = NULL,
  x_hat = NULL
)
```

## Arguments

- nMx:

  numeric. Vector of mortality rates in abridged or single age classes.

- Age:

  integer. Vector with lower bound for each age class (could be integer
  or abridged). Last age is assumed as lower age from open age group.

- Sex:

  character. Either male `"m"`, or female `"f"`.

- method:

  character. This indicates how to calculate `e_{OAG}` with methods:
  `classical` (default), `H-C` or `Mitra` (Ediev, 2014).

- extrapLaw:

  character. Available options: `Gompertz` or `Kannisto`.

- eOAG:

  numeric. An estimate of life expectancy in the input OAG age. If this
  has a value, then replace `method`.

- OAnew:

  integer. After extrapolating, pick a new OAG.

- alpha:

  numeric. Parameter for `H-C` method. Default values are automatically
  selected for input OAG.

- beta:

  numeric. Parameter for `H-C` method. Default values are automatically
  selected for input OAG.

- r:

  numeric. Parameter for `H-C` and `Mitra` methods. Annual growth rate
  of the population in the open age interval.

- x_hat:

  numeric. Parameter for `Mitra` method. Mean age of the population in
  the open age interval.

## Examples

``` r
if (FALSE) { # \dontrun{
# Pakistan - 1968-1971 - males - OAG=65 - HLD (2020)
nMx <- c(0.13328, 0.01539, 0.0031, 0.00155, 0.00169, 0.00185, 0.00201,
         0.0024, 0.00289, 0.00367, 0.00498, 0.00736, 0.01214, 0.02301, 0.0879)
Age <- c(0L, 1L, 5L, 10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L, 65L)
fit_extrap_constr <- lt_extrap_constrained(nMx, Age, OAnew = 100, extrapLaw = "Gompertz")
fit_extrap_constr$ex[fit_extrap_constr$Age == 65] - (1/nMx[Age == 65])
# if some other value is required on e(65), like 10. The function forces extrapolation to fit that value.
nMx <- c(0.13328, 0.01539, 0.0031, 0.00155, 0.00169, 0.00185, 0.00201,
         0.0024, 0.00289, 0.00367, 0.00498, 0.00736, 0.01214, 0.02301, 0.0879)
Age <- c(0L, 1L, 5L, 10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L, 65L)
fit_extrap_constr <- lt_extrap_constrained(nMx, Age, OAnew = 100, extrapLaw = "Gompertz", eOAG = 10)$lt
fit_extrap_constr$ex[fit_extrap_constr$Age == 65] - 10
} # }
```
