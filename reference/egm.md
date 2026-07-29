# Extinct Generation Method (EGM) from Lexis-square deaths

Compute cohort population counts (the extinct generation method) by
accumulating deaths along cohort lines in a Lexis diagram, locating
estimates at at the beginning of each year. Inputs are calendar year
`y`, age at the Lexis square `x`, and deaths `d` in each square. Cohort
(year of birth, `yb`) can be supplied directly, or derived by splitting
each Lexis square into upper/lower triangles using `alpha`.

## Usage

``` r
egm(y, x, d, yb = NULL, alpha = NULL)
```

## Arguments

- y:

  Integer or numeric vector. Calendar year (Lexis square).

- x:

  Integer or numeric vector. Age (Lexis square).

- d:

  Numeric vector. Deaths in the Lexis square \\(y, x)\\.

- yb:

  Integer or numeric vector (optional). Year of birth (cohort). If
  supplied, it must be the same length as `y`/`x`/`d` and will be used
  directly (no triangle split). If `NULL`, `yb` is derived from `y` and
  `x` using `alpha`.

- alpha:

  Numeric scalar in \\\[0, 1\]\\ (optional). Upper triangle share used
  to split Lexis squares when `yb` is `NULL`. Default is `0.5`. Ignored
  when `yb` is provided.

## Value

A data frame with one row per cohort–age–year of population stock at the
beggining of the year:

- yb:

  Year of birth (cohort).

- x:

  Age (exact age on vertical segment).

- y:

  Calendar year corresponding to age \\x\\.

- AgeDeath105:

  Inform if max age with data is 105 or more. A flag for extinct
  interpretation.

- pop:

  Cohort population at exact age \\x\\, reconstructed by reverse
  cumulative deaths.

## Details

The extinct generation method reconstructs cohort population size at age
\\x\\ as the (reverse) cumulative sum of cohort deaths from age \\x\\
onward. Implementation note: this uses `dplyr::mutate(..., .by = yb)`,
which requires **dplyr \>= 1.1.0**.

## Assumptions

- Death registration is complete for extinct (or nearly extinct)
  cohorts.

- When `yb` is not provided, the within-square timing of deaths is
  approximated by a constant split `alpha` between triangles.
