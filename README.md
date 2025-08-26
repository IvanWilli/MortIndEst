
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MortIndEst

This package includes several functions for indirect estimates,
specially about adult mortality.

Main development was done during World Population Prospects revision
2024, with the goal to include as inputs in the B3 model (United
Nations: Department of Economic and Social Affairs 2024).

Functions included are related to intercensal survival estimates,
orphanhood, siblings survival, widowhood survival, census method for
${}_{15}q_{60}$ and old ages extension of mortality rates. Some of the
contents consists in R implementations of materials in Moultrie et al.
(2012) and United Nations (1983).

The idea of the package is to continue expanding its content, so
suggestions are welcome. For other indirect methods not included here
you can see [fertsestr](https://github.com/josehcms/fertestr), …

## Installation

You can install the development version of MortIndEst from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("IvanWilli/MortIndEst")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MortIndEst)
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-iussp2012" class="csl-entry">

Moultrie, Tom A., Rob E. Dorrington, Allan G. Hill, Kenneth Hill, Ian M.
Timæus, and Basia Zaba, eds. 2012. *Tools for Demographic Estimation*.
Paris: International Union for the Scientific Study of Population
(IUSSP). <https://demographicestimation.iussp.org/>.

</div>

<div id="ref-un_manualX" class="csl-entry">

United Nations. 1983. *Manual x: Indirect Techniques for Demographic
Estimation*. Population Studies 81. New York: United Nations.
<https://www.un.org/development/desa/pd/sites/www.un.org.development.desa.pd/files/files/documents/2020/Jan/un_1983_manual_x.pdf>.

</div>

<div id="ref-wpp2024" class="csl-entry">

United Nations: Department of Economic and Social Affairs. 2024. “World
Population Prospects 2024.” United Nations.
<https://population.un.org/wpp/>.

</div>

</div>
