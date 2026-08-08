#' Orphanhood database
#'
#' Census and survey records of parental survival by respondent age, sex,
#' location, source, period, and related metadata.
#'
#' @format A `data.table` and `data.frame` with 59,937 rows and 93 columns.
#' Important fields include `LocName`, `LocID`, `DataCatalogName`,
#' `ReferenceYearStart`, `ReferenceYearEnd`, `SexName`, `AgeStart`, `AgeEnd`,
#' `PeriodStart`, `PeriodEnd`, and `DataValue`.
#' @source United Nations structured demographic data extract.
"Orphanhood_DB"

#' Census population database
#'
#' Census population records by location, sex, age, source, period, and related
#' metadata.
#'
#' @format A `data.frame` with 59,524 rows and 93 columns.
#' Important fields include `LocName`, `LocID`, `DataCatalogName`,
#' `ReferenceYearStart`, `ReferenceYearEnd`, `SexName`, `AgeStart`, `AgeEnd`,
#' `PeriodStart`, `PeriodEnd`, and `DataValue`.
#' @source United Nations structured demographic data extract.
"all_censuses"

#' WPP 2024 five-year life tables for Brazil
#'
#' Brazil Country-level five-year life table values from the 2024 World Population
#' Prospects extract.
#'
#' @format A `data.frame` with 29 columns:
#' \describe{
#'   \item{SortOrder, LocID, VarID, Time, SexID, AgeGrpStart, AgeGrpSpan}{Integer identifiers and time or age fields.}
#'   \item{Notes, ISO3_code, ISO2_code, LocTypeName, Location, Variant, Sex, AgeGrp}{Character metadata fields.}
#'   \item{MidPeriod, mx, qx, px, lx, dx, Lx, Sx, Tx, ex, ax}{Numeric life table values.}
#' }
#' @source United Nations World Population Prospects 2024.
"lt5_wpp24_brazil"

#' WPP 2022 mean age of childbearing
#'
#' Mean age of childbearing by location, year, and parent sex.
#'
#' @format A tibble with 37,488 rows and 4 columns:
#' \describe{
#'   \item{LocID}{Integer location identifier.}
#'   \item{Year}{Integer calendar year.}
#'   \item{parent}{Parent sex category.}
#'   \item{mac}{Mean age of childbearing.}
#' }
#' @source United Nations World Population Prospects 2022.
"mac_wpp22"

#' Spectrum mortality models
#'
#' Fitted Spectrum mortality model components for female and male mortality
#' estimation.
#'
#' @format A named list of length 2 with elements `female` and `male`.
#' @source Spectrum mortality model files.
"mods.r"

#' Argentina 2001 census population
#'
#' Argentina 2001 census population counts by year and single year of age.
#'
#' @format A `data.frame` with 303 rows and 3 columns:
#' \describe{
#'   \item{year}{Integer census year.}
#'   \item{age}{Integer age in completed years.}
#'   \item{pop}{Numeric population count for both sexes combined.}
#' }
#' @source Instituto Nacional de Estadistica y Censos (INDEC), Argentina.
"arg_pop"

#' Argentina deaths
#'
#' Argentina priod death counts by year and single year of age, years 2001-2023.
#'
#' @format A `data.frame` with 2,553 rows and 3 columns:
#' \describe{
#'   \item{year}{Integer calendar year.}
#'   \item{age}{Integer age in completed years.}
#'   \item{deaths}{Numeric death count for both sexes combined.}
#' }
#' @source Direccion de Estadisticas e Informacion de la Salud (DEIS),
#' Ministerio de Salud, Argentina.
"arg_deaths"
