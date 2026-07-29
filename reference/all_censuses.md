# Census population database

Census population records by location, sex, age, source, period, and
related metadata.

## Usage

``` r
all_censuses
```

## Format

A `data.table` and `data.frame` with 59,524 rows and 93 columns.
Important fields include `LocName`, `LocID`, `DataCatalogName`,
`ReferenceYearStart`, `ReferenceYearEnd`, `SexName`, `AgeStart`,
`AgeEnd`, `PeriodStart`, `PeriodEnd`, and `DataValue`.

## Source

United Nations structured demographic data extract.
