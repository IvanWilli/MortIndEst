# get data

# census
library(tidyverse)
load("D:/UN_WPP/WPP_indirect_methods/intercensal_survival/all_censuses.rda")
all_censuses %>%
  filter(LocName == "Guatemala") %>%
  distinct(TimeMid)
