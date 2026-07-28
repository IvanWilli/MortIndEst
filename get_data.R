# get data

# census
library(tidyverse)
load("D:/UN_WPP/WPP_indirect_methods/intercensal_survival/all_censuses.rda")

load("D:/UN_WPP/WPP_indirect_methods/orphanhood/Orphanhood_DB_full.rda")
guatemala_orphanhood <- Orphanhood_DB %>%
  filter(LocName == "Guatemala")
load("D:/UN_WPP/WPP_indirect_methods/widowhood/Widowhood_DB_full.rda")
Widowhood_DB %>%
  filter(LocName == "Guatemala")



census_data <- all_censuses %>%
  filter(LocName == "Guatemala")
guatemala_census %>%
  distinct(TimeMid)
load("D:/UN_WPP/WPP_indirect_methods/intercensal_survival/wpp22_e0.rda")
load("D:/UN_WPP/WPP_indirect_methods/intercensal_survival/wpp22_lx.rda")


census_variant <- expand_grid(date1 = unique(census_data$TimeMid),
            date2 = unique(census_data$TimeMid),
            sex = unique(census_data$SexName),
            method = c("match", "bproj", "fproj", "var-r", "feeney", "logit")) %>%
  filter(date2 > date1) %>%
  mutate(distance = abs(date2 - date1)) %>%
  group_by(date1, sex, method) %>% arrange(distance) %>% slice(1) %>%
  mutate(census_pair = paste0(round(date1,1),"-", round(date2,1))) %>%
  ungroup() %>%
  mutate(t_middle = trunc((date1 + date2)/2)) %>%
  filter(t_middle>=1950) %>%
  left_join(data.frame(li =     c(15,      0,       0,       10,       5,        0),
                       ui =     c(60,      40,      40,      30,       30,       60),
                       method = c("match", "bproj", "fproj", "var-r", "feeney", "logit")))


map_df(1:nrow(census_variant),
               function(Y){
                 Y <- 2
                 X <- census_variant[Y,]
                 # X <- census_variant %>% filter(sex == "Female", method == "fproj") %>% slice(5)
                 # print(paste(X$LocID, X$census_pair, X$sex, X$method, collapse = " - "))
                 census_selection_sex <- guatemala_census %>% filter(SexName == X$sex, TimeMid %in% c(X$date1, X$date2)) %>% arrange(TimeMid, AgeStart)
                 e0 <- wpp22_e0 %>% filter(LocationId == 320, Sex == X$sex, trunc(TimeMid) == trunc((X$date1+X$date2)/2)) %>% pull(Value)
                 q0_1 <- 1 - wpp22_lx %>% filter(AgeStart ==1, LocationId == 320, Sex == X$sex, trunc(TimeMid) == trunc(X$t_middle)) %>% pull(Value)/1e5
                 q0_5 <- 1 - wpp22_lx %>% filter(AgeStart ==5, LocationId == 320, Sex == X$sex, trunc(TimeMid) == trunc(X$t_middle)) %>% pull(Value)/1e5
                   out <- intercensal_survival(c1 = census_selection_sex$DataValue[census_selection_sex$TimeMid == X$date1],
                                               c2 = census_selection_sex$DataValue[census_selection_sex$TimeMid == X$date2],
                                               date1 = X$date1,
                                               date2 = X$date2,
                                               age1 = census_selection_sex$AgeStart[census_selection_sex$TimeMid == X$date1],
                                               age2 = census_selection_sex$AgeStart[census_selection_sex$TimeMid == X$date2],
                                               sex = ifelse(X$sex == "Female", "f", "m"),
                                               mlt_family = NULL,
                                               method = X$method,
                                               ages_fit = seq(X$li, X$ui, 5),
                                               mlt_e0_logit_feeney = e0,
                                               q01_q05 = c(q0_1, q0_5))
                  out$lt_fit %>%
                    inner_join(out$surv_data %>%
                                 select(Age = age, nSx, n), by = "Age")
                })

