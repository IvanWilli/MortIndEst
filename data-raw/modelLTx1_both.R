# Run from the package root to regenerate the model life tables.
library(dplyr)

mlt <- DemoToolsData::modelLTx1
srb <- 1.04

mlt_both <- mlt %>%
  summarise(ax1 = weighted.mean(ax1, dx1 * ifelse(sex == "male", srb, 1)),
            sx1 = weighted.mean(sx1, Lx1 * ifelse(sex == "male", srb, 1)),
            across(c(lx1, dx1, Lx1, Tx1),
                   ~ weighted.mean(.x, ifelse(sex == "male", srb, 1))),
            .by = c(source, family, type_mlt, e0, age)) %>%
  mutate(sex = "both", mx1 = dx1/Lx1, qx1 = dx1/lx1, ex1 = Tx1/lx1)

modelLTx1_both <- bind_rows(mlt, mlt_both) %>%
  arrange(source, family, sex, e0, age) %>%
  mutate(lx = lx1/lx1[age == 0], .by = c(source, family, sex, e0))

save(modelLTx1_both, file = "data/modelLTx1_both.rda", compress = "xz", version = 2)
