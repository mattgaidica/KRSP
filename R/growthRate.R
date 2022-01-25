library(RMySQL)
library(dplyr)
library(lubridate)
con = dbConnect(MySQL(),
                dbname = "KRSP v20", 
                host = "localhost", 
                port = 3306,
                user = "root",
                password = "root")

juvenile <- tbl(con, "juvenile") %>%
  select(squirrel_id, litter_id, sex, weight, tagwt = tagWT)
# litter level information
litter <- tbl(con, "litter") %>%
  select(litter_id = id, year = yr, field_bdate = fieldBDate, n1_date = date1, tag_date = tagDt,
         mother_id = squirrel_id, grid, locx, locy, ln, food)

results = inner_join(juvenile, litter, by = "litter_id") %>% 
  collect()

# clean up variables
results <- results %>%
  # convert dates from character to date class
  mutate(field_bdate = lubridate::ymd(field_bdate),
         n1_date = lubridate::ymd(n1_date),
         tag_date = lubridate::ymd(tag_date))

growthTable <- results %>% 
  group_by(litter_id) %>% 
  mutate(litter_size = n()) %>%
  ungroup() %>% 
  mutate(
    # growth rate
    nest_days = as.numeric(difftime(tag_date, n1_date, units = "days")),
    growth = (tagwt - weight) / nest_days,
    # calculate birth date as julian day
    bdate = lubridate::yday(field_bdate),
    growth = case_when(
      (is.na(weight) | !between(weight, 1, 50)) ~ NA_real_,
      (is.na(weight) | !between(weight, 1, 100)) ~ NA_real_,
      nest_days < 5 ~NA_real_,
      (tagwt-weight)==0 ~ NA_real_,
      TRUE ~ growth
    )) %>% 
  select(-nest_days) %>% 
  select(squirrel_id, growth)

write.csv(growthTable,'/Users/matt/Documents/MATLAB/KRSP/R/krsp_growth.csv')
