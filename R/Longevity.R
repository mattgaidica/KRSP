library(RMySQL)
library(dplyr)
library(lubridate)

con = dbConnect(MySQL(),
                dbname = "KRSP v20", 
                host = "localhost", 
                port = 3306,
                user = "root",
                password = "root")
litter = tbl(con, "litter") %>% 
  select(squirrel_id, litter_id = id, part = fieldBDate)

flastall = tbl(con, "flastall") %>% 
  select(squirrel_id, byear, bcert, litter_id, dates, datee) %>% 
  left_join(., litter) %>% 
  collect() %>% 
  mutate(longevity = as.integer(difftime(datee, dates, units = "days")))

write.csv(flastall,'/Users/matt/Downloads/longevity.csv')

#To calculate age add this line:
# mutate(your_data_set, age = if_else(bcert == "Y", year - byear, NA_real_))