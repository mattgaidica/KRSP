# cones -------------------------------------------------------------------
library(RMySQL)
library(dplyr)
library(lubridate)

con = dbConnect(MySQL(),
                dbname = "KRSP v20", 
                host = "localhost", 
                port = 3306,
                user = "root",
                password = "root")

###############
# Import Data #
###############
historic_ids<- tbl(con, "historic_squirrel_ids") %>% 
  select (old_squirrel_id, new_squirrel_id) %>% 
  collect()
# Note that supplementary tables are not updated in the annual data cleanup so squirrel_id values must be updated from 
# the historic_squirrel_ids table

flastall <- tbl(con, "flastall2") %>% 
  #flastall2 contains juveniles that were not tagged
  # exclusions
  filter(gr %in% c("SU", "KL", "CH", "AG", "LL", "JO", "RR", "SX")) %>% 
  select(squirrel_id, sex, byear) %>% 
  collect()


####################
# Midden Cone Data #
####################
sup = dbConnect(MySQL(),
                dbname = "krsp_suppLite", 
                host = "localhost", 
                port = 3306,
                user = "root",
                password = "root")
midden_cones<-tbl(sup, "midden_cones") %>% 
  filter(squirrel_id !="UTS") %>% #remove UTS squirrels from data
  collect() %>% 
  left_join (historic_ids, by=c("squirrel_id" = "old_squirrel_id")) %>% 
  mutate (squirrel_id = ifelse(is.na(new_squirrel_id), squirrel_id, new_squirrel_id),
          date = as.POSIXct(Date, format = '%Y-%m-%d %H:%M:%S')) # date looks to be in datetime format

# There are some problems with the squirrel_id numbers for some squirrels in this table.  
# They were checked based on tag numbers and locations and are manually changed below.  

midden_cones<- midden_cones %>% 
  mutate (squirrel_id = ifelse(squirrel_id == 19851, 19537, squirrel_id),
          squirrel_id = ifelse(squirrel_id == 19911, 11895, squirrel_id))

#########################
# Calculate Cache Sizes #
#########################
midden_cones<-midden_cones %>%   
  mutate(total_new_2019 = total_newclosed + total_newopen,
         total_new = coalesce(total_new, total_new_2019),
         total_cones = total_old+total_new,
         total_cones2 = total_old+total_newclosed,
         cache_size_total = (total_cones/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         cache_size_total2 = (total_cones2/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         cache_size_new = (total_new/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         cache_size_new_closed = (total_newclosed/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         cache_size_new_open = (total_newopen/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         cache_size_old = (total_old/no_quad) * ((pi * (width/2) * (length/2))/area_quad),
         squirrel_id = as.numeric(as.character(squirrel_id))) %>%  #Needed to match variable types.  The as.character is needed otherwise the squirrel_ids get recoded as 1...n
  group_by(squirrel_id, year) %>% # This line of code and the one below removes cases where a squirrel owns more than one midden
  slice(which.max(cache_size_total)) %>% # keeps the midden with more cones
  select(year, grid, midden, date, locx, locy, squirrel_id, cache_size_total, cache_size_total2, cache_size_new_open, cache_size_new_closed, cache_size_new, cache_size_old) 
# Haines recommended that total cache size be log10(x+1) transformed.

#################
# Link flastall #
#################

midden_cones<-midden_cones %>% 
  left_join(flastall, by="squirrel_id")

midden_cones = midden_cones %>% 
  select(year, grid, squirrel_id, contains("cache_size"))

write.csv(midden_cones,'/Users/matt/Downloads/midden_cones.csv')