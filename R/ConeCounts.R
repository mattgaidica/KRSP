ConeCounts = function(grids = c('KL', 'SU', 'BT', 'JO','AG', 'LL'), years = seq(1988,2020,1)){
  #has not been optimized by MRW. Clunky as hell
  cone_counts<-tbl(con, "cones") %>%
    filter(Year>=1988) %>%
    collect() %>%
    mutate(Year = as.numeric(Year), 
           LocX = as.numeric(LocX), 
           LocY = as.numeric(LocY), 
           DBH = as.numeric(DBH), 
           Per = as.numeric(Per), 
           NumNew = as.numeric(NumNew),
           cone_index = log(NumNew + 1),
           total_cones = 1.11568 * exp(0.1681 + 1.1891 * log(NumNew + 0.01)) # according to Krebs et al. 2012
    )
  
  ##################################
  # Means calculated per Grid Year #
  ##################################
  cones_grids_years <- cone_counts %>% group_by(Grid, Year) %>% 
    summarize(num_trees = sum(!is.na(NumNew)),
              cone_counts = mean(NumNew, na.rm = TRUE),
              cone_index = mean(cone_index, na.rm = TRUE)) %>% 
    mutate(Year_tp1 = Year+1,
           cone_index_t = ifelse(is.finite(cone_index), cone_index, NA))
  
  #link in cones from the previous year
  cone_temp<-cones_grids_years %>% 
    select(Grid, Year, Year_tp1, cone_index_tm1=cone_index_t)
  
  cones_grids_years<-left_join(cones_grids_years, cone_temp, by=c("Grid", "Year" = "Year_tp1")) %>% 
    select(-Year_tp1, -Year.y)
  
  # Manually code mast years
  # Modified to be cleaner and faster
  
  MAST_YEARS = c(1993, 1998, 2005, 2010, 2014, 2019)
  
  cones_grids_years<-cones_grids_years %>%
    mutate (mast = case_when(
      Year %in% MAST_YEARS ~ "y",
      TRUE ~ "n")
    ) %>% 
    mutate (Exp = "c") %>% 
    mutate (Exp = ifelse(Grid=="AG"&Year>2004&Year<2018, "f", Exp),
            Exp = ifelse(Grid=="JO"&Year>2006&Year<2013, "f", Exp),
            Exp = ifelse(Grid=="LL"&Year>2005&Year<2012, "f", Exp)) %>% 
    mutate (EXP_label = 1) %>% 
    mutate (EXP_label = ifelse(Exp=="f", 19, EXP_label)) %>% filter(Grid %in% c("KL", "SU", "JO", "BT", "SUX", "AG", "LL", "CH", "RR"))
  cones_grids_years = cones_grids_years %>% 
    mutate(cone_index_tm1 = ifelse(Year == 2005 & Grid == "AG", 1.045008523, cone_index_tm1))
  colnames(cones_grids_years) = tolower(colnames(cones_grids_years))
  cones_grids_years = cones_grids_years %>% filter(grid %in% grids, year %in% years) %>% ungroup()
  return(cones_grids_years)
}

