library(tidyverse)
library(ggthemes)

events <- read_csv("DebrisOnDebris_with_names.csv") %>% 
  mutate(pc = if_else(is.na(pc_best), pc_nom, pc_best))

event_summary <- events %>% 
  group_by(event_number) %>% 
  summarize(days_tracked = n(), 
            last_notice = min(days_to_tca), 
            first_notice = (max(days_to_tca)),
            pc_at_last_notice = pc[which.min(days_to_tca)],
            pc_max = max(pc),
            pc_min = min(pc),
            lead = first_notice - last_notice)


pc_threshold = .0001

sim <- event_summary %>% 
  mutate(notice_issued = if_else(pc_max > pc_threshold,1,0),
         collision_happend = if_else(pc_at_last_notice > runif(n = 1, min = 0, max = 1),1,0))

summary(sim) #look at the collisions that happend

table(sim$notice_issued, sim$collision_happend)



