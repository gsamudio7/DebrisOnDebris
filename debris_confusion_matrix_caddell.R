library(tidyverse)
library(dtplyr)
library(ggthemes)
library(rlang)

events <- read_csv("data/events.csv") %>% 
  mutate(PcFrag10000 = as.numeric(PcFrag10000))

tca_notice_threshold  = 5   #earliest day to get notifed
Pc_threshold = .00001        #this is the level the decision maker sets for notification
collision_prob = .00015     #this will end up being a probability pull
Frag_considered = "100"      #which category we are looking at

event_summary <- events %>% 
  lazy_dt() %>% 
  filter(time2TCA < tca_notice_threshold) %>% 
  select(eventNumber, time2TCA, matches(paste0(Frag_considered,"$"))) %>% #only return cols with Frag Considered
  group_by(eventNumber) %>% 
  summarize(observations= n(), 
            Frag = Frag_considered,
            last_notice = min(time2TCA),
            first_notice = max(time2TCA),
            Pc_max = max(.data[[paste0("PcFrag",Frag_considered)]]),
            Pc_first = .data[[paste0("PcFrag",Frag_considered)]][which.max(time2TCA)],      #This funky call allows the frag to be dynamic based off user input
            Pc_last= .data[[paste0("PcFrag",Frag_considered)]][which.min(time2TCA)]) %>% 
  as_tibble()

final <- event_summary %>% 
  na.omit() %>% 
  mutate(collision = Pc_last > collision_prob, #collision_prob should be calculated every time in the future
         warning_issued = Pc_max >= Pc_threshold, #warning issued if Pc ever breaks threshold
         true_pos = warning_issued & collision,
         true_neg = warning_issued & !collision,
         false_pos = !warning_issued & collision,
         false_neg = !warning_issued & !collision) %>% 
  summarise(total_warnings = sum(warning_issued),
            total_collsions = sum(collision),
            true_pos_rate = sum(true_pos, na.rm = TRUE)/total_warnings,
            true_neg_rate = sum(true_neg, na.rm = TRUE)/total_warnings,
            false_pos_rate = sum(false_pos, na.rm = TRUE)/(n()-total_warnings),
            false_neg_rate = sum(false_neg, na.rm = TRUE)/(n()-total_warnings))
