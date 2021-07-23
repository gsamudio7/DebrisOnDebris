library(tidyverse)
library(dtplyr)
library(ggthemes)
library(rlang)

events <- read_csv("data/events.csv") %>% 
  mutate(PcFrag10000 = as.numeric(PcFrag10000))

tca_notice_threshold  = 5 #threshold for days to notification

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
  as_tibble() %>% 
  na.omit()

#this function allows us to get a confusion matrix with static data, but with uncertain collision if desired.
#We could add bootstrapping to it in the future if desired.
debris_cm <- function(sim_num = 1, Pc_threshold = .000001, collision_prob = .00015 , Frag_considered = "100" ) {
# Pc_threshold = .000001        #this is the level the decision maker sets for notification
# collision_prob = .00015     #this will end up being a probability pull
# Frag_considered = "100"      #which category we are looking at

collision_sim <- runif(n = nrow(event_summary), min = 0, max = 1)

final <- event_summary %>% 
  mutate(warning_issued = Pc_max >= Pc_threshold, #warning issued if Pc ever breaks threshold
         collision_prob = collision_sim[row_number()], #you can comment this on/off
         collision = Pc_last >= collision_prob, #collision_prob should be calculated every time in the future
         true_pos = warning_issued & collision,
         true_neg = !warning_issued & !collision,
         false_pos = warning_issued & !collision,
         false_neg = !warning_issued & collision) %>% 
  summarise(total_warnings = sum(warning_issued),
            total_collisions = sum(collision),
            true_pos_rate = sum(true_pos, na.rm = TRUE)/total_collisions,
            true_neg_rate = sum(true_neg, na.rm = TRUE)/(n()-total_collisions),
            false_pos_rate = sum(false_pos, na.rm = TRUE)/total_collisions,
            false_neg_rate = sum(false_neg, na.rm = TRUE)/(n()-total_collisions))

return(final)
}

#run the model several times
sim <- vector(mode = 'list', length = 1000) %>% 
  lapply(debris_cm)%>% 
  bind_rows()

sim

ggplot(sim, aes(total_warnings, total_collisions))+
  geom_point()

ggplot(sim, aes(total_collisions))+
  geom_histogram()+
  theme_minimal()
