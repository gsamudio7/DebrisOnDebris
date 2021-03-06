library(tidyverse)
library(dtplyr)
library(ggthemes)
library(rlang)

events <- read_csv("data/events.csv") %>% 
  mutate(PcFrag10000 = as.numeric(PcFrag10000))

tca_notice_threshold  = 5 #threshold for days to notification
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
  as_tibble() %>% 
  na.omit()

#this function allows us to get a confusion matrix with static data, but with uncertain of collision if desired.
#We could add bootstrapping to it in the future if desired.
debris_cm <- function(Pc_threshold = .000001, collision_prob = .00015, sim_num = 1) {
# Pc_threshold = .000001        #this is the level the decision maker sets for notification
# collision_prob = .00015     #this will end up being a probability pull

collision_sim <- runif(n = nrow(event_summary), min = 0, max = 1)

final <- event_summary %>% 
  na.omit() %>% 
  mutate(warning_issued = Pc_max >= Pc_threshold, #warning issued if Pc ever breaks threshold
         #collision_prob = collision_sim[row_number()], #you can comment this on/off
         collision = Pc_last >= collision_prob, #collision_prob should be calculated every time in the future
         true_pos = warning_issued & collision,
         true_neg = !warning_issued & !collision,
         false_pos = warning_issued & !collision,
         false_neg = !warning_issued & collision) %>% 
  summarise(total_warnings = sum(warning_issued),
            total_collisions = sum(collision),
            Pc_threshold = Pc_threshold,
            true_pos_rate = sum(true_pos)/total_collisions,
            true_neg_rate = sum(true_neg)/(n()-total_collisions),
            false_pos_rate = sum(false_pos)/(n()-total_collisions),
            false_neg_rate = sum(false_neg)/(total_collisions))

return(final)
}

#run the model several times
sim <- vector(mode = 'list', length = 1000) %>% 
  lapply(debris_cm)%>% 
  bind_rows()

sim

ggplot(sim, aes(total_warnings, total_collisions))+
  geom_jitter(width = .1, height = 0)+
  theme_minimal()

ggplot(sim, aes(total_collisions))+
  geom_histogram(binwidth = .5)+
  theme_minimal()

ggplot(sim, aes(true_neg_rate, false_neg_rate))+
  geom_point()+
  theme_minimal()

ggplot(sim, aes(true_pos_rate, false_pos_rate))+
  geom_point()+
  theme_minimal()

#ROC development
#creates a range of thresholds to test
#Pc_thresholds <- data.frame(Pcthreshold = seq(from = .00000000000000000001, to = .0000001, by = .000000001))
Pc_thresholds <- data.frame(Pcthreshold = c(.0000001, .0000001, .000001, .00001, .0001, .001, .01))


#runs the model with the different thresholds
roc_data <- pmap_dfr(list(Pc_thresholds[1:100,]), debris_cm)

#draws ROC curve
roc_data %>% 
  add_row(false_pos_rate = 0, true_pos_rate = 0) %>% 
  add_row(false_pos_rate = 1, true_pos_rate = 1) %>% #so the line starts at 0,0 and goes to 1,1
  arrange(false_pos_rate, true_pos_rate) %>% 
  ggplot(aes(x = false_pos_rate, y = true_pos_rate))+
  geom_line()+
  theme_minimal()

#only works if there are collisions. 
roc_data %>% 
  summarise(sum(total_collisions))
