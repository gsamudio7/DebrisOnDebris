library(tidyverse)
library(ggthemes)
library(dtplyr)

events <- read_csv("data/events.csv") %>% 
  mutate(PcFrag10000 = as.numeric(PcFrag10000))


ggplot(events, aes(time2TCA))+
  geom_histogram()+
  theme_minimal()

events %>% 
ggplot(aes(PcBest, time2TCA))+
  geom_point()+
  theme_minimal()

events %>% 
  filter(eventNumber < 5000,
         time2TCA < 5) %>%
  arrange(eventNumber, desc(TCA)) %>% 
  ggplot(aes(time2TCA, PcBest))+
  geom_point(alpha = .5)+
  geom_line(aes(group = eventNumber), alpha = .5, arrow = arrow(angle = 25))+
  scale_x_reverse()+
  theme_minimal()

tca_of_concern <- 5

event_summary <- events %>% 
  filter(time2TCA < 5) %>% 
  group_by(eventNumber) %>% 
  summarize(observations= n(), 
            last_notice = min(time2TCA), 
            first_notice = max(time2TCA),
            PcBest_max = max(PcBest),
            PcBest_min= min(PcBest),
            PcBest_at_last_notice = PcBest[which.min(time2TCA)],
            PcBest_at_tca_of_concern = PcBest[which.max(time2TCA)],
            PcBest_range = PcBest_at_tca_of_concern - PcBest_at_last_notice,
            PcFrag10_max = max(PcFrag10),
            PcFrag10_min= min(PcFrag10),
            PcFrag10_at_last_notice = PcFrag10[which.min(time2TCA)],
            PcFrag100_max = max(PcFrag100),
            PcFrag100_min= min(PcFrag100),
            PcFrag100_at_last_notice = PcFrag100[which.min(time2TCA)],
            PcFrag1000_max = max(PcFrag1000),
            PcFrag1000_min= min(PcFrag1000),
            PcFrag1000_at_last_notice = PcFrag1000[which.min(time2TCA)],
            PcFrag10000_max = max(as.numeric(PcFrag10000)),
            PcFrag10000_min= min(as.numeric(PcFrag10000)),
            PcFrag10000_at_last_notice = as.numeric(PcFrag10000[which.min(time2TCA)]),
            days_tracked = first_notice - last_notice,
            TCA = min(TCA),
            single_notice = as.factor(if_else(first_notice == last_notice, "Single Notice", "Multiple-Notices"))) %>% 
  unique()

event_summary %>% 
  filter(observations > 1) %>% 
ggplot(aes(x = PcBest_range))+
  geom_histogram(binwidth = .00001)+
  theme_minimal()

event_summary %>% 
  filter(observations > 1) %>% 
  ggplot(aes(x = PcBest_range))+
  geom_density()+
  scale_x_continuous(limits = c(-.00001, .00001))+
  theme_minimal()


ggplot(event_summary, aes(x = days_tracked))+
  geom_density(adjust = 2)+
  theme_minimal()

ggplot(event_summary, aes(last_notice))+
  geom_density(adjust = 2)+
  theme_minimal()

#60% of events appear to stop being of concern before the 2.5 day mark
event_summary %>% 
  summarise(stopped_tracking = sum(last_notice > 2.5, na.rm = TRUE), stopped_tracking_PcBest = stopped_tracking/n())


#only 11 events contain a PcBest_at_last > .0001 within the 5 day mark ... that's .005% of the dataset
event_summary %>% 
  filter(last_notice < 5,
         PcBest_at_last_notice > .0001) %>% 
  count()

event_summary %>% 
  filter(PcBest_at_last_notice > .0001) %>% 
  count()

event_summary %>% 
  arrange(desc(first_notice,last_notice)) %>% 
  mutate(order = row_number()) %>% 
  pivot_longer(cols = c(last_notice, first_notice), names_to = "notice", values_to = "days") %>%
  select(days, order, notice, single_notice) %>% 
  drop_na() %>% 
  ggplot(aes(days, order))+
  geom_line(aes(group = order), color = 'lightgrey')+
  geom_point(aes(color = notice), size = 1.5)+
  scale_color_colorblind()+
  scale_x_continuous(breaks = seq(from = 0 , to = 10, by = 1))+
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(y = "Event",
       x = "Days",
       color = "Notice")+
  facet_grid(rows = vars(single_notice))

#suggested analysis from Matt

events %>% 
ggplot(aes(log10(PcBest)))+
  stat_ecdf()+
  #scale_x_continuous(limits = c(0,.0000025))+
  facet_wrap(~floor(time2TCA), drop = TRUE, ncol = 1)+
  theme_minimal()


events %>% 
  filter(!is.na(TCA)) %>% 
  ggplot(aes(PcBest, color= as.factor(floor(time2TCA))))+
  stat_ecdf()+
  scale_x_continuous(limits = c(0,.00001))+
  theme_minimal()+
  labs(color = "TCA",
       y = "CDF")

#only 27 events with PcBest > .0001 or 
events %>% 
  filter(PcBest > .0001) %>% 
  count()

#that's .004% of all reports
27/nrow(events)*100

#there is prob_cat that makes the number of fragments harder to use
#here I'll use the .5 as yes to make it easy
# events <- events %>% 
#   mutate(frag = if_else(ProbCatIfColl > .5, NumFragIfCatColl, NumFragIfNonCatColl),       #this doesn't work anymore no split for numfrag
#          frag_breaks = cut(frag, breaks = c(0,10, 50, 100, 200, Inf)))

events <- events %>%
  mutate(frag_breaks = cut(NumFrag, breaks = c(0,10, 50, 100, 200, Inf)))

events %>% 
  filter(!is.na(TCA)) %>% 
  select(PcBest, frag_breaks, eventNumber) %>% 
  unique() %>% 
  ggplot(aes(PcBest, color = frag_breaks))+
  stat_ecdf()+
  theme_minimal()+
  labs(color = "Fragment",
       y = "CDF")

events %>% 
  filter(!is.na(TCA)) %>% 
  select(frag_breaks, eventNumber) %>% 
  unique() %>% 
  ggplot(aes(frag_breaks))+
  geom_histogram(stat = 'count')+
  theme_minimal()

#with a PcBest set of .0001 and only less than 5 days to collislion, there would have been 392 events
events %>% 
  filter(PcBest > .0001,
         TCA < 5)

#If I'm reading this right, there was never a PcBest above.0009
events %>% 
  arrange(desc(PcBest))

ggplot(events, aes(PcBest))+
  geom_histogram()+
  theme_minimal()

#its because all PcBest_nom above that threshold had small values for PcBest_best

event_summary %>% 
  filter(last_notice < 5) %>% 
  count(PcBest_max > .001)


event_summary %>% 
  filter(last_notice < 5,
         frag > 10) %>% 
  summarise(total_PcBest = sum(PcBest_at_last_notice, na.rm = TRUE))

event_summary %>% 
  filter(first_notice >= 5,
         PcBest_at_last_notice > .0001)
