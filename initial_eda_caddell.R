library(tidyverse)
library(ggthemes)

events <- read_csv("events.csv")

ggplot(events, aes(`Days to TCA`))+
  geom_histogram()+
  theme_minimal()

events %>% 
ggplot(aes(Pc, `Days to TCA`))+
  geom_point()+
  theme_minimal()

events %>% 
  filter(`Event Number` < 500,
         `Days to TCA` < 5) %>%
  arrange(`Event Number`, desc(`Days to TCA`)) %>% 
  ggplot(aes(`Days to TCA`, Pc))+
  geom_point(alpha = .5)+
  geom_line(aes(group = `Event Number`), alpha = .5, arrow = arrow(angle = 25))+
  scale_x_reverse()+
  theme_minimal()

event_summary <- events %>% 
  group_by(`Event Number`) %>% 
  summarize(observations= n(), 
            last_notice = min(`Days to TCA`), 
            first_notice = (max(`Days to TCA`)),
            Pc_max = max(Pc),
            Pc_min= min(Pc),
            Pc_at_last_notice = Pc[which.min(`Days to TCA`)],
            days_tracked = first_notice - last_notice,
            single_notice = as.factor(if_else(first_notice == last_notice, "Single Notice", "Multiple-Notices")))
  
ggplot(event_summary, aes(x = days_tracked))+
  geom_density(adjust = 2)+
  theme_minimal()

ggplot(event_summary, aes(last_notice))+
  geom_density(adjust = 2)+
  theme_minimal()

#57% of events appear to stop being of concern before the 2.5 day mark
event_summary %>% 
  summarise(stopped_tracking = sum(last_notice > 2.5, na.rm = TRUE), stopped_tracking_pct = stopped_tracking/n())

ggplot(event_summary, aes(first_notice))+
  geom_density(adjust = 2)+
  theme_minimal()

#only 80 events contain a Pc > .0001 within the 2 day mark ... that's .07% of the dataset
event_summary %>% 
  filter(last_notice < 2,
         Pc_at_last_notice > .0001) %>% 
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
events <- events %>% 
  mutate(pc = if_else(is.na(pc_best), pc_nom, pc_best))


events %>% 
  filter(!is.na(`Days to TCA`)) %>% 
ggplot(aes(log10(Pc)))+
  stat_ecdf()+
  #scale_x_continuous(limits = c(0,.0000025))+
  facet_wrap(~floor(`Days to TCA`), drop = TRUE, ncol = 1)+
  theme_minimal()


events %>% 
  filter(!is.na(`Days to TCA`)) %>% 
  ggplot(aes(Pc, color= as.factor(floor(`Days to TCA`))))+
  stat_ecdf()+
  scale_x_continuous(limits = c(0,.00001))+
  theme_minimal()+
  labs(color = "`Days to TCA`",
       y = "CDF")

#only 392 events with pc > .0001 or 
events %>% 
  filter(Pc > .0001) %>% 
  count()

#that's .09% of all reports
392/410889

#there is prob_cat that makes the number of fragments harder to use
#here I'll use the .5 as yes to make it easy
events <- events %>% 
  mutate(frag = if_else(ProbCatIfColl > .5, NumFragIfCatColl, NumFragIfNonCatColl),
         frag_breaks = cut(frag, breaks = c(0,10, 50, 100, 200, Inf)))

events %>% 
  filter(!is.na(`Days to TCA`)) %>% 
  select(Pc, frag_breaks, `Event Number`) %>% 
  unique() %>% 
  ggplot(aes(Pc, color = frag_breaks))+
  stat_ecdf()+
  theme_minimal()+
  labs(color = "Fragment",
       y = "CDF")

events %>% 
  filter(!is.na(`Days to TCA`)) %>% 
  select(frag_breaks, `Event Number`) %>% 
  unique() %>% 
  ggplot(aes(frag_breaks))+
  geom_histogram(stat = 'count')+
  theme_minimal()

#with a pc set of .0001 and only less than 5 days to collislion, there would have been 392 events
events %>% 
  filter(Pc > .0001,
         `Days to TCA` < 5)

#If I'm reading this right, there was never a PC above.0044
events %>% 
  arrange(desc(Pc))

ggplot(events, aes(Pc))+
  geom_histogram()+
  theme_minimal()

#its because all pc_nom above that threshold had small values for pc_best

event_summary %>% 
  count(Pc_at_last_notice > .001)
