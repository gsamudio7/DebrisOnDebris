# Imports ####

library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

load("data/concernData.RData")

# How does quantity change over time ####
concernSummary_at_TOI[fragNum=="PcBest" & Pc_at_TOI >= 1e-5,.N,by=TCA_Bin] %>% 
  plot_ly(
    type="bar",
    x=~TCA_Bin,
    y=~N,
    color=I("royalblue4")
  )




# How does quantity of NEW concern events change over time
concernEvents <- concernSummary_at_TOI[fragNum=="PcBest" & Pc_at_TOI >= 1e-5]

newConcernEvents <- list("10"=concernEvents[TCA_Bin==10,unique(eventNumber)])
for (i in (9:1)) {
  oldEvents <- concernEvents[TCA_Bin==i + 1,unique(eventNumber)] 
  newEvents <- concernEvents[TCA_Bin==i,unique(eventNumber)] 
  newConcernEvents[[as.character(i)]] = setdiff(newEvents,oldEvents)
}
newConcernSummary <- merge(
  data.table("TCA_Bin"=as.factor(10:1),
             "newCount"=unlist(lapply(newConcernEvents,length))),
  concernEvents[,.(concernCount=.N),by=TCA_Bin],
  by="TCA_Bin"
) 

newConcernSummary_plot <- newConcernSummary %>% 
  plot_ly(
    type="bar",
    x=~TCA_Bin,
    y=~concernCount,
    color=I("#222222"),
    name="Total",
    hoverinfo="none"
  ) %>%
  
  add_trace(
    y=~newCount,
    name="New",
    color=I("#2359c4")
  ) %>%
  
  layout(
    xaxis=list(title="<b>Days to TCA Bin</b>",
               gridcolor="#333333"),
    yaxis=list(title="<b>Concern Events</b>",
               gridcolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )



save(newConcernSummary_plot,file="products/newConcernSummary_plot.RData")


# Count the number of new events plus events that are carried through
'%ni%' <- Negate('%in%')
events_at_5 <- concernSummary_at_TOI[TCA_Bin==5,unique(eventNumber)]
length(events_at_5)
events_at_5_and_1 <- concernSummary_at_TOI[eventNumber %in% events_at_5 & TCA_Bin==1,unique(eventNumber)]
length(events_at_5_and_1)
events_at_1_and_not_at_5 <- concernSummary_at_TOI[eventNumber %ni% events_at_5 & TCA_Bin==1,
                                                  unique(eventNumber)]
length(events_at_1_and_not_at_5)
events_at_5_and_not_at_1 <- setdiff(events_at_5,concernSummary[,unique(eventNumber)])

# Verify
length(events_at_5_and_1) + length(events_at_1_and_not_at_5) == concernSummary_at_TOI[TCA_Bin==1,uniqueN(eventNumber)]

# Wrangle data to show
## event PC at 5 / PC at 1 / makes it through (Y/N)
deltaSummary <- rbindlist(list(
  concernSummary_at_TOI[eventNumber %in% events_at_5_and_1 &
                          fragNum=="PcBest" &
                          TCA_Bin==5, .(eventNumber,TCA_Bin,Pc=Pc_at_TOI,Censored=FALSE)],
  concernSummary[eventNumber %in% events_at_5_and_1 &
                   fragNum=="PcBest", .(eventNumber,TCA_Bin=1,Pc=Pc_min,Censored=FALSE)],
  concernSummary_at_TOI[eventNumber %in% events_at_5_and_not_at_1 &
                          fragNum=="PcBest" &
                          TCA_Bin==5, .(eventNumber,TCA_Bin,Pc=Pc_at_TOI,Censored=TRUE)],
  data.table("eventNumber"=events_at_5_and_not_at_1,
             "TCA_Bin"=1,
             "Pc"=0,
             "Censored"=TRUE)
))
setorder(deltaSummary,-"TCA_Bin")
delta_plot_summary <- deltaSummary[,.(delta=diff(Pc)),by=c("eventNumber","Censored")] 



delta_plot_summary[,"Delta" := ifelse(delta > 0, "> 0",
                                      ifelse(delta == 0, "= 0","< 0"))]


deltaCountSummary <- merge(
  data.table("Total"=c(delta_plot_summary[Censored==TRUE,.N],
                       delta_plot_summary[Censored==FALSE,.N]),
             "Censored"=c(TRUE,FALSE)
  ),
  
  delta_plot_summary[,.(Count=.N),by=c("Censored","Delta")],
  by="Censored"
)
deltaCountSummary[,"Proportion" := round(Count/Total,3)]
deltaCountSummary[,!"Total"]

delta_plot_summary %>%
  plot_ly(
    type="histogram",
    x=~delta,
    nbinsx=10,
    color=~Censored
  ) %>%
  layout(
    xaxis = list(title="<b>Collision Probability Delta\nat 5 and 1 days to TCA</b>",
                 tickvals = seq(-700,0,100),
                 ticktext = seq(-700,0,100) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10),
                 gridcolor="#333333"))


# What is the distribution of Pc at 5 for events that make it to 1
concernSummary_at_TOI[fragNum=="PcBest" & 
                        eventNumber %in% events_at_5_and_1 & 
                        TCA_Bin == 5 |
                        
                        fragNum=="PcBest" & 
                        eventNumber %in% events_at_5_and_1 & 
                        TCA_Bin == 1] %>%
  plot_ly(
    type="box",
    x=~TCA_Bin,
    y=~log(Pc_at_TOI)
  ) 

# What is the distribution of the delta?
concernSummary_Pc_at_1_from_5 <- concernSummary_at_TOI[
  fragNum=="PcBest" & 
    eventNumber %in% events_at_5_and_1 & 
    TCA_Bin == 1 |
    
    fragNum=="PcBest" & 
    eventNumber %in% events_at_5_and_1 & 
    TCA_Bin == 5] %>% setorder(eventNumber,-TCA_Bin)

concernSummary_diff_at_1_from_5 <- concernSummary_Pc_at_1_from_5[
  ,.(Delta=diff(Pc_at_TOI)),by=eventNumber] 

# How many 0
concernSummary_diff_at_1_from_5[Delta==0,.N]/concernSummary_diff_at_1_from_5[,.N] # 0.64

# How many increased
concernSummary_diff_at_1_from_5[Delta > 0,.N]/concernSummary_diff_at_1_from_5[,.N] # 0.085

# How many decreased
concernSummary_diff_at_1_from_5[Delta < 0,.N]/concernSummary_diff_at_1_from_5[,.N] 

# Box Plot
concernSummary_Pc_at_1_from_5 %>%
  plot_ly(
    type="box",
    y=~log(Pc_at_TOI),
    x=~as.integer(TCA_Bin) - 1
  ) %>%
  layout(
    yaxis = list(title="<b>Probability of Collision</b>",
                 tickvals = seq(-700,0,100),
                 ticktext = seq(-700,0,100) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10)),
    xaxis = list(title="<b>Days to TCA</b>")
  )

concernSummary_Pc_at_1_from_5 %>%
  plot_ly(
    type="scatter",
    mode="lines",
    y=~log(Pc_at_TOI),
    x=~as.integer(TCA_Bin) - 1
  ) 

# Conclusion
## No obvious advantage in having a higher Pc_warn threshold at 5 days than the Pc_concern rate

# Question remains, what Pc_warn value at 5 days is the ideal to use to 
## minimize False negatives and minimize False Positives




# Plot
Pc_by_TCA_plot <- 
  
  concernSummary_at_TOI[fragNum=="PcBest"] %>%
  plot_ly(
    type="box",
    y=~log(Pc_at_TOI),
    x=~TCA_Bin,
    hoverinfo="none",
    color=I("#2359c4")
  ) %>%
  layout(
    yaxis = list(title="<b>Collision Probability</b>",
                 tickvals = seq(-30,-6,2),
                 ticktext = seq(-30,-6,2) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10),
                 gridcolor="#333333"),
    xaxis = list(title="<b>Days to TCA</b>",
                 gridcolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )

save(Pc_by_TCA_plot,file="products/Pc_by_TCA_plot.RData")

# Count total positives (events of concern) that we have visibility on at 5 days out
Pc_concern <- 1e-5
(P_at_5 <- concernSummary_at_TOI[TCA_Bin==5 & PcBest_at_TOI >= Pc_concern,.N]) # 21
(N_at_5 <- concernSummary_at_TOI[TCA_Bin==5 & PcBest_at_TOI < Pc_concern,.N]) # 29323

(P <- concernSummary[fragNum=="PcBest" & Pc_min >= Pc_concern,.N]) # 58
(N <- concernSummary[fragNum=="PcBest" & Pc_min < Pc_concern,.N]) # 48193