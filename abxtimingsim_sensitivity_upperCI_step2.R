##########################################################################
# Project:      Simulations - Using upper confidence interval of 
#                   time to abx for probabilities
#
# Author:       Sarah Seelye
#
# Date Updated: 2023 July 19
##########################################################################


#install.packages("haven")
#install.packages("dplyr")
#install.packages("crosstable")
#install.packages("sampling")
#install.packages("table1")
library(sampling)
library(dplyr)
library(haven)
library(crosstable)
library(data.table)
library(parallel)
library(tidyverse)
library(table1)

n_percent = function(x, value = 1) {
  return(paste0(
    format(sum(x %in% value, na.rm = T), big.mark = ","), " (",
    format(round(sum(x %in% value, na.rm = T)/n()*100, digits = 1), nsmall = 1), "%)"))
}


memory.limit(size=56000)

#Start with harmonized KP and VA dataset
happi<-read_dta("kp_va_happi_aim3_upperCI.dta") #1559523 obs, 184 vars

##Study Sub-populations
# septic_shock
# severe_sepsis
# infection_cohort
# early_abx_stop
# never_treated

#Indicator for having received abx
happi$receivedabx<-ifelse(happi$never_treated==1, 0, 1)

#Hospital Rates of Septic Shock
septic_shock_freq<-happi%>%
  count(hospid, septic_shock)%>%
  group_by(hospid)%>%
  mutate(septic_shock_freq=n/sum(n))

septic_shock_freq<-septic_shock_freq[which(septic_shock_freq$septic_shock==1),c(1,4)] 

uniq_hosps<-as.data.frame(table(happi$hospid)) 
colnames(uniq_hosps)<-c("hospid", "N")
uniq_hosps<-left_join(uniq_hosps, septic_shock_freq, by="hospid")
uniq_hosps[is.na(uniq_hosps)]<-0 #A couple hospitals had an event rate of 0%

#152 hospitals. Bottom 26 hospitals are low rate, Top 26 are high rate (rank 127-152), remaining 100 are middle rate

uniq_hosps = uniq_hosps%>%
  arrange(septic_shock_freq)  %>%
  mutate(group_septic_shock = if_else(order(septic_shock_freq)>=127,"high",
                                      if_else(order(septic_shock_freq)<=26,"low",
                                              "middle")))  

#Join indicators of hospital type to the happi dataset
happi<-left_join(happi, uniq_hosps[,c(1,4)], by="hospid")

rates<-crosstable(happi, c(septic_shock,severe_sepsis, infection_cohort, early_abx_stop, never_treated), by=group_septic_shock)
rates[which(rates$variable==1),] 

#Create categorical variable of subpopulations from the 5 dummy variables
happi$type<-names(happi[c("septic_shock", "severe_sepsis", "infection_cohort", "early_abx_stop", "never_treated")])[max.col(happi[c("septic_shock", "severe_sepsis", "infection_cohort", "early_abx_stop", "never_treated")])]

#Antibiotic-associated outcome rates

table<-happi%>%
  group_by(type)%>%
  summarise(
    allergy_abx_ind=n_percent(allergy),
    thrombocytopenia_hosp=n_percent(ttp),
    leukopenia_hosp=n_percent(leukopenia),
    acute_renal_injury=n_percent(acute_renal_injury),
    acute_liver_injury=n_percent(acute_liver_injury),
    cdiff_infection=n_percent(cdiff),
    any_mdro_or_escr=n_percent(any_mdro_or_escr),
    mort90_ed=n_percent(mort90_ed),
    any_adverse_nomort=n_percent(any_adverse_nomort),
    any_adverse_wmort=n_percent(any_adverse_wmort),
    
    aod_lung=n_percent(aod_lung),
    aod_kidney=n_percent(aod_kidney),
    aod_liver=n_percent(aod_liver),
    aod_heme=n_percent(aod_heme),
    aod_lactate=n_percent(aod_lactate),
    renal=n_percent(renal),
    liver=n_percent(liver),
    chf=n_percent(chf)
  )

table<-as.data.frame(t(table)) 
colnames(table)<-table[1,]
table<-table[-1,]

#write.csv(table, "table with raw rates by subgroup_upperCI.csv", row.names=TRUE)

#Rates of each of the 5 subgroups for the low/median/high 
out<-happi%>%  
  count(group_septic_shock, type)%>%
  group_by(group_septic_shock)%>%
  mutate(freq=n/sum(n))  

out$sample_counts<-round(20000*out$freq, 0)
high<-out[which(out$group_septic_shock=="high"),]
middle<-out[which(out$group_septic_shock=="middle"),]
low<-out[which(out$group_septic_shock=="low"),]

#write.csv(out, "subgroup_rates_upperCI.csv")

#30-day mortality

#Save original and 50% faster TTA probabilities of 30-day mortality for sepsis/septic shock 
#into sepsis dataset (probabilities created in Stata do file) 
sepsis<-happi[which(happi$septic_shock==1 | happi$severe_sepsis==1),] #n=314654, 
sepsis<-as.data.frame(cbind(uniqid=sepsis$uniqid, pmort30_orig_ul=sepsis$pmort30_orig_ul, pmort30_fast_ul=sepsis$pmort30_fast_ul))

#Expected number of deaths
mean(sepsis$pmort30_orig_ul)*nrow(sepsis) #n=45155
sum(sepsis$pmort30_orig_ul)

### Simulation
all.results<-list()

sim<-function(hospital.type){ 
  
  ##Stratified sampling of 20,000 hospitalizations 
  counts<-t(hospital.type$sample_counts) 
  
  #Sort happi so that "counts" and "happi" datasets are in the same order
  happi<-happi[order(happi$type),] 
  
  #simple random sampling w/o replacement. stratified by type
  sample<-strata(happi, stratanames='type', counts, method = "srswor") 
  sample<-getdata(happi, sample) 
  
  #Output the median time to abx before and after speeding it up
  check<-sample[which(sample$type %in% c("septic_shock", "severe_sepsis" )),]
  
  #median time for sepsis patients
  orig.time.sepsis<-median(check$time_to_abx_min, na.rm = TRUE)
  new.time.sepsis<-median(check$tta_fast_ul, na.rm = TRUE)
  
  #median time for overall sample, which is impacted by the number of sepsis patients
  orig.time<-median(sample$time_to_abx_min, na.rm = TRUE)
  new.time<-median(sample$tta_fast_ul, na.rm = TRUE)
  
  median_difference_time_hr<-(orig.time.sepsis-new.time.sepsis)/60 
  
  #Absolute percent increase in the number treated (never treated to treated given spillover)
  additional_treated_percent<-median_difference_time_hr*c(.002, .004, .008, .016, .01, .05, .10 ) 
  
  #Number of newly treated (those who go from never treated to treated)
  number_newly_treated<-counts[3]*additional_treated_percent 
  
  #Report 30-day sepsis mortality using original and sped up TTA for septic shock and severe sepsis
  sepsis.sample<-sample[which(sample$septic_shock==1 | sample$severe_sepsis==1),] 
  
  #Expected number of sepsis deaths before and after increasing the time to abx
  new.mort30<-mean(sepsis.sample$pmort30_fast_ul)*nrow(sepsis.sample)
  old.mort30<-mean(sepsis.sample$pmort30_orig_ul)*nrow(sepsis.sample)
  
  #Expected probability of sepsis death before and after increasing the time to abx
  new.mort30prob<-mean(sepsis.sample$pmort30_fast_ul)*100
  old.mort30prob<-mean(sepsis.sample$pmort30_orig_ul)*100
  
  ### Compute the number of new abx events due to increasing the number treated
  never.treated.sample<-sample[which(sample$never_treated==1),] 
  
  #Randomly sample n=number_newly_treated of the never treated group to get converted to stopped early 
  newly.treated1<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[1]),] 
  newly.treated2<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[2]),]
  newly.treated3<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[3]),]
  newly.treated4<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[4]),]
  newly.treated5<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[5]),]
  newly.treated6<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[6]),]
  newly.treated7<-never.treated.sample[sample(nrow(sample[which(sample$never_treated==1),]), number_newly_treated[7]),]
  
  #treated probabilities  
  new.probs1<-newly.treated1[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat", 
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs2<-newly.treated2[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs3<-newly.treated3[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs4<-newly.treated4[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs5<-newly.treated5[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs6<-newly.treated6[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  
  new.probs7<-newly.treated7[,c("pallergy_treat", "pttp_treat", "pleukopenia_treat", 
                                "prenal_treat", "pliver_treat", "pcdiff_treat",  
                                "pmdro_treat", "pmdroblood_treat", "pmort90_treat",
                                "pany_adverse_nomort_treat", "pany_adverse_wmort_treat")]
  #original probabilities
  old.probs1<-newly.treated1[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs2<-newly.treated2[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs3<-newly.treated3[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs4<-newly.treated4[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs5<-newly.treated5[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs6<-newly.treated6[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  old.probs7<-newly.treated7[,c("pallergy_orig", "pttp_orig", "pleukopenia_orig", 
                                "prenal_orig", "pliver_orig", "pcdiff_orig",  
                                "pmdro_orig", "pmdroblood_orig", "pmort90_orig",
                                "pany_adverse_nomort_orig", "pany_adverse_wmort_orig")]
  
  #Output the number of expected events for the "new" and "old" predicted probabilities for the never treated
  new1<-as.data.frame(t(apply(new.probs1, 2, sum)))
  old1<-as.data.frame(t(apply(old.probs1, 2, sum)))
  
  new2<-as.data.frame(t(apply(new.probs2, 2, sum)))
  old2<-as.data.frame(t(apply(old.probs2, 2, sum)))
  
  new3<-as.data.frame(t(apply(new.probs3, 2, sum)))
  old3<-as.data.frame(t(apply(old.probs3, 2, sum)))
  
  new4<-as.data.frame(t(apply(new.probs4, 2, sum)))
  old4<-as.data.frame(t(apply(old.probs4, 2, sum)))
  
  new5<-as.data.frame(t(apply(new.probs5, 2, sum)))
  old5<-as.data.frame(t(apply(old.probs5, 2, sum)))
  
  new6<-as.data.frame(t(apply(new.probs6, 2, sum)))
  old6<-as.data.frame(t(apply(old.probs6, 2, sum)))
  
  new7<-as.data.frame(t(apply(new.probs7, 2, sum)))
  old7<-as.data.frame(t(apply(old.probs7, 2, sum)))
  
  #Output the expected number of additional abx events due to increasing the number treated
  abx.events1<-t(new1-old1) 
  abx.events2<-t(new2-old2)
  abx.events3<-t(new3-old3)
  abx.events4<-t(new4-old4)
  abx.events5<-t(new5-old5)
  abx.events6<-t(new6-old6)
  abx.events7<-t(new7-old7)
  
  #Number of new events due to stopping early (total events under stopped early-total events for never treated)
  abx.events<-as.data.frame(cbind(abx.events1, abx.events2, abx.events3, abx.events4,
                                  abx.events5, abx.events6, abx.events7))
  
  rownames(abx.events)<-c("allergy",
                          "ttp",
                          "leukopenia",
                          "acute_renal_injury",
                          "acute_liver_injury",
                          "cdiff",
                          "mdro",
                          "mdroblood",
                          "mort90_ed",
                          "any_adverse_nomort",
                          "any_adverse_wmort")
  
  colnames(abx.events)<-c("newly.treated.2", "newly.treated.4", "newly.treated.8", "newly.treated.16", 
                          "newly.treated1", "newly.treated5", "newly.treated10")
  
  setDT(abx.events, keep.rownames = "outcome")
  
  data1<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated.2, rate=.002))
  data2<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated.4, rate=.004))
  data3<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated.8, rate=.008))
  data4<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated.16, rate=.016))
  data5<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated1, rate=.01))
  data6<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated5, rate=.05))
  data7<-as.data.frame(cbind(outcome=abx.events$outcome, newly.treated=abx.events$newly.treated10, rate=.1))
  
  data1<-reshape(data1, idvar="rate", timevar = "outcome", direction = "wide")
  data2<-reshape(data2, idvar="rate", timevar = "outcome", direction = "wide")
  data3<-reshape(data3, idvar="rate", timevar = "outcome", direction = "wide")
  data4<-reshape(data4, idvar="rate", timevar = "outcome", direction = "wide")
  data5<-reshape(data5, idvar="rate", timevar = "outcome", direction = "wide")
  data6<-reshape(data6, idvar="rate", timevar = "outcome", direction = "wide")
  data7<-reshape(data7, idvar="rate", timevar = "outcome", direction = "wide")
  
  data<-rbind(data1, data2, data3, data4, data5, data6, data7)
  
  #Output the results
  results<-as.data.frame(cbind(
    number.deaths.original=old.mort30, number.deaths.faster=new.mort30, lives.saved=old.mort30-new.mort30,
    percent.deaths.orginal=old.mort30prob, percent.deaths.faster=new.mort30prob,
    newly.treated.2=number_newly_treated[1], newly.treated.4=number_newly_treated[2], 
    newly.treated.8=number_newly_treated[3], newly.treated.16=number_newly_treated[4],
    newly.treated1=number_newly_treated[5], newly.treated5=number_newly_treated[6],
    newly.treated10=number_newly_treated[7],
    orig.time=orig.time, new.time=new.time, 
    orig.time.sepsis=orig.time.sepsis, new.time.sepsis=new.time.sepsis
  ))
  
  all.results[[1]]<-results
  all.results[[2]]<-data
  names(all.results)<-c("results", "outcomes")
  
  return(all.results)
}


RNGkind("L'Ecuyer-CMRG")
iterations<-1000

set.seed(123)
sim.high<-mclapply(seq_len(iterations), function(hospital.type) 
  sim(hospital.type = high))

sim.high.results<-lapply(sim.high, "[[",1) 
sim.high.results<-do.call(rbind, lapply(sim.high.results, as.data.frame))
sim.high.outcomes<-lapply(sim.high, "[[",2)
sim.high.outcomes<-do.call(rbind, lapply(sim.high.outcomes, as.data.frame))
#write.csv(sim.high.results, "sim.high.results_upperCI.csv", row.names=TRUE)
#write.csv(sim.high.outcomes, "sim.high.outcomes_upperCI.csv", row.names=TRUE)


set.seed(123)
sim.middle<-mclapply(seq_len(iterations), function(hospital.type) 
  sim(hospital.type = middle))

sim.middle.results<-lapply(sim.middle, "[[",1)
sim.middle.results<-do.call(rbind, lapply(sim.middle.results, as.data.frame))
sim.middle.outcomes<-lapply(sim.middle, "[[",2)
sim.middle.outcomes<-do.call(rbind, lapply(sim.middle.outcomes, as.data.frame))
#write.csv(sim.middle.results, "sim.middle.results_upperCI.csv", row.names=TRUE)
#write.csv(sim.middle.outcomes, "sim.middle.outcomes_upperCI.csv", row.names=TRUE)


set.seed(123)
sim.low<-mclapply(seq_len(iterations), function(hospital.type) 
  sim(hospital.type = low))

sim.low.results<-lapply(sim.low, "[[",1)
sim.low.results<-do.call(rbind, lapply(sim.low.results, as.data.frame))
sim.low.outcomes<-lapply(sim.low, "[[",2)
sim.low.outcomes<-do.call(rbind, lapply(sim.low.outcomes, as.data.frame))
#write.csv(sim.low.results, "sim.low.results_upperCI.csv", row.names=TRUE)
#write.csv(sim.low.outcomes, "sim.low.outcomes_upperCI.csv", row.names=TRUE)


out1<-as.data.frame(apply(sim.high.results, 2, quantile,probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
out2<-as.data.frame(apply(sim.middle.results, 2, quantile,probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))
out3<-as.data.frame(apply(sim.low.results, 2, quantile,probs=c(0.1, 0.25, 0.5, 0.75, 0.9)))


#write.csv(out1, "Simulation results of high_upperCI.csv", row.names=TRUE)
#write.csv(out2, "Simulation results of middle_upperCI.csv", row.names=TRUE)
#write.csv(out3, "Simulation results of low_upperCI.csv", row.names=TRUE)


#Summarize the results of the abx-associated adverse outcomes

# high septic shock hospitals
allergy<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.allergy), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
allergy$condition<-"allergy"

thromb<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.ttp), probs=c(.1, .25, .5, .75, .9)))) %>%  
  unnest_wider(q)
thromb$condition<-"ttp"

leuk<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.leukopenia), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
leuk$condition<-"leuk"

renal<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_renal_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
renal$condition<-"renal"

liver<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_liver_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
liver$condition<-"liver"

cdiff<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.cdiff), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
cdiff$condition<-"cdiff"

mdro<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdro), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mdro$condition<-"mdro"

mdroblood<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdroblood), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mdroblood$condition<-"mdroblood"

mort90<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mort90_ed), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mort90$condition<-"mort90"

adversenomort<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_nomort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversenomort$condition<-"adversenomort"

adversewmort<-sim.high.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_wmort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversewmort$condition<-"adversewmort"

results.out<-cbind(type="high", rbind(allergy, thromb, leuk, renal, liver, cdiff, mdro, mdroblood, mort90, 
                                      adversenomort, adversewmort))


# middle septic shock hospitals

allergy<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.allergy), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
allergy$condition<-"allergy"

thromb<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.ttp), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
thromb$condition<-"ttp"

leuk<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.leukopenia), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
leuk$condition<-"leuk"

renal<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_renal_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
renal$condition<-"renal"

liver<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_liver_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
liver$condition<-"liver"

cdiff<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.cdiff), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
cdiff$condition<-"cdiff"

mdro<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdro), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mdro$condition<-"mdro"

mdroblood<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdroblood), probs=c(.1, .25, .5, .75, .9)))) %>%  
  unnest_wider(q)
mdroblood$condition<-"mdroblood"

mort90<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mort90_ed), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mort90$condition<-"mort90"

adversenomort<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_nomort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversenomort$condition<-"adversenomort"

adversewmort<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_wmort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversewmort$condition<-"adversewmort"


results.out2<-cbind(type="middle", rbind(allergy, thromb, leuk, renal, liver, cdiff, mdro, mdroblood, mort90,
                                         adversenomort, adversewmort))


#low septic shock hospitals

allergy<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.allergy), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
allergy$condition<-"allergy"

thromb<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.ttp), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
thromb$condition<-"ttp"

leuk<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.leukopenia), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
leuk$condition<-"leuk"

renal<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_renal_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
renal$condition<-"renal"

liver<-sim.middle.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.acute_liver_injury), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
liver$condition<-"liver"

cdiff<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.cdiff), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
cdiff$condition<-"cdiff"

mdro<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdro), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mdro$condition<-"mdro"

mdroblood<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mdroblood), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
mdroblood$condition<-"mdroblood"

mort90<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.mort90_ed), probs=c(.1, .25, .5, .75, .9)))) %>%  
  unnest_wider(q)
mort90$condition<-"mort90"

adversenomort<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_nomort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversenomort$condition<-"adversenomort"

adversewmort<-sim.low.outcomes %>% 
  group_by(rate) %>% 
  summarise(q = list(quantile(as.numeric(newly.treated.any_adverse_wmort), probs=c(.1, .25, .5, .75, .9)))) %>% 
  unnest_wider(q)
adversewmort$condition<-"adversewmort"

results.out3<-cbind(type="low", rbind(allergy, thromb, leuk, renal, liver, cdiff, mdro, mdroblood, mort90,
                                      adversenomort, adversewmort))

results<-rbind(results.out, results.out2, results.out3)


#write.csv(results, "Simulation results of abx outcomes_upperCI.csv", row.names=FALSE)

