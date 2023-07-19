*------------------------------------
* Program Setup and Import
*------------------------------------

version 16
cap log close
set more off
clear all
set linesize 80

cd ""

local c_date = c(current_date)
local c_time = c(current_time)
local c_time_date = "`c_date'"+"_" +"`c_time'"
local time_string = subinstr("`c_time_date'", ":", "_", .)
local time_string = subinstr("`time_string'", " ", "_", .)
display "`time_string'"

log using "abxtimingsim_lower_step1_`time_string'.log", replace


* ------------------------------------------------------------------------------
*
*  	Preparing Dataset for Simulations using Lower CI for Time to Antibiotics
*
* 	Created: 		2023 July 13
* 	Last updated: 	2023 July 19
* 	Author: 		S. Seelye
*
*-------------------------------------------------------------------------------

* open dataset
use aim1_analytic_dataset_20211102, clear		

* count
count  // 1,560,126
	
* drop sample variable
drop sample 

* exclude hospitalizations in hospitals with fewer than 300 hospitalizations    
bysort hospid: gen counthosp = _N
drop if counthosp<300
drop counthosp 
count // 1,559,523
 
* rename outcome variables 
rename allergy_abx_ind allergy 
rename thrombocytopenia_hosp ttp
rename leukopenia_hosp leukopenia 
rename cdiff_infection cdiff 

* examine outcomes by subpopulation
gen subpop = .
replace subpop = 1 if never_treated
replace subpop = 2 if early_abx_stop
replace subpop = 3 if infection_cohort
replace subpop = 4 if severe_sepsis
replace subpop = 5 if septic_shock 

lab def subpop 1 "never treated" 2 "early abx stop" 3 "infection" 4 "severe sepsis" 5 "septic shock"
lab val subpop subpop

* any adverse event (without 90-day mortality)
gen any_adverse_nomort = inlist(1, allergy, ttp, leukopenia, ///
								  acute_renal_injury, acute_liver_injury, ///
								  cdiff, any_mdro_or_escr)
tab any_adverse_nomort								  
tab subpop any_adverse_nomort , ro


local aod 																		///
		 aod_lactate aod_kidney pressor_in_72hr aod_liver		///
		 aod_heme aod_lung 			
		
local sirs																		///
		sirs_temp sirs_rr sirs_pulse sirs_wbc

local comorbid 															///
		cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
		renal htn cardic_arrhym valvular_d2 pulm_circ pvd paralysis 	///
		pud hypothyroid ah lymphoma ra coag obesity wtloss fen			///
		anemia_cbl anemia_def etoh drug psychoses depression 
		
logit any_adverse_nomort i.abx_in_48hr c.age male `sirs' `aod' `comorbid' va , or
margins abx_in_48hr	

* any adverse event (with 90-day mortality)
gen any_adverse_wmort = inlist(1, allergy, ttp, leukopenia, ///
								  acute_renal_injury, acute_liver_injury, ///
								  cdiff, any_mdro_or_escr, mort90_ed)
tab any_adverse_wmort								  
tab subpop any_adverse_wmort , ro


local aod 																		///
		 aod_lactate aod_kidney pressor_in_72hr aod_liver		///
		 aod_heme aod_lung 			
		
local sirs																		///
		sirs_temp sirs_rr sirs_pulse sirs_wbc

local comorbid 															///
		cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
		renal htn cardic_arrhym valvular_d2 pulm_circ pvd paralysis 	///
		pud hypothyroid ah lymphoma ra coag obesity wtloss fen			///
		anemia_cbl anemia_def etoh drug psychoses depression 
		
logit any_adverse_wmort i.abx_in_48hr c.age male `sirs' `aod' `comorbid' va , or
margins abx_in_48hr	

*-------------------------------------------------------------------------------
* Predict mortality under baseline and 50% faster TTA
*-------------------------------------------------------------------------------
* expand the dataset with a new observation for each patient
expand 2 , gen(dupindicator)

* for the duplicates, reduce TTA by 50%
replace time_to_abx_min = time_to_abx_min*0.5 if dupindicator==1

sort uniqid dupindicator
	
*-------------------------------------------------------------------------------
* Use lower CI of time-to-antibiotics to create new adjusted probabilities
*-------------------------------------------------------------------------------
local aod 																		///
		 aod_lactate aod_kidney pressor_in_72hr aod_liver		///
		 aod_heme aod_lung 			
		
local sirs																		///
		sirs_temp sirs_rr sirs_pulse sirs_wbc

local comorbid 															///
		cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
		renal htn cardic_arrhym valvular_d2 pulm_circ pvd paralysis 	///
		pud hypothyroid ah lymphoma ra coag obesity wtloss fen			///
		anemia_cbl anemia_def etoh drug psychoses depression 
		
logit mort30_ed time_to_abx_min c.age male `sirs' `aod' `comorbid' va if dupindicator==0 & (subpop==4 | subpop==5)

* convert time_to_abx coef to lower CI to use to obtain predicted probabilities
mat c = r(table)
mat list c
mat d = c[5,1]
mat list d
matrix c[1,1] = d
matrix list c

forval i=1/45 {
    local betall`i' = c[1,`i']   
}
display `betall1' 
display `betall2' 
display `betall3'
display `betall4'
display `betall5'
display `betall44'
display `betall45'

gen manual_logodds_ll = time_to_abx_min*`betall1' + age*`betall2' + male*`betall3' +    ///
                  sirs_temp*`betall4' + sirs_rr*`betall5' + sirs_pulse*`betall6' +   ///
                  sirs_wbc*`betall7' +                                         ///
                  aod_lactate*`betall8' + aod_kidney*`betall9' +                  ///
                  pressor_in_72hr*`betall10' + aod_liver*`betall11' + 	            ///
                  aod_heme*`betall12' + aod_lung*`betall13' +                       ///
                  cancer_nonmet*`betall14' + cancer_met*`betall15' + pulm*`betall16' +  ///
                  chf*`betall17' + dm_uncomp*`betall18' + dm_comp*`betall19' +	        ///
                  liver*`betall20' + neuro*`betall21' + renal*`betall22' + htn*`betall23' + ///
                  cardic_arrhym*`betall24' + valvular_d2*`betall25' +               ///
                  pulm_circ*`betall26' + pvd *`betall27' + paralysis*`betall28' + 	    ///
                  pud*`betall29' + hypothyroid*`betall30' + ah*`betall31' +             ///
                  lymphoma*`betall32' + ra*`betall33' + coag*`betall34' +               ///
                  obesity*`betall35' + wtloss*`betall36' + fen*`betall37' +			    ///
                  anemia_cbl*`betall38' + anemia_def*`betall39' + etoh*`betall40' +     ///
                  drug*`betall41' + psychoses*`betall42' + depression*`betall43' +      ///
                  va*`betall44' +  `betall45' 

gen manual_prob_ll = (exp(manual_logodds_ll))/(1+exp(manual_logodds_ll))                 
drop manual_logodds_ll

rename manual_prob_ll pmort30_orig_ll

* create new variable for p-hat faster 
gen pmort30_fast_ll = pmort30_orig_ll if dupindicator==1

* create new variable for TTA 50% faster 
gen tta_fast_ll = time_to_abx_min if dupindicator==1
gen tta_fast_hr_ll = tta_fast_ll/60

* recode the original observation so that it fills in missing values for 
* pmort30_fast and tta_fast
bysort uniqid: replace pmort30_fast_ll=pmort30_fast_ll[_n+1]  
bysort uniqid: replace tta_fast_ll=tta_fast_ll[_n+1]  

*-------------------------------------------------------------------------------
* Predict adverse abx-associated events under baseline and abx-treated scenarios
*-------------------------------------------------------------------------------
* For the duplicates, recode abx_in_48hr to 1 so that everyone is considered 
* abx treated 
replace abx_in_48hr=1 if dupindicator==1

sort uniqid dupindicator

* estimate abx associated outcomes on original observations 
foreach var in allergy ttp leukopenia acute_renal_injury acute_liver_injury ///
			cdiff any_mdro_or_escr any_mdro_blood_or_escr mort90_ed ///
			mort30_ed any_adverse_nomort any_adverse_wmort {		

	
	* unadjusted rates 
	tab subpop `var' if dupindicator==0 , ro
	
	
	* adjusted rates
	local aod 																		///
			 aod_lactate aod_kidney pressor_in_72hr aod_liver		///
			 aod_heme aod_lung 			
			
	local sirs																		///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

	local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_circ pvd paralysis 	///
			pud hypothyroid ah lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def etoh drug psychoses depression 
			
	logit `var' i.abx_in_48hr c.age male `sirs' `aod' `comorbid' va if dupindicator==0 
	margins abx_in_48hr
	margins abx_in_48hr if dupindicator==0 

	predict p`var'_orig
	
	* create new variable for p-hat treated 
	gen p`var'_treat = p`var'_orig if dupindicator==1 
	
	* recode the orginal observation so that it fills in missing values for 
	* pvar_treated 
	bysort uniqid: replace p`var'_treat=p`var'_treat[_n+1]
}	


* drop the duplicate indicator 
drop if dupindicator==1
drop dupindicator 

* rename variables 
rename pacute_renal_injury_* prenal_*
rename pacute_liver_injury_* pliver_*
rename pany_mdro_or_escr_* pmdro_*
rename pany_mdro_blood_or_escr_* pmdroblood_*
rename pmort90_ed_* pmort90_*
     

* save dataset 
save kp_va_happi_aim3_lowerCI, replace

log close