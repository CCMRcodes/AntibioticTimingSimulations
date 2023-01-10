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

log using "Aim3_analysis_step1_`time_string'.log", replace


* ----------------------------------------------------------------------
*
*  	Tables 1 and 2; Preparing Dataset for Simulations
*
* 	Created: 		2022 Apr 27
* 	Last updated: 	2022 Dec 8
* 	Author: 		S. Seelye
*
*-----------------------------------------------------------------------


* open dataset 

use Data\aim1_analytic_dataset_20211102, clear		

* count
count 
	* 1,560,126
	
* drop sample variable
drop sample 

* rename outcome variables 
rename allergy_abx_ind allergy 
rename thrombocytopenia_hosp ttp
rename leukopenia_hosp leukopenia 
rename cdiff_infection cdiff 

* check counts 
tab allergy admityear if va==1
tab allergy admityear if va==0

tab ttp admityear if va==1
tab ttp admityear if va==0

tab leukopenia admityear if va==1
tab leukopenia admityear if va==0

tab acute_renal_injury admityear if va==1
tab acute_renal_injury admityear if va==0

tab acute_liver_injury admityear if va==1
tab acute_liver_injury admityear if va==0

tab cdiff admityear if va==1
tab cdiff admityear if va==0

* examine outcomes by subpopulation
gen subpop = .
replace subpop = 1 if never_treated
replace subpop = 2 if early_abx_stop
replace subpop = 3 if infection_cohort
replace subpop = 4 if severe_sepsis
replace subpop = 5 if septic_shock 

lab def subpop 1 "never treated" 2 "early abx stop" 3 "infection" 4 "severe sepsis" 5 "septic shock"
lab val subpop subpop

tab subpop never_treated, m
tab subpop early_abx_stop, m
tab subpop infection_cohort, m
tab subpop severe_sepsis, m
tab subpop septic_shock, m

tab allergy subpop, chi co
tab ttp subpop, chi co
tab leukopenia subpop, chi co
tab acute_renal_injury subpop, chi co
tab acute_liver_injury subpop, chi co
tab cdiff subpop, chi co
tab any_mdro_or_escr subpop, chi co
tab any_mdro_blood_or_escr subpop, chi co
tab mort90_ed subpop, chi co


*-----------
* Table 1
*-----------

tab subpop 

univar age
version 16: table subpop, c(n age median age p25 age p75 age)

tab subpop sirs_wbc, ro
tab subpop sirs_temp, ro 
tab subpop sirs_pulse, ro
tab subpop sirs_rr, ro

tab subpop aod_lung, ro
tab subpop aod_kidney, ro
tab subpop aod_liver, ro
tab subpop aod_heme, ro
tab subpop aod_lactate, ro
tab subpop pressor_in_72hr, ro

tab subpop renal, ro
tab subpop liver, ro
tab subpop chf, ro
tab subpop pulm, ro 
tab subpop htn, ro 
tab subpop obesity, ro 
tab subpop pvd, ro
tab subpop cardic_arrhym, ro


*------------
* Table 2
*------------

* estimate abx associated outcomes on original observations 
foreach var in allergy ttp leukopenia acute_renal_injury acute_liver_injury ///
			cdiff any_mdro_or_escr any_mdro_blood_or_escr mort90_ed mort30_ed {		

	
	* unadjusted rates 
	describe `var'
	tab subpop `var' , ro
	
	
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
			
	logit `var' i.abx_in_48hr c.age male `sirs' `aod' `comorbid' va , or
	margins abx_in_48hr

}	


	
	
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
* Prepare data for  simulations
*-------------------------------------------------------------------------------

*-------------------------------------------------------------------------------
* First, predict mortality under baseline and 50% faster TTA
*-------------------------------------------------------------------------------

* expand the dataset with a new observation for each patient
expand 2 , gen(dupindicator)

* for the duplicates, reduce TTA by 50%
replace time_to_abx_min = time_to_abx_min*0.5 if dupindicator==1

sort uniqid dupindicator
	*br uniqid subpop dupindicator time_to_abx_min


* estimate mortality on original observations for severe sepsis / septic shock cohort
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

* predict mortality for all observations, including duplicates (i.e. 50% faster)
predict pmort30_orig, pr

* create new variable for p-hat faster 
gen pmort30_fast=pmort30_orig if dupindicator==1

* create new variable for TTA 50% faster 
gen tta_fast = time_to_abx_min if dupindicator==1

* recode the original observation so that it fills in missing values for 
* pmort30_fast and tta_fast
bysort uniqid: replace pmort30_fast=pmort30_fast[_n+1]  
bysort uniqid: replace tta_fast=tta_fast[_n+1]  


*-------------------------------------------------------------------------------
* Next, predict adverse abx-associated events under baseline and abx-treated 
* scenarios
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


* save dataset for simulations to run in R
*save kp_va_happi_aim3_20221103, replace

log close