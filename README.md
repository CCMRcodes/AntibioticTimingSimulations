# AntibioticTimingSimulations

Overview:

This repository contains code needed to perform the simulations and data analysis described in the paper "Impact of reducing time-to-antibiotics on sepsis mortality, excess antibiotic use, and adverse events: a cohort and simulation study." The paper used simulations to estimate the benefits and harms of shortening time-to-antibiotic receipt for sepsis.

The code is structured in two steps:

* Step 1: abxtimingsim_step1.do - dataset preparation for step 2 simulations; tables 1-2
* Step 2: abxtimingsim_step2.R - simulations; output of results for tables 3-5 (R script)

Sensitivity analyses:

Sensitivity analyses were conducted to estimate a range of sepsis deaths averted and newly antibiotic-treated patients and adverse events per sepsis death averted. The code follows the same structure as above, separately for the lower and upper bounds of the time-to-antibiotics estimate.

Lower Bounds
* Step 1: abxtimingsim_step1_lower.do 
* Step 2: abxtimingsim_step2_lower.R

Upper Bounds
* Step 1: abxtimingsim_step1_upper.do 
* Step 2: abxtimingsim_step2_upper.R 
