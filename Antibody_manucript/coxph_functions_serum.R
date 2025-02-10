# Cox regression under a case-cohort sampling design
library(CaseCohortCoxSurvival)
library(survival)
library(dplyr)

# The first case-cohort Cox regression methods
# based on the R pacakge CaseCohortCoxSurvival.
# The package assumes the phase 2 samples, e.g.,
# subcohort | case have all their phase 2 covariates
# measured, and if not, then the missingness is 
# at random.
case_cohort_Cox_1 <- function(data, caseInd, subcohortInd, baseline_strata,
                              time_to_event, phase1_var_list, phase2_var_list){
  md <- caseCohortCoxSurvival(data = data, status = caseInd, 
                              subcohort = subcohortInd, strata = baseline_strata, 
                              time = time_to_event, 
                              cox.phase1 = phase1_var_list, 
                              cox.phase2 = phase2_var_list,
                              print = 0)
  return(md)
}

# Function processing the return from the cccs package
process_cccs_output <- function(md){
  tbl_output = broom::tidy(md$coxph.fit, exp = TRUE, conf.int = TRUE)
  #CI = as.data.frame(exp(confint(md$coxph.fit)))
  
  tbl_output = tbl_output %>%
    dplyr::select(-robust.se, -statistic, -std.error)
  #tbl_output = cbind(tbl_output, CI)
  return(tbl_output)
}

# Use the survey design package + survival package
# W is the straitification variable: baseline strata | case
case_cohort_Cox_2 <- function(data, adj_vars, mk){
  
  # Prepare the design weights
  CaseCohortDesign<-twophase(id=list(~Subjectid,~Subjectid), 
                             strata=list(NULL,~W),
                             subset=~ph2.m1.sera_logical, data = data)
  
  # Prepare the formula
  RHS = paste(c(adj_vars, mk), collapse = '+')
  f = as.formula(paste0('Surv(EventTimePrimary, EventIndPrimary1) ~', RHS))
  md = svycoxph(formula = f, design = CaseCohortDesign, data = data)
  
  return(md)
}