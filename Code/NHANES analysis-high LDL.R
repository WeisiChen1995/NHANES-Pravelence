################################################################################
# Program purpose: Analyzing NHANES data to get US prevalence for high LDL-c for midlife and late-life
# Input data: 
# Demographic  data: "DEMO_L.xpt"
# LDL-C data: "TRIGLY_L.xpt"
# Output files: 
# Date started: 23 Oct 2026
# Author: Weisi Chen
# Resources:
# Source report (2015-18 NHANES): https://www.cdc.gov/nchs/data/nhsr/nhsr220.pdf
# * Source report not yet available for LDL-C from 2021-23 NHANES
# Population age group weights (2000): https://wwwn.cdc.gov/nchs/nhanes/tutorials/samplecode.aspx
# https://www.cdc.gov/nchs/data/statnt/statnt20.pdf
################################################################################
################################################################################
################################################################################
# Load libraries
library(haven)
library(dplyr)
library(survey)
library(surveytable)
library(tidyr)

# Input data
demo <- read_xpt("Y:/Research Projects/CARE/6. Cohorts and datasets/2. Datasets/Working datasets/NHANES 2021-23/Demographic/DEMO_L.xpt")
ldl <- read_xpt("Y:/Research Projects/CARE/6. Cohorts and datasets/2. Datasets/Working datasets/NHANES 2021-23/laboratory data/TRIGLY_L.xpt")

# 1. Merge LDL-cholesterol data with the demographic data
dat <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR, SDMVSTRA, SDMVPSU) %>%
  left_join(
    ldl %>% select(SEQN, WTSAF2YR, LBXTLG, LBDLDL, LBDLDLM, LBDLDLN),
    by = "SEQN"
  ) %>%
  mutate(
    # CARE binary LDL cut-offs
    # LDL Friedewald >= 160mg/gL
    ldl160_friedewald = if_else(
      is.na(LBDLDL), NA_real_,
      as.numeric(LBDLDL >= 160)
    ),
    
    # LDL Martin-Hopkins >= 160mg/gL
    ldl160_martin = if_else(
      is.na(LBDLDLM), NA_real_,
      as.numeric(LBDLDLM >= 160)
    ),
    
    # LDL Sampson-NIH equation 2 >= 160mg/gL
    ldl160_nih = if_else(
      is.na(LBDLDLN), NA_real_,
      as.numeric(LBDLDLN >= 160)
    ),
    
    # CDC source report cut-offs - 5 categories (for comparison): <70, 70-<100, 100-<160, 160-<190, >=190
    # LDL Friedewald
    ldl_fried_cat = cut(
      LBDLDL, breaks = c(-Inf, 70, 100, 160, 190, Inf),
      labels = c("<70", "70-<100", "100-<160", "160->190", ">=190"), 
      right = FALSE
    ),
    
    # LDL Martin-Hopkins
    ldl_martin_cat = cut(
      LBDLDLM, breaks = c(-Inf, 70, 100, 160, 190, Inf),
      labels = c("<70", "70-<100", "100-<160", "160->190", ">=190"), 
      right = FALSE
    ),
    
    # LDL Sampson-NIH equation 2
    ldl_nih_cat = cut(
      LBDLDLN, breaks = c(-Inf, 70, 100, 160, 190, Inf),
      labels = c("<70", "70-<100", "100-<160", "160->190", ">=190"), 
      right = FALSE
    ),
    
    
    # Official CDC age-standardization categories
    age_CDC = factor(
      case_when(
        RIDAGEYR >= 20 & RIDAGEYR < 40 ~ "20-39",
        RIDAGEYR >= 40 & RIDAGEYR < 60 ~ "40-59",
        RIDAGEYR >= 60 ~ "60+",
        TRUE ~ NA_character_
      ),
      levels = c("20-39", "40-59", "60+")
    )
    
  )


# 2. NHANES survey design
design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTSAF2YR,
  data = filter(dat, !is.na(WTSAF2YR),
                WTSAF2YR > 0),
  nest = TRUE
)

# Adult population with observed LDL-cholesterol (20+)
adults <- subset(
  design,
  RIDAGEYR >= 20 & !is.na(LBDLDL)
)

# 3. Age-standardization weights:
# 2000 US standard population (published): 
# Source: https://wwwn.cdc.gov/nchs/nhanes/tutorials/samplecode.aspx
# https://www.cdc.gov/nchs/data/statnt/statnt20.pdf
standard_pop_CDC <- c(
  0.396579,  # age group: 20-39
  0.371795,  # age group: 40-59
  0.231626   # age group: 60+
)

standard_pop_CARE <- c(
  0.396579,  # age group: 20-39
  0.426194,  # age group: 40-64, calculated by 83470/195850
  0.177227   # age group: 65+, calculated by 34710/195850
)


##############
# Function to calculate crude and age-adjusted prevalence (%) with Korn-Graubard 95% confidence intervals 
# Edit so that it works for both binary and categorical variables
adjusted_prev_ci_korn <- function(
    survey_design,
    variable,
    age_adjusted = TRUE
) {
  
  x <- survey_design$variables[[variable]]
  
  # For categorical variables, calculate each category separately
  if (is.factor(x)) {
    return(t(sapply(levels(x), function(cat) {
      d <- survey_design
      d$variables$.indicator <- as.numeric(x == cat)
      adjusted_prev_ci_korn(d, ".indicator", age_adjusted)
    })))
  }
  
  # Exclude missing outcomes
  survey_design <- survey_design[!is.na(x), ]
  
  formula_var <- as.formula(paste0("~", variable))
  
  if (age_adjusted) {
    
    analysis_design <- svystandardize(
      survey_design,
      by = ~age_CDC,
      over = ~1,
      population = standard_pop_CDC
    )
    
  } else {
    
    analysis_design <- survey_design
    
  }
  
  result <- surveytable::svyciprop_adjusted(
    formula = formula_var,
    design = analysis_design,
    level = 0.95,
    adj = "NCHS",
    aa_vr = if (age_adjusted) "age_CDC" else NULL,
    aa_pop = if (age_adjusted) data.frame(
      Level = levels(survey_design$variables$age_CDC),
      Population = standard_pop_CDC
    ) else NULL
  )
  
  
  ci <- attr(result, "ci")
  
  c(
    Prevalence = as.numeric(result)[1] * 100,
    KG_Lower = as.numeric(ci[1]) * 100,
    KG_Upper = as.numeric(ci[2]) * 100
  )
}

# 4. Compare against CDC published estimates
groups <- list(
  Overall = adults,
  Male = subset(adults, RIAGENDR == 1),
  Female = subset(adults, RIAGENDR == 2)
)

################
# Calculate prevalence for age 20+ males and females
# Crude prevalence
# LDL-C Friedewald
lapply(
  groups,
  adjusted_prev_ci_korn,
  variable = "ldl_fried_cat",
  age_adjusted = FALSE
)

# LDL-C Martin-Hopkins
lapply(
  groups,
  adjusted_prev_ci_korn,
  variable = "ldl_martin_cat",
  age_adjusted = FALSE
)

# LDL-C Sampson-NIH equation 2
lapply(
  groups,
  adjusted_prev_ci_korn,
  variable = "ldl_nih_cat",
  age_adjusted = FALSE
)


################
# Age 20-39 males and females

# LDL-C Friedewald
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "20-39")),
  adjusted_prev_ci_korn,
  variable = "ldl_fried_cat",
  age_adjusted = FALSE
)

# LDL-C Martin-Hopkins
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "20-39")),
  adjusted_prev_ci_korn,
  variable = "ldl_martin_cat",
  age_adjusted = FALSE
)

# LDL-C Sampson-NIH equation 2
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "20-39")),
  adjusted_prev_ci_korn,
  variable = "ldl_nih_cat",
  age_adjusted = FALSE
)

################
# Age 40-59 males and females

# LDL-C Friedewald
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "40-59")),
  adjusted_prev_ci_korn,
  variable = "ldl_fried_cat",
  age_adjusted = FALSE
)

# LDL-C Martin-Hopkins
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "40-59")),
  adjusted_prev_ci_korn,
  variable = "ldl_martin_cat",
  age_adjusted = FALSE
)

# LDL-C Sampson-NIH equation 2
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "40-59")),
  adjusted_prev_ci_korn,
  variable = "ldl_nih_cat",
  age_adjusted = FALSE
)



################
# Age 60+ males and females

# LDL-C Friedewald
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "60+")),
  adjusted_prev_ci_korn,
  variable = "ldl_fried_cat",
  age_adjusted = FALSE
)

# LDL-C Martin-Hopkins
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "60+")),
  adjusted_prev_ci_korn,
  variable = "ldl_martin_cat",
  age_adjusted = FALSE
)

# LDL-C Sampson-NIH equation 2
lapply(
  lapply(groups, function(d) subset(d, age_CDC == "60+")),
  adjusted_prev_ci_korn,
  variable = "ldl_nih_cat",
  age_adjusted = FALSE
)


# Cross-check with the ref: https://www.cdc.gov/nchs/data/nhsr/nhsr220.pdf
# Note Table 1 and Table 2 reported LDL prevalence using NHANES 2015-18 not NHANES 2021-23, so some differences are expected

