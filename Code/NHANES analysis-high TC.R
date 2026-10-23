################################################################################
# Program purpose: Analyzing NHANES data to get US prevalence for high total cholesterol for midlife and late-life
# Input data: 
# Demographic  data:
# "Y:\Research Projects\CARE\6. Cohorts and datasets\2. Datasets\Working datasets\NHANES 2021-23\Demographic\DEMO_L.xpt"
# Cholesterol data:
# "Y:\Research Projects\CARE\6. Cohorts and datasets\2. Datasets\Working datasets\NHANES 2021-23\laboratory data\TCHOL_L.xpt"
# Output files: 
# Date started: 13 Oct 2026
# Author: Weisi Chen
# Resources:
# Source report: Prevalence of total cholesterol in US (CDC): https://www.cdc.gov/nchs/data/databriefs/db515.pdf
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

# Import and merge NHANES files
demo <- read_xpt("Demographic/DEMO_L.xpt")
tchol <- read_xpt("laboratory data/TCHOL_L.xpt")

# 1. Merge cholesterol data with the demographic data
dat <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR,
         SDMVSTRA, SDMVPSU) %>%
  left_join(
    tchol %>% select(SEQN, LBXTC, WTPH2YR),
    by = "SEQN"
  ) %>%
  mutate(
    
    # CARE cut-off
    tc200 = if_else(
      is.na(LBXTC), NA_real_,
      as.numeric(LBXTC >= 200)
    ),
    
    # CDC cut-off
    tc240 = if_else(
      is.na(LBXTC), NA_real_,
      as.numeric(LBXTC >= 240)
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
    ),
    
  )


# 2. NHANES survey design
design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTPH2YR,
  data = filter(dat, !is.na(WTPH2YR),
                WTPH2YR > 0),
  nest = TRUE
)

# Adult population with observed cholesterol (20+)
adults <- subset(
  design,
  RIDAGEYR >= 20 & !is.na(LBXTC)
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


##############
# Function to calculate crude and age-adjusted prevalence (%) with Korn-Graubard 95% confidence intervals 
adjusted_prev_ci_korn <- function(
    survey_design,
    variable,
    age_adjusted = TRUE
) {
  
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
    adj = "NCHS"
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


# Calculate prevalence for age 20+ males and females
# Crude prevalence
t(sapply(
  groups,
  adjusted_prev_ci_korn,
  variable = "tc240",
  age_adjusted = FALSE
))


# Age adjusted prevalence
t(sapply(
  groups,
  adjusted_prev_ci_korn,
  variable = "tc240",
  age_adjusted = TRUE
))

# Age 40-59 males and females
# Crude prevalence
t(sapply(
  lapply(groups, function(d) subset(d, age_CDC == "40-59")),
  adjusted_prev_ci_korn,
  variable = "tc240",
  age_adjusted = FALSE
))


# Age 60+ males and females
# Crude prevalence
t(sapply(
  lapply(groups, function(d) subset(d, age_CDC == "60+")),
  adjusted_prev_ci_korn,
  variable = "tc240",
  age_adjusted = FALSE
))

# Cross-check with the ref: https://www.cdc.gov/nchs/data/databriefs/db515.pdf (Data table for Figure 1)

