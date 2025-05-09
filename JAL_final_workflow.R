
# =====================================================
# Judgements of Loyalty – A Bayesian Workflow
# 
# Final models after prior sensitivity diagnostics
# =====================================================

########### Load packages ###########

library(brms)
library(tidyverse)
library(tidybayes)
library(posterior)
library(bayesplot)
library(marginaleffects)
library(priorsense)
theme_set(theme_minimal())

# HPDI utility for 89% intervals (used by tidybayes)
hpdi_level <- 0.89

#Coding of the data
# A = gender of the assessor / judge (0 = Man, 1 = Women)
# G = gender of the candidate (0 = Man, 1 = Women) 
# C = target country
# CI = Sweden
# CJ = Enlgand
# CK = Germany
# CL = Brazil
# CM = Thailand
# T = Type of relationship
# TA = Parents from target country
# TB = Married to a person from target country
# TC = Born in target country

#Q1 = I vilken utsträckning instämmer du i påståendet att personen kan antas 
#vara lojal mot de intressen som Säkerhetsskyddslag (2018:585) ska skydda på en 
#skala från 1 till 7

#Q1 in Enlgish
#To what extent do you agree with the statement that the person can be assumed 
#to be loyal to the interests that the Security Protection Act (2018:585) aims 
#to protect, on a scale from 1 to 7? 

#Q2 = I vilken utsträckning instämmer du i påståendet att där föreligger en 
#risk för en lojalitetskonflikt på en skala från 1 till 7?

#Q2 in Enlgish
#To what extent do you agree with the statement that there is a risk of a 
#conflict of loyalty, on a scale from 1 to 7?

#Q3 = I vilken utsträckning instämmer du i påståendet att personen kan antas 
#värdera lojaliteten till Sverige högre än lojaliteten mot annan part om 
#personen tvingas välja mellan Sverige och den andra parten på en skala från 
#1 till 7?

#Q3 in Enlish
#To what extent do you agree with the statement that the person can be assumed 
#to value loyalty to Sweden higher than loyalty to another party if forced to 
#choose between Sweden and the other party, on a scale from 1 to 7?

#The response scale
#1.	Do not agree at all
#2.	Do not agree
#3.	Partially do not agree
#4.	Neither agree nor disagree
#5.	Partially agree
#6.	Agree
#7.	Fully agree

# DATA LOADING AND PREPROCESSING
library(readxl)
IRM4R <- read_excel("IRM4R.xlsx")
View(IRM4R)

d <- reshape2::melt(IRM4R, id.vars=c("ID", "A", "G", "C", "T", "D"), measure.vars=c("L1", "L2", "L3" ),
                    variable.name="item",value.name="rating")

d <- d[ order(d$ID, d$C, d$item), ]

d$C <- as.factor(d$C)
d$T <- as.factor(d$T)
d$G <- as.factor(d$G)
d$A <- as.factor(d$A)

view(d)

nCores = parallel::detectCores()

########### DESCRIPTIVES ###########

#Descriptives
with(IRM4R, tapply(L1, list(C, T), median))
#     TA TB TC
# CI 6.5  7  6
# CJ 6.0  6  6
# CK 6.0  6  5
# CL 6.0  6  5
# CM 6.0  6  5

with(IRM4R, tapply(L2, list(C, T), median))
#     TA TB TC
# CI 6.0  6  6
# CJ 6.0  5  3
# CK 5.5  5  4
# CL 4.5  3  3
# CM 5.5  3  3

with(IRM4R, tapply(L3, list(C, T), median))
#     TA TB TC
# CI 7.0  7  7
# CJ 6.0  5  5
# CK 5.5  6  4
# CL 5.5  5  3
# CM 5.5  5  5

with(IRM4R, tapply(D == 1, list(C, T), sum)) #Counts of Yes given condition
#    TA TB TC
# CI  5  4  2
# CJ  6 11  9
# CK  6 12 10
# CL  6  9 14
# CM  7 10 12

#Number and percentage of yes in experimental conditions

# Ensure D is numeric (0 = no, 1 = yes)
d$D <- as.numeric(d$D)

# One row per judgment: each ID's rating of a country
d_judgment <- d[!duplicated(d[c("ID", "C")]), ]

# Exclude control condition
d_exp <- subset(d_judgment, C != "CI")

# Count "yes" responses
n_yes <- sum(d_exp$D == 1, na.rm = TRUE)
n_total <- sum(!is.na(d_exp$D))

# Compute percentage
yes_percent <- round(100 * n_yes / n_total, 1)

# Output
cat("Yes responses:", n_yes, "\n")
cat("Total experimental observations:", n_total, "\n")
cat("Percentage of 'yes' responses:", yes_percent, "%\n")

########### MODEL FITTING ###########

d$D <- as.integer(d$D)

###### 3.1 Binary probit model (yes/no judgments) ######
prob_prior_scaled <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 2), class = "b"),
  prior(exponential(0.5), class = "sd"),
  prior(lkj(2), class = "cor")
)

fitProb_final_scaled <- brm(
  formula = D ~ 1 + C + G + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 4000, iter = 16000,
  chains = 4, cores = parallel::detectCores(),
  seed = 5476,
  control = list(adapt_delta = 0.99),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fitProb_final_scaled"
)

#Trace and density plots
plot(fitProb_final_scaled) #Looks healthy

#Posterior predictive check
pp_check(fitProb_final_scaled) #Looks good

fitProb_final_scaled = add_criterion(fitProb_final_scaled, criterion = "loo")
# Automatically saving the model object in 'fitProb_final_scaled.rds'
# Warning message:
# Found 27 observations with a pareto_k > 0.7 in model 'fitProb_final_scaled'. 
# We recommend to set 'moment_match = TRUE' in order to perform moment matching 
# for problematic observations. 

loo_fitProb <- loo(fitProb_final_scaled)
loo_fitProb_mm <- loo_moment_match(fitProb_final_scaled, loo = loo_fitProb)

loo_fitProb_mm
# Computed from 48000 by 870 log-likelihood matrix.
# 
#          Estimate  SE
# elpd_loo    -21.6 0.8
# p_loo         7.0 0.3
# looic        43.1 1.5
# 
# MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.2]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

summary(fitProb_final_scaled)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          5.91      1.39     3.62     9.04 1.00    19477    32602
# sd(CCJ)                8.05      3.03     3.61    15.30 1.00    41686    34426
# sd(CCK)                5.68      3.30     0.45    13.45 1.00     8004     5969
# sd(CCL)               11.58      3.53     6.19    19.89 1.00    46175    35058
# sd(CCM)               12.71      3.53     7.20    20.89 1.00    56809    37168
# cor(Intercept,CCJ)     0.31      0.23    -0.19     0.70 1.00    36354    32834
# cor(Intercept,CCK)     0.35      0.28    -0.29     0.79 1.00    34276    29173
# cor(CCJ,CCK)           0.48      0.23    -0.09     0.80 1.00    14709    14769
# cor(Intercept,CCL)     0.02      0.21    -0.40     0.41 1.00    25888    31972
# cor(CCJ,CCL)           0.35      0.18    -0.03     0.67 1.00    18025    30471
# cor(CCK,CCL)           0.21      0.23    -0.28     0.61 1.00     9511    11638
# cor(Intercept,CCM)    -0.14      0.19    -0.51     0.24 1.00    27727    36077
# cor(CCJ,CCM)           0.43      0.17     0.06     0.74 1.00    21331    31520
# cor(CCK,CCM)           0.31      0.22    -0.18     0.70 1.00    11591    10072
# cor(CCL,CCM)           0.57      0.13     0.28     0.79 1.00    34101    40239
# 
# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -4.09      1.09    -6.36    -2.09 1.00    21943    34500
# CCJ           1.91      1.38    -1.01     4.45 1.00    31945    38471
# CCK           3.40      1.38     0.48     5.88 1.00    13868    23798
# CCL           2.30      1.48    -0.73     5.07 1.00    43885    39578
# CCM           1.99      1.51    -1.10     4.82 1.00    46805    39279
# G1           -0.82      0.88    -2.59     0.90 1.00    35390    37016
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).



###### 3.2 Ordinal model with discrimination ######
fit_dt_priors_adjusted <- c(
  # Threshold priors (per item)
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L1),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L1),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L1),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L1),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L1),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L1),
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L2),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L2),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L2),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L2),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L2),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L2),
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L3),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L3),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L3),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L3),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L3),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L3),
  # Effects and disc
  prior(normal(0, 2), class = b),
  prior(normal(0, 1), class = b, dpar = disc),
  prior(exponential(0.5), class = sd),
  prior(lkj(2), class = cor)
)

fit_dt_reduced_adjusted_thresholds <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + (1 + C | ID) + (1 | item)
  ) +
    lf(disc ~ 0 + C + G + D + (1 | ID), cmc = FALSE),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_adjusted,
  warmup = 4000, iter = 16000,
  chains = 4, cores = parallel::detectCores(),
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = "fit_dt_reduced_adjusted_thresholds"
)

#Trace and density plots
plot(fit_dt_reduced_adjusted_thresholds) #Looks healthy

#Posterior predictive check
pp_check(fit_dt_reduced_adjusted_thresholds) #Looks good

fit_dt_reduced_adjusted_thresholds = add_criterion(fit_dt_reduced_adjusted_thresholds, criterion = "loo")
# Automatically saving the model object in 'fit_dt_reduced_adjusted_thresholds.rds'
# Warning message:
# Found 31 observations with a pareto_k > 0.7 in model 
# 'fit_dt_reduced_adjusted_thresholds'. We recommend to set 
# 'moment_match = TRUE' in order to perform moment matching for problematic 
# observations. 

loo_fitdt <- loo(fit_dt_reduced_adjusted_thresholds)
loo_fitdt_mm <- loo_moment_match(fit_dt_reduced_adjusted_thresholds, loo = loo_fitdt)
loo_fitdt_mm

# Computed from 48000 by 870 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo  -1085.2 28.3
# p_loo       170.7  9.7
# looic      2170.4 56.5
# 
# MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.2, 1.5]).
# 
# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     852   97.9%   988     
# (0.7, 1]      (bad)       18    2.1%   <NA>    
# (1, Inf)      (very bad)   0    0.0%   <NA>    
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_reduced_adjusted_thresholds)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + C + G + D + (1 + C | ID) + (1 | item) 
# disc ~ 0 + C + G + D + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.47      0.08     0.32     0.65 1.00     8116    15986
# sd(CCJ)                0.15      0.10     0.01     0.37 1.00     5205     8369
# sd(CCK)                0.21      0.10     0.03     0.43 1.00     4387     6801
# sd(CCL)                0.31      0.10     0.12     0.51 1.00     4850     9435
# sd(CCM)                0.11      0.08     0.00     0.30 1.00     4341     5286
# sd(disc_Intercept)     0.65      0.10     0.49     0.86 1.00    10894    21052
# cor(Intercept,CCJ)     0.20      0.30    -0.44     0.73 1.00    28145    25973
# cor(Intercept,CCK)     0.32      0.26    -0.23     0.76 1.00    24258    27736
# cor(CCJ,CCK)           0.25      0.34    -0.50     0.80 1.00     8187    15341
# cor(Intercept,CCL)     0.04      0.23    -0.38     0.50 1.00    17032    20486
# cor(CCJ,CCL)          -0.07      0.32    -0.68     0.56 1.00     5551    11686
# cor(CCK,CCL)           0.25      0.29    -0.38     0.73 1.00     7706    16990
# cor(Intercept,CCM)    -0.02      0.33    -0.62     0.62 1.00    32354    33474
# cor(CCJ,CCM)           0.10      0.36    -0.61     0.74 1.00    10833    15572
# cor(CCK,CCM)           0.16      0.36    -0.58     0.78 1.00     8405    12292
# cor(CCL,CCM)           0.20      0.34    -0.52     0.77 1.00     9609    21152
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.07      1.01     0.84     4.67 1.00    14548    24045
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.12      0.67    -3.55    -0.92 1.00    44956    34164
# Intercept[L1,2]    -1.41      0.53    -2.49    -0.40 1.00    43856    38918
# Intercept[L1,3]    -0.30      0.45    -1.17     0.61 1.00    54053    37336
# Intercept[L1,4]     0.18      0.45    -0.68     1.07 1.00    62276    37758
# Intercept[L1,5]     0.89      0.46     0.01     1.79 1.00    45247    34706
# Intercept[L1,6]     1.96      0.50     1.00     2.95 1.00    19228    28326
# Intercept[L2,1]    -1.35      0.44    -2.21    -0.50 1.00    22884    29408
# Intercept[L2,2]    -0.70      0.41    -1.50     0.10 1.00    42497    35407
# Intercept[L2,3]    -0.16      0.40    -0.95     0.63 1.00    58727    37222
# Intercept[L2,4]    -0.04      0.40    -0.83     0.76 1.00    59672    37053
# Intercept[L2,5]     0.31      0.41    -0.48     1.11 1.00    54895    36826
# Intercept[L2,6]     1.47      0.44     0.60     2.35 1.00    23390    31962
# Intercept[L3,1]    -1.30      0.44    -2.18    -0.44 1.00    29791    32375
# Intercept[L3,2]    -0.75      0.41    -1.57     0.06 1.00    48124    35101
# Intercept[L3,3]    -0.34      0.41    -1.14     0.45 1.00    61841    36823
# Intercept[L3,4]    -0.07      0.41    -0.87     0.72 1.00    65173    36748
# Intercept[L3,5]     0.52      0.41    -0.29     1.34 1.00    51977    36385
# Intercept[L3,6]     1.40      0.44     0.54     2.28 1.00    23932    31809
# CCJ                -0.80      0.14    -1.08    -0.55 1.00     7947    14180
# CCK                -0.77      0.14    -1.06    -0.52 1.00     8169    14845
# CCL                -0.99      0.16    -1.32    -0.71 1.00     7214    13117
# CCM                -0.93      0.15    -1.24    -0.67 1.00     6948    12640
# G1                  0.03      0.05    -0.06     0.13 1.00    33896    31187
# D                  -0.49      0.09    -0.69    -0.34 1.00     9217    18035
# disc_CCJ            0.82      0.13     0.58     1.07 1.00    16395    29334
# disc_CCK            0.92      0.13     0.67     1.18 1.00    17106    28208
# disc_CCL            0.83      0.13     0.57     1.09 1.00    16805    29258
# disc_CCM            1.08      0.13     0.83     1.34 1.00    16045    28187
# disc_G1            -0.03      0.09    -0.20     0.14 1.00    30061    34094
# disc_D              0.12      0.12    -0.11     0.35 1.00    20625    29581
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

###### Exploratory subset models for participant gender ######
fit_dt_priorsAG_scaled <- c(
  # Threshold priors
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L1),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L1),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L1),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L1),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L1),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L1),
  
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L2),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L2),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L2),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L2),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L2),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L2),
  
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L3),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L3),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L3),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L3),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L3),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L3),
  
  # Main effects and interaction (we can keep this prior as-is or slightly relaxed if needed)
  prior(normal(0, 1), class = b),
  
  # Revised discrimination prior (weaker, consistent with power-scaling)
  prior(normal(0, 1), class = b, dpar = disc),
  
  # Random effect SD prior
  prior(exponential(1), class = sd)
)

fit_dt_AG_scaled <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + A*G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + A*G + (1 | ID),
      cmc = FALSE  
    ),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priorsAG_scaled,
  warmup = 4000, iter = 16000,
  chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = 'fit_dt_AG_scaled')

#Trace and density plots
plot(fit_dt_AG_scaled) #Looks healthy

#Posterior predictive check
pp_check(fit_dt_AG_scaled) #Looks good

fit_dt_AG_scaled = add_criterion(fit_dt_AG_scaled, criterion = "loo")
# Warning message:
# Found 11 observations with a pareto_k > 0.7 in model 'fit_dt_AG_scaled'. 
# We recommend to set 'moment_match = TRUE' in order to perform moment matching 
# for problematic observations.  

loo_fitdt_AG <- loo(fit_dt_AG_scaled)
loo_fitdt_AG_mm <- loo_moment_match(fit_dt_AG_scaled, loo = loo_fitdt_AG)
loo_fitdt_AG_mm
# Computed from 48000 by 870 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo  -1288.2 26.3
# p_loo       116.3  8.9
# looic      2576.3 52.7
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.5]).
# 
# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     869   99.9%   821     
# (0.7, 1]      (bad)        0    0.0%   <NA>    
# (1, Inf)      (very bad)   1    0.1%   <NA>    
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_AG_scaled)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + A * G + (1 | ID) + (1 | item) 
# disc ~ 0 + A * G + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.76      0.12     0.55     1.01 1.00     9007    17008
# sd(disc_Intercept)     0.45      0.07     0.33     0.61 1.00    13171    23942
# 
# ~item (Number of levels: 3) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.00      0.59     0.15     2.47 1.00    12769     9622
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.78      0.61    -4.04    -1.64 1.00    32190    32354
# Intercept[L1,2]    -2.27      0.53    -3.32    -1.24 1.00    31652    31646
# Intercept[L1,3]    -0.96      0.49    -1.91     0.00 1.00    25877    21573
# Intercept[L1,4]    -0.01      0.49    -0.96     0.94 1.00    24714    16321
# Intercept[L1,5]     1.24      0.50     0.26     2.22 1.00    20473    13915
# Intercept[L1,6]     2.68      0.55     1.60     3.76 1.00    14762    14505
# Intercept[L2,1]    -1.96      0.43    -2.83    -1.12 1.00    17787    25496
# Intercept[L2,2]    -0.97      0.38    -1.74    -0.23 1.00    28439    33454
# Intercept[L2,3]    -0.03      0.37    -0.76     0.69 1.00    40590    36005
# Intercept[L2,4]     0.17      0.37    -0.56     0.89 1.00    41204    37356
# Intercept[L2,5]     0.68      0.37    -0.06     1.42 1.00    37302    36883
# Intercept[L2,6]     2.00      0.43     1.17     2.84 1.00    20942    32786
# Intercept[L3,1]    -1.83      0.44    -2.70    -0.98 1.00    21032    30935
# Intercept[L3,2]    -1.17      0.39    -1.94    -0.38 1.00    29319    32895
# Intercept[L3,3]    -0.50      0.38    -1.22     0.26 1.00    39571    35946
# Intercept[L3,4]     0.00      0.37    -0.71     0.76 1.00    43206    38022
# Intercept[L3,5]     0.93      0.39     0.19     1.72 1.00    33144    36923
# Intercept[L3,6]     1.91      0.43     1.11     2.78 1.00    20147    31406
# A1                  0.45      0.26    -0.03     0.97 1.00     8103    14604
# G1                  0.40      0.15     0.12     0.69 1.00    19761    26183
# A1:G1              -0.43      0.18    -0.80    -0.09 1.00    20270    26821
# disc_A1            -0.14      0.14    -0.40     0.13 1.00     8399    16121
# disc_G1             0.01      0.14    -0.26     0.28 1.00    17594    26590
# disc_A1:G1         -0.00      0.16    -0.32     0.32 1.00    18626    28313
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

d_women <- d[d$A == "1", ]

#Women participants
fit_dt_women_scaled <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + G + (1 | ID),
      cmc = FALSE
    ),
  data = d_women,
  family = cumulative(probit),
  prior = fit_dt_priorsAG_scaled,  
  warmup = 4000, iter = 16000, chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = "fit_dt_women_scaled"
)

#Trace and density plots
plot(fit_dt_women_scaled) #Looks healthy

#Posterior predictive check
pp_check(fit_dt_women_scaled) #Looks good

fit_dt_women_scaled = add_criterion(fit_dt_women_scaled, criterion = "loo")
# Warning message:
# Found 5 observations with a pareto_k > 0.7 in model 'fit_dt_women_scaled'. 
# We recommend to set 'moment_match = TRUE' in order to perform moment matching 
# for problematic observations.  

loo_fitdt_women <- loo(fit_dt_women_scaled)
loo_fitdt_women_mm <- loo_moment_match(fit_dt_women_scaled, loo = loo_fitdt_women)
loo_fitdt_women_mm
# Computed from 48000 by 630 log-likelihood matrix.
# 
# Estimate   SE
# elpd_loo   -956.9 22.5
# p_loo        88.3  7.9
# looic      1913.9 45.1
# ------
# MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.6]).
# 
# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     629   99.8%   636     
# (0.7, 1]      (bad)        0    0.0%   <NA>    
# (1, Inf)      (very bad)   1    0.2%   <NA>    
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_women_scaled)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item) 
# disc ~ 0 + G + (1 | ID)
# Data: d_women (Number of observations: 630) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 42) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.70      0.11     0.51     0.95 1.00    11736    22249
# sd(disc_Intercept)     0.40      0.08     0.26     0.57 1.00    13765    24634
# 
# ~item (Number of levels: 3) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.23      0.60     0.42     2.71 1.00    19074    20464
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.50      0.63    -3.82    -1.35 1.00    42373    34058
# Intercept[L1,2]    -1.92      0.52    -2.97    -0.91 1.00    50267    38201
# Intercept[L1,3]    -0.81      0.47    -1.74     0.12 1.00    47059    28787
# Intercept[L1,4]     0.02      0.47    -0.90     0.94 1.00    47660    26494
# Intercept[L1,5]     1.20      0.48     0.25     2.14 1.00    43496    25442
# Intercept[L1,6]     2.43      0.51     1.43     3.41 1.00    34306    24871
# Intercept[L2,1]    -1.70      0.42    -2.53    -0.86 1.00    29854    32565
# Intercept[L2,2]    -0.92      0.39    -1.69    -0.15 1.00    38058    36152
# Intercept[L2,3]    -0.11      0.38    -0.85     0.65 1.00    42859    35869
# Intercept[L2,4]     0.10      0.38    -0.64     0.86 1.00    43606    37012
# Intercept[L2,5]     0.57      0.39    -0.19     1.33 1.00    43271    37195
# Intercept[L2,6]     1.66      0.41     0.86     2.49 1.00    36869    35845
# Intercept[L3,1]    -1.61      0.44    -2.46    -0.75 1.00    32492    31293
# Intercept[L3,2]    -1.01      0.40    -1.80    -0.21 1.00    39270    33534
# Intercept[L3,3]    -0.52      0.40    -1.29     0.26 1.00    43423    35313
# Intercept[L3,4]    -0.04      0.39    -0.79     0.75 1.00    45366    35280
# Intercept[L3,5]     0.81      0.40     0.05     1.61 1.00    42374    35705
# Intercept[L3,6]     1.61      0.42     0.81     2.44 1.00    36799    33476
# G1                 -0.02      0.10    -0.21     0.17 1.00    41092    37620
# disc_G1            -0.01      0.09    -0.19     0.16 1.00    24194    31896
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


#Male participants

d_men <- d[d$A == "0", ]
# Same setup as fit_women, but on d_men

fit_dt_men_scaled <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + G + (1 | ID),
      cmc = FALSE
    ),
  data = d_men,
  family = cumulative(probit),
  prior = fit_dt_priorsAG_scaled,
  warmup = 4000, iter = 16000, chains = 4, cores = nCores,
  seed = 5476,
  init_r = 0.2,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  file = "fit_dt_men_scaled"
)

#Trace and density plots
plot(fit_dt_men_scaled) #Looks healthy

#Posterior predictive check
pp_check(fit_dt_men_scaled) #Looks good

fit_dt_men_scaled = add_criterion(fit_dt_men_scaled, criterion = "loo")
# Warning message:
# Found 3 observations with a pareto_k > 0.7 in model 'fit_dt_men_scaled'. 
# We recommend to set 'moment_match = TRUE' in order to perform moment matching
# for problematic observations.

loo_fitdt_men <- loo(fit_dt_men_scaled)
loo_fitdt_men_mm <- loo_moment_match(fit_dt_men_scaled, loo = loo_fitdt_men)
loo_fitdt_men_mm
# Computed from 48000 by 240 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -343.4 13.7
# p_loo        40.1  3.8
# looic       686.7 27.4
# ------
# MCSE of elpd_loo is 0.1.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.3]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_men)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item) 
# disc ~ 0 + G + (1 | ID)
# Data: d_men (Number of observations: 240) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 16) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.71      0.19     0.42     1.14 1.00    12244    21374
# sd(disc_Intercept)     0.55      0.17     0.29     0.94 1.00    12680    23285
# 
# ~item (Number of levels: 3) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.82      0.56     0.06     2.22 1.00    14717    14386
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.52      0.63    -3.81    -1.36 1.00    33964    33212
# Intercept[L1,2]    -1.96      0.53    -3.00    -0.93 1.00    36699    37942
# Intercept[L1,3]    -0.98      0.49    -1.90    -0.00 1.00    30612    33552
# Intercept[L1,4]    -0.14      0.49    -1.03     0.84 1.00    27612    29415
# Intercept[L1,5]     0.89      0.51    -0.01     1.91 1.00    24790    25232
# Intercept[L1,6]     2.40      0.56     1.36     3.53 1.00    20938    24213
# Intercept[L2,1]    -2.18      0.51    -3.21    -1.23 1.00    25339    30614
# Intercept[L2,2]    -0.94      0.39    -1.73    -0.17 1.00    41415    38793
# Intercept[L2,3]     0.01      0.36    -0.71     0.74 1.00    53632    40916
# Intercept[L2,4]     0.13      0.36    -0.57     0.87 1.00    54120    41253
# Intercept[L2,5]     0.61      0.37    -0.11     1.37 1.00    49505    40322
# Intercept[L2,6]     2.07      0.46     1.20     3.00 1.00    28381    33831
# Intercept[L3,1]    -1.98      0.53    -3.07    -0.99 1.00    28470    31954
# Intercept[L3,2]    -1.31      0.43    -2.16    -0.46 1.00    38270    37933
# Intercept[L3,3]    -0.36      0.38    -1.09     0.41 1.00    47965    38807
# Intercept[L3,4]     0.03      0.37    -0.68     0.80 1.00    50583    38791
# Intercept[L3,5]     0.85      0.39     0.13     1.65 1.00    43612    38964
# Intercept[L3,6]     1.99      0.45     1.15     2.91 1.00    27913    34906
# G1                  0.43      0.16     0.14     0.75 1.00    28908    30752
# disc_G1             0.00      0.15    -0.29     0.29 1.00    26235    32122
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
# There were 12 divergent transitions after warmup. Increasing adapt_delta 
# above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

# Figures 1–3 

library(tidybayes)
library(ggplot2)
library(dplyr)

################# Figure 1 - Between Judge Variability #################

# Extract standard deviations for country slopes (random slopes by ID for C)
country_sds_scaled <- spread_draws(fitProb_final_scaled,
                                   sd_ID__CCJ,
                                   sd_ID__CCK,
                                   sd_ID__CCL,
                                   sd_ID__CCM)

# Reshape to long format
country_sds_long_scaled <- country_sds_scaled %>%
  pivot_longer(cols = starts_with("sd_ID__"),
               names_to = "country",
               values_to = "sd") %>%
  mutate(country = recode(country,
                          sd_ID__CCJ = "UK",
                          sd_ID__CCK = "Germany",
                          sd_ID__CCL = "Brazil",
                          sd_ID__CCM = "Thailand"))

# Violin + boxplot of standard deviations (noise)
ggplot(country_sds_long_scaled, aes(x = country, y = sd, fill = country)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "white") +
  labs(
    #title = "Between-Judge Variability by Country",
    y = "Standard Deviation of Random Slopes",
    x = "Country"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

################# Figure 2 - Noise and effect plots #################

library(tidyverse)
library(brms)
library(posterior)
library(forcats)
library(patchwork)

# 1. FIXED EFFECTS — pull posterior draws for b_CC*
fixed_draws <- as_draws_df(fit_dt_reduced_adjusted_thresholds) %>%
  select(.draw, b_CCJ, b_CCK, b_CCL, b_CCM) %>%
  pivot_longer(cols = starts_with("b_"), names_to = "coef", values_to = "estimate") %>%
  mutate(
    country = recode(coef,
                     b_CCJ = "UK",
                     b_CCK = "Germany",
                     b_CCL = "Brazil",
                     b_CCM = "Thailand")
  )

# 2. RANDOM SLOPES — participant-level variation in CCJ–CCM
random_slope_vars <- names(as_draws_df(fit_dt_reduced_adjusted_thresholds)) %>%
  str_subset("^r_ID\\[.*,(CCJ|CCK|CCL|CCM)\\]")

random_slopes <- as_draws_df(fit_dt_reduced_adjusted_thresholds) %>%
  select(all_of(random_slope_vars)) %>%
  pivot_longer(cols = everything(), names_to = "term", values_to = "slope") %>%
  mutate(
    country_code = str_extract(term, "CCJ|CCK|CCL|CCM"),
    country = recode(country_code,
                     CCJ = "UK",
                     CCK = "Germany",
                     CCL = "Brazil",
                     CCM = "Thailand")
  )

# 3. PLOT — Posterior distributions of fixed effects
p_fixed <- ggplot(fixed_draws, aes(x = estimate, y = fct_rev(country), fill = country)) +
  geom_violin(alpha = 0.6, scale = "width", color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Posterior Distributions: Fixed Effects (relative to Sweden)",
       x = "Effect on Latent Rating Scale", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# 4. PLOT — Random slopes per participant
p_random <- ggplot(random_slopes, aes(x = slope, y = fct_rev(country), fill = country)) +
  geom_violin(alpha = 0.6, scale = "width", color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Participant-Level Variation: Random Slopes (relative to Sweden)",
       x = "Participant-Specific Effect", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

# 5. Combine vertically
p_fixed / p_random + plot_layout(ncol = 1)

############## Figure 3 - General tendency ###############

############## General Tendency — Figure 3 ###############

library(tidyverse)
library(brms)
library(posterior)
library(ggplot2)

# ------- Option 1: Density of Intercepts ---------

# Extract random intercepts for ID
re_id <- ranef(fit_dt_reduced_adjusted_thresholds)$ID[, , "Intercept"]

# Turn into a tidy tibble
intercepts_df <- as_tibble(re_id) %>%
  rename_with(~ paste0("ID", seq_along(.)), everything()) %>%
  pivot_longer(everything(), names_to = "participant", values_to = "intercept")

# Density plot
ggplot(intercepts_df, aes(x = intercept)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Participant-Level Intercepts",
    subtitle = "Variability in general rating tendency (harshness vs. leniency)",
    x = "Intercept on latent scale",
    y = "Density"
  ) +
  theme_minimal()


# ------- Option 2: Ranked Posterior Means + 89% CI ---------

# Step 1: Extract posterior draws for participant-level intercepts
draws <- as_draws_df(fit_dt_reduced_adjusted_thresholds)

# Get all intercept terms from r_ID[ , Intercept]
intercepts_df <- draws %>%
  select(starts_with("r_ID[")) %>%
  select(matches(",Intercept]")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "intercept") %>%
  mutate(
    ID = str_extract(param, "(?<=r_ID\\[)\\d+"),
    ID = as.integer(ID)
  )

# Step 2: Summarize posterior means and 89% CI
summary_df <- intercepts_df %>%
  group_by(ID) %>%
  summarise(
    mean = mean(intercept),
    lower = quantile(intercept, 0.055),
    upper = quantile(intercept, 0.945),
    .groups = "drop"
  ) %>%
  arrange(mean) %>%
  mutate(ID = factor(ID, levels = ID))  # to preserve order in plot

# Step 3: Plot
ggplot(summary_df, aes(x = ID, y = mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Participant-Level Intercepts (General Leniency vs. Harshness)",
    subtitle = "Higher = harsher ratings; 89% credible intervals shown",
    x = "Participant (ranked)",
    y = "Estimated Intercept (Latent Rating Scale)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# ------- Option 3: Color by participant gender Ranked Posterior Means + 89% CI ---------

library(tidyverse)
library(posterior)
library(ggplot2)

# Step 1: Extract intercept draws
draws <- as_draws_df(fit_dt_reduced_adjusted_thresholds)

intercepts_df <- draws %>%
  select(starts_with("r_ID[")) %>%
  select(matches(",Intercept]")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "intercept") %>%
  mutate(
    ID = str_extract(param, "(?<=r_ID\\[)\\d+"),
    ID = as.integer(ID)
  )

# Step 2: Add gender info from data (distinct per ID)
id_gender <- d %>%
  select(ID, A) %>%
  distinct() %>%
  mutate(A = recode(as.character(A), "0" = "Man", "1" = "Woman"))

# Step 3: Summarize posterior draws
summary_df <- intercepts_df %>%
  group_by(ID) %>%
  summarise(
    mean = mean(intercept),
    lower = quantile(intercept, 0.055),
    upper = quantile(intercept, 0.945),
    .groups = "drop"
  ) %>%
  left_join(id_gender, by = "ID") %>%
  arrange(mean) %>%
  mutate(ID = factor(ID, levels = ID))

# Step 4: Plot with gender color
ggplot(summary_df, aes(x = ID, y = mean, color = A)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Participant-Level Intercepts (General Leniency vs. Harshness)",
    subtitle = "Lower = harsher ratings; color = gender; 89% credible intervals shown",
    x = "Participant (ranked)",
    y = "Estimated Intercept (Latent Rating Scale)",
    color = "Gender"
  ) +
  scale_color_manual(values = c("Man" = "blue", "Woman" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Step 4.2 Alternative without X-axis labels and titles
ggplot(summary_df, aes(x = ID, y = mean, color = A)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Participant (ordered)",
    y = "Estimated Intercept (Latent Rating Scale)",
    color = "Gender"
  ) +
  scale_color_manual(values = c("Man" = "blue", "Woman" = "red")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # remove x-axis text
    axis.ticks.x = element_blank(), # remove x-axis ticks
    legend.position = "right"
  )

# ⑤ SENSITIVITY ANALYSIS ---------------------------------------

# Probit model
sens_prob <- powerscale_sensitivity(fitProb_final_scaled)
print(sens_prob)
# Sensitivity based on cjs_dist
# Prior selection: all priors
# Likelihood selection: all data
# 
# variable prior likelihood
# b_Intercept 0.228      0.113
# b_CCJ 0.109      0.027
# b_CCK 0.143      0.045
# b_CCL 0.115      0.016
# b_CCM 0.111      0.029
# b_G1 0.109      0.047
# sd_ID__Intercept 0.311      0.196
# sd_ID__CCJ 0.419      0.208
# sd_ID__CCK 0.343      0.271
# sd_ID__CCL 0.424      0.184
# sd_ID__CCM 0.448      0.204
# cor_ID__Intercept__CCJ 0.073      0.028
# cor_ID__Intercept__CCK 0.071      0.080
# cor_ID__CCJ__CCK 0.092      0.143
# cor_ID__Intercept__CCL 0.031      0.020
# cor_ID__CCJ__CCL 0.057      0.027
# cor_ID__CCK__CCL 0.048      0.071
# cor_ID__Intercept__CCM 0.024      0.050
# cor_ID__CCJ__CCM 0.046      0.030
# cor_ID__CCK__CCM 0.038      0.086
# cor_ID__CCL__CCM 0.050      0.029
# Intercept 0.192      0.133

# Ordinal IRT model
sens_dt <- powerscale_sensitivity(fit_dt_reduced_adjusted_thresholds)
print(sens_dt)
# Sensitivity based on cjs_dist
# Prior selection: all priors
# Likelihood selection: all data
# 
# variable prior likelihood
# b_Intercept[L1,1] 0.188      0.033
# b_Intercept[L1,2] 0.131      0.049
# b_Intercept[L1,3] 0.078      0.022
# b_Intercept[L1,4] 0.091      0.023
# b_Intercept[L1,5] 0.109      0.042
# b_Intercept[L1,6] 0.124      0.092
# b_Intercept[L2,1] 0.109      0.108
# b_Intercept[L2,2] 0.096      0.044
# b_Intercept[L2,3] 0.084      0.023
# b_Intercept[L2,4] 0.082      0.023
# b_Intercept[L2,5] 0.074      0.018
# b_Intercept[L2,6] 0.091      0.119
# b_Intercept[L3,1] 0.106      0.052
# b_Intercept[L3,2] 0.101      0.020
# b_Intercept[L3,3] 0.094      0.021
# b_Intercept[L3,4] 0.087      0.019
# b_Intercept[L3,5] 0.080      0.049
# b_Intercept[L3,6] 0.093      0.102
# b_CCJ 0.101      0.277
# b_CCK 0.094      0.307
# b_CCL 0.110      0.244
# b_CCM 0.119      0.308
# b_G1 0.019      0.053
# b_D 0.090      0.100
# b_disc_CCJ 0.018      0.322
# b_disc_CCK 0.015      0.344
# b_disc_CCL 0.008      0.440
# b_disc_CCM 0.012      0.333
# b_disc_G1 0.038      0.233
# b_disc_D 0.042      0.042
# sd_ID__Intercept 0.062      0.138
# sd_ID__CCJ 0.154      1.201
# sd_ID__CCK 0.160      1.172
# sd_ID__CCL 0.124      0.902
# sd_ID__CCM 0.196      1.487
# sd_item__Intercept 0.195      0.067
# sd_ID__disc_Intercept 0.054      0.702
# cor_ID__Intercept__CCJ 0.088      0.357
# cor_ID__Intercept__CCK 0.061      0.195
# cor_ID__CCJ__CCK 0.099      0.528
# cor_ID__Intercept__CCL 0.035      0.160
# cor_ID__CCJ__CCL 0.052      0.206
# cor_ID__CCK__CCL 0.071      0.492
# cor_ID__Intercept__CCM 0.045      0.181
# cor_ID__CCJ__CCM 0.069      0.288
# cor_ID__CCK__CCM 0.083      0.452
# cor_ID__CCL__CCM 0.068      0.478

# Optional: Plot ECDFs
# powerscale_plot_ecdf(powerscale_sequence(fitProb_final_scaled, parameter = "b_CGermany"), parameter = "b_CGermany")
# powerscale_plot_ecdf(powerscale_sequence(fit_dt_reduced_adjusted_thresholds, parameter = "b_CCJ"), parameter = "b_CCJ")

############# Prior vs. Posterior #############

# For probit model: fixed effects ~ Normal(0, 2)
prior_sd_probit <- 2

# For ordinal IRT model:
# location effects ~ Normal(0, 2)
# discrimination effects ~ Normal(0, 1)
location_prior_sd <- 2
disc_prior_sd <- 1

library(ggplot2)
library(tibble)
library(dplyr)
library(posterior)

##Probit model
# Get posterior draws
draws_probit <- as_draws_df(fitProb_final_scaled)

# Extract fixed effect names
param_names_probit <- names(draws_probit)[grepl("^b_", names(draws_probit))]

for (param in param_names_probit) {
  posterior_vals <- draws_probit[[param]]
  prior_vals <- rnorm(10000, 0, prior_sd_probit)
  
  plot_data <- bind_rows(
    tibble(value = posterior_vals, distribution = "Posterior"),
    tibble(value = prior_vals, distribution = "Prior")
  )
  
  p <- ggplot(plot_data, aes(x = value, fill = distribution, color = distribution)) +
    geom_density(alpha = 0.5, linewidth = 0.7) +
    labs(
      #title = paste("Prior vs. Posterior for", param),
      x = "Estimate",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_manual(values = c("Prior" = "gray", "Posterior" = "steelblue")) +
    scale_color_manual(values = c("Prior" = "gray40", "Posterior" = "steelblue"))
  
  ggsave(filename = paste0("prior_posterior_", param, ".png"),
         plot = p, width = 6, height = 4, dpi = 300)
}

##IRT model
library(tibble)
library(dplyr)
library(ggplot2)
library(posterior)

# Prior SDs
location_prior_sd <- 2
disc_prior_sd <- 1

# Extract posterior draws
draws_irt <- as_draws_df(fit_dt_reduced_adjusted_thresholds)

# Get parameter names
location_params <- names(draws_irt)[grepl("^b_", names(draws_irt))]
disc_params <- names(draws_irt)[grepl("^b_disc_", names(draws_irt))]

# Function to make and save plot
make_prior_posterior_plot <- function(param_name, posterior_vals, prior_vals, folder = "IRT_prior_post_plots") {
  plot_data <- bind_rows(
    tibble(value = posterior_vals, distribution = "Posterior"),
    tibble(value = prior_vals, distribution = "Prior")
  )
  
  p <- ggplot(plot_data, aes(x = value, fill = distribution, color = distribution)) +
    geom_density(alpha = 0.5, linewidth = 0.7) +
    labs(
      #title = paste("Prior vs. Posterior for", param_name),
      x = "Estimate",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_manual(values = c("Prior" = "gray", "Posterior" = "steelblue")) +
    scale_color_manual(values = c("Prior" = "gray40", "Posterior" = "steelblue"))
  
  # Create folder if needed
  if (!dir.exists(folder)) dir.create(folder)
  
  ggsave(
    filename = file.path(folder, paste0("prior_posterior_", param_name, ".png")),
    plot = p,
    width = 6, height = 4, dpi = 300, bg = "white"
  )
}

# Generate plots for location parameters
for (param in location_params) {
  posterior_vals <- draws_irt[[param]]
  prior_vals <- rnorm(10000, 0, location_prior_sd)
  make_prior_posterior_plot(param, posterior_vals, prior_vals)
}

# Generate plots for discrimination parameters
for (param in disc_params) {
  posterior_vals <- draws_irt[[param]]
  prior_vals <- rnorm(10000, 0, disc_prior_sd)
  make_prior_posterior_plot(param, posterior_vals, prior_vals)
}

sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.4.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Stockholm
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.0        marginaleffects_0.25.0 bayesplot_1.11.1      
# [4] posterior_1.6.1        tidybayes_3.0.7        lubridate_1.9.4       
# [7] forcats_1.0.0          stringr_1.5.1          dplyr_1.1.4           
# [10] purrr_1.0.4            readr_2.1.5            tidyr_1.3.1           
# [13] tibble_3.2.1           ggplot2_3.5.2          tidyverse_2.0.0       
# [16] readxl_1.4.5           priorsense_1.1.1.9000  brms_2.22.0           
# [19] Rcpp_1.0.14           
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6         tensorA_0.36.2.1     QuickJSR_1.6.0      
# [4] processx_3.8.6       inline_0.3.21        lattice_0.22-6      
# [7] tzdb_0.4.0           callr_3.7.6          vctrs_0.6.5         
# [10] tools_4.4.1          ps_1.9.0             generics_0.1.3      
# [13] stats4_4.4.1         curl_6.2.1           parallel_4.4.1      
# [16] pkgconfig_2.0.3      Matrix_1.7-2         data.table_1.17.0   
# [19] checkmate_2.3.2      RColorBrewer_1.1-3   distributional_0.5.0
# [22] RcppParallel_5.1.10  lifecycle_1.0.4      compiler_4.4.1      
# [25] farver_2.1.2         Brobdingnag_1.2-9    codetools_0.2-20    
# [28] pillar_1.10.2        arrayhelpers_1.1-0   StanHeaders_2.32.10 
# [31] bridgesampling_1.1-2 abind_1.4-8          nlme_3.1-167        
# [34] rstan_2.32.6         svUnit_1.0.6         tidyselect_1.2.1    
# [37] mvtnorm_1.3-3        stringi_1.8.4        reshape2_1.4.4      
# [40] labeling_0.4.3       grid_4.4.1           colorspace_2.1-1    
# [43] cli_3.6.5            magrittr_2.0.3       loo_2.8.0           
# [46] dichromat_2.0-0.1    pkgbuild_1.4.6       withr_3.0.2         
# [49] scales_1.4.0         backports_1.5.0      timechange_0.3.0    
# [52] estimability_1.5.1   matrixStats_1.5.0    emmeans_1.10.7      
# [55] gridExtra_2.3        cellranger_1.1.0     hms_1.1.3           
# [58] coda_0.19-4.1        ggdist_3.3.3         V8_6.0.1            
# [61] rstantools_2.4.0     rlang_1.1.6          xtable_1.8-4        
# [64] glue_1.8.0           rstudioapi_0.17.1    jsonlite_1.9.1      
# [67] R6_2.6.1             plyr_1.8.9   