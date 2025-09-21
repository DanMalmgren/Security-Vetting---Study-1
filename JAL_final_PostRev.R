# =====================================================
# Judgements of Loyalty – A Bayesian Workflow
# 
# Final models after revision
# =====================================================

########### Packages ###########

library(brms) #Modelling
library(tidyverse)
library(future) #Multisession
library(dplyr)
library(marginaleffects)
library(priorsense)
library(bayesplot)
library(readxl) #Data import
library(tidyr)
library(psych) #For Cronbach

plan(multisession)
options(future.globals.maxSize = 8 * 1024^3)  # Increase memory cap

########### General information ###########

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

#Question D = Har information som kan ha betydelse för säkerhetsprövningen 
#framkommit?

#Question D in English
#Has information that may be relevant to the security assessment emerged?

#L1 = I vilken utsträckning instämmer du i påståendet att personen kan antas 
#vara lojal mot de intressen som Säkerhetsskyddslag (2018:585) ska skydda på en 
#skala från 1 till 7

#L1 in Enlgish
#To what extent do you agree with the statement that the person can be assumed 
#to be loyal to the interests that the Security Protection Act (2018:585) aims 
#to protect, on a scale from 1 to 7? 

#L2 = I vilken utsträckning instämmer du i påståendet att där föreligger en 
#risk för en lojalitetskonflikt på en skala från 1 till 7?

#L2 in Enlgish
#To what extent do you agree with the statement that there is a risk of a 
#conflict of loyalty, on a scale from 1 to 7?

#L3 = I vilken utsträckning instämmer du i påståendet att personen kan antas 
#värdera lojaliteten till Sverige högre än lojaliteten mot annan part om 
#personen tvingas välja mellan Sverige och den andra parten på en skala från 
#1 till 7?

#L3 in Enlish
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

########### Loading data and preprocessing ###########


IRM4R <- read_excel("IRM4R.xlsx")
View(IRM4R)

d <- reshape2::melt(IRM4R, id.vars=c("ID", "A", "G", "C", "T", "D"), measure.vars=c("L1", "L2", "L3" ),
                    variable.name="item",value.name="rating")

d <- d[ order(d$ID, d$C, d$item), ]

d$C <- as.factor(d$C)
d$T <- as.factor(d$T)
d$G <- as.factor(d$G)
d$A <- as.factor(d$A)
d$ID <- as.factor(d$ID)

nCores = parallel::detectCores()

# Read in the Excel file with group info
groups <- read_excel("Groups.xlsx")   

# Check that the IDs match your main dataset
head(groups)
head(d$ID)

# Both ID variables as character
d <- d %>% mutate(ID = as.character(ID))
groups <- groups %>% mutate(ID = as.character(ID))

# Join
d <- d %>%
  left_join(groups, by = "ID")

# Check 
table(d$Group, useNA = "ifany")   

view(d)

d %>% 
  select(rating, C, G, D, Group, ID, item) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) 

anyNA(d)  #Check for missing values
#[1] FALSE

########### DESCRIPTIVES ###########

with(IRM4R, tapply(L1, list(C, T), median)) #Median scores for L1
#     TA TB TC
# CI 6.5  7  6
# CJ 6.0  6  6
# CK 6.0  6  5
# CL 6.0  6  5
# CM 6.0  6  5

with(IRM4R, tapply(L2, list(C, T), median)) #Median scores for L2
#     TA TB TC
# CI 6.0  6  6
# CJ 6.0  5  3
# CK 5.5  5  4
# CL 4.5  3  3
# CM 5.5  3  3

with(IRM4R, tapply(L3, list(C, T), median)) #Median scores for L3
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

########### Cronbach and friends ###########

# Check for duplicates per participant × country × stimulus gender × item
d %>%
  filter(item %in% c("L1","L2","L3")) %>%
  count(ID, C, G, item) %>%
  filter(n > 1)
# [1] ID   C    G    item n   
# <0 rows> (or 0-length row.names)

rm(loyalty_wide)
# One row per protocol 
loyalty_wide <- d %>%
  filter(item %in% c("L1","L2","L3")) %>%
  select(ID, C, G, item, rating) %>%
  pivot_wider(names_from = item, values_from = rating) %>%
  # keep only complete cases for L1–L3 to avoid NA warnings in alpha
  filter(!if_any(c(L1, L2, L3), is.na))

# Quick check
head(loyalty_wide)
# A tibble: 6 × 6
# ID    C     G        L1    L2    L3
# <chr> <fct> <fct> <dbl> <dbl> <dbl>
# 1 1     CI    1         5     3     6
# 2 1     CJ    0         4     4     4
# 3 1     CK    0         5     1     7
# 4 1     CL    0         4     2     6
# 5 1     CM    0         4     3     5
# 6 2     CI    1         6     6     7


# Cronbach's alpha
alpha_classic <- psych::alpha(loyalty_wide[, c("L1","L2","L3")])

alpha_classic$total$raw_alpha      # Cronbach’s alpha
# [1] 0.5736464

alpha_classic$total$std.alpha      # Standardized alpha
# [1] 0.6193364

# Ordinal alpha
poly_r <- psych::polychoric(loyalty_wide[, c("L1","L2","L3")])$rho

alpha_ordinal <- psych::alpha(poly_r)

alpha_ordinal$total$raw_alpha     
# [1] 0.6834393


########### Probit Pre-Scaling ###########

# Priors
prob_priors_full <- c(
  prior(normal(0, 2), class = "b"),                 # fixed effects: C, G, T, A, Group
  prior(exponential(0.5), class = "sd"),            # random-effect SDs
  prior(lkj(2), class = "cor")                      # random-effect correlation
)

# Full probit (pre-scaling)
fit_prob_full <- brm(
  D ~ 1 + C + G + T + A + Group + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_priors_full,
  seed = 5476,
  chains = 4, iter = 8000, warmup = 2000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_full_preScaling"
)

# 828.199 seconds in total

summary(fit_prob_full)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + T + A + Group + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 8000; warmup = 2000; thin = 1;
# total post-warmup draws = 24000
# 
# Multilevel Hyperparameters:
# ~ID (Number of levels: 58) 
#                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          5.53      1.32     3.35     8.50 1.00     8927    14626
# sd(CCJ)                7.77      2.88     3.60    14.73 1.00    16184    16970
# sd(CCK)                6.46      3.34     0.89    14.23 1.00     4931     3010
# sd(CCL)               11.87      3.56     6.33    20.16 1.00    20214    18920
# sd(CCM)               12.94      3.50     7.40    20.89 1.00    23414    19304
# cor(Intercept,CCJ)     0.33      0.23    -0.16     0.71 1.00    12745    15199
# cor(Intercept,CCK)     0.43      0.25    -0.16     0.83 1.00    12562    15106
# cor(CCJ,CCK)           0.50      0.21    -0.01     0.81 1.00     8323     9151
# cor(Intercept,CCL)     0.04      0.22    -0.39     0.45 1.00     7659    11845
# cor(CCJ,CCL)           0.35      0.18    -0.03     0.67 1.00     7413    12748
# cor(CCK,CCL)           0.22      0.22    -0.25     0.60 1.00     5291     7736
# cor(Intercept,CCM)    -0.11      0.19    -0.48     0.27 1.00     9848    15992
# cor(CCJ,CCM)           0.44      0.17     0.07     0.74 1.00     8350    14354
# cor(CCK,CCM)           0.30      0.21    -0.14     0.68 1.00     7181     8358
# cor(CCL,CCM)           0.56      0.13     0.27     0.79 1.00    15002    19679
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          -4.51      1.82    -8.23    -1.08 1.00    15508    17355
# CCJ                 1.85      1.38    -1.15     4.34 1.00    13914    17310
# CCK                 3.13      1.42     0.15     5.69 1.00     7859    10810
# CCL                 2.16      1.52    -0.94     5.00 1.00    20064    19605
# CCM                 1.86      1.54    -1.29     4.77 1.00    20816    19270
# G1                 -0.84      0.90    -2.65     0.90 1.00    17148    17777
# TTB                -0.02      1.46    -2.85     2.89 1.00    15585    17196
# TTC                -0.17      1.51    -3.07     2.82 1.00    15294    17516
# A1                 -0.78      1.45    -3.65     2.05 1.00    18347    18460
# GroupSpecialist     2.77      1.44    -0.15     5.55 1.00    18350    16395
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


########### IRM Pre-Scaling ###########

# Item-specific thresholds
thresh_priors_items <- c(
  # L1
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L1),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L1),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L1),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L1),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L1),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L1),
  
  # L2
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L2),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L2),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L2),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L2),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L2),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L2),
  
  # L3
  prior(normal(-1.07, 1), class = Intercept, coef = 1, group = L3),
  prior(normal(-0.57, 1), class = Intercept, coef = 2, group = L3),
  prior(normal(-0.18, 1), class = Intercept, coef = 3, group = L3),
  prior(normal( 0.18, 1), class = Intercept, coef = 4, group = L3),
  prior(normal( 0.57, 1), class = Intercept, coef = 5, group = L3),
  prior(normal( 1.07, 1), class = Intercept, coef = 6, group = L3)
)

# Original priors
irm_priors_original <- c(
  thresh_priors_items,
  prior(normal(0, 1), class = b),                  # location effects
  prior(normal(0, 0.25), class = b, dpar = disc),  # discrimination effects
  prior(exponential(1), class = sd),               # RE sds
  prior(lkj(2), class = cor)                       # RE correlations
)

# Full IRM (pre-scaling): location uses C,G,D,T,A,Group; discrimination uses same
fit_irm_full <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + T + A + Group + (1 + C | ID) + (1 | item),
    disc ~ 0 + C + G + D + T + A + Group + (1 | ID),
    cmc = FALSE
  ),
  data = d,
  family = cumulative(probit),
  prior = irm_priors_original,
  seed = 5476,
  chains = 4, iter = 8000, warmup = 2000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_irm_full_preScaling"
)

summary(fit_irm_full)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + C + G + D + T + A + Group + (1 + C | ID) + (1 | item) 
# disc ~ 0 + C + G + D + T + A + Group + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 8000; warmup = 2000; thin = 1;
# total post-warmup draws = 24000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.68      0.12     0.48     0.94 1.00     5136    10350
# sd(CCJ)                0.19      0.12     0.01     0.46 1.00     4303     8153
# sd(CCK)                0.21      0.13     0.01     0.48 1.00     4261     8249
# sd(CCL)                0.41      0.14     0.11     0.69 1.00     3712     4128
# sd(CCM)                0.13      0.09     0.01     0.35 1.00     5243     8328
# sd(disc_Intercept)     0.62      0.09     0.46     0.81 1.00     6667    13167
# cor(Intercept,CCJ)     0.13      0.31    -0.50     0.68 1.00    20211    16547
# cor(Intercept,CCK)     0.17      0.30    -0.46     0.71 1.00    18434    16170
# cor(CCJ,CCK)           0.14      0.35    -0.58     0.75 1.00     8691    11603
# cor(Intercept,CCL)     0.01      0.23    -0.42     0.49 1.00    11171    13109
# cor(CCJ,CCL)          -0.14      0.33    -0.72     0.54 1.00     4316     8148
# cor(CCK,CCL)           0.14      0.32    -0.53     0.71 1.00     4828    10381
# cor(Intercept,CCM)    -0.07      0.33    -0.67     0.59 1.00    21652    17713
# cor(CCJ,CCM)           0.04      0.35    -0.64     0.69 1.00    13625    15044
# cor(CCK,CCM)           0.06      0.35    -0.62     0.70 1.00    12182    12626
# cor(CCL,CCM)           0.13      0.35    -0.58     0.74 1.00    11534    16612
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.49      0.88     1.24     4.62 1.00     7609    11816
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]         -2.43      0.63    -3.74    -1.28 1.00    22635    17620
# Intercept[L1,2]         -1.83      0.52    -2.87    -0.83 1.00    21230    20262
# Intercept[L1,3]         -0.61      0.46    -1.50     0.29 1.00    22978    17669
# Intercept[L1,4]          0.12      0.45    -0.75     1.00 1.00    30943    17521
# Intercept[L1,5]          1.20      0.46     0.32     2.10 1.00    26359    18805
# Intercept[L1,6]          2.85      0.52     1.84     3.87 1.00    11126    17187
# Intercept[L2,1]         -2.05      0.48    -3.00    -1.11 1.00    12315    15299
# Intercept[L2,2]         -1.04      0.43    -1.88    -0.21 1.00    22255    16461
# Intercept[L2,3]         -0.19      0.41    -1.00     0.61 1.00    33571    15710
# Intercept[L2,4]         -0.00      0.41    -0.82     0.81 1.00    34981    15993
# Intercept[L2,5]          0.53      0.42    -0.30     1.36 1.00    30312    16384
# Intercept[L2,6]          2.30      0.49     1.35     3.28 1.00    11019    15788
# Intercept[L3,1]         -1.98      0.49    -2.96    -1.04 1.00    12276    15852
# Intercept[L3,2]         -1.13      0.43    -1.98    -0.29 1.00    21729    17015
# Intercept[L3,3]         -0.48      0.42    -1.29     0.34 1.00    35426    18769
# Intercept[L3,4]         -0.05      0.41    -0.86     0.77 1.00    40121    18107
# Intercept[L3,5]          0.86      0.43     0.03     1.71 1.00    28276    17356
# Intercept[L3,6]          2.23      0.48     1.30     3.19 1.00    11947    14864
# CCJ                     -1.12      0.19    -1.52    -0.76 1.00     5682     9530
# CCK                     -1.08      0.19    -1.47    -0.73 1.00     5753     9594
# CCL                     -1.41      0.21    -1.85    -1.01 1.00     5086     9046
# CCM                     -1.33      0.20    -1.74    -0.96 1.00     4857     9223
# G1                       0.04      0.07    -0.10     0.19 1.00    21492    17525
# D                       -0.76      0.13    -1.03    -0.52 1.00     6734    11728
# TTB                     -0.01      0.24    -0.48     0.47 1.00     5675    10097
# TTC                     -0.72      0.25    -1.24    -0.24 1.00     5477     9755
# A1                      -0.01      0.26    -0.50     0.50 1.00     6086    10163
# GroupSpecialist         -0.10      0.25    -0.59     0.40 1.00     6322    10630
# disc_CCI                -0.71      0.13    -0.96    -0.46 1.00    10386    16141
# disc_CCJ                 0.17      0.12    -0.07     0.41 1.00     9501    15894
# disc_CCK                 0.25      0.12     0.01     0.49 1.00     9639    16199
# disc_CCL                 0.15      0.13    -0.10     0.39 1.00     8867    14807
# disc_CCM                 0.39      0.13     0.15     0.64 1.00    10094    15723
# disc_G1                  0.04      0.08    -0.11     0.20 1.00    20557    18406
# disc_D                   0.20      0.10    -0.01     0.40 1.00    15238    16316
# disc_TTB                -0.02      0.16    -0.33     0.29 1.00    10012    14510
# disc_TTC                -0.12      0.16    -0.43     0.19 1.00     8572    12077
# disc_A1                  0.03      0.14    -0.25     0.31 1.00     7402    12389
# disc_GroupSpecialist     0.01      0.16    -0.30     0.32 1.00     9398    13602
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

########### Sensitivity Analysis Script ###########
# Packages
library(brms)
library(priorsense)     # Kallioinen et al. power-scaling
library(loo)
library(tidyverse)

### Helpers ### 
pp_overlay <- function(fit, par_regex = "^b_") {
  # Quick prior vs posterior overlay for fixed effects
  # par_regex examples: "^b_" (location), "^b_disc_" (discrimination)
  pp <- posterior::as_draws_df(fit) %>% as_tibble()
  pars <- names(pp)[grepl(par_regex, names(pp))]
  if (length(pars) == 0) return(invisible(NULL))
  pri <- prior_samples(fit) %>% as_tibble()
  for (p in pars) {
    if (!p %in% names(pri)) next
    df <- bind_rows(
      tibble(val = pri[[p]],   src = "Prior"),
      tibble(val = pp[[p]],    src = "Posterior")
    )
    g <- ggplot(df, aes(val, after_stat(density), linetype = src)) +
      geom_density() + theme_minimal(base_size = 12) +
      labs(title = paste("Prior vs Posterior:", p), x = p, y = NULL)
    print(g)
  }
}

### Power-scaling sensitivity (priorsense) ###

run_powerscale <- function(fit, label, grid = c(0, 0.25, 0.5, 0.75, 1)) {
  # powerscale_sensitivity() quantifies prior vs likelihood influence
  # For larger models, keep grid short to manage runtime.
  ps <- powerscale_sensitivity(fit, scalings = grid)
  print(ps)
  # Nice summary plot per parameter:
  plot(ps) + ggtitle(paste("Power-scaling sensitivity:", label))
  # Table of influence indices:
  out <- summary(ps)
  out$label <- label
  out
}

ps_prob <- run_powerscale(fit_prob_full, "Probit full (pre-scaling)")
# Sensitivity based on cjs_dist
# Prior selection: all priors
# Likelihood selection: all data
# 
# variable          prior likelihood                                diagnosis
# b_Intercept       0.201      0.082           potential prior-data conflict
# b_CCJ             0.099      0.036           potential strong prior / weak likelihood
# b_CCK             0.144      0.044           potential strong prior / weak likelihood
# b_CCL             0.122      0.026           potential strong prior / weak likelihood
# b_CCM             0.105      0.037           potential strong prior / weak likelihood
# b_G1              0.108      0.043           potential strong prior / weak likelihood
# b_TTB             0.085      0.047           potential strong prior / weak likelihood
# b_TTC             0.089      0.043           potential strong prior / weak likelihood
# b_A1              0.082      0.022           potential strong prior / weak likelihood
# b_GroupSpecialist 0.168      0.042           potential strong prior / weak likelihood
# sd_ID__Intercept  0.324      0.179           potential prior-data conflict
# sd_ID__CCJ        0.409      0.236           potential prior-data conflict
# sd_ID__CCK        0.375      0.259           potential prior-data conflict
# sd_ID__CCL        0.446      0.194           potential prior-data conflict
# sd_ID__CCM        0.467      0.206           potential prior-data conflict
# cor_ID__Intercept__CCJ 0.080      0.031      potential strong prior / weak likelihood
# cor_ID__Intercept__CCK 0.102      0.087      potential prior-data conflict
# cor_ID__CCJ__CCK 0.103            0.121      potential prior-data conflict
# cor_ID__Intercept__CCL 0.028      0.022                                        -
# cor_ID__CCJ__CCL 0.051      0.031            potential strong prior / weak likelihood
# cor_ID__CCK__CCL 0.043      0.057                                        -
# cor_ID__Intercept__CCM 0.025      0.044                                        -
# cor_ID__CCJ__CCM 0.046      0.024                                        -
# cor_ID__CCK__CCM 0.025      0.054                                        -
# cor_ID__CCL__CCM 0.055      0.027            potential strong prior / weak likelihood

ps_irm <- run_powerscale(fit_irm_full, "IRM original priors (pre-scaling)")
# Sensitivity based on cjs_dist
# Prior selection: all priors
# Likelihood selection: all data
# variable            prior likelihood                                diagnosis
# b_Intercept[L1,1]   0.230      0.049 potential strong prior / weak likelihood
# b_Intercept[L1,2]   0.210      0.055            potential prior-data conflict
# b_Intercept[L1,3]   0.106      0.048 potential strong prior / weak likelihood
# b_Intercept[L1,4]   0.085      0.036 potential strong prior / weak likelihood
# b_Intercept[L1,5]   0.177      0.031 potential strong prior / weak likelihood
# b_Intercept[L1,6]   0.294      0.039 potential strong prior / weak likelihood
# b_Intercept[L2,1]   0.249      0.040 potential strong prior / weak likelihood
# b_Intercept[L2,2]   0.166      0.044 potential strong prior / weak likelihood
# b_Intercept[L2,3]   0.098      0.023 potential strong prior / weak likelihood
# b_Intercept[L2,4]   0.085      0.021 potential strong prior / weak likelihood
# b_Intercept[L2,5]   0.114      0.024 potential strong prior / weak likelihood
# b_Intercept[L2,6]   0.249      0.053            potential prior-data conflict
# b_Intercept[L3,1]   0.231      0.049 potential strong prior / weak likelihood
# b_Intercept[L3,2]   0.168      0.043 potential strong prior / weak likelihood
# b_Intercept[L3,3]   0.117      0.032 potential strong prior / weak likelihood
# b_Intercept[L3,4]   0.083      0.024 potential strong prior / weak likelihood
# b_Intercept[L3,5]   0.139      0.035 potential strong prior / weak likelihood
# b_Intercept[L3,6]   0.241      0.029 potential strong prior / weak likelihood
# b_CCJ               0.320      0.245            potential prior-data conflict
# b_CCK               0.315      0.257            potential prior-data conflict
# b_CCL               0.355      0.152            potential prior-data conflict
# b_CCM               0.380      0.227            potential prior-data conflict
# b_G1                0.060      0.060            potential prior-data conflict
# b_D                 0.304      0.076            potential prior-data conflict
# b_TTB               0.048      0.038                                        -
# b_TTC               0.157      0.075            potential prior-data conflict
# b_A1                0.050      0.028 potential strong prior / weak likelihood
# b_GroupSpecialist   0.060      0.041 potential strong prior / weak likelihood

# Create output folder if it doesn't exist
if (!dir.exists("replication_outputs")) dir.create("replication_outputs")

# ---------------- IRM model ----------------
ps_irm <- powerscale_sensitivity(fit_irm_full, scalings = c(0, 0.25, 0.5, 0.75, 1))
print(ps_irm)

# Save summary table
ps_summary_irm <- summary(ps_irm)
write.csv(ps_summary_irm,
          file = "replication_outputs/ps_irm_original_summary.csv",
          row.names = FALSE)

# Save plot
pdf("replication_outputs/ps_irm_original_plot.pdf", width = 8, height = 6)
plot(ps_irm, main = "Power-scaling sensitivity: IRM original priors")
dev.off()

# ---------------- Probit model ----------------
ps_prob <- powerscale_sensitivity(fit_prob_full, scalings = c(0, 0.25, 0.5, 0.75, 1))
print(ps_prob)

# Save summary table
ps_summary_prob <- summary(ps_prob)
write.csv(ps_summary_prob,
          file = "replication_outputs/ps_prob_original_summary.csv",
          row.names = FALSE)

# Save plot
pdf("replication_outputs/ps_prob_full_plot.pdf", width = 8, height = 6)
plot(ps_prob, main = "Power-scaling sensitivity: Probit full (original priors)")
dev.off()


########### Ordinal (IRM) Model Scaled Priors ###########

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

#For kfold comparison of included predictors
#Group included

fit_dt_withGroup_light <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + Group + (1 + C | ID) + (1 | item),
    disc ~ 0 + C + G + D + Group + (1 | ID),
    cmc = FALSE
  ),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_adjusted,
  warmup = 2000, iter = 6000,
  chains = 2, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE)
)
# 1134.69 seconds in total

pp_check(fit_dt_withGroup_light) #Looks good

plot(fit_dt_withGroup_light) #Looks healthy

summary(fit_dt_withGroup_light)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + C + G + D + Group + (1 + C | ID) + (1 | item) 
# disc ~ 0 + C + G + D + Group + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 2 chains, each with iter = 6000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.77      0.14     0.52     1.06 1.00     1992     3454
# sd(CCJ)                0.21      0.13     0.01     0.50 1.00     1544     2703
# sd(CCK)                0.27      0.14     0.03     0.58 1.00     1662     2798
# sd(CCL)                0.44      0.15     0.13     0.74 1.00     1475     1531
# sd(CCM)                0.14      0.10     0.01     0.38 1.00     1845     2919
# sd(disc_Intercept)     0.63      0.09     0.48     0.82 1.00     3136     5137
# cor(Intercept,CCJ)     0.16      0.30    -0.46     0.71 1.00     9110     5392
# cor(Intercept,CCK)     0.31      0.27    -0.29     0.76 1.00     7551     5625
# cor(CCJ,CCK)           0.20      0.35    -0.54     0.78 1.00     3128     5030
# cor(Intercept,CCL)     0.03      0.23    -0.41     0.49 1.00     5196     4611
# cor(CCJ,CCL)          -0.12      0.32    -0.70     0.54 1.00     1862     2724
# cor(CCK,CCL)           0.17      0.30    -0.48     0.69 1.00     2097     4472
# cor(Intercept,CCM)    -0.08      0.33    -0.68     0.58 1.00     8417     5713
# cor(CCJ,CCM)           0.05      0.35    -0.64     0.70 1.00     5191     5768
# cor(CCK,CCM)           0.08      0.36    -0.61     0.72 1.00     3989     3734
# cor(CCL,CCM)           0.14      0.34    -0.56     0.74 1.00     3805     4846
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.89      1.26     1.29     6.10 1.00     3424     4293
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]         -2.40      0.64    -3.76    -1.22 1.00    10490     6571
# Intercept[L1,2]         -1.79      0.53    -2.85    -0.79 1.00     8185     7062
# Intercept[L1,3]         -0.58      0.46    -1.48     0.33 1.00     6799     6545
# Intercept[L1,4]          0.13      0.45    -0.75     1.02 1.00     8385     6145
# Intercept[L1,5]          1.19      0.46     0.30     2.12 1.00     7477     6372
# Intercept[L1,6]          2.87      0.55     1.80     3.96 1.00     3550     5471
# Intercept[L2,1]         -2.04      0.49    -3.00    -1.08 1.00     4039     4890
# Intercept[L2,2]         -1.03      0.42    -1.86    -0.19 1.00     7521     6027
# Intercept[L2,3]         -0.19      0.41    -0.99     0.60 1.00    15517     6248
# Intercept[L2,4]         -0.00      0.41    -0.80     0.80 1.00    16316     6152
# Intercept[L2,5]          0.53      0.41    -0.28     1.33 1.00    13784     6432
# Intercept[L2,6]          2.38      0.52     1.38     3.42 1.00     3884     5320
# Intercept[L3,1]         -1.96      0.50    -2.98    -1.01 1.00     4658     5499
# Intercept[L3,2]         -1.10      0.43    -1.97    -0.27 1.00     7243     6016
# Intercept[L3,3]         -0.47      0.41    -1.27     0.35 1.00    12079     6545
# Intercept[L3,4]         -0.05      0.41    -0.85     0.76 1.00    14291     6876
# Intercept[L3,5]          0.86      0.43     0.02     1.69 1.00     9627     6492
# Intercept[L3,6]          2.28      0.51     1.27     3.29 1.00     3878     5441
# CCJ                     -1.28      0.24    -1.79    -0.86 1.00     2070     3597
# CCK                     -1.23      0.24    -1.73    -0.81 1.00     2176     3554
# CCL                     -1.57      0.27    -2.14    -1.09 1.00     1946     3617
# CCM                     -1.49      0.25    -2.02    -1.04 1.00     1863     3468
# G1                       0.04      0.08    -0.10     0.19 1.00     8799     6492
# D                       -0.77      0.14    -1.06    -0.51 1.00     2328     4239
# GroupSpecialist          0.01      0.25    -0.49     0.51 1.00     2524     4087
# disc_CCI                -0.80      0.19    -1.17    -0.42 1.00     1744     3442
# disc_CCJ                 0.20      0.19    -0.16     0.57 1.00     1709     3297
# disc_CCK                 0.30      0.19    -0.07     0.69 1.00     1702     2911
# disc_CCL                 0.19      0.20    -0.20     0.57 1.00     1737     3414
# disc_CCM                 0.46      0.19     0.09     0.84 1.00     1728     3419
# disc_G1                  0.03      0.09    -0.14     0.21 1.00     9215     6486
# disc_D                   0.18      0.12    -0.05     0.42 1.00     6045     6138
# disc_GroupSpecialist    -0.01      0.20    -0.42     0.40 1.00     3236     4812
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

kfold_withGroup_light <- kfold(
  fit_dt_withGroup_light,
  K = 10,
  chains = 2,
  iter = 6000,
  save_fits = TRUE,
  seed = 5476,
  file = "kfold_withGroup_light"
)

print(kfold_withGroup_light)
# Based on 10-fold cross-validation.
# 
# Estimate   SE
# elpd_kfold  -1067.8 26.9
# p_kfold       148.9 10.0
# kfoldic      2135.7 53.8

# Main model with matched settings. Group excluded as predictor
fit_dt_reduced_light <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + (1 + C | ID) + (1 | item)
  ) +
    lf(disc ~ 0 + C + G + D + (1 | ID), cmc = FALSE),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_adjusted,
  warmup = 2000, iter = 6000,  
  chains = 2, cores = 2,
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_dt_reduced_light"
)

# 11867.17 seconds

summary(fit_dt_reduced_light)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + C + G + D + (1 + C | ID) + (1 | item) 
# disc ~ 0 + C + G + D + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 2 chains, each with iter = 6000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.47      0.08     0.32     0.64 1.00     1412     2430
# sd(CCJ)                0.15      0.09     0.01     0.35 1.00     1071     2543
# sd(CCK)                0.20      0.10     0.03     0.41 1.00      882     1792
# sd(CCL)                0.30      0.09     0.12     0.50 1.01     1074     1724
# sd(CCM)                0.10      0.07     0.00     0.27 1.00      962     1148
# sd(disc_Intercept)     0.65      0.10     0.49     0.86 1.00     1847     3872
# cor(Intercept,CCJ)     0.20      0.30    -0.42     0.72 1.00     4312     4321
# cor(Intercept,CCK)     0.33      0.26    -0.23     0.77 1.00     4191     4122
# cor(CCJ,CCK)           0.24      0.34    -0.49     0.79 1.00     1618     2631
# cor(Intercept,CCL)     0.05      0.22    -0.37     0.51 1.00     3316     3488
# cor(CCJ,CCL)          -0.08      0.33    -0.67     0.57 1.00     1064     1600
# cor(CCK,CCL)           0.24      0.29    -0.39     0.73 1.00     1553     2687
# cor(Intercept,CCM)    -0.02      0.33    -0.63     0.62 1.00     5113     5697
# cor(CCJ,CCM)           0.09      0.36    -0.61     0.72 1.00     2536     3685
# cor(CCK,CCM)           0.14      0.36    -0.58     0.76 1.00     1877     3227
# cor(CCL,CCM)           0.19      0.34    -0.53     0.77 1.00     2481     4462
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.06      1.00     0.83     4.70 1.00     2567     4108
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.12      0.67    -3.56    -0.90 1.00     6762     5739
# Intercept[L1,2]    -1.41      0.53    -2.46    -0.41 1.00     7497     7084
# Intercept[L1,3]    -0.29      0.44    -1.17     0.58 1.00     9385     6451
# Intercept[L1,4]     0.19      0.43    -0.66     1.04 1.00    10439     6543
# Intercept[L1,5]     0.89      0.44     0.02     1.76 1.00     8054     6507
# Intercept[L1,6]     1.96      0.48     1.02     2.92 1.00     3635     5542
# Intercept[L2,1]    -1.34      0.43    -2.19    -0.49 1.00     4585     4836
# Intercept[L2,2]    -0.69      0.40    -1.47     0.12 1.00     7903     5953
# Intercept[L2,3]    -0.16      0.40    -0.93     0.64 1.00     9875     6022
# Intercept[L2,4]    -0.04      0.40    -0.81     0.76 1.00     9872     6195
# Intercept[L2,5]     0.31      0.40    -0.48     1.10 1.00     9301     6052
# Intercept[L2,6]     1.46      0.44     0.61     2.35 1.00     4361     5423
# Intercept[L3,1]    -1.29      0.45    -2.16    -0.44 1.00     5874     5196
# Intercept[L3,2]    -0.75      0.42    -1.56     0.06 1.00     8526     5772
# Intercept[L3,3]    -0.34      0.41    -1.13     0.46 1.00    10035     5856
# Intercept[L3,4]    -0.07      0.41    -0.86     0.74 1.00    10415     5755
# Intercept[L3,5]     0.52      0.42    -0.31     1.34 1.00     8320     5560
# Intercept[L3,6]     1.40      0.45     0.52     2.28 1.00     4686     4766
# CCJ                -0.80      0.13    -1.07    -0.55 1.00     1660     3082
# CCK                -0.77      0.13    -1.04    -0.53 1.00     1731     2982
# CCL                -0.99      0.15    -1.31    -0.72 1.00     1439     2712
# CCM                -0.93      0.14    -1.22    -0.67 1.00     1430     2864
# G1                  0.03      0.05    -0.06     0.13 1.00     6900     6436
# D                  -0.49      0.09    -0.68    -0.34 1.00     1678     3144
# disc_CCJ            0.83      0.13     0.58     1.08 1.00     3220     5155
# disc_CCK            0.92      0.13     0.67     1.17 1.00     2890     4364
# disc_CCL            0.83      0.13     0.59     1.08 1.00     2758     4657
# disc_CCM            1.08      0.13     0.83     1.34 1.00     3436     5199
# disc_G1            -0.03      0.09    -0.20     0.14 1.00     4892     5805
# disc_D              0.12      0.12    -0.11     0.35 1.00     3844     5035
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# kfold
kfold_reduced_light <- kfold(
  fit_dt_reduced_light,
  K = 10,
  chains = 2,
  iter = 6000,
  save_fits = TRUE,
  seed = 5476,
  file = "kfold_reduced_light"
)

print(kfold_reduced_light)
# Based on 10-fold cross-validation.
# 
#            Estimate   SE
# elpd_kfold  -1052.3 27.0
# p_kfold       137.0  9.3
# kfoldic      2104.6 54.0

loo_compare(kfold_reduced_light, kfold_withGroup_light)
#                         elpd_diff se_diff
# fit_dt_reduced_light     0.0       0.0  
# fit_dt_withGroup_light -15.5       8.6  

#Same model with type of relationship and participant gender as predictors
fit_dt_with_TA_light <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + T + A + (1 + C | ID) + (1 | item),
    disc ~ 0 + C + G + D + T + A + (1 | ID),
    cmc = FALSE
  ),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_adjusted,
  warmup = 2000, iter = 6000,
  chains = 2, cores = 2,
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_dt_with_TA_light"
)

summary(fit_dt_with_TA_light)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: rating | thres(gr = item) ~ 1 + C + G + D + T + A + (1 + C | ID) + (1 | item) 
# disc ~ 0 + C + G + D + T + A + (1 | ID)
# Data: d (Number of observations: 870) 
# Draws: 2 chains, each with iter = 6000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.72      0.13     0.50     1.00 1.00     2101     3851
# sd(CCJ)                0.21      0.13     0.01     0.50 1.00     1713     3509
# sd(CCK)                0.23      0.14     0.01     0.53 1.00     1818     3177
# sd(CCL)                0.45      0.15     0.15     0.75 1.00     1774     2282
# sd(CCM)                0.14      0.10     0.01     0.38 1.00     1944     3947
# sd(disc_Intercept)     0.64      0.09     0.48     0.84 1.00     2318     4720
# cor(Intercept,CCJ)     0.12      0.30    -0.49     0.68 1.00     8129     5207
# cor(Intercept,CCK)     0.16      0.30    -0.46     0.70 1.00     7739     5846
# cor(CCJ,CCK)           0.14      0.35    -0.58     0.74 1.00     3745     5378
# cor(Intercept,CCL)    -0.01      0.23    -0.44     0.45 1.00     5019     5420
# cor(CCJ,CCL)          -0.15      0.32    -0.72     0.53 1.00     1801     2866
# cor(CCK,CCL)           0.13      0.32    -0.53     0.70 1.00     2054     3950
# cor(Intercept,CCM)    -0.10      0.32    -0.67     0.56 1.00     7209     5388
# cor(CCJ,CCM)           0.03      0.35    -0.65     0.69 1.00     5185     5721
# cor(CCK,CCM)           0.05      0.35    -0.64     0.68 1.00     4577     5685
# cor(CCL,CCM)           0.13      0.34    -0.57     0.73 1.00     5027     6391
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     3.19      1.35     1.44     6.65 1.00     3074     3966
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.43      0.62    -3.74    -1.29 1.00     7651     6092
# Intercept[L1,2]    -1.84      0.53    -2.89    -0.84 1.00     7448     6214
# Intercept[L1,3]    -0.60      0.45    -1.48     0.30 1.00     7350     6514
# Intercept[L1,4]     0.14      0.44    -0.70     1.01 1.00     9735     6745
# Intercept[L1,5]     1.26      0.45     0.39     2.15 1.00     8835     6200
# Intercept[L1,6]     2.99      0.54     1.95     4.08 1.00     4240     5280
# Intercept[L2,1]    -2.12      0.48    -3.09    -1.18 1.00     4265     5214
# Intercept[L2,2]    -1.06      0.42    -1.88    -0.23 1.00     7378     5774
# Intercept[L2,3]    -0.19      0.40    -0.97     0.60 1.00    12137     5928
# Intercept[L2,4]     0.01      0.40    -0.78     0.81 1.00    12544     5770
# Intercept[L2,5]     0.57      0.41    -0.25     1.39 1.00    11079     6164
# Intercept[L2,6]     2.45      0.51     1.47     3.48 1.00     4047     5242
# Intercept[L3,1]    -2.05      0.51    -3.10    -1.06 1.00     5001     5383
# Intercept[L3,2]    -1.15      0.44    -2.00    -0.30 1.00     6803     6218
# Intercept[L3,3]    -0.47      0.42    -1.28     0.34 1.00     9887     6383
# Intercept[L3,4]    -0.03      0.41    -0.83     0.78 1.00    11224     6327
# Intercept[L3,5]     0.92      0.43     0.07     1.76 1.00     8628     5960
# Intercept[L3,6]     2.40      0.52     1.40     3.44 1.00     3965     5354
# CCJ                -1.30      0.24    -1.81    -0.87 1.00     2417     3676
# CCK                -1.25      0.24    -1.76    -0.82 1.00     2453     3826
# CCL                -1.60      0.27    -2.17    -1.12 1.00     2151     3697
# CCM                -1.51      0.26    -2.05    -1.06 1.00     2123     3502
# G1                  0.04      0.08    -0.11     0.19 1.00     7568     5788
# D                  -0.79      0.14    -1.10    -0.53 1.00     2485     3630
# TTB                -0.06      0.26    -0.57     0.46 1.00     2476     3577
# TTC                -0.80      0.28    -1.38    -0.27 1.00     2236     3535
# A1                  0.03      0.24    -0.44     0.52 1.00     2436     3693
# disc_CCI           -0.68      0.24    -1.15    -0.22 1.00     1699     3434
# disc_CCJ            0.33      0.24    -0.14     0.79 1.00     1676     3273
# disc_CCK            0.43      0.24    -0.04     0.90 1.00     1665     3199
# disc_CCL            0.31      0.24    -0.17     0.79 1.00     1675     3272
# disc_CCM            0.59      0.24     0.10     1.06 1.00     1672     3443
# disc_G1             0.04      0.09    -0.12     0.21 1.00     8202     6479
# disc_D              0.19      0.12    -0.04     0.42 1.00     5758     6224
# disc_TTB           -0.15      0.22    -0.58     0.29 1.00     2990     4217
# disc_TTC           -0.28      0.23    -0.71     0.18 1.00     2687     3543
# disc_A1            -0.04      0.20    -0.44     0.36 1.00     2406     3354
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

kfold_with_TA_light <- kfold(
  fit_dt_with_TA_light,
  K = 10,
  chains = 2,
  iter = 6000,
  save_fits = TRUE,
  seed = 5476,
  file = "kfold_with_TA_light"
)

print(kfold_with_TA_light)
# Based on 10-fold cross-validation.
# 
# Estimate   SE
# elpd_kfold  -1080.0 27.2
# p_kfold       158.3  9.9
# kfoldic      2160.0 54.5

loo_compare(kfold_reduced_light, kfold_withGroup_light, kfold_with_TA_light)
# elpd_diff se_diff
# fit_dt_reduced_light     0.0       0.0  
# fit_dt_withGroup_light -15.5       8.6  
# fit_dt_with_TA_light   -27.7       7.8  


#### Final version of IRM ####
# 4 chains for inference

fit_dt_reduced_final <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + C + G + D + (1 + C | ID) + (1 | item),
    disc ~ 0 + C + G + D + (1 | ID),
    cmc = FALSE
  ),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_adjusted,
  warmup = 4000, iter = 16000,
  chains = 4, cores = parallel::detectCores(),
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_dt_reduced_final"
)

#Trace and density plots
plot(fit_dt_reduced_final) #Looks healthy
plots <- plot(fit_dt_reduced_final)
pdf("fit_dt_reduced_final_all_plots.pdf", width = 8, height = 6)
for (p in plots) print(p)
dev.off()

pp_check(fit_dt_reduced_final) #Looks good

summary(fit_dt_reduced_final) 
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.75      0.14     0.51     1.05 1.00     7426    14081
# sd(CCJ)                0.21      0.13     0.01     0.49 1.00     7326    16341
# sd(CCK)                0.27      0.14     0.02     0.56 1.00     7238    11993
# sd(CCL)                0.43      0.15     0.13     0.72 1.00     5837     5578
# sd(CCM)                0.14      0.10     0.01     0.37 1.00     7875    15856
# sd(disc_Intercept)     0.62      0.09     0.47     0.81 1.00    13582    25494
# cor(Intercept,CCJ)     0.16      0.30    -0.47     0.70 1.00    39926    30793
# cor(Intercept,CCK)     0.31      0.27    -0.29     0.78 1.00    36830    28917
# cor(CCJ,CCK)           0.20      0.34    -0.53     0.77 1.00    13238    25828
# cor(Intercept,CCL)     0.03      0.23    -0.40     0.51 1.00    20645    22633
# cor(CCJ,CCL)          -0.12      0.32    -0.70     0.52 1.00     8421    16143
# cor(CCK,CCL)           0.18      0.30    -0.47     0.70 1.00     9671    19748
# cor(Intercept,CCM)    -0.08      0.34    -0.68     0.60 1.00    33783    33980
# cor(CCJ,CCM)           0.04      0.36    -0.64     0.69 1.00    21226    31042
# cor(CCK,CCM)           0.08      0.35    -0.62     0.71 1.00    19333    27357
# cor(CCL,CCM)           0.14      0.34    -0.56     0.74 1.00    18539    31304
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.87      1.23     1.29     6.01 1.00    15400    21135
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.39      0.64    -3.74    -1.22 1.00    41352    35581
# Intercept[L1,2]    -1.78      0.53    -2.84    -0.76 1.00    35938    39410
# Intercept[L1,3]    -0.57      0.46    -1.47     0.32 1.00    40339    36919
# Intercept[L1,4]     0.13      0.44    -0.73     1.00 1.00    60775    36784
# Intercept[L1,5]     1.19      0.46     0.30     2.09 1.00    43374    36791
# Intercept[L1,6]     2.85      0.54     1.80     3.94 1.00    13182    23540
# Intercept[L2,1]    -2.02      0.49    -3.00    -1.07 1.00    16302    27271
# Intercept[L2,2]    -1.02      0.43    -1.86    -0.18 1.00    34415    33206
# Intercept[L2,3]    -0.19      0.41    -0.99     0.61 1.00    75980    36206
# Intercept[L2,4]    -0.01      0.41    -0.80     0.80 1.00    80070    36788
# Intercept[L2,5]     0.53      0.42    -0.29     1.34 1.00    63412    37439
# Intercept[L2,6]     2.36      0.52     1.36     3.39 1.00    14231    25781
# Intercept[L3,1]    -1.94      0.50    -2.95    -0.98 1.00    19376    29110
# Intercept[L3,2]    -1.09      0.44    -1.95    -0.25 1.00    33815    33228
# Intercept[L3,3]    -0.46      0.42    -1.27     0.36 1.00    59682    37990
# Intercept[L3,4]    -0.05      0.41    -0.86     0.77 1.00    77532    38450
# Intercept[L3,5]     0.86      0.43     0.02     1.71 1.00    41211    36313
# Intercept[L3,6]     2.27      0.51     1.28     3.30 1.00    13376    26926
# CCJ                -1.27      0.24    -1.77    -0.84 1.00     7722    15893
# CCK                -1.23      0.24    -1.72    -0.80 1.00     7990    16209
# CCL                -1.56      0.27    -2.12    -1.08 1.00     7053    14492
# CCM                -1.48      0.25    -2.01    -1.02 1.00     6879    14474
# G1                  0.04      0.07    -0.10     0.19 1.00    37251    35130
# D                  -0.76      0.14    -1.06    -0.51 1.00     8687    19213
# disc_CCI           -0.79      0.18    -1.15    -0.43 1.00     7048    14067
# disc_CCJ            0.21      0.18    -0.14     0.57 1.00     6886    13125
# disc_CCK            0.30      0.18    -0.05     0.67 1.00     7168    13508
# disc_CCL            0.19      0.19    -0.18     0.57 1.00     6901    14357
# disc_CCM            0.46      0.18     0.11     0.83 1.00     6999    13161
# disc_G1             0.03      0.09    -0.14     0.21 1.00    38382    37617
# disc_D              0.19      0.12    -0.05     0.42 1.00    26186    33092
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).




#### Subset models for gender ####
fit_dt_priors_subset <- c(
  # Thresholds
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
  
  # Effects
  prior(normal(0, 2), class = b),
  prior(normal(0, 1), class = b, dpar = disc),
  
  # Random effects
  prior(exponential(1), class = sd)
)

fit_dt_AG_final <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + A*G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + A*G + (1 | ID),
      cmc = FALSE
    ),
  data = d,
  family = cumulative(probit),
  prior = fit_dt_priors_subset,
  warmup = 4000, iter = 16000, chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = "fit_dt_AG_final"
)

pp_check(fit_dt_AG_final) #Looks good

plot(fit_dt_AG_final) #Looks healthy
plots <- plot(fit_dt_AG_final)
pdf("fit_dt_AG_final_all_plots.pdf", width = 8, height = 6)
for (p in plots) print(p)
dev.off()


add_criterion(fit_dt_AG_final, criterion = "loo")
# Warning message:
# Found 5 observations with a pareto_k > 0.7 in model 'fit_dt_AG_final'. 
# We recommend to set 'moment_match = TRUE' in order to perform moment matching 
# for problematic observations.

loo_dt_AG <- loo(fit_dt_AG_final)
loo_dt_AG_mm <- loo_moment_match(fit_dt_AG_final, loo = loo_dt_AG)
fit_dt_AG_final$criteria$loo <- loo_dt_AG_mm
print(fit_dt_AG_final$criteria$loo)
# Computed from 48000 by 870 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo  -1286.9 26.0
# p_loo       115.2  8.4
# looic      2573.8 52.1
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.7]).
# 
# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     865   99.4%   883     
# (0.7, 1]      (bad)        3    0.3%   <NA>    
# (1, Inf)      (very bad)   2    0.2%   <NA>    
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_AG_final)
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.77      0.12     0.56     1.03 1.00     7986    16653
# sd(disc_Intercept)     0.45      0.07     0.33     0.61 1.00    14375    25230
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.98      0.58     0.14     2.41 1.00    14248    11592
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.78      0.61    -4.05    -1.64 1.00    32710    33239
# Intercept[L1,2]    -2.28      0.53    -3.34    -1.25 1.00    31568    34023
# Intercept[L1,3]    -0.97      0.50    -1.93     0.00 1.00    26324    22789
# Intercept[L1,4]    -0.02      0.49    -0.97     0.94 1.00    26892    18266
# Intercept[L1,5]     1.25      0.51     0.26     2.24 1.00    22328    15415
# Intercept[L1,6]     2.70      0.56     1.61     3.80 1.00    15610    14774
# Intercept[L2,1]    -1.96      0.44    -2.84    -1.12 1.00    17967    30348
# Intercept[L2,2]    -0.98      0.38    -1.74    -0.23 1.00    29894    36399
# Intercept[L2,3]    -0.03      0.37    -0.76     0.70 1.00    44179    38008
# Intercept[L2,4]     0.18      0.37    -0.56     0.90 1.00    45312    38618
# Intercept[L2,5]     0.68      0.38    -0.06     1.43 1.00    40438    37859
# Intercept[L2,6]     2.01      0.43     1.18     2.87 1.00    20977    31560
# Intercept[L3,1]    -1.84      0.44    -2.70    -0.99 1.00    21124    28903
# Intercept[L3,2]    -1.18      0.39    -1.93    -0.40 1.00    30267    34534
# Intercept[L3,3]    -0.49      0.37    -1.21     0.26 1.00    42387    38144
# Intercept[L3,4]     0.00      0.37    -0.70     0.76 1.00    45391    38329
# Intercept[L3,5]     0.94      0.39     0.21     1.73 1.00    34508    36215
# Intercept[L3,6]     1.93      0.43     1.12     2.80 1.00    20234    30721
# A1                  0.49      0.26    -0.01     1.03 1.00     8313    16476
# G1                  0.42      0.15     0.13     0.72 1.00    20855    26255
# A1:G1              -0.46      0.19    -0.83    -0.10 1.00    21628    27097
# disc_A1            -0.14      0.14    -0.41     0.13 1.00     8301    15273
# disc_G1            -0.00      0.14    -0.27     0.27 1.00    17408    27795
# disc_A1:G1          0.00      0.17    -0.32     0.33 1.00    17926    26358
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


## Refit women only ##
d_women <- d[d$A == "1", ]

fit_dt_women_final <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + G + (1 | ID),
      cmc = FALSE
    ),
  data = d_women,
  family = cumulative(probit),
  prior = fit_dt_priors_subset,
  warmup = 4000, iter = 16000, chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = "fit_dt_women_final"
)

pp_check(fit_dt_women_final) #Looks good

plot(fit_dt_women_final) #Looks healthy

plots <- plot(fit_dt_women_final)
pdf("fit_dt_women_final_all_plots.pdf", width = 8, height = 6)
for (p in plots) print(p)
dev.off()

add_criterion(fit_dt_women_final, criterion = "loo")
loo_women <- loo(fit_dt_women_final)
print(loo_women)
loo_women_mm <- loo_moment_match(fit_dt_women_final, loo = loo_women)
print(loo_women_mm)
# Computed from 48000 by 630 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -956.5 22.4
# p_loo        88.0  7.7
# looic      1913.1 44.8
# ------
# MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.7]).
# 
# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     628   99.7%   651     
# (0.7, 1]      (bad)        1    0.2%   <NA>    
# (1, Inf)      (very bad)   1    0.2%   <NA>    
# See help('pareto-k-diagnostic') for details.

summary(fit_dt_women_final)
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.70      0.11     0.51     0.95 1.00    12839    24850
# sd(disc_Intercept)     0.40      0.08     0.26     0.58 1.00    15062    25995
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.24      0.61     0.42     2.76 1.00    21026    22143
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.49      0.63    -3.81    -1.35 1.00    49311    35097
# Intercept[L1,2]    -1.92      0.53    -2.97    -0.90 1.00    53867    37994
# Intercept[L1,3]    -0.81      0.47    -1.73     0.12 1.00    46392    29958
# Intercept[L1,4]     0.02      0.47    -0.90     0.94 1.00    46798    27295
# Intercept[L1,5]     1.20      0.48     0.25     2.14 1.00    42529    25912
# Intercept[L1,6]     2.43      0.51     1.43     3.42 1.00    35579    26151
# Intercept[L2,1]    -1.70      0.43    -2.53    -0.86 1.00    34439    33774
# Intercept[L2,2]    -0.92      0.39    -1.69    -0.14 1.00    43477    35392
# Intercept[L2,3]    -0.11      0.38    -0.85     0.65 1.00    50623    36161
# Intercept[L2,4]     0.10      0.39    -0.64     0.87 1.00    51421    36914
# Intercept[L2,5]     0.57      0.39    -0.18     1.34 1.00    50888    38390
# Intercept[L2,6]     1.67      0.42     0.87     2.49 1.00    42434    37879
# Intercept[L3,1]    -1.61      0.43    -2.46    -0.76 1.00    38473    34629
# Intercept[L3,2]    -1.01      0.40    -1.79    -0.21 1.00    46287    35961
# Intercept[L3,3]    -0.52      0.40    -1.29     0.27 1.00    50536    34989
# Intercept[L3,4]    -0.04      0.40    -0.80     0.76 1.00    51937    35778
# Intercept[L3,5]     0.81      0.40     0.04     1.62 1.00    49886    35382
# Intercept[L3,6]     1.61      0.42     0.81     2.45 1.00    42969    34278
# G1                 -0.02      0.10    -0.21     0.17 1.00    47788    39288
# disc_G1            -0.01      0.09    -0.19     0.16 1.00    29927    36061
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

## Refit men only ##
d_men <- d[d$A == "0", ]

fit_dt_men_final <- brm(
  bf(
    rating | thres(gr = item) ~ 1 + G + (1 | ID) + (1 | item)
  ) +
    lf(
      disc ~ 0 + G + (1 | ID),
      cmc = FALSE
    ),
  data = d_men,
  family = cumulative(probit),
  prior = fit_dt_priors_subset,
  warmup = 4000, iter = 16000, chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2,
  file = "fit_dt_men_final"
)

add_criterion(fit_dt_men_final, criterion = "loo")
loo_men <- loo(fit_dt_men_final)
loo_men_mm <- loo_moment_match(fit_dt_men_final, loo = loo_men)
print(loo_men_mm)

summary(fit_dt_men_final)
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.72      0.19     0.42     1.15 1.00    13532    24279
# sd(disc_Intercept)     0.55      0.17     0.29     0.94 1.00    13590    23699
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.82      0.55     0.06     2.20 1.00    16461    15331
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.52      0.63    -3.85    -1.36 1.00    44363    34607
# Intercept[L1,2]    -1.95      0.52    -3.00    -0.94 1.00    41578    37440
# Intercept[L1,3]    -0.98      0.48    -1.89    -0.01 1.00    31585    31670
# Intercept[L1,4]    -0.14      0.48    -1.02     0.83 1.00    29015    30371
# Intercept[L1,5]     0.89      0.50    -0.00     1.90 1.00    26955    27503
# Intercept[L1,6]     2.41      0.56     1.37     3.53 1.00    24205    25693
# Intercept[L2,1]    -2.18      0.50    -3.21    -1.23 1.00    30480    32854
# Intercept[L2,2]    -0.95      0.39    -1.71    -0.18 1.00    45749    38339
# Intercept[L2,3]     0.01      0.36    -0.70     0.73 1.00    61362    41091
# Intercept[L2,4]     0.14      0.36    -0.57     0.86 1.00    61691    41206
# Intercept[L2,5]     0.61      0.37    -0.11     1.35 1.00    57869    41686
# Intercept[L2,6]     2.08      0.46     1.20     3.00 1.00    33728    37760
# Intercept[L3,1]    -1.98      0.53    -3.10    -0.99 1.00    33648    34199
# Intercept[L3,2]    -1.30      0.43    -2.16    -0.45 1.00    45968    39795
# Intercept[L3,3]    -0.35      0.38    -1.08     0.42 1.00    56229    42007
# Intercept[L3,4]     0.04      0.38    -0.68     0.81 1.00    59335    41262
# Intercept[L3,5]     0.85      0.39     0.12     1.66 1.00    50809    39575
# Intercept[L3,6]     1.99      0.45     1.14     2.92 1.00    32729    36811
# G1                  0.44      0.16     0.15     0.77 1.00    32386    33031
# disc_G1            -0.00      0.15    -0.29     0.29 1.00    31130    34275
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
#   There were 14 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

fit_dt_men_final_refit <- update(fit_dt_men_final, control = list(adapt_delta = 0.999))

summary(fit_dt_men_final_refit)
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.72      0.18     0.43     1.14 1.00    14690    25709
# sd(disc_Intercept)     0.55      0.17     0.29     0.94 1.00    13796    23786
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.82      0.57     0.05     2.21 1.00    16244    15396
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.53      0.63    -3.84    -1.37 1.00    46536    36991
# Intercept[L1,2]    -1.96      0.53    -3.02    -0.93 1.00    44881    38540
# Intercept[L1,3]    -0.98      0.49    -1.90     0.00 1.00    33986    34637
# Intercept[L1,4]    -0.14      0.49    -1.03     0.84 1.00    30538    32890
# Intercept[L1,5]     0.89      0.51    -0.01     1.92 1.00    27334    27773
# Intercept[L1,6]     2.41      0.57     1.37     3.55 1.00    24028    25923
# Intercept[L2,1]    -2.18      0.51    -3.23    -1.23 1.00    32557    33565
# Intercept[L2,2]    -0.95      0.39    -1.72    -0.18 1.00    52008    39348
# Intercept[L2,3]     0.01      0.36    -0.71     0.73 1.00    67864    40615
# Intercept[L2,4]     0.14      0.36    -0.58     0.86 1.00    66313    41136
# Intercept[L2,5]     0.61      0.37    -0.11     1.36 1.00    60419    41464
# Intercept[L2,6]     2.08      0.46     1.20     3.01 1.00    33866    36639
# Intercept[L3,1]    -1.98      0.53    -3.07    -1.00 1.00    34915    34824
# Intercept[L3,2]    -1.30      0.43    -2.14    -0.47 1.00    45564    37913
# Intercept[L3,3]    -0.35      0.38    -1.07     0.42 1.00    57327    41171
# Intercept[L3,4]     0.04      0.37    -0.67     0.81 1.00    61589    40591
# Intercept[L3,5]     0.86      0.38     0.14     1.65 1.00    52378    40238
# Intercept[L3,6]     2.00      0.45     1.16     2.92 1.00    33906    35740
# G1                  0.44      0.16     0.15     0.77 1.00    33734    33504
# disc_G1            -0.00      0.15    -0.29     0.29 1.00    32550    37073
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
#   There were 5 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

fit_dt_men_final_refit2 <- update(fit_dt_men_final, control = list(adapt_delta = 0.9999))

summary(fit_dt_men_final_refit2)
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
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.71      0.19     0.42     1.15 1.00    10926    20645
# sd(disc_Intercept)     0.55      0.16     0.29     0.94 1.00    12460    23087
# 
# ~item (Number of levels: 3) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.82      0.56     0.06     2.22 1.00    12414    11993
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[L1,1]    -2.52      0.63    -3.84    -1.36 1.00    30556    31358
# Intercept[L1,2]    -1.95      0.53    -3.00    -0.92 1.00    34231    37886
# Intercept[L1,3]    -0.98      0.49    -1.90     0.00 1.00    26821    30828
# Intercept[L1,4]    -0.14      0.49    -1.03     0.84 1.00    24324    26294
# Intercept[L1,5]     0.89      0.51    -0.00     1.91 1.00    21545    22776
# Intercept[L1,6]     2.41      0.56     1.38     3.55 1.00    18567    21611
# Intercept[L2,1]    -2.18      0.51    -3.22    -1.23 1.00    21856    29484
# Intercept[L2,2]    -0.95      0.39    -1.72    -0.17 1.00    37917    38344
# Intercept[L2,3]     0.01      0.36    -0.70     0.73 1.00    53649    39877
# Intercept[L2,4]     0.14      0.36    -0.57     0.86 1.00    53798    40410
# Intercept[L2,5]     0.62      0.37    -0.10     1.36 1.00    48699    40695
# Intercept[L2,6]     2.08      0.46     1.21     3.01 1.00    26054    35388
# Intercept[L3,1]    -1.98      0.53    -3.07    -1.00 1.00    24555    31290
# Intercept[L3,2]    -1.31      0.43    -2.16    -0.46 1.00    36142    37942
# Intercept[L3,3]    -0.36      0.38    -1.08     0.42 1.00    48948    37920
# Intercept[L3,4]     0.03      0.37    -0.67     0.81 1.00    52111    39331
# Intercept[L3,5]     0.85      0.39     0.13     1.66 1.00    44215    37264
# Intercept[L3,6]     1.99      0.45     1.16     2.92 1.00    25998    33985
# G1                  0.44      0.16     0.14     0.77 1.00    25383    31093
# disc_G1            -0.00      0.15    -0.29     0.29 1.00    23638    33194
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Warning message: There were 2 divergent transitions after warmup. Increasing 
# adapt_delta above 0.9999 may help. 
# See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

saveRDS(fit_dt_men_final_refit2, "fit_dt_men_final_refit2.rds")

pp_check(fit_dt_men_final_refit2) #Looks good

plot(fit_dt_men_final_refit2) #Looks healthy
plots <- plot(fit_dt_men_final_refit2)
pdf("fit_dt_men_final_refit2_all_plots.pdf", width = 8, height = 6)
for (p in plots) print(p)
dev.off()

# LOO for final men’s subset model
loo_men_final <- loo(fit_dt_men_final_refit2)
loo_men_final_mm <- loo_moment_match(fit_dt_men_final_refit2, loo = loo_men_final)

# Attach loo values to the model
fit_dt_men_final_refit2$criteria$loo <- loo_men_final_mm
print(fit_dt_men_final_refit2$criteria$loo)
# Computed from 48000 by 240 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -343.2 13.7
# p_loo        39.9  3.8
# looic       686.4 27.4
# ------
# MCSE of elpd_loo is 0.1.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.5]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

# Save the model
saveRDS(fit_dt_men_final_refit2, "fit_dt_men_final_refit2.rds")


########### Probit models ###########

prob_prior_scaled <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 2), class = "b"),
  prior(exponential(0.5), class = "sd"),
  prior(lkj(2), class = "cor")
)

fit_prob_reduced_light <- brm(
  D ~ 1 + C + G + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 2000, iter = 6000,
  chains = 2, cores = 2,
  seed = 5476,
  control = list(adapt_delta = 0.99),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_reduced_light"
)

kfold_prob_reduced_light <- kfold(
  fit_prob_reduced_light,
  K = 10,
  chains = 2,
  iter = 6000,
  seed = 5476,
  save_fits = TRUE,
  file = "kfold_prob_reduced_light"
)

print(kfold_prob_reduced_light)
# Based on 10-fold cross-validation.
# 
# Estimate  SE
# elpd_kfold    -25.2 1.0
# p_kfold        10.7 0.7
# kfoldic        50.3 2.0

fit_prob_withGroup_light <- brm(
  D ~ 1 + C + G + Group + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 2000, iter = 6000,
  chains = 2, cores = 2,
  seed = 5476,
  control = list(adapt_delta = 0.99),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_withGroup_light"
)

kfold_prob_withGroup_light <- kfold(
  fit_prob_withGroup_light,
  K = 10,
  chains = 2,
  iter = 6000,
  seed = 5476,
  save_fits = TRUE,
  file = "kfold_prob_withGroup_light"
)

summary(fit_prob_withGroup_light)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + Group + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 2 chains, each with iter = 6000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          5.05      1.21     3.06     7.78 1.00     2115     3720
# sd(CCJ)                7.80      2.90     3.59    14.62 1.00     4233     4880
# sd(CCK)                6.58      3.21     1.53    14.37 1.00     1752     1127
# sd(CCL)               11.47      3.54     6.06    19.76 1.00     4335     5011
# sd(CCM)               12.38      3.49     7.03    20.62 1.00     5393     5479
# cor(Intercept,CCJ)     0.34      0.22    -0.13     0.72 1.00     2579     3943
# cor(Intercept,CCK)     0.44      0.24    -0.10     0.82 1.00     2905     4142
# cor(CCJ,CCK)           0.52      0.19     0.06     0.81 1.00     2303     2495
# cor(Intercept,CCL)     0.04      0.20    -0.36     0.43 1.00     2018     3720
# cor(CCJ,CCL)           0.35      0.18    -0.04     0.66 1.00     1752     3710
# cor(CCK,CCL)           0.22      0.21    -0.21     0.59 1.00     1504     3011
# cor(Intercept,CCM)    -0.14      0.19    -0.50     0.25 1.00     2147     4013
# cor(CCJ,CCM)           0.43      0.17     0.07     0.74 1.00     2429     4176
# cor(CCK,CCM)           0.30      0.20    -0.13     0.66 1.00     1837     2540
# cor(CCL,CCM)           0.56      0.13     0.27     0.79 1.00     4094     5744
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          -4.73      1.10    -7.05    -2.73 1.00     2532     4444
# CCJ                 1.78      1.33    -0.98     4.25 1.00     4072     4874
# CCK                 3.01      1.35     0.19     5.47 1.00     2388     3783
# CCL                 2.18      1.46    -0.83     4.84 1.00     4215     4816
# CCM                 1.94      1.49    -1.09     4.76 1.00     4593     4671
# G1                 -0.83      0.87    -2.55     0.86 1.00     3626     5043
# GroupSpecialist     3.00      1.31     0.36     5.55 1.00     3419     4719
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

#Follow up model with group and 4 chains 
fit_prob_withGroup_4 <- brm(
  D ~ 1 + C + G + Group + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 2000, iter = 8000,
  chains = 4, cores = 4,
  seed = 5476,
  control = list(adapt_delta = 0.99),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_withGroup_4"
)

pp_check(fit_prob_withGroup_4) #Looks good
plot(fit_prob_withGroup_4) #Looks healthy

loo_fit_prob_withGroup_4 <- loo(fit_prob_withGroup_4)
loo_fit_prob_withGroup_4_mm <- loo_moment_match(fit_prob_withGroup_4, loo = loo_fit_prob_withGroup_4)
print(loo_fit_prob_withGroup_4_mm)
# Computed from 24000 by 870 log-likelihood matrix.
# 
#          Estimate  SE
# elpd_loo    -21.5 0.8
# p_loo         7.0 0.3
# looic        43.1 1.5
# ------
# MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.1]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

summary(fit_prob_withGroup_4)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + Group + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 8000; warmup = 2000; thin = 1;
# total post-warmup draws = 24000
# 
# Multilevel Hyperparameters:
# ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          5.05      1.24     3.04     7.87 1.00     7065    11598
# sd(CCJ)                7.77      2.88     3.57    14.62 1.00    14474    15498
# sd(CCK)                6.54      3.25     1.22    14.25 1.00     4299     2548
# sd(CCL)               11.44      3.45     6.14    19.55 1.00    15499    15602
# sd(CCM)               12.37      3.45     7.01    20.39 1.00    19411    17798
# cor(Intercept,CCJ)     0.35      0.22    -0.13     0.72 1.00     9496    13205
# cor(Intercept,CCK)     0.44      0.25    -0.13     0.83 1.00     9117    12241
# cor(CCJ,CCK)           0.51      0.20     0.04     0.81 1.00     6404     6847
# cor(Intercept,CCL)     0.04      0.21    -0.37     0.43 1.00     6686    11113
# cor(CCJ,CCL)           0.35      0.18    -0.03     0.67 1.00     6249    12447
# cor(CCK,CCL)           0.22      0.21    -0.22     0.60 1.00     4626     7023
# cor(Intercept,CCM)    -0.14      0.19    -0.50     0.24 1.00     7251    12194
# cor(CCJ,CCM)           0.43      0.17     0.07     0.73 1.00     7752    13007
# cor(CCK,CCM)           0.30      0.21    -0.14     0.66 1.00     5623     5934
# cor(CCL,CCM)           0.56      0.13     0.27     0.79 1.00    13115    18145
# 
# Regression Coefficients:
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept          -4.74      1.10    -7.03    -2.74 1.00     8967    14867
# CCJ                 1.81      1.32    -0.97     4.28 1.00    11543    15755
# CCK                 3.02      1.36     0.18     5.54 1.00     6497     8857
# CCL                 2.20      1.46    -0.80     4.95 1.00    15150    16646
# CCM                 1.93      1.48    -1.10     4.72 1.00    15870    16534
# G1                 -0.81      0.87    -2.56     0.88 1.00    12071    16409
# GroupSpecialist     3.01      1.31     0.39     5.56 1.00    10914    14114
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

fit_prob_with_TA_light <- brm(
  D ~ 1 + C + G + T + A + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 2000, iter = 6000,
  chains = 2, cores = 2,
  seed = 5476,
  control = list(adapt_delta = 0.99),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_with_TA_light"
)

kfold_prob_with_TA_light <- kfold(
  fit_prob_with_TA_light,
  K = 10,
  chains = 2,
  iter = 6000,
  seed = 5476,
  save_fits = TRUE,
  file = "kfold_prob_with_TA_light"
)

summary(fit_prob_with_TA_light)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + T + A + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 2 chains, each with iter = 6000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
# ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          6.01      1.42     3.66     9.24 1.00     2520     4420
# sd(CCJ)                7.59      2.86     3.46    14.47 1.00     5015     6247
# sd(CCK)                5.51      3.30     0.33    13.10 1.00      991     1102
# sd(CCL)               11.86      3.59     6.29    20.37 1.00     5735     5450
# sd(CCM)               13.00      3.57     7.39    21.30 1.00     7791     5621
# cor(Intercept,CCJ)     0.29      0.23    -0.20     0.69 1.00     4967     5779
# cor(Intercept,CCK)     0.34      0.28    -0.29     0.80 1.00     4536     4764
# cor(CCJ,CCK)           0.45      0.24    -0.17     0.80 1.00     1746     2121
# cor(Intercept,CCL)     0.03      0.22    -0.40     0.45 1.00     2857     4178
# cor(CCJ,CCL)           0.34      0.18    -0.05     0.67 1.00     2176     3735
# cor(CCK,CCL)           0.21      0.23    -0.27     0.62 1.00     1301     1895
# cor(Intercept,CCM)    -0.11      0.19    -0.48     0.27 1.00     3438     4805
# cor(CCJ,CCM)           0.43      0.18     0.05     0.74 1.00     2725     4666
# cor(CCK,CCM)           0.30      0.23    -0.20     0.70 1.00     1603     1836
# cor(CCL,CCM)           0.56      0.14     0.27     0.79 1.00     4784     5946
# 
# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -3.03      1.73    -6.58     0.15 1.00     3173     5846
# CCJ           2.08      1.36    -0.79     4.55 1.00     4317     5871
# CCK           3.53      1.37     0.68     6.02 1.00     1642     4156
# CCL           2.31      1.53    -0.83     5.18 1.00     5913     6374
# CCM           2.01      1.52    -1.04     4.89 1.00     7327     6344
# G1           -0.84      0.89    -2.63     0.89 1.00     5394     6268
# TTB           0.22      1.51    -2.73     3.21 1.00     3590     5465
# TTC          -0.37      1.56    -3.42     2.71 1.00     3603     5084
# A1           -1.59      1.40    -4.27     1.20 1.00     4476     5273
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

loo_compare(kfold_prob_reduced_light, kfold_prob_withGroup_light, kfold_prob_with_TA_light)
# elpd_diff se_diff
# fit_prob_withGroup_light  0.0       0.0   
# fit_prob_reduced_light   -0.7       0.8   
# fit_prob_with_TA_light   -0.9       1.0   

#Same predictors as in the IRM for consistency
fit_prob_reduced_final <- brm(
  D ~ 1 + C + G + (1 + C | ID),
  data = d,
  family = bernoulli(link = "probit"),
  prior = prob_prior_scaled,
  warmup = 4000, iter = 16000,  
  chains = 4, cores = parallel::detectCores(),
  seed = 5476,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  init_r = 0.2,
  save_pars = save_pars(all = TRUE),
  file = "fit_prob_reduced_final"
)

pp_check(fit_prob_reduced_final) #Looks good

#Trace and density plots
plot(fit_prob_reduced_final) #Looks healthy
plots <- plot(fit_prob_reduced_final)
pdf("fit_prob_reduced_final_all_plots.pdf", width = 8, height = 6)
for (p in plots) print(p)
dev.off()

loo_fit_prob_reduced_final <- loo(fit_prob_reduced_final)
loo_fit_prob_reduced_final_mm <- loo_moment_match(fit_prob_reduced_final, loo = loo_fit_prob_reduced_final)

# Attach the criterion to the model for future reference
loo_fit_prob_reduced_final$criteria$loo <- loo_fit_prob_reduced_final_mm
print(loo_fit_prob_reduced_final$criteria$loo)
# Computed from 48000 by 870 log-likelihood matrix.
# 
# Estimate  SE
# elpd_loo    -21.6 0.8
# p_loo         7.1 0.3
# looic        43.2 1.5
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.2]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

summary(fit_prob_reduced_final)
# Family: bernoulli 
# Links: mu = probit 
# Formula: D ~ 1 + C + G + (1 + C | ID) 
# Data: d (Number of observations: 870) 
# Draws: 4 chains, each with iter = 16000; warmup = 4000; thin = 1;
# total post-warmup draws = 48000
# 
# Multilevel Hyperparameters:
#   ~ID (Number of levels: 58) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          5.91      1.39     3.61     9.03 1.00    16559    26579
# sd(CCJ)                8.06      3.05     3.57    15.34 1.00    33644    31020
# sd(CCK)                5.71      3.31     0.51    13.49 1.00     7403     5827
# sd(CCL)               11.60      3.52     6.18    19.81 1.00    36187    33843
# sd(CCM)               12.72      3.58     7.14    21.05 1.00    51401    36058
# cor(Intercept,CCJ)     0.31      0.23    -0.19     0.70 1.00    30197    31853
# cor(Intercept,CCK)     0.35      0.28    -0.29     0.80 1.00    30500    26435
# cor(CCJ,CCK)           0.48      0.23    -0.09     0.81 1.00    12676    14928
# cor(Intercept,CCL)     0.02      0.21    -0.39     0.41 1.00    21223    28856
# cor(CCJ,CCL)           0.34      0.18    -0.04     0.67 1.00    14727    24732
# cor(CCK,CCL)           0.21      0.23    -0.27     0.62 1.00     9061    13615
# cor(Intercept,CCM)    -0.14      0.19    -0.51     0.24 1.00    23423    34563
# cor(CCJ,CCM)           0.43      0.18     0.06     0.75 1.00    18044    30444
# cor(CCK,CCM)           0.31      0.22    -0.17     0.70 1.00    11194    11013
# cor(CCL,CCM)           0.57      0.13     0.27     0.79 1.00    28293    37429
# 
# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -4.08      1.09    -6.35    -2.08 1.00    18384    29586
# CCJ           1.91      1.38    -1.00     4.43 1.00    26007    33895
# CCK           3.39      1.37     0.46     5.86 1.00    12381    22733
# CCL           2.30      1.49    -0.75     5.08 1.00    38023    36494
# CCM           1.97      1.52    -1.15     4.83 1.00    39432    35495
# G1           -0.81      0.88    -2.57     0.91 1.00    32159    35473
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

########### Prior vs. Posterior plots ###########
# Load required packages
library(ggplot2)
library(tibble)
library(dplyr)
library(posterior)

# Set prior standard deviations
location_prior_sd <- 2   # for b_ parameters
disc_prior_sd <- 1       # for disc_ parameters

# Function for plotting
make_prior_posterior_plot <- function(param_name, posterior_vals, prior_vals, folder) {
  plot_data <- bind_rows(
    tibble(value = posterior_vals, distribution = "Posterior"),
    tibble(value = prior_vals, distribution = "Prior")
  )
  
  p <- ggplot(plot_data, aes(x = value, fill = distribution, color = distribution)) +
    geom_density(alpha = 0.5, linewidth = 0.7) +
    labs(x = "Estimate", y = "Density") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_fill_manual(values = c("Prior" = "gray", "Posterior" = "steelblue")) +
    scale_color_manual(values = c("Prior" = "gray40", "Posterior" = "steelblue"))
  
  # Create output folder if needed
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  ggsave(
    filename = file.path(folder, paste0("prior_posterior_", param_name, ".png")),
    plot = p, width = 6, height = 4, dpi = 300, bg = "white"
  )
}

# Wrap
plot_model_priors <- function(fit, model_name, location_pattern = "^b_", disc_pattern = "^disc_") {
  draws <- as_draws_df(fit)
  
  # Plot location parameters (e.g., b_)
  location_params <- names(draws)[grepl(location_pattern, names(draws))]
  for (param in location_params) {
    posterior_vals <- draws[[param]]
    prior_vals <- rnorm(10000, 0, location_prior_sd)
    make_prior_posterior_plot(param, posterior_vals, prior_vals,
                              folder = paste0("priors_", model_name, "_location"))
  }
  
  # Plot discrimination parameters (e.g., disc_)
  disc_params <- names(draws)[grepl(disc_pattern, names(draws))]
  for (param in disc_params) {
    posterior_vals <- draws[[param]]
    prior_vals <- rnorm(10000, 0, disc_prior_sd)
    make_prior_posterior_plot(param, posterior_vals, prior_vals,
                              folder = paste0("priors_", model_name, "_discrimination"))
  }
}

plot_model_priors(fit_prob, "fit_prob")
plot_model_priors(fit_irm, "fit_irm")


################ Bayes Factors ################
# --- RQ1: Is each foreign country > Sweden in "Yes" probability? ---
bf_sd <- function(fit, term, side = c("two", "greater", "less"), prior_sd) {
  side <- match.arg(side)
  draws <- as_draws_df(fit)[[paste0("b_", term)]]
  if (is.null(draws)) stop("Parameter not found: b_", term)
  
  # Prior from N(0, prior_sd^2)
  prior <- rnorm(200000, 0, prior_sd)
  
  if (side == "two") {
    post_h0 <- mean(abs(draws) < 1e-6)
    prior_h0 <- mean(abs(prior) < 1e-6)
  } else if (side == "greater") {
    post_h0 <- mean(draws <= 0)
    prior_h0 <- mean(prior <= 0)
  } else if (side == "less") {
    post_h0 <- mean(draws >= 0)
    prior_h0 <- mean(prior >= 0)
  }
  
  bf10 <- (1 - post_h0) / post_h0 * (prior_h0 / (1 - prior_h0))
  
  data.frame(
    term = term,
    side = side,
    BF10 = bf10,
    prior_prob_H0 = prior_h0,
    post_prob_H0 = post_h0
  )
}

# UK > Sweden
bf_sd(fit_prob_reduced_final, "CCJ", side = "greater", prior_sd = 2)
#   term    side     BF10 prior_prob_H0 post_prob_H0
# 1  CCJ greater 10.32887      0.499095   0.08797917

# Germany > Sweden
bf_sd(fit_prob_reduced_final, "CCK", side = "greater", prior_sd = 2)
#   term    side     BF10 prior_prob_H0 post_prob_H0
# 1  CCK greater 77.13923      0.499055      0.01275

# Brazil > Sweden
bf_sd(fit_prob_reduced_final, "CCL", side = "greater", prior_sd = 2)
#   term    side     BF10 prior_prob_H0 post_prob_H0
# 1  CCL greater 14.01584       0.50045   0.06670833

# Thailand > Sweden
bf_sd(fit_prob_reduced_final, "CCM", side = "greater", prior_sd = 2)
#   term    side     BF10 prior_prob_H0 post_prob_H0
# 1  CCM greater 8.929058      0.499465    0.1005208

#RQ1 v2
library(dplyr)
library(tidybayes)

# List of foreign countries
countries <- c("CCJ", "CCK", "CCL", "CCM")  # UK, Germany, Brazil, Thailand

# Function for one-sided BF: H1 = β > 0
bf_one_sided <- function(term) {
  bf_contrast(fit_prob_reduced_final,
              term = term,
              prior_sd_b = 2,
              side = "greater")
}

# Loop over countries
library(brms)

# Loop over foreign countries for H1: β_country > 0
rq1_results <- purrr::map_dfr(
  c("CCJ", "CCK", "CCL", "CCM"),
  ~ {
    h <- hypothesis(fit_prob_reduced_final,
                    paste0(.x, " > 0"),
                    class = "b")
    tibble(
      country = .x,
      estimate = h$hypothesis$Estimate,
      l95 = h$hypothesis$CI.Lower,
      u95 = h$hypothesis$CI.Upper,
      BF10 = h$hypothesis$Evid.Ratio
    )
  }
)

rq1_results
# # A tibble: 4 × 5
#   country estimate    l95   u95  BF10
# <chr>      <dbl>  <dbl> <dbl> <dbl>
# 1 CCJ         1.91 -0.468  4.06 10.4 
# 2 CCK         3.39  1.01   5.49 77.4 
# 3 CCL         2.30 -0.240  4.63 14.0 
# 4 CCM         1.97 -0.629  4.38  8.95

# Extract between-judge SDs
rq1_sd <- spread_draws(fit_prob_reduced_final,
                       sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM) %>%
  pivot_longer(cols = everything(), names_to = "term", values_to = "sd") %>%
  mutate(country = recode(term,
                          sd_ID__CCJ = "UK",
                          sd_ID__CCK = "Germany",
                          sd_ID__CCL = "Brazil",
                          sd_ID__CCM = "Thailand")) %>%
  group_by(country) %>%
  summarise(med_sd = median(sd),
            l95 = quantile(sd, 0.025),
            u95 = quantile(sd, 0.975))

rq1_results
rq1_sd

# --- RQ2: Contrast for cultural distance v.1---
bf_contrast_one_sided <- function(fit, terms, weights, direction = c("greater", "less"), prior_sd_b) {
  direction <- match.arg(direction)
  draws <- as_draws_df(fit)[, paste0("b_", terms)]
  if (any(is.na(draws))) stop("Some terms not found.")
  
  contrast_draws <- as.matrix(draws) %*% weights
  prior <- rnorm(200000, 0, sqrt(sum((weights * prior_sd_b)^2)))
  
  if (direction == "greater") {
    post_tail <- mean(contrast_draws > 0)
    prior_tail <- mean(prior > 0)
  } else {
    post_tail <- mean(contrast_draws < 0)
    prior_tail <- mean(prior < 0)
  }
  
  bf10 <- post_tail / (1 - post_tail) * ((1 - prior_tail) / prior_tail)
  
  data.frame(
    contrast = paste(terms, collapse = " + "),
    direction = direction,
    BF10 = bf10,
    post_tail_prob = post_tail,
    post_l95 = quantile(contrast_draws, 0.025),
    post_u95 = quantile(contrast_draws, 0.975)
  )
}

# Higher weights for Brazil & Thailand
w <- c(UK = 1, Germany = 1, Brazil = 2, Thailand = 2)
bf_contrast_one_sided(
  fit_prob_reduced_final,
  terms = c("CCJ","CCK","CCL","CCM"),
  weights = as.numeric(w[c("UK","Germany","Brazil","Thailand")]),
  direction = "greater",
  prior_sd_b = 2
)

#                   contrast direction     BF10 post_tail_prob post_l95 post_u95
# 2.5% CCJ + CCK + CCL + CCM   greater 167.8159      0.9940833 3.052935 24.44812

# --- RQ3: Does 'Yes' (D) lower loyalty ratings? ---
bf_sd(fit_dt_reduced_final, "D", side = "less", prior_sd = 2)  # Expecting negative effect
#   term side BF10 prior_prob_H0 post_prob_H0
# 1    D less  Inf      0.498955            0

#RQ2 v.2
# One-sided BF for a linear contrast (Savage–Dickey), where b_* are fixed effects
bf_contrast_diff <- function(fit, terms_pos, terms_neg, prior_sd_b = 2, side = c("greater","less","two")) {
  side <- match.arg(side)
  d <- as_draws_df(fit)
  pos <- rowSums(as.matrix(d[, paste0("b_", terms_pos), drop = FALSE]))
  neg <- rowSums(as.matrix(d[, paste0("b_", terms_neg), drop = FALSE]))
  contr <- pos - neg
  
  k <- length(terms_pos) + length(terms_neg)
  prior_sd_contrast <- sqrt(k) * prior_sd_b
  prior <- rnorm(200000, 0, prior_sd_contrast)
  
  if (side %in% c("greater","less")) {
    post_tail <- if (side == "greater") mean(contr > 0) else mean(contr < 0)
    prior_tail <- if (side == "greater") mean(prior > 0) else mean(prior < 0)
    bf10 <- post_tail/(1-post_tail) * ((1-prior_tail)/prior_tail)
    post_l95 <- quantile(contr, 0.025)
    post_u95 <- quantile(contr, 0.975)
  } else if (side == "two") {
    # Savage–Dickey for point null at 0
    eps <- 1e-6
    post_h0 <- mean(abs(contr) < eps)
    prior_h0 <- mean(abs(prior) < eps)
    bf10 <- (1 - post_h0)/post_h0 * (prior_h0/(1 - prior_h0))
    post_tail <- NA
    post_l95 <- quantile(contr, 0.025)
    post_u95 <- quantile(contr, 0.975)
  }
  
  data.frame(
    contrast = paste(paste(terms_pos, collapse="+"), "- (", paste(terms_neg, collapse="+"), ")"),
    side = side,
    BF10 = bf10,
    post_tail = post_tail,
    post_l95 = post_l95,
    post_u95 = post_u95
  )
}

# Regions
# Europe = mean(UK, Germany)  -> terms_pos/neg just compare sums; scaling cancels in one-sided tests.
# Asia   = Thailand
# S. America = Brazil

# Europe vs Asia (is Asia > Europe?)
bf_contrast_diff(fit_prob_reduced_final,
                 terms_pos = c("CCM"),                 # Thailand
                 terms_neg = c("CCJ","CCK"),           # UK + Germany
                 prior_sd_b = 2, side = "greater")
#               contrast    side       BF10  post_tail  post_l95 post_u95
# 2.5% CCM - ( CCJ+CCK ) greater 0.09740295 0.08879167 -8.088486 1.551573

# Europe vs S. America (is Brazil > Europe?)
bf_contrast_diff(fit_prob_reduced_final,
                 terms_pos = c("CCL"),                 # Brazil
                 terms_neg = c("CCJ","CCK"),
                 prior_sd_b = 2, side = "greater")
#               contrast    side      BF10 post_tail  post_l95 post_u95
# 2.5% CCL - ( CCJ+CCK ) greater 0.1259165   0.11125 -7.849373 1.918253

# Asia vs S. America (two-sided difference Thailand vs Brazil)
bf_contrast_diff(fit_prob_reduced_final,
                 terms_pos = c("CCM"),
                 terms_neg = c("CCL"),
                 prior_sd_b = 2, side = "two")
#           contrast side BF10 post_tail  post_l95 post_u95
# 2.5% CCM - ( CCL )  two  Inf        NA -4.153413 3.515645

region_sd_contrast <- function(fit) {
  s <- tidybayes::spread_draws(fit, sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM) %>%
    dplyr::mutate(
      sd_Europe = (sd_ID__CCJ + sd_ID__CCK) / 2,  # UK + Germany
      sd_Asia   =  sd_ID__CCM,                    # Thailand
      sd_SA     =  sd_ID__CCL,                    # Brazil
      diff_Asia_vs_Eur = sd_Asia - sd_Europe,
      diff_SA_vs_Eur   = sd_SA  - sd_Europe,
      diff_SA_vs_Asia  = sd_SA  - sd_Asia
    )
  
  data.frame(
    contrast = c("Asia > Europe","S. America > Europe","S. America ≠ Asia"),
    post_prob = c(mean(s$diff_Asia_vs_Eur > 0),
                  mean(s$diff_SA_vs_Eur > 0),
                  NA),
    l95 = c(quantile(s$diff_Asia_vs_Eur, 0.025),
            quantile(s$diff_SA_vs_Eur, 0.025),
            quantile(s$diff_SA_vs_Asia, 0.025)),
    u95 = c(quantile(s$diff_Asia_vs_Eur, 0.975),
            quantile(s$diff_SA_vs_Eur, 0.975),
            quantile(s$diff_SA_vs_Asia, 0.975)),
    two_sided_note = c("", "", "report 95% CrI; no one-sided prob here")
  )
}

# Run for Model 1 (probit D) and Model 2 (IRM ratings):
region_sd_contrast(fit_prob_reduced_final)
contrast post_prob        l95       u95
# 1       Asia > Europe 0.9314583  -1.849763 14.885871
# 2 S. America > Europe 0.8778333  -3.016074 13.777031
# 3   S. America ≠ Asia        NA -10.796789  8.498971
# two_sided_note
# 1                                       
# 2                                       
# 3 report 95% CrI; no one-sided prob here
region_sd_contrast(fit_dt_reduced_final)
# contrast post_prob         l95       u95
# 1       Asia > Europe 0.2307917 -0.36111257 0.1835809
# 2 S. America > Europe 0.8708750 -0.16316683 0.5191072
# 3   S. America ≠ Asia        NA -0.03285241 0.5939974
# two_sided_note
# 1                                       
# 2                                       
# 3 report 95% CrI; no one-sided prob here

#RQ2 v3
foreign_countries <- c("CCJ", "CCK", "CCL", "CCM")
contrast_pairs <- combn(foreign_countries, 2, simplify = FALSE)

rq2_results <- purrr::map_dfr(
  contrast_pairs,
  ~ {
    h <- hypothesis(fit_prob_reduced_final,
                    paste0(.x[1], " - ", .x[2], " = 0"),
                    class = "b")
    tibble(
      contrast = paste(.x, collapse = " - "),
      estimate = h$hypothesis$Estimate,
      l95 = h$hypothesis$CI.Lower,
      u95 = h$hypothesis$CI.Upper,
      BF10 = h$hypothesis$Evid.Ratio
    )
  }
)

rq2_results
# # A tibble: 6 × 5
# contrast  estimate   l95   u95  BF10
# <chr>        <dbl> <dbl> <dbl> <dbl>
# 1 CCJ - CCK  -1.48   -4.81  1.80    NA
# 2 CCJ - CCL  -0.385  -4.16  3.39    NA
# 3 CCJ - CCM  -0.0580 -3.87  3.73    NA
# 4 CCK - CCL   1.09   -2.69  4.85    NA
# 5 CCK - CCM   1.42   -2.36  5.20    NA
# 6 CCL - CCM   0.327  -3.52  4.15    NA

#RQ5
library(tidybayes)
library(dplyr)
library(tidyr)
library(ggplot2)

# Extract available SDs: intercept and country slopes
sds_long <- spread_draws(
  fit_dt_reduced_final,
  sd_ID__Intercept, sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM
) %>%
  pivot_longer(everything(), names_to = "param", values_to = "sd") %>%
  mutate(effect = recode(param,
                         sd_ID__Intercept = "Intercept (level noise)",
                         sd_ID__CCJ = "UK (pattern noise)",
                         sd_ID__CCK = "Germany (pattern noise)",
                         sd_ID__CCL = "Brazil (pattern noise)",
                         sd_ID__CCM = "Thailand (pattern noise)"
  ))

# Posterior summaries
sds_summary <- sds_long %>%
  group_by(effect) %>%
  summarise(
    mean = mean(sd),
    l95 = quantile(sd, 0.025),
    u95 = quantile(sd, 0.975),
    pr_gt_0p5 = mean(sd > 0.5),
    .groups = "drop"
  )

print(sds_summary)
# A tibble: 8 × 5
# effect                        mean        l95       u95 pr_gt_0p5
# <chr>                        <dbl>      <dbl>     <dbl>     <dbl>
# 1 .chain                       2.5      1           4       1      
# 2 .draw                    24000.    1201.      46800.      1      
# 3 .iteration                6000.     301.      11700.      1      
# 4 Brazil (pattern noise)       0.431    0.129       0.723   0.305  
# 5 Germany (pattern noise)      0.268    0.0236      0.557   0.0524 
# 6 Intercept (level noise)      0.754    0.510       1.05    0.980  
# 7 Thailand (pattern noise)     0.139    0.00598     0.374   0.00277
# 8 UK (pattern noise)           0.208    0.0107      0.494   0.0231 

#RQ4
# Posterior draws of the "distance effect" on variability
sds_draws <- spread_draws(
  fit_dt_reduced_final,
  sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM
) %>%
  mutate(
    high_dist = (sd_ID__CCL + sd_ID__CCM) / 2,  # Brazil, Thailand
    low_dist  = (sd_ID__CCJ + sd_ID__CCK) / 2,  # UK, Germany
    diff = high_dist - low_dist
  )

# Posterior probability high_dist > low_dist
mean(sds_draws$diff > 0)
# [1] 0.6425
quantile(sds_draws$diff, c(.025, .975))
#         2.5%      97.5% 
#   -0.2211671  0.3078155 

library(tidybayes)
library(dplyr)
library(tidyr)

# Compare pattern-noise SDs: (Brazil, Thailand) vs (UK, Germany)
rq4_variability_contrast <- function(fit, label) {
  draws <- spread_draws(
    fit,
    sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM
  ) %>%
    mutate(
      high_dist = (sd_ID__CCL + sd_ID__CCM) / 2,  # Brazil + Thailand
      low_dist  = (sd_ID__CCJ + sd_ID__CCK) / 2,  # UK + Germany
      diff = high_dist - low_dist
    )
  
  tibble::tibble(
    Model = label,
    post_prob_high_gt_low = mean(draws$diff > 0),
    diff_l95 = quantile(draws$diff, 0.025),
    diff_u95 = quantile(draws$diff, 0.975),
    # Optional: country-wise SD medians for context
    med_UK       = median(draws$sd_ID__CCJ),
    med_Germany  = median(draws$sd_ID__CCK),
    med_Brazil   = median(draws$sd_ID__CCL),
    med_Thailand = median(draws$sd_ID__CCM)
  )
}
res_rq4_irm   <- rq4_variability_contrast(fit_dt_reduced_final,   "Model 2: IRM (ratings)")
res_rq4_probit<- rq4_variability_contrast(fit_prob_reduced_final, "Model 1: Probit (D)")

print(rbind(res_rq4_probit, res_rq4_irm), n = Inf, width = Inf)

print(res_rq4_irm)

library(ggplot2)

plot_rq4_contrast <- function(fit, title = NULL) {
  d <- spread_draws(
    fit, sd_ID__CCJ, sd_ID__CCK, sd_ID__CCL, sd_ID__CCM
  ) %>%
    mutate(
      high_dist = (sd_ID__CCL + sd_ID__CCM) / 2,
      low_dist  = (sd_ID__CCJ + sd_ID__CCK) / 2,
      diff = high_dist - low_dist
    )
  
  pprob <- mean(d$diff > 0)
  ci    <- quantile(d$diff, c(.025, .975))
  
  ggplot(d, aes(x = diff)) +
    geom_density() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_segment(aes(x = ci[1], xend = ci[2], y = 0, yend = 0), linewidth = 1) +
    labs(
      x = "Difference in SD (High distance − Low distance)",
      y = "Posterior density",
      title = title %||% "RQ4: Distance vs. Disagreement (Pattern Noise)",
      subtitle = sprintf("Pr(diff > 0) = %.3f; 95%% CrI = [%.3f, %.3f]", pprob, ci[1], ci[2])
    ) +
    theme_minimal()
}

# Examples:
plot_rq4_contrast(fit_prob_reduced_final, "Model 1: Probit (D)")
plot_rq4_contrast(fit_dt_reduced_final,   "Model 2: IRM (ratings)")

#Plots
library(tidybayes)
library(dplyr)
library(tidyr)
library(ggplot2)

################# Figure 2 - Between-Judge Variability (95% CI) #################

make_fig2 <- function(fit_prob_reduced_final,
                      outfile = "figure2_between_judge_variability_95CI.png",
                      width = 7, height = 5, dpi = 300) {
  library(tidyverse)
  library(tidybayes)
  
  # Extract draws for the relevant SD parameters
  country_sds <- spread_draws(
    fit_prob_reduced_final,
    sd_ID__CCJ,
    sd_ID__CCK,
    sd_ID__CCL,
    sd_ID__CCM
  )
  
  # Reshape to long format
  country_sds_long <- country_sds %>%
    pivot_longer(cols = starts_with("sd_ID__"),
                 names_to = "country_raw", values_to = "sd") %>%
    mutate(country = recode(country_raw,
                            "sd_ID__CCJ" = "UK",
                            "sd_ID__CCK" = "Germany",
                            "sd_ID__CCL" = "Brazil",
                            "sd_ID__CCM" = "Thailand"),
           country = factor(country, levels = c("UK","Germany","Brazil","Thailand")))
  
  # Compute posterior medians + 95% CrI
  country_sds_ci <- country_sds_long %>%
    group_by(country) %>%
    summarise(
      median = median(sd),
      lower  = quantile(sd, 0.025),
      upper  = quantile(sd, 0.975),
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(country_sds_long, aes(x = country, y = sd, fill = country)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_pointrange(
      data = country_sds_ci,
      aes(y = median, ymin = lower, ymax = upper),
      color = "black", size = 0.5
    ) +
    labs(
      x = "Country (ordered by cultural distance from Sweden)",
      y = "Between-judge variability: SD of random slopes (95% CrI)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  invisible(p)
}

# Run
make_fig2(fit_prob_reduced_final)

################# Figure 3 — Effects (fixed) and variability (random slopes) #################

make_fig3_combined <- function(fit_irm,
                               outfile = "figure3_combined_AB_C.png",
                               width = 7, height = 10, dpi = 300) {
  library(tidyverse)
  library(patchwork)
  library(posterior)
  library(stringr)
  library(ggplot2)
  
  # ---------- helpers ----------
  summ_ci95 <- function(x) tibble(
    median = median(x),
    l95 = quantile(x, 0.025),
    u95 = quantile(x, 0.975)
  )
  
  expected_rating <- function(fit, newdata) {
    P <- posterior_epred(fit, newdata = newdata)  # draws x N x K
    K <- dim(P)[3]
    # expected category (1..K) for each draw/obs
    apply(sweep(P, 3, 1:K, `*`), c(1, 2), sum)   # draws x N
  }
  
  # ---------- posterior draws ----------
  draws <- as_draws_df(fit_irm)
  
  # ---------- (A) Fixed effects: b_CC* (relative to Sweden), violin + 95% CrI ----------
  fixed_draws <- draws %>%
    select(.draw, b_CCJ, b_CCK, b_CCL, b_CCM) %>%
    pivot_longer(cols = starts_with("b_"), names_to = "coef", values_to = "estimate") %>%
    mutate(
      country = recode(coef,
                       b_CCJ = "UK",
                       b_CCK = "Germany",
                       b_CCL = "Brazil",
                       b_CCM = "Thailand"),
      # order by theoretical/cultural distance: UK, Germany, Brazil, Thailand
      country = factor(country, levels = c("UK","Germany","Brazil","Thailand"))
    )
  
  fixed_ci <- fixed_draws %>%
    group_by(country) %>%
    summarise(summ_ci95(estimate), .groups = "drop")
  
  p_fixed <- ggplot(fixed_draws, aes(x = estimate, y = fct_rev(country), fill = country)) +
    geom_violin(alpha = 0.6, scale = "width", color = NA, trim = FALSE) +
    geom_pointrange(
      data = fixed_ci,
      aes(y = fct_rev(country), x = median, xmin = l95, xmax = u95),
      inherit.aes = FALSE,
      color = "black", size = 0.5
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Fixed effects (relative to Sweden)",
      x = "Effect on latent loyalty (probit scale, 95% CrI)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      # add left margin so the panel tag "A" doesn't hug the y-axis
      plot.margin = margin(t = 5, r = 5, b = 5, l = 18)
    )
  
  # ---------- (B) Between-judge variability: SD across judges of random slopes ----------
  rslope_vars <- names(draws) %>%
    str_subset("^r_ID\\[\\d+,(CCJ|CCK|CCL|CCM)\\]")
  
  random_slopes <- draws %>%
    select(.draw, all_of(rslope_vars)) %>%
    pivot_longer(-.draw, names_to = "term", values_to = "slope") %>%
    mutate(
      country_code = str_extract(term, "CCJ|CCK|CCL|CCM"),
      country = recode(country_code,
                       CCJ = "UK",
                       CCK = "Germany",
                       CCL = "Brazil",
                       CCM = "Thailand"),
      country = factor(country, levels = c("UK","Germany","Brazil","Thailand"))
    )
  
  sd_by_draw <- random_slopes %>%
    group_by(country, .draw) %>%
    summarise(sd_across_judges = sd(slope), .groups = "drop")
  
  sd_ci <- sd_by_draw %>%
    group_by(country) %>%
    summarise(summ_ci95(sd_across_judges), .groups = "drop")
  
  p_var <- ggplot(sd_by_draw, aes(x = country, y = sd_across_judges, fill = country)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_pointrange(
      data = sd_ci,
      aes(y = median, ymin = l95, ymax = u95),
      color = "black", size = 0.5
    ) +
    labs(
      title = "Between-judge variability (random slopes SD)",
      x = "Country (ordered by cultural distance from Sweden)",
      y = "SD across judges of country-specific random slopes (95% CrI)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 18)
    )
  
  # ---------- (C) Expected ratings by D (levels) ----------
  base_dat <- model.frame(fit_irm) %>% as.data.frame()
  # preserve D's type
  if (is.factor(base_dat$D)) {
    new0 <- base_dat; new0$D <- factor(0, levels = levels(base_dat$D))
    new1 <- base_dat; new1$D <- factor(1, levels = levels(base_dat$D))
  } else {
    new0 <- base_dat; new0$D <- 0
    new1 <- base_dat; new1$D <- 1
  }
  
  E0 <- expected_rating(fit_irm, new0)   # draws x N
  E1 <- expected_rating(fit_irm, new1)   # draws x N
  
  m0 <- rowMeans(E0)
  m1 <- rowMeans(E1)
  
  levels_df <- tibble(
    rating = c(as.numeric(m0), as.numeric(m1)),
    D = factor(
      rep(c("D = 0 (not judged relevant)", "D = 1 (judged relevant)"),
          each = length(m0)),
      levels = c("D = 0 (not judged relevant)", "D = 1 (judged relevant)")
    )
  )
  
  levels_summ <- levels_df %>%
    group_by(D) %>%
    summarise(
      median = median(rating),
      l95 = quantile(rating, 0.025),
      u95 = quantile(rating, 0.975),
      .groups = "drop"
    )
  
  p_levels <- ggplot(levels_df, aes(x = D, y = rating, fill = D)) +
    geom_violin(alpha = 0.6, trim = FALSE, color = NA) +
    geom_pointrange(
      data = levels_summ,
      aes(x = D, y = median, ymin = l95, ymax = u95),
      color = "black", size = 0.6, inherit.aes = FALSE
    ) +
    labs(
      title = "Expected loyalty rating by relevance judgment",
      x = NULL,
      y = "Expected rating (1–7), 95% CrI"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 18)
    )
  
  # ---------- combine (A/B stacked; C full width) ----------
  p <- (p_fixed / p_var) / p_levels +
    plot_annotation(tag_levels = "A") &
    theme(
      # add a bit of left margin to *all* tags
      plot.tag = element_text(margin = margin(l = 8))
    )
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi, units = "in")
  p  # return so it prints
}

# Usage:
# make_fig3_combined(fit_dt_reduced_final)

################# Figure 4 — Participant-level intercepts (95% CrI) #################

library(tidyverse)
library(posterior)
library(stringr)

# 1) Posterior draws (final model)
draws <- as_draws_df(fit_dt_reduced_final)

# 2) Extract participant-level intercept draws: r_ID[<id>,Intercept]
intercepts_df <- draws %>%
  select(starts_with("r_ID[")) %>%
  select(matches(",Intercept]")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "param",
    values_to = "intercept"
  ) %>%
  mutate(
    ID = str_extract(param, "(?<=r_ID\\[)\\d+"),
    ID = as.integer(ID)
  )

# 3) Summarize posterior mean and 95% CrI per participant
summary_df <- intercepts_df %>%
  group_by(ID) %>%
  summarise(
    mean  = mean(intercept),
    lower = quantile(intercept, 0.025),
    upper = quantile(intercept, 0.975),
    .groups = "drop"
  ) %>%
  arrange(mean) %>%                      # order by harshness/leniency
  mutate(ID = factor(ID, levels = ID))  # preserve order on x-axis

# 4) Plot
ggplot(summary_df, aes(x = ID, y = mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Participant-Level Intercepts (General Leniency vs. Harshness)",
    subtitle = "Higher = harsher ratings; 95% credible intervals",
    x = "Participant (ranked by posterior mean intercept)",
    y = "Estimated intercept (latent probit scale)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.6
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
#   [1] priorsense_1.1.1.9000 posterior_1.6.1       bayesplot_1.11.1     
# [4] waldo_0.6.2           compare_0.2-6         readxl_1.4.5         
# [7] future_1.34.0         lubridate_1.9.4       forcats_1.0.0        
# [10] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.4          
# [13] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
# [16] ggplot2_3.5.2         tidyverse_2.0.0       brms_2.22.0          
# [19] Rcpp_1.0.14          
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6         tensorA_0.36.2.1     QuickJSR_1.6.0      
# [4] processx_3.8.6       inline_0.3.21        lattice_0.22-6      
# [7] callr_3.7.6          tzdb_0.4.0           ps_1.9.0            
# [10] vctrs_0.6.5          tools_4.4.1          generics_0.1.3      
# [13] curl_6.2.1           stats4_4.4.1         parallel_4.4.1      
# [16] pkgconfig_2.0.3      Matrix_1.7-2         checkmate_2.3.2     
# [19] RColorBrewer_1.1-3   distributional_0.5.0 RcppParallel_5.1.10 
# [22] lifecycle_1.0.4      compiler_4.4.1       farver_2.1.2        
# [25] textshaping_1.0.0    Brobdingnag_1.2-9    codetools_0.2-20    
# [28] crayon_1.5.3         pillar_1.10.2        StanHeaders_2.32.10 
# [31] bridgesampling_1.1-2 abind_1.4-8          nlme_3.1-167        
# [34] parallelly_1.42.0    rstan_2.32.6         tidyselect_1.2.1    
# [37] digest_0.6.37        mvtnorm_1.3-3        stringi_1.8.4       
# [40] diffobj_0.3.5        reshape2_1.4.4       listenv_0.9.1       
# [43] labeling_0.4.3       grid_4.4.1           colorspace_2.1-1    
# [46] cli_3.6.5            magrittr_2.0.3       loo_2.8.0           
# [49] utf8_1.2.4           dichromat_2.0-0.1    pkgbuild_1.4.6      
# [52] future.apply_1.11.3  withr_3.0.2          scales_1.4.0        
# [55] backports_1.5.0      estimability_1.5.1   timechange_0.3.0    
# [58] emmeans_1.10.7       matrixStats_1.5.0    globals_0.16.3      
# [61] gridExtra_2.3        cellranger_1.1.0     ragg_1.3.3          
# [64] hms_1.1.3            coda_0.19-4.1        V8_6.0.1            
# [67] rstantools_2.4.0     rlang_1.1.6          xtable_1.8-4        
# [70] glue_1.8.0           jsonlite_1.9.1       rstudioapi_0.17.1   
# [73] R6_2.6.1             plyr_1.8.9           systemfonts_1.2.1   