#Survey 1, version 9
library(tidyr) 
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(tidyverse)
library(tidybayes)
library(performance)
library(Rcpp)
library(forcats)
library(stringr)
library(purrr)
library(brms)
library(bayesplot)
library(marginaleffects)
library(collapse)


#Building the model
#Priors on the thresholds from Solomon Kurz blog
#Credit due to Solomon Kurz
#See https://solomonkurz.netlify.app/blog/2021-12-29-notes-on-the-bayesian-cumulative-probit/

#See S Kurz blog on cumulative probit; 
#https://solomonkurz.netlify.app/blog/2021-12-29-notes-on-the-bayesian-cumulative-probit/
tibble(rating = 1:7) %>% 
  mutate(proportion = 1/7) %>% 
  mutate(cumulative_proportion = cumsum(proportion)) %>% 
  mutate(right_hand_threshold = qnorm(cumulative_proportion))


# A tibble: 7 × 4
# rating proportion cumulative_proportion right_hand_threshold
#      <int>      <dbl>                 <dbl>                <dbl>
# 1      1      0.143                 0.143               -1.07 
# 2      2      0.143                 0.286               -0.566
# 3      3      0.143                 0.429               -0.180
# 4      4      0.143                 0.571                0.180
# 5      5      0.143                 0.714                0.566
# 6      6      0.143                 0.857                1.07 
# 7      7      0.143                 1                    8.13
# > 

#Visualising the latent scale
tibble(z = seq(from = -3.75, to = 3.75, length.out = 1e3)) %>% 
  mutate(d = dnorm(x = z, mean = 0, sd = 1)) %>% 
  
  ggplot(aes(x = z, y = d)) +
  geom_line(color = "blue") +
  geom_area(aes(fill = z >= qnorm(p = 0 / 7)), alpha = 1/4) + 
  geom_area(aes(fill = z >= qnorm(p = 1 / 7)), alpha = 1/4) + 
  geom_area(aes(fill = z >= qnorm(p = 2 / 7)), alpha = 1/4) + 
  geom_area(aes(fill = z >= qnorm(p = 3 / 7)), alpha = 1/4) + 
  geom_area(aes(fill = z >= qnorm(p = 4 / 7)), alpha = 1/4) + 
  geom_area(aes(fill = z >= qnorm(p = 5 / 7)), alpha = 1/4) +
  geom_area(aes(fill = z >= qnorm(p = 6 / 7)), alpha = 1/4) +
  geom_vline(xintercept = qnorm(p = 1:6 / 7), linetype = 3) +
  scale_fill_manual(values = c("transparent", "blue"), breaks = NULL) +
  scale_x_continuous(expression(Phi), breaks = -3:3,
                     sec.axis = dup_axis(
                       name = NULL,
                       breaks = qnorm(p = 1:6 / 7),
                       labels = parse(text = str_c("tau[", 1:6, "]"))
                     )) +
  scale_y_continuous(NULL, breaks = NULL) +
  coord_cartesian(xlim = c(-3.25, 3.25)) +
  labs(title = expression(Phi*" for an evenly-distributed 7-point ordinal variable"),
       subtitle = "Each shaded section, defined by the thresholds, has the same probability mass.")


#Treatment contrasts and manually setting the reference category
# Set 'CM' as the reference category for C
result$C <- relevel(factor(result$C), ref = "CM")

# Set 'TC' as the reference category for T
result$T <- relevel(factor(result$T), ref = "TC")

# Set male as the reference category for G
result$G <- relevel(factor(result$G), ref = "M")

#Set male as the reference category for A
result$A <- relevel(factor(result$A), ref = "M")


#Check contrasts
contrasts(result$C)
contrasts(result$T)

levels(result$T)

#Set cores
nCores = parallel::detectCores()

#Model 1 - Exluding the discrimination parameter
fit1 <- brm(
  formula = bf(Y ~ 1 + C*T + G*A + (1 | id) + (1 | Item)),
  data = result, 
  family = cumulative("probit"),
  prior = c(prior(normal(-1.07, 1), class = Intercept, coef = 1),
            prior(normal(-0.57, 1), class = Intercept, coef = 2),
            prior(normal(-0.18, 1), class = Intercept, coef = 3),
            prior(normal( 0.18, 1), class = Intercept, coef = 4),
            prior(normal( 0.57, 1), class = Intercept, coef = 5),
            prior(normal( 1.07, 1), class = Intercept, coef = 6),
            prior(normal( 0, 0.25 ), class = b),
            prior(exponential(1) , class = sd)),
  chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99),
  init_r = 0.2
)

#Trace and density plots
plot(fit1) #Looks healthy

#Posterior predictive check
pp_check(fit1) #Looks healthy

fit1 = add_criterion(fit1, criterion = "loo")

#Posterior summary
summary(fit1)
# Family: cumulative 
# Links: mu = probit; disc = identity 
# Formula: Y ~ 1 + C * T + G * A + (1 | id) + (1 | Item) 
# Data: result (Number of observations: 450) 
# Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
# total post-warmup draws = 4000
# 
# Multilevel Hyperparameters:
#   ~id (Number of levels: 90) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.08      0.06     0.00     0.22 1.00     1778     1424
# 
# ~Item (Number of levels: 15) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.83      0.23     0.45     1.34 1.00     1725     2279
# 
# Regression Coefficients:
#              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[1]    -2.36      0.24    -2.83    -1.88 1.00     2608     3033
# Intercept[2]    -1.38      0.21    -1.78    -0.94 1.00     2540     2869
# Intercept[3]    -0.64      0.21    -1.03    -0.20 1.00     2495     2806
# Intercept[4]     0.15      0.20    -0.24     0.59 1.00     2463     3244
# Intercept[5]     0.87      0.21     0.48     1.32 1.00     2509     3036
# Intercept[6]     1.65      0.22     1.23     2.10 1.00     2592     2985
# C1               0.21      0.23    -0.26     0.66 1.00     2517     2490
# C2              -0.36      0.25    -0.85     0.11 1.00     2650     2879
# C3              -0.17      0.21    -0.57     0.24 1.00     2806     2992
# C4              -0.01      0.21    -0.42     0.38 1.00     3129     2598
# T1              -0.00      0.19    -0.37     0.36 1.00     2731     2742
# T2              -0.00      0.19    -0.37     0.37 1.00     2669     2631
# G1              -0.03      0.05    -0.13     0.07 1.00     6848     3051
# A1               0.01      0.05    -0.09     0.10 1.00     6435     3083
# C1:T1            0.01      0.22    -0.41     0.43 1.00     3850     2912
# C2:T1            0.00      0.22    -0.42     0.44 1.00     3941     2941
# C3:T1            0.01      0.22    -0.41     0.43 1.00     4041     3480
# C4:T1           -0.01      0.22    -0.45     0.41 1.00     3848     2802
# C1:T2           -0.00      0.22    -0.43     0.42 1.00     3922     3038
# C2:T2           -0.03      0.22    -0.46     0.38 1.00     3943     3302
# C3:T2           -0.02      0.22    -0.44     0.42 1.00     4287     2654
# C4:T2           -0.02      0.22    -0.44     0.41 1.00     4200     2350
# G1:A1           -0.02      0.05    -0.12     0.08 1.00     6742     2985
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# disc     1.00      0.00     1.00     1.00   NA       NA       NA
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

#Priors that include the discrimination parameter
#The variance component of the latent variable Y

fit2priors <- c(prior(normal(-1.07, 1), class = Intercept, coef = 1),
                prior(normal(-0.57, 1), class = Intercept, coef = 2),
                prior(normal(-0.18, 1), class = Intercept, coef = 3),
                prior(normal( 0.18, 1), class = Intercept, coef = 4),
                prior(normal( 0.57, 1), class = Intercept, coef = 5),
                prior(normal( 1.07, 1), class = Intercept, coef = 6),
                prior(normal( 0, 0.25 ), class = b),
                prior(exponential(1) , class = sd),
                prior(normal(0, log(2) / 2), class = b, dpar = disc))

#Fit 2 worked as well, all chains converged

fit2 <- brm(
  formula = bf(Y ~ 1 + C*T + G*A + (1 | id) + (1 | Item)) +
            lf(disc ~ 0 + C + T, cmc = FALSE),
  data = result, 
  family = cumulative("probit"),
  prior = fit2priors,
  warmup = 3000, iter = 9000,
  chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99),
  init_r = 0.2
)

fit2 = add_criterion(fit2, criterion = "loo")

loo(fit2)

# Computed from 52000 by 450 log-likelihood matrix.
# 
#         Estimate   SE
# elpd_loo   -684.0 16.1
# p_loo        31.4  2.6
# looic      1368.1 32.2
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.9, 1.6]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

#Trace and density plots
plot(fit2) #Looks healthy

#Posterior predictive check
pp_check(fit2) #Looks healthy

summary(fit2)
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: Y ~ 1 + C * T + G * A + (1 | id) + (1 | Item) 
# disc ~ 0 + C + T
# Data: result (Number of observations: 450) 
# Draws: 4 chains, each with iter = 9000; warmup = 3000; thin = 1;
# total post-warmup draws = 24000
# 
# Multilevel Hyperparameters:
#   ~id (Number of levels: 90) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.06      0.04     0.00     0.16 1.00    11869    12484
# 
# ~Item (Number of levels: 15) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.86      0.23     0.50     1.38 1.00     9048    12827
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[1]    -1.97      0.32    -2.61    -1.36 1.00    12369    16188
# Intercept[2]    -1.14      0.29    -1.72    -0.59 1.00    12984    16302
# Intercept[3]    -0.54      0.27    -1.09    -0.02 1.00    14284    17080
# Intercept[4]     0.11      0.26    -0.41     0.63 1.00    15731    18096
# Intercept[5]     0.80      0.27     0.28     1.35 1.00    15338    17785
# Intercept[6]     1.71      0.32     1.11     2.36 1.00    12784    16620
# CCJ             -0.09      0.23    -0.54     0.35 1.00    27059    18931
# CCK              0.04      0.23    -0.40     0.49 1.00    28254    18633
# CCL              0.13      0.23    -0.33     0.59 1.00    29799    18612
# CCM              0.26      0.24    -0.22     0.74 1.00    26212    18593
# TTB             -0.01      0.22    -0.45     0.42 1.00    26898    19530
# TTC             -0.01      0.22    -0.44     0.43 1.00    23581    18968
# GF               0.11      0.09    -0.06     0.29 1.00    28356    19128
# AF              -0.07      0.09    -0.26     0.11 1.00    27466    19719
# CCJ:TTB         -0.03      0.24    -0.50     0.45 1.00    42952    18498
# CCK:TTB          0.02      0.24    -0.45     0.50 1.00    41585    18423
# CCL:TTB          0.03      0.24    -0.44     0.50 1.00    40424    18245
# CCM:TTB          0.07      0.24    -0.41     0.54 1.00    43764    17966
# CCJ:TTC         -0.03      0.24    -0.51     0.45 1.00    38077    18302
# CCK:TTC          0.00      0.24    -0.47     0.48 1.00    40432    17479
# CCL:TTC          0.04      0.24    -0.44     0.51 1.00    40664    18333
# CCM:TTC          0.09      0.24    -0.39     0.56 1.00    42858    17483
# GF:AF           -0.04      0.12    -0.28     0.20 1.00    23981    20065
# disc_CCJ         0.12      0.12    -0.11     0.36 1.00    13916    17147
# disc_CCK        -0.03      0.13    -0.28     0.22 1.00    12866    17768
# disc_CCL        -0.44      0.13    -0.69    -0.19 1.00    13534    17529
# disc_CCM        -0.69      0.14    -0.97    -0.42 1.00    15126    17966
# disc_TTB         0.54      0.10     0.35     0.72 1.00    23960    20189
# disc_TTC         0.56      0.10     0.36     0.76 1.00    25194    20873
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

#Model with random effects for repeated measures of C in id
#I.e. estimate an "intercept" given id
#An effect specific for each participant for Country
#Adaptive pooling of item and id


fit3priors = c(
  prior(normal(-1.07, 1), class = Intercept, coef = 1),
  prior(normal(-0.57, 1), class = Intercept, coef = 2),
  prior(normal(-0.18, 1), class = Intercept, coef = 3),
  prior(normal(0.18, 1), class = Intercept, coef = 4),
  prior(normal(0.57, 1), class = Intercept, coef = 5),
  prior(normal(1.07, 1), class = Intercept, coef = 6),
  
  prior(normal(0, 0.1), class = b),  # Shrinking fixed effects towards zero
  
  prior(exponential(2), class = sd, group = id),  # Stronger shrinkage for random effects
  prior(exponential(2), class = sd, group = Item),
  
  prior(lkj_corr_cholesky(4), class = cor),  # Stronger shrinkage for correlations
  
  prior(normal(0, 0.1), class = b, dpar = disc)  # Shrinking discrimination parameter
)

fit3 <- brm(
  formula = bf(Y ~ 1 + C*T + G*A + (1 + C | id) + (1 | Item)) +
    lf(disc ~ 0 + C + T, cmc = FALSE),
  data = result, 
  family = cumulative("probit"),
  prior = fit3priors,
  warmup = 3000, iter = 9000,
  chains = 4, cores = 8,
  seed = 5476,
  control = list(adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  init_r = 0.2
)

#Fit3 worked, but see loo values

fit3 = add_criterion(fit3, criterion = "loo")

loo(fit3)
# Computed from 24000 by 450 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -677.1 14.0
# p_loo       124.5  7.8
# looic      1354.2 27.9
# ------
#   MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.3]).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     367   81.6%   482     
# (0.7, 1]   (bad)       82   18.2%   <NA>    
#   (1, Inf)   (very bad)   1    0.2%   <NA>    
#   See help('pareto-k-diagnostic') for details.


# Run LOO with moment matching
# Run reloo for the problematic observations
loo_result <- loo(fit3)

# Increase the maximum size for globals to 2 GB
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB in bytes

# Now rerun reloo
reloo_result <- reloo(loo_result, fit3, cores = 6) #Resolved the issue

print(reloo_result)
# Computed from 24000 by 450 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -684.8 14.7
# p_loo       132.2  9.0
# looic      1369.6 29.3
# ------
#   MCSE of elpd_loo is 0.6.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.3]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

co <- loo_compare(fit1, fit2, fit3)

print(co)
#    elpd_diff se_diff
# fit3   0.0       0.0  
# fit2  -7.0       6.0  
# fit1 -49.3       8.3  


fit3
# Family: cumulative 
# Links: mu = probit; disc = log 
# Formula: Y ~ 1 + C * T + G * A + (1 + C | id) + (1 | Item) 
# disc ~ 0 + C + T
# Data: result (Number of observations: 450) 
# Draws: 4 chains, each with iter = 9000; warmup = 3000; thin = 1;
# total post-warmup draws = 24000
# 
# Multilevel Hyperparameters:
#   ~id (Number of levels: 90) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)          0.07      0.06     0.00     0.20 1.00    12123    12139
# sd(CCJ)                0.13      0.11     0.00     0.40 1.00     8773    11631
# sd(CCK)                0.21      0.16     0.01     0.58 1.00     5025     8758
# sd(CCL)                0.92      0.24     0.38     1.37 1.00     2803     2405
# sd(CCM)                1.43      0.28     0.91     2.00 1.00     3705     5124
# cor(Intercept,CCJ)    -0.04      0.29    -0.59     0.53 1.00    31801    18676
# cor(Intercept,CCK)    -0.03      0.29    -0.57     0.53 1.00    20523    17115
# cor(CCJ,CCK)          -0.04      0.29    -0.59     0.53 1.00    17193    16040
# cor(Intercept,CCL)    -0.07      0.28    -0.60     0.49 1.00     3575     6682
# cor(CCJ,CCL)          -0.07      0.28    -0.58     0.48 1.00     3634     8487
# cor(CCK,CCL)           0.15      0.27    -0.42     0.63 1.00     4231     8718
# cor(Intercept,CCM)    -0.11      0.29    -0.63     0.47 1.00     1806     4821
# cor(CCJ,CCM)          -0.04      0.28    -0.56     0.51 1.00     2264     5527
# cor(CCK,CCM)          -0.07      0.27    -0.57     0.48 1.00     2892     6434
# cor(CCL,CCM)           0.22      0.16    -0.10     0.53 1.00     7877    11782
# 
# ~Item (Number of levels: 15) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.02      0.22     0.68     1.52 1.00     6950    10871
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept[1]    -2.40      0.29    -2.97    -1.84 1.00     6782    11179
# Intercept[2]    -1.44      0.26    -1.96    -0.93 1.00     6516    11246
# Intercept[3]    -0.71      0.25    -1.21    -0.22 1.00     6717    11975
# Intercept[4]     0.10      0.25    -0.39     0.59 1.00     7079    12165
# Intercept[5]     0.92      0.26     0.43     1.44 1.00     6983    11723
# Intercept[6]     1.97      0.30     1.40     2.59 1.00     6309    11585
# CCJ             -0.02      0.10    -0.21     0.17 1.00    35774    18821
# CCK              0.01      0.10    -0.19     0.20 1.00    37425    18919
# CCL              0.02      0.10    -0.17     0.21 1.00    37784    18102
# CCM              0.04      0.10    -0.16     0.24 1.00    41318    16592
# TTB             -0.00      0.10    -0.19     0.19 1.00    36585    18347
# TTC             -0.00      0.10    -0.19     0.20 1.00    36343    18601
# GF               0.06      0.07    -0.08     0.20 1.00    35740    19464
# AF              -0.05      0.07    -0.20     0.09 1.00    38092    18779
# CCJ:TTB         -0.01      0.10    -0.20     0.19 1.00    44422    16690
# CCK:TTB          0.00      0.10    -0.19     0.20 1.00    44493    18596
# CCL:TTB          0.00      0.10    -0.19     0.20 1.00    45034    16484
# CCM:TTB          0.01      0.10    -0.18     0.21 1.00    42018    18669
# CCJ:TTC         -0.01      0.10    -0.20     0.19 1.00    40570    18740
# CCK:TTC          0.00      0.10    -0.19     0.19 1.00    46236    17922
# CCL:TTC          0.01      0.10    -0.19     0.20 1.00    42837    17646
# CCM:TTC          0.01      0.10    -0.18     0.21 1.00    42978    17580
# GF:AF           -0.00      0.08    -0.16     0.15 1.00    34966    19989
# disc_CCJ         0.06      0.08    -0.09     0.21 1.00    18051    18618
# disc_CCK        -0.01      0.08    -0.16     0.15 1.00    13982    17144
# disc_CCL        -0.02      0.10    -0.23     0.18 1.00     6270     8399
# disc_CCM        -0.01      0.10    -0.22     0.19 1.00    10258    10403
# disc_TTB         0.19      0.07     0.04     0.33 1.00    19552    18744
# disc_TTC         0.37      0.09     0.20     0.53 1.00    10165    13538
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

#Produce bars plots posterior predictive check
yrep_charUL <- posterior_predict(fit3)
yrep_int <- sapply(data.frame(yrep_charUL, stringsAsFactors = TRUE), as.integer)
y_int <- as.integer(result$Y)

ppc_bars_grouped(
  y = y_int,
  yrep = yrep_int,
  group = interaction(result$C, result$T),
  freq=FALSE,
  prob = 0.5,
  fatten = 1,
  size = 1.5
)


#Preparing conditional effects for type of relationship
cond_TA <- data.frame(T = "TA",
                      cond__ = "TA")

cond_TB <- data.frame(T = "TB",
                      cond__ = "TB")

cond_TC <- data.frame(T = "TC",
                      cond__ = "TC")

#Plotting conditional effects
plot(conditional_effects(fit3, effects = "C",
                         conditions = cond_TA,
                         categorical = TRUE))

plot(conditional_effects(fit3, effects = "C",
                         conditions = cond_TB,
                         categorical = TRUE))

plot(conditional_effects(fit3, effects = "C",
                         conditions = cond_TC,
                         categorical = TRUE))

#Plotting predicitions
plot_predictions(fit3, condition = list(
  "C",
  "T"),
  type = "link")


############################Logistic regression################################

logprior =  c(prior(normal(0, 1), class = "Intercept"),
              
              prior(normal(0, 1), class = "b"),
              
              prior(exponential(1) , class = sd))

fitLog1 <- brm(
  formula = D ~ 1 + C*T + G*A + (1 | id) + (1 | Item),
  data = Dresult, 
  family = bernoulli(link = "logit"),
  prior = logprior,
  warmup = 3000, iter = 16000,
  chains = 4, cores = nCores,
  seed = 5476,
  control = list(adapt_delta = .99),
  save_pars = save_pars(all = TRUE),
  file = 'fitLog1'
)

fitLog1 = add_criterion(fitLog1, criterion = "loo")

loo(fitLog1)
# Computed from 52000 by 450 log-likelihood matrix.
# 
#          Estimate   SE
# elpd_loo   -202.7 12.2
# p_loo        23.8  2.3
# looic       405.5 24.4
# ------
#   MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.4]).
# 
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.


fitLog1
# Family: bernoulli 
# Links: mu = logit 
# Formula: D ~ 1 + C * T + G * A + (1 | id) + (1 | Item) 
# Data: Dresult (Number of observations: 450) 
# Draws: 4 chains, each with iter = 16000; warmup = 3000; thin = 1;
# total post-warmup draws = 52000
# 
# Multilevel Hyperparameters:
#   ~id (Number of levels: 90) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.34      0.24     0.01     0.86 1.00    11751    19545
# 
# ~Item (Number of levels: 15) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.54      0.45     0.02     1.67 1.00     7182    12143
# 
# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept    -1.30      0.33    -1.95    -0.63 1.00    28598    25695
# C1            2.02      0.45     0.95     2.74 1.00    15645    15128
# C2            1.34      0.40     0.41     2.04 1.00    20952    15895
# C3           -0.77      0.47    -1.63     0.24 1.00    25191    19936
# C4           -0.85      0.46    -1.70     0.16 1.00    24015    20915
# T1            0.52      0.31    -0.12     1.13 1.00    28496    22802
# T2           -0.17      0.35    -0.87     0.52 1.00    28264    27235
# GM           -0.48      0.35    -1.17     0.19 1.00    38710    38959
# AM           -0.00      0.33    -0.65     0.65 1.00    37995    38980
# C1:T1        -0.55      0.47    -1.43     0.46 1.00    30775    25249
# C2:T1        -0.33      0.46    -1.25     0.59 1.00    31636    26586
# C3:T1         0.13      0.51    -0.93     1.14 1.00    36397    28686
# C4:T1        -0.43      0.54    -1.47     0.67 1.00    38304    28618
# C1:T2         0.46      0.49    -0.57     1.42 1.00    33387    26401
# C2:T2        -0.06      0.48    -0.99     0.95 1.00    34414    26738
# C3:T2         0.19      0.56    -0.89     1.32 1.00    41384    32868
# C4:T2         0.68      0.55    -0.45     1.73 1.00    34924    28784
# GM:AM        -0.02      0.45    -0.91     0.86 1.00    32238    37008
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

plot(fitLog1) #Looks healthy

pp_check(fitLog1) #Looks healthy

#Plotting conditional effects
conditional_effects(fitLog1, effects = "C",
                    conditions = cond_TA)

conditional_effects(fitLog1, effects = "C",
                    condtions = cond_TB)

conditional_effects(fitLog1, effects = "C",
                    conditions = cond_TC)