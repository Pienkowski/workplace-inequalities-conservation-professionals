#################################
######### Data analysis #########
#################################

# The following explores the data used in the manuscript: Psychological distress and workplace risk inequalities among conservation professionals

### This script contains the following steps ###
# 1) Set the environment and load data 
# 2) Estimate K6 latent variable
# 3) Creating informative priors from previous study
# 4) Perform the model diagnostics from the WAMBS check list (Depaoli & van de Schoot 2017):
# 4.1) 'Do you understand the priors?'
# 4.2) 'Does the trace-plot exhibit convergence?' (After running the model.)
# 4.3) 'Does convergence remain after doubling the number of iterations?'
# 4.4) 'Does the histogram have enough information?'
# 4.5) 'Do the chains exhibit a strong degree of autocorrelation?'
# 4.6) 'Does the posterior distribution make substantive sense?'
# 4.7) 'Do different specifications of the multivariate variance priors influence the results?'
# 4.8) 'Is there a notable effect of the prior when compared with noninformative priors?'
# 4.9) 'Are the results stable from a sensitivity analysis?'
# 4.10) 'Is the Bayesian way of interpreting and reporting model results used?' 
# 5) Run the model on the ten imputed datasets
# 6) Extract and pool the model results 
# 7) Run the model on the ten imputed datasets with default weakly informative priors 
# 8) Run the model on the ten imputed datasets using the original ERI
# 9) Run the model on the ten imputed datasets replacing position with organisation

######### 1) Set the environment and load data #########

# Load the data
load("mice.imp.DEEP.Rdata")

# Load packages 
library(ggplot2)
library(mirt)
library(blavaan)
library(psych)
library(coda)
library(plyr)

# Display up to 10 decimal places 
options("scipen" = 10) 

######### 2) Estimate K6 latent variable ######### 
K6_col_name_num <- paste0(c("K10_2","K10_4", "K10_5", "K10_8", "K10_9","K10_10"),"_num")
SS_col_name_num <- paste0(c("SS1", "SS2","SS3"),"_num")

### Kessler-6 ###
# Convert to numeric 
for (i in seq_along(1:length(mice.imp.DEEP))){
  mice.imp.DEEP[[i]][K6_col_name_num] <- apply(mice.imp.DEEP[[i]][K6_col_name_num], 2, as.numeric)
}

# Parallel analysis (with training data) - suggests the extraction of three factors. This is probably because some of the items are paired. 
fa.parallel(mice.imp.DEEP[[1]][K6_col_name_num] , cor = "poly", fm="wls", fa="fa",   main = "Parallel analysis")

### Graded response model ###
K6_GRM <- mirt(mice.imp.DEEP[[1]][K6_col_name_num], model = 1, "graded")

# Inspect fit statistics 
M2(K6_GRM, type = "C2", calcNULL = FALSE)

# Item information curve (or 'information and trace lines') and information and SE
plot(K6_GRM, type = 'infotrace', facet_items= F)
plot(K6_GRM, type = 'infoSE', facet_items= F) 

# Item characteristic curve (or 'item scoring traceline plots')
plot(K6_GRM, type = 'trace') 

# Extract ten sets of plausible values and inspect them 
K6_PV <- fscores(K6_GRM, plausible.draws = 10)
K6_PV_combined <- data.frame(Iteration = "1", K6 = K6_PV[[1]] )
for (i in seq_along(2:length(K6_PV))){
  K6_PV_combined <- rbind(K6_PV_combined, data.frame(Iteration = as.factor(i+1), K6 = K6_PV[[i+1]] ))
}

# Plot
ggplot(K6_PV_combined, aes(x=K6, color=Iteration)) +
  geom_density()

# Extract 1 set of plausible values for each of the ten imputed datasets and scale and center
for (i in seq_along(1:length(mice.imp.DEEP))){
  
  # Graded response model 
  K6_GRM <- mirt(mice.imp.DEEP[[i]][K6_col_name_num], model = 1, "graded")
  
  # Extract 1 draw of plausible values
  K6_PV <- fscores(K6_GRM, plausible.draws = 1)
  
  # Append it to the imputed dataset 
  mice.imp.DEEP[[i]]$K6_est <- K6_PV
  
  # Scale and centrer 
  mice.imp.DEEP[[i]]$K6_est <- scale(mice.imp.DEEP[[i]]$K6_est, center = T, scale = T)
  mice.imp.DEEP[[i]]$K6_est <- as.numeric(mice.imp.DEEP[[i]]$K6_est)
}


### Create composite variable for social support ### 
for (i in seq_along(1:length( mice.imp.DEEP))){
  mice.imp.DEEP[[i]]$SS_comp <- rowMeans( mice.imp.DEEP[[i]][SS_col_name_num]) 
  mice.imp.DEEP[[i]]$SS_comp_scaled <- scale(mice.imp.DEEP[[i]]$SS_comp, scale = T, center = T)
  mice.imp.DEEP[[i]]$SS_comp_scaled <- as.numeric(mice.imp.DEEP[[i]]$SS_comp_scaled)
}

### Scale and center personal security ###
for (i in seq_along(1:length( mice.imp.DEEP))){
  mice.imp.DEEP[[i]]$PS_2_num_scaled <- scale(mice.imp.DEEP[[i]]$PS_2_num, scale = T, center = T)
}


######### 3) Creating informative priors from previous study ######### 

# Load the results of the prior study 
load("results_table_1_regre.RData") 
load("results_table_1_regre_sub_3.b.RData") 


### Extract the key variables
# Adapted ERI imbalance score 
ERI_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "ERI-score"),][,c(12,17, 18, 19)] 
ERI_est$variance <- abs(ERI_est$SE)^2
ERI_est$SD <- abs(ERI_est$SE)
format(round(ERI_est,3),3)

# Satisfied with “your personal relationships?”
SS1_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Personal relationships"),][,c(12,17, 18, 19)] 
SS1_est$variance <- abs(SS1_est$SE)^2
SS1_est$SD <- abs(SS1_est$SE)
format(round(SS1_est,3),3)

# Satisfied with “the support you get from your friends and family?”
SS2_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Friends and family support"),][,c(12,17, 18, 19)] 
SS2_est$variance <- abs(SS2_est$SE)^2
SS2_est$SD <- abs(SS2_est$SE)
format(round(SS2_est,3),3)

# Satisfied with “the amount of time you are able to spend with friends and family”
SS3_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Friends and family time"),][,c(12,17, 18, 19)] 
SS3_est$variance <- abs(SS3_est$SE)^2
SS3_est$SD <- abs(SS3_est$SE)
format(round(SS3_est,3),3)

# Take means accross the composite social support variables 
SS_comp <- rbind(SS1_est,SS2_est,SS3_est)
SS_comp_mean <- as.data.frame(t(colMeans(SS_comp)))
format(round(SS_comp_mean,3),3)

# Dangerous situations
PS_2_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Dangerous situations"),][,c(12,17, 18, 19)] 
PS_2_est$variance <- abs(PS_2_est$SE)^2
PS_2_est$SD <- abs(PS_2_est$SE)
format(round(PS_2_est,3),3)

# Dispositional optimism 
DO_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Dispositional optimism"),][,c(12,17, 18, 19)]
DO_est$variance <- abs(DO_est$SE)^2 # Ignore the minus sign - relic from calculation in companion manuscript
DO_est$SD <- abs(DO_est$SE)
format(round(DO_est,3),3)

# Age 
age_est <- results_table_1_regre_sub_3.b[ which(results_table_1_regre_sub_3.b$lhs=='K10' & results_table_1_regre_sub_3.b$op == "~" & results_table_1_regre_sub_3.b$rhs == "Age"),][,c(12,17, 18, 19)] 
age_est$variance <- abs(age_est$SE)^2
age_est$SD <- abs(age_est$SE)
format(round(age_est,3),3)

# Years in conservation 
years_est <- results_table_1_regre_sub_3.b[ which(results_table_1_regre_sub_3.b$lhs=='K10' & results_table_1_regre_sub_3.b$op == "~" & results_table_1_regre_sub_3.b$rhs == "Years in conservation"),][,c(12,17, 18, 19)] 
years_est$variance <- abs(years_est$SE)^2
years_est$SD <- abs(years_est$SE)
format(round(years_est,3),3)

# Gender = male 
male_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Men"),][,c(12,17, 18, 19)]
male_est$variance <- abs(male_est$SE)^2
male_est$SD <- abs(male_est$SE)
format(round(male_est,3),3)

# Health
health_est <- results_table_1_regre[ which(results_table_1_regre$lhs=='K10' & results_table_1_regre$op == "~" & results_table_1_regre$rhs == "Physical health"),][,c(12,17, 18, 19)]
health_est$variance <- abs(health_est$SE)^2
health_est$SD <- abs(health_est$SE)
format(round(health_est,3),3)

######### 4) Perform the model diagnostics from the WAMBS check list (Depaoli & van de Schoot 2017): ######### 
######  4.1) 'Do you understand the priors?' ###### 

# Function for plotting strongly informative priors 
prior.plot_strong <- function(df_est, mean= 0, variance=1, sec.min=-0.5, sec.max=0.5, step=.001, label=label) {
  x <- seq(sec.min, sec.max, by = step)
  
  # For a normally distributed prior 
  prior.d <- dnorm(x, mean = mean, sd = sqrt(variance))
  df <- data.frame(cbind(x, prior.d))
  
  # Plot 
  ggplot(data=df, aes(x=x, y=prior.d, group=1)) +
    geom_line(size = 0.5) +
    xlab("Estimate") +
    ylab("Prob. den.") + 
   # scale_x_continuous(
   #   labels = scales::number_format(accuracy = 0.05))  + theme_minimal() +
    geom_ribbon(aes(xmin = df_est$CI.lower.lv, xmax = df_est$CI.upper.lv), fill = "blue", alpha = .30) +
    geom_vline(xintercept = df_est$std.lv, color = "blue", size= 0.5) 
  
  # Save last plot 
  ggsave(filename =  paste0("C:/Users/wolf5246/Dropbox/Oxford/PhD/Chapter_7/Manuscript/Preperation/Priors/", label, ".png"), plot = last_plot(), width = 50, height = 50, units = "mm", dpi = 400,)
  
}

# Function for plotting weakly informative priors 
prior.plot_weak <- function(mean=0, variance=1, sec.min=-6, sec.max=6, step=.001, label=label) {
  x <- seq(sec.min, sec.max, by = step)
  
  # For a normally distributed prior 
  prior.d <- dnorm(x,mean = mean, sd = sqrt(variance))
 
  
  # Plot 
  df <- data.frame(x = x, prior.d = prior.d)  
  ggplot(data=df, aes(x=x, y=prior.d, group=1)) +
    geom_line(size = 0.5) +
    xlab("Estimate") +
    ylab("Prob. den.") # + 
   # scale_x_continuous(
   #   labels = scales::number_format(accuracy = 0.05))
  
  # Save last plot 
  ggsave(filename =  paste0("C:/Users/wolf5246/Dropbox/Oxford/PhD/Chapter_7/Manuscript/Preperation/Priors/", label, ".png"), plot = last_plot(), width = 50, height = 50, units = "mm", dpi = 400,)
  
}


# ERI ~ Position
prior.plot_weak(mean = 0, variance = 100,  sec.min=-30, sec.max=30,
           label="1. ERI ~ position (all)")

# Social support ~ Position
prior.plot_weak(mean = 0, variance = 100,  sec.min=-30, sec.max=30,
           label="2. Social support ~ position (all)")

# Dangerous situation ~ Position
prior.plot_weak(mean = 0, variance = 100,  sec.min=-30, sec.max=30,
           label="3. Dangerous situation ~ position (all)")

# Distress ~ Position
prior.plot_weak(type = "normal", mean = 0, variance = 100,  sec.min=-30, sec.max=30,
           label="4. Distress ~ position (all)")

# Distress ~ ERI
prior.plot_strong(df_est = DO_est, mean = DO_est$std.lv, 
                  variance = DO_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="5. Distress ~ Dispositional optimism")

# Distress ~ Social support
prior.plot_strong(df_est = SS_comp_mean, 
                  mean = SS_comp_mean$std.lv, 
                  variance = SS_comp_mean$variance, sec.min=-0.5, sec.max=0.5,
                  label="6. Distress ~ Social support")

# Distress ~ Dangerous situations
prior.plot_strong(df_est = PS_2_est, 
                  mean = PS_2_est$std.lv, 
                  variance = PS_2_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="7. Distress ~ Dangerous situations")

# Distress ~ Dispositional optimism
prior.plot_strong(df_est = DO_est, 
                  mean = DO_est$std.lv, 
                  variance = DO_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="8. Distress ~ Dispositional optimism")

# Distress ~ Age
prior.plot_strong(df_est = age_est, 
                  mean = age_est$std.lv, 
                  variance = age_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="9. Distress ~ Age")

# Distress ~ Years conservation
prior.plot_strong(df_est = years_est, 
                  mean = years_est$std.lv, 
                  variance = years_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="10. Distress ~ Years conservation")

# Distress ~ Male
prior.plot_strong(df_est = male_est, 
                  mean = male_est$std.lv, 
                  variance = male_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="11. Distress ~ Male")


# Distress ~ Health
prior.plot_strong(df_est = health_est, 
                  mean = health_est$std.lv, 
                  variance = health_est$variance, sec.min=-0.5, sec.max=0.5,
                  label="12. Distress ~ Health")

# Organisation ~ Position
prior.plot_weak( mean = 0, variance = 100,  sec.min=-30, sec.max=30,
           label="13. Distress ~ organisation (all)")


######  2.2) 'Does the trace-plot exhibit convergence?' ###### 
### Define the model and set the priors ###
# Stan uses standard deviations as the priors (rather than variance, like in Jags). See here for useful details: http://ecmerkle.github.io/blavaan/articles/prior.html 
# This means we have to take the square root of the variance to get the standard deviation, which is supplied as the prior.
# Diagnostics are performed using the first of the ten imputed datasets.

# Chains
chains = 4 

# Burn-in iterations 
burnin <- 5000 # 5000

# Post-burn-in iterations
sample <- 5000 # 5000

# Set seed
seed = 123

# Check second argument of normal distribution function 
# variance to SD 
sqrt(0.04) 
sqrt(0.16)
sqrt(0.36)
sqrt(1) 
sqrt(100) 

# SD to variance 
0.2^2
0.4^2
0.6^2
1^2
10^2


# The model structure  - with informative priors
model_1 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.2941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.2110364, 0.04008927)")*gender_simple_Male + prior("normal(0, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(-0.1363494, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(-0.08829758, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(-0.1799596, 0.02504651)")*health_num_scaled +

# Position 
j*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + k*prior("normal(0, 10)")*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*prior("normal(0.2738432, 0.02124379)")*ERI_n_scaled + 

# Social support # SS_comp_mean
h*prior("normal(-0.1129176, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(-0.005403912, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_n_scaled ~ a*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + d*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + e*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + f*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ prior("normal(0, 10)")*Organisation_Organisationtwo + prior("normal(0, 10)")*Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'

# Bayesian SEM
fit_main_1 <- bsem(model_1, data = mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_1, "fit_main_1.rds")

# Load the model (if needed)
fit_main_1 <- readRDS("fit_main_1.rds")

# Traceplots for key parameters - first ten 
key_params <- 22
plot(fit_main_1, pars = 1:10, plot.type = "trace")

# Traceplots for key parameters - remaining 
plot(fit_main_1, pars = 11:key_params, plot.type = "trace")

# Geweke diagnostic
fit_main_1_mcmc.list <- blavInspect(fit_main_1, what = "mcmc")
geweke.plot(fit_main_1_mcmc.list) # We expect about 5% to be out more than +/- 1.96.

# Gelman and Rubin diagnostic (rule of thumb is that everything below 1.1 is OK): https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
gelman.diag(fit_main_1_mcmc.list)
gelman.plot(fit_main_1_mcmc.list)


######  4.3) 'Does convergence remain after doubling the number of iterations?' ###### 

# Double the burn-in 
burnin.d <- burnin*2

# Double the post burn-in
sample.d <- sample*2

# Bayesian SEM  
fit_main_2 <- bsem(model_1, data=mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin.d, sample = sample.d, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_2, "fit_main_2.rds")

# Load the model (if needed)
fit_main_2 <- readRDS("fit_main_2.rds")

# Traceplots for key parameters - first ten 
plot(fit_main_2, pars = 1:10, plot.type = "trace")

# Traceplots for key parameters - remaining 
plot(fit_main_2, pars = 11:key_params, plot.type = "trace")

# Geweke diagnostic
fit_main_2_mcmc.list <- blavInspect(fit_main_2, what = "mcmc")
geweke.plot(fit_main_2_mcmc.list) # We expect about 5% to be out more than +/- 1.96.

# Gelman and Rubin diagnostic (rule of thumb is that everything below 1.1 is OK): https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
gelman.diag(fit_main_2_mcmc.list)
gelman.plot(fit_main_2_mcmc.list)


######  4.4) 'Does the histogram have enough information?' ###### 
plot(fit_main_1 , pars = 1:key_params, plot.type = "hist")


######  4.5) 'Do the chains exhibit a strong degree of autocorrelation?'###### 
# Params 1 - 4
plot(fit_main_1, pars = 1:4, plot.type = "acf")

# Params 5 - 8
plot(fit_main_1, pars = 5:8, plot.type = "acf")

# Params 9 - 12
plot(fit_main_1, pars = 9:12, plot.type = "acf")

# Params 13 - 16
plot(fit_main_1, pars = 14:16, plot.type = "acf")

# Params 17 - 20
plot(fit_main_1, pars = 17:20, plot.type = "acf")

# Params 21 - 22
plot(fit_main_1, pars = 21:22, plot.type = "acf")


######  4.6) 'Does the posterior distribution make substantive sense?' ###### 
# Combing the results of the three chains together 
fit_1_MCMCbinded <- as.matrix(fit_main_1_mcmc.list)

# Examine the posterior for each key variable  
par(mfrow = c(4,7))
for (i in seq_along(1:key_params)) {
  plot(density(fit_1_MCMCbinded[,i]))
}
dev.off()

######  4.7) 'Do different specifications of the multivariate variance priors influence the results?' ###### 

# Examine the observed variable precision parameter (observed, because we are using plausible valued) - theta (which appears to be called 'itheta' in the manual?) - which has a default prior of gamma(1, 0.5).
dpriors(target = "stan")


# The model structure 
model_2 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.2941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.2110364, 0.04008927)")*gender_simple_Male + prior("normal(0, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(-0.1363494, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(-0.08829758, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(-0.1799596, 0.02504651)")*health_num_scaled +

# Position 
j*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + k*prior("normal(0, 10)")*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*prior("normal(0.2738432, 0.02124379)")*ERI_n_scaled + 

# Social support # SS_comp_mean
h*prior("normal(-0.1129176, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(-0.005403912, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_n_scaled ~ a*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + d*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + e*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + f*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ prior("normal(0, 10)")*Organisation_Organisationtwo + prior("normal(0, 10)")*Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k

### Using a alternative to the default 
K6_est ~~ gamma(1, 0.05)[sd]*K6_est
'

# Bayesian SEM  
fit_main_3 <- bsem(model_2, data=mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_3, "fit_main_3.rds")

# Load the model (if needed)
fit_main_3 <- readRDS("fit_main_3.rds")

# Calculate the bias associated with using this different specification
fit_main_1_sum <- data.frame(summary(fit_main_1))
fit_main_2_sum <- data.frame(summary(fit_main_3))
round(100*(as.numeric(fit_main_1_sum$PE.est)-as.numeric(fit_main_2_sum$PE.est))/as.numeric(fit_main_2_sum$PE.est),2) # Fine, if it is only small (either small percentage differences, or where there is small absolute variability)

#####  2.8) 'Is there a notable effect of the prior when compared with noninformative priors?' ###### 
# blavaan has default very weakly informative priors, so the analysis is repeated without specified priors

# The model structure - with default weakly informative priors
model_3 <- '
K6_est ~ 

# Dispositional optimism # DO_est
LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
gender_simple_Male + gender_simple_Unknown +

# Age # age_est
age_year_scaled +

# Years in conservation # years_est
years_cons_scaled +

# Physical health # health_est
health_num_scaled +

# Position 
j*position_simple_Frontlinepractitioner + position_simple_Otherunknown + k*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*ERI_n_scaled + 

# Social support # SS_comp_mean
h*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_n_scaled ~ a*position_simple_Frontlinepractitioner + position_simple_Otherunknown + d*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*position_simple_Frontlinepractitioner + position_simple_Otherunknown + e*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*position_simple_Frontlinepractitioner + position_simple_Otherunknown + f*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ Organisation_Organisationtwo + Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'

# Bayesian SEM
fit_main_4 <- bsem(model_3, data = mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_4, "fit_main_4.rds")

# Load the model (if needed)
fit_main_4 <- readRDS("fit_main_4.rds")

# Calculate the bias associated with using weakly informative priors
fit_main_1_sum <- data.frame(summary(fit_main_1))
fit_main_4_sum <- data.frame(summary(fit_main_4))
round(100*(as.numeric(fit_main_1_sum$PE.est)-as.numeric(fit_main_4_sum$PE.est))/as.numeric(fit_main_4_sum$PE.est),2) 

######  2.9) 'Are the results stable from a sensitivity analysis?' ###### 
# Shifting the hyperparameters 'up' - all mean values by 0.2


# The model structure  - with informative priors
model_4 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.0941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.0110364, 0.04008927)")*gender_simple_Male + prior("normal(0.2, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(0.0636506, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(0.11170242, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(0.0200404, 0.02504651)")*health_num_scaled +

# Position 
j*prior("normal(0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(0.2, 10)")*position_simple_Otherunknown + k*prior("normal(0.2, 10)")*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*prior("normal(0.4738432, 0.02124379)")*ERI_n_scaled + 

# Social support # SS_comp_mean
h*prior("normal(0.0870824, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(0.194596088, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_n_scaled ~ a*prior("normal(0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(0.2, 10)")*position_simple_Otherunknown + d*prior("normal(0.2, 10)")*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(0.2, 10)")*position_simple_Otherunknown + e*prior("normal(0.2, 10)")*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(0.2, 10)")*position_simple_Otherunknown + f*prior("normal(0.2, 10)")*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ prior("normal(0.2, 10)")*Organisation_Organisationtwo + prior("normal(0.2, 10)")*Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'

# Bayesian SEM
fit_main_5 <- bsem(model_4, data = mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_5, "fit_main_5.rds")

# # Load the model (if needed)
# fit_main_5 <- readRDS("fit_main_5.rds")

# Calculate the bias associated with using weakly informative priors
fit_main_1_sum <- data.frame(summary(fit_main_1))
fit_main_5_sum <- data.frame(summary(fit_main_5))
round(100*(as.numeric(fit_main_1_sum$PE.est)-as.numeric(fit_main_5_sum$PE.est))/as.numeric(fit_main_5_sum$PE.est),2)

# Shifting the hyperparameters 'down' - by 0.5

model_5 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.4941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.4110364, 0.04008927)")*gender_simple_Male + prior("normal(-0.2, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(-0.3363494, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(-0.28829758, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(-0.3799596, 0.02504651)")*health_num_scaled +

# Position 
j*prior("normal(-0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(-0.2, 10)")*position_simple_Otherunknown + k*prior("normal(-0.2, 10)")*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*prior("normal(0.0738432, 0.02124379)")*ERI_n_scaled + 

# Social support # SS_comp_mean
h*prior("normal(-0.3129176, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(-0.205403912, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_n_scaled ~ a*prior("normal(-0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(-0.2, 10)")*position_simple_Otherunknown + d*prior("normal(-0.2, 10)")*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(-0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(-0.2, 10)")*position_simple_Otherunknown + e*prior("normal(-0.2, 10)")*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(-0.2, 10)")*position_simple_Frontlinepractitioner + prior("normal(-0.2, 10)")*position_simple_Otherunknown + f*prior("normal(-0.2, 10)")*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ prior("normal(-0.2, 10)")*Organisation_Organisationtwo + prior("normal(-0.2, 10)")*Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'


# Bayesian SEM
fit_main_6 <- bsem(model_5, data = mice.imp.DEEP[[1]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))

# Save the model 
saveRDS(fit_main_6, "fit_main_6.rds")

# # Load the model (if needed)
# fit_main_5 <- readRDS("fit_main_5.rds")

# Calculate the bias associated with using weakly informative priors
fit_main_1_sum <- data.frame(summary(fit_main_1))
fit_main_6_sum <- data.frame(summary(fit_main_6))
round(100*(as.numeric(fit_main_1_sum$PE.est)-as.numeric(fit_main_6_sum$PE.est))/as.numeric(fit_main_6_sum$PE.est),2) 

######  2.10) 'Is the Bayesian way of interpreting and reporting model results used?' ###### 
# Refer to the BARG steps: https://www.nature.com/articles/s41562-021-01177-7#Sec3
# Credibility interval - there is a 95% probability that the true coefficient exists between the credibility interval. 
summary(fit_main_1)


#########  5) Run the model on the ten imputed datasets ######### 
# Model list
fit_sub_list <- list()

# Number of iterations 
iters <- 10 # (change back to run full model)

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){
  
  # Bayesian SEM
  fit_sub_list[[i]] <- bsem(model_1, data=mice.imp.DEEP[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_sub_list, "fit_sub_list.rds")

#########  6) Extract and pool the model results ######### 
fit_sub_list <- readRDS("fit_sub_list.rds")

# Extract the variable names 
summary_DF_names_sub <- paste0(fit_sub_list[[1]]@ParTable$lhs, fit_sub_list[[1]]@ParTable$op , fit_sub_list[[1]]@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_sub_list_sample <-list()
fit_2_MCMCbinded_DF <- list()
for (i in seq_along(1:iters)) {
  fit_sub_list_sample[[i]] <- blavInspect(fit_sub_list[[i]], what="mcmc")
  fit_sub_list_sample[[i]] <- as.matrix(fit_sub_list_sample[[i]])
  fit_2_MCMCbinded_DF[[i]] <- data.frame(fit_sub_list_sample[[i]])
}

# Pool posterior 
fit_2_MCMCbinded_DF <- do.call("rbind", fit_2_MCMCbinded_DF)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF) <- make.unique(summary_DF_names_sub)[1:length(colnames(fit_2_MCMCbinded_DF))]

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF) <- make.unique(summary_DF_names_sub)[1:length(colnames(fit_2_MCMCbinded_DF))]

# Even nicer names for key variables
colnames(fit_2_MCMCbinded_DF) <- revalue(colnames(fit_2_MCMCbinded_DF), c("K6_est~LOTR_total_scaled" = "Distress ~ Dispositional optimism", #
                                                                          "K6_est~gender_simple_Male" = "Distress ~ Male", #
                                                                          "K6_est~age_year_scaled" = "Distress ~ Age", #
                                                                          "K6_est~years_cons_scaled" = "Distress ~ Years in conservation", #
                                                                          "K6_est~health_num_scaled" = "Distress ~ Physical health", #
                                                                          "K6_est~position_simple_Frontlinepractitioner" = "Distress ~ Field-based", #
                                                                          "K6_est~position_simple_Otherunknown" = "Distress ~ Unknown position", #
                                                                          "K6_est~position_simple_Office_basedpractitioner" = "Distress ~ Office-based", #
                                                                          "K6_est~ERI_n_scaled" = "Distress ~ Effort-reward imbalance", #
                                                                          "K6_est~SS_comp_scaled" = "Distress ~ Social support", #
                                                                          "K6_est~PS_2_num_scaled" = "Distress ~ Dangerous situations", #
                                                                          "ERI_n_scaled~position_simple_Frontlinepractitioner" = "Effort-reward imbalance ~ Field-based", #
                                                                          "ERI_n_scaled~position_simple_Otherunknown" = "Effort-reward imbalance ~ Unknown position", #
                                                                          "ERI_n_scaled~position_simple_Office_basedpractitioner" = "Effort-reward imbalance ~ Office-based", #
                                                                          "SS_comp_scaled~position_simple_Frontlinepractitioner" = "Social support ~ Field-based", #
                                                                          "SS_comp_scaled~position_simple_Otherunknown" = "Social support ~ Unknown position", #
                                                                          "SS_comp_scaled~position_simple_Office_basedpractitioner" = "Social support ~ Office-based", #
                                                                          "PS_2_num_scaled~position_simple_Frontlinepractitioner" = "Dangerous situations ~ Field-based", #
                                                                          "PS_2_num_scaled~position_simple_Otherunknown" = "Dangerous situations ~ Unknown position", #
                                                                          "PS_2_num_scaled~position_simple_Office_basedpractitioner" = "Dangerous situations ~ Office-based", #
                                                                          "K6_est~Organisation_Organisationtwo" = "Distress ~ Organisation two", #
                                                                          "K6_est~Organisation_Organisationthree" = "Distress ~ Organisation three" #

                                                                          
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_2_MCMCbinded_DF <- fit_2_MCMCbinded_DF[c("Effort-reward imbalance ~ Field-based", "Effort-reward imbalance ~ Office-based",
                                             "Social support ~ Field-based" , "Social support ~ Office-based", 
                                             "Dangerous situations ~ Field-based", "Dangerous situations ~ Office-based",
                                             "Distress ~ Field-based", "Distress ~ Office-based",
                                             "Distress ~ Effort-reward imbalance", "Distress ~ Social support", "Distress ~ Dangerous situations",
                                             "Distress ~ Dispositional optimism",
                                             "Distress ~ Age",
                                             "Distress ~ Years in conservation",
                                             "Distress ~ Male",
                                             "Distress ~ Physical health",
                                             "Distress ~ Organisation two",
                                             "Distress ~ Organisation three", 
                                             colnames(fit_2_MCMCbinded_DF)[13:length(colnames(fit_2_MCMCbinded_DF))] )]

# Save pooled samples
saveRDS(fit_2_MCMCbinded_DF, "fit_2_MCMCbinded_DF.rds")


#########  7) Run the model on the ten imputed datasets with default weakly informative priors ######### 
# Model list
fit_sub_list_weak <- list()

# Number of iterations 
iters <- 10 # (change back to run full model)

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){
  
  # Bayesian SEM
  fit_sub_list_weak[[i]] <- bsem(model_3, data=mice.imp.DEEP[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_sub_list_weak, "fit_sub_list_weak.rds")
fit_sub_list_weak <- readRDS("fit_sub_list_weak.rds")

# Extract the variable names 
summary_DF_names_sub_weak <- paste0(fit_sub_list_weak[[1]]@ParTable$lhs, fit_sub_list_weak[[1]]@ParTable$op , fit_sub_list_weak[[1]]@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_sub_list_sample_weak <-list()
fit_2_MCMCbinded_DF_weak <- list()
for (i in seq_along(1:iters)) {
  fit_sub_list_sample_weak[[i]] <- blavInspect(fit_sub_list_weak[[i]], what="mcmc")
  fit_sub_list_sample_weak[[i]] <- as.matrix(fit_sub_list_sample_weak[[i]])
  fit_2_MCMCbinded_DF_weak[[i]] <- data.frame(fit_sub_list_sample_weak[[i]])
}

# Pool posterior 
fit_2_MCMCbinded_DF_weak <- do.call("rbind", fit_2_MCMCbinded_DF_weak)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF_weak) <- make.unique(summary_DF_names_sub_weak)[1:length(colnames(fit_2_MCMCbinded_DF_weak))]

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF_weak) <- make.unique(summary_DF_names_sub_weak)[1:length(colnames(fit_2_MCMCbinded_DF_weak))]

# Even nicer names for key variables
colnames(fit_2_MCMCbinded_DF_weak) <- revalue(colnames(fit_2_MCMCbinded_DF_weak), c("K6_est~LOTR_total_scaled" = "Distress ~ Dispositional optimism", #
                                                                          "K6_est~gender_simple_Male" = "Distress ~ Male", #
                                                                          "K6_est~age_year_scaled" = "Distress ~ Age", #
                                                                          "K6_est~years_cons_scaled" = "Distress ~ Years in conservation", #
                                                                          "K6_est~health_num_scaled" = "Distress ~ Physical health", #
                                                                          "K6_est~position_simple_Frontlinepractitioner" = "Distress ~ Field-based", #
                                                                          "K6_est~position_simple_Otherunknown" = "Distress ~ Unknown position", #
                                                                          "K6_est~position_simple_Office_basedpractitioner" = "Distress ~ Office-based", #
                                                                          "K6_est~ERI_n_scaled" = "Distress ~ Effort-reward imbalance", #
                                                                          "K6_est~SS_comp_scaled" = "Distress ~ Social support", #
                                                                          "K6_est~PS_2_num_scaled" = "Distress ~ Dangerous situations", #
                                                                          "ERI_n_scaled~position_simple_Frontlinepractitioner" = "Effort-reward imbalance ~ Field-based", #
                                                                          "ERI_n_scaled~position_simple_Otherunknown" = "Effort-reward imbalance ~ Unknown position", #
                                                                          "ERI_n_scaled~position_simple_Office_basedpractitioner" = "Effort-reward imbalance ~ Office-based", #
                                                                          "SS_comp_scaled~position_simple_Frontlinepractitioner" = "Social support ~ Field-based", #
                                                                          "SS_comp_scaled~position_simple_Otherunknown" = "Social support ~ Unknown position", #
                                                                          "SS_comp_scaled~position_simple_Office_basedpractitioner" = "Social support ~ Office-based", #
                                                                          "PS_2_num_scaled~position_simple_Frontlinepractitioner" = "Dangerous situations ~ Field-based", #
                                                                          "PS_2_num_scaled~position_simple_Otherunknown" = "Dangerous situations ~ Unknown position", #
                                                                          "PS_2_num_scaled~position_simple_Office_basedpractitioner" = "Dangerous situations ~ Office-based", #
                                                                          "K6_est~Organisation_Organisationtwo" = "Distress ~ Organisation two", #
                                                                          "K6_est~Organisation_Organisationthree" = "Distress ~ Organisation three" #
                                                                          
                                                                          
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_2_MCMCbinded_DF_weak <- fit_2_MCMCbinded_DF_weak[c("Effort-reward imbalance ~ Field-based", "Effort-reward imbalance ~ Office-based",
                                             "Social support ~ Field-based" , "Social support ~ Office-based", 
                                             "Dangerous situations ~ Field-based", "Dangerous situations ~ Office-based",
                                             "Distress ~ Field-based", "Distress ~ Office-based",
                                             "Distress ~ Effort-reward imbalance", "Distress ~ Social support", "Distress ~ Dangerous situations",
                                             "Distress ~ Dispositional optimism",
                                             "Distress ~ Age",
                                             "Distress ~ Years in conservation",
                                             "Distress ~ Male",
                                             "Distress ~ Physical health",
                                             "Distress ~ Organisation two",
                                             "Distress ~ Organisation three", 
                                             colnames(fit_2_MCMCbinded_DF_weak)[13:length(colnames(fit_2_MCMCbinded_DF_weak))] )]

# Save pooled samples
saveRDS(fit_2_MCMCbinded_DF_weak, "fit_2_MCMCbinded_DF_weak.rds")


### Comparing models 
compare_models <- list()
for (i in seq_along(1:length(fit_sub_list))){
  results <- blavCompare(fit_sub_list[[i]], fit_sub_list_weak[[i]])
  compare_models[[i]] <- results$waic[[1]][3] - results$waic[[2]][3]
}

# Mean difference in WAIC - the lower the better (and > 2 delta WAIC implies notably better fir)
mean(do.call("rbind",compare_models)) 
# Weakly informative prior model performs better. 
# However, I'd expect models with weakly informative priors to fit the data better than those with 
# strongly informative priors (if the effect size implied by the prior differs to that implied by the data, 
# which could easily happen with a small dataset.)


#########  8) Run the model on the ten imputed datasets using the original ERI ######### 

# Drop two observations which don't have original ERI scores 
mice.imp.DEEP_OERI <- list()
for (i in seq_along(1:iters)){
  mice.imp.DEEP_OERI[[i]] <- mice.imp.DEEP[[i]][!is.na(mice.imp.DEEP[[i]]$ERI_o_scaled), ]
}


# The model structure  - with informative priors
model_6 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.2941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.2110364, 0.04008927)")*gender_simple_Male + prior("normal(0, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(-0.1363494, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(-0.08829758, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(-0.1799596, 0.02504651)")*health_num_scaled +

# Position 
j*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + k*prior("normal(0, 10)")*position_simple_Office_basedpractitioner + 

# ERI # ERI_est
g*prior("normal(0.2738432, 0.02124379)")*ERI_o_scaled + 

# Social support # SS_comp_mean
h*prior("normal(-0.1129176, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(-0.005403912, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = position_simple_Researchacademia)
ERI_o_scaled ~ a*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + d*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + e*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(0, 10)")*position_simple_Frontlinepractitioner + prior("normal(0, 10)")*position_simple_Otherunknown + f*prior("normal(0, 10)")*position_simple_Office_basedpractitioner 

### Organisation (RL = Organisation_Organisationone)
K6_est ~ prior("normal(0, 10)")*Organisation_Organisationtwo + prior("normal(0, 10)")*Organisation_Organisationthree

###### Total effects 
#### Frontline and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Frontline -> ERI -> distress 
bh := b*h # Frontline -> Social support -> distress 
ci := c*i # Frontline -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'

# Model list
fit_sub_list_OERI <- list()

# Number of iterations 
iters <- 10 # (change back to run full model)

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){
  
  # Bayesian SEM
  fit_sub_list_OERI[[i]] <- bsem(model_6, data=mice.imp.DEEP_OERI[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_sub_list_OERI, "fit_sub_list_OERI.rds")
fit_sub_list_OERI <- readRDS("fit_sub_list_OERI.rds")

# Extract the variable names 
summary_DF_names_sub_OERI <- paste0(fit_sub_list_OERI[[1]]@ParTable$lhs, fit_sub_list_OERI[[1]]@ParTable$op , fit_sub_list_OERI[[1]]@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_sub_list_sample_OERI <-list()
fit_2_MCMCbinded_DF_OERI <- list()
for (i in seq_along(1:iters)) {
  fit_sub_list_sample_OERI[[i]] <- blavInspect(fit_sub_list_OERI[[i]], what="mcmc")
  fit_sub_list_sample_OERI[[i]] <- as.matrix(fit_sub_list_sample_OERI[[i]])
  fit_2_MCMCbinded_DF_OERI[[i]] <- data.frame(fit_sub_list_sample_OERI[[i]])
}

# Pool posterior 
fit_2_MCMCbinded_DF_OERI <- do.call("rbind", fit_2_MCMCbinded_DF_OERI)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF_OERI) <- make.unique(summary_DF_names_sub_OERI)[1:length(colnames(fit_2_MCMCbinded_DF_OERI))]

# Even nicer names for key variables
colnames(fit_2_MCMCbinded_DF_OERI) <- revalue(colnames(fit_2_MCMCbinded_DF_OERI), c("K6_est~LOTR_total_scaled" = "Distress ~ Dispositional optimism", #
                                                                                    "K6_est~gender_simple_Male" = "Distress ~ Male", #
                                                                                    "K6_est~age_year_scaled" = "Distress ~ Age", #
                                                                                    "K6_est~years_cons_scaled" = "Distress ~ Years in conservation", #
                                                                                    "K6_est~health_num_scaled" = "Distress ~ Physical health", #
                                                                                    "K6_est~position_simple_Frontlinepractitioner" = "Distress ~ Field-based", #
                                                                                    "K6_est~position_simple_Otherunknown" = "Distress ~ Unknown position", #
                                                                                    "K6_est~position_simple_Office_basedpractitioner" = "Distress ~ Office-based", #
                                                                                    "K6_est~ERI_o_scaled" = "Distress ~ Effort-reward imbalance (original)", #
                                                                                    "K6_est~SS_comp_scaled" = "Distress ~ Social support", #
                                                                                    "K6_est~PS_2_num_scaled" = "Distress ~ Dangerous situations", #
                                                                                    "ERI_o_scaled~position_simple_Frontlinepractitioner" = "Effort-reward imbalance (original) ~ Field-based", #
                                                                                    "ERI_o_scaled~position_simple_Otherunknown" = "Effort-reward imbalance (original) ~ Unknown position", #
                                                                                    "ERI_o_scaled~position_simple_Office_basedpractitioner" = "Effort-reward imbalance (original) ~ Office-based", #
                                                                                    "SS_comp_scaled~position_simple_Frontlinepractitioner" = "Social support ~ Field-based", #
                                                                                    "SS_comp_scaled~position_simple_Otherunknown" = "Social support ~ Unknown position", #
                                                                                    "SS_comp_scaled~position_simple_Office_basedpractitioner" = "Social support ~ Office-based", #
                                                                                    "PS_2_num_scaled~position_simple_Frontlinepractitioner" = "Dangerous situations ~ Field-based", #
                                                                                    "PS_2_num_scaled~position_simple_Otherunknown" = "Dangerous situations ~ Unknown position", #
                                                                                    "PS_2_num_scaled~position_simple_Office_basedpractitioner" = "Dangerous situations ~ Office-based", #
                                                                                    "K6_est~Organisation_Organisationtwo" = "Distress ~ Organisation two", #
                                                                                    "K6_est~Organisation_Organisationthree" = "Distress ~ Organisation three" #
                                                                                    
                                                                                    
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_2_MCMCbinded_DF_OERI <- fit_2_MCMCbinded_DF_OERI[c("Effort-reward imbalance (original) ~ Field-based", "Effort-reward imbalance (original) ~ Office-based",
                                                       "Social support ~ Field-based" , "Social support ~ Office-based", 
                                                       "Dangerous situations ~ Field-based", "Dangerous situations ~ Office-based",
                                                       "Distress ~ Field-based", "Distress ~ Office-based",
                                                       "Distress ~ Effort-reward imbalance (original)", "Distress ~ Social support", "Distress ~ Dangerous situations",
                                                       "Distress ~ Dispositional optimism",
                                                       "Distress ~ Age",
                                                       "Distress ~ Years in conservation",
                                                       "Distress ~ Male",
                                                       "Distress ~ Physical health",
                                                       "Distress ~ Organisation two",
                                                       "Distress ~ Organisation three", 
                                                       colnames(fit_2_MCMCbinded_DF_OERI)[13:length(colnames(fit_2_MCMCbinded_DF_OERI))] )]

# Save pooled samples
saveRDS(fit_2_MCMCbinded_DF_OERI, "fit_2_MCMCbinded_DF_OERI.rds")


#########  9) Run the model on the ten imputed datasets replacing position with organisation ######### 

# The model structure  - with informative priors
model_7 <- '
K6_est ~ 

# Dispositional optimism # DO_est
prior("normal(-0.2941245, 0.02650791)")*LOTR_total_scaled +

# Gender (RL = gender_Female) # male_est
prior("normal(-0.2110364, 0.04008927)")*gender_simple_Male + prior("normal(0, 10)")*gender_simple_Unknown +

# Age # age_est
prior("normal(-0.1363494, 0.03282156)")*age_year_scaled +

# Years in conservation # years_est
prior("normal(-0.08829758, 0.03433606)")*years_cons_scaled +

# Physical health # health_est
prior("normal(-0.1799596, 0.02504651)")*health_num_scaled +

# Organisation (RL = Organisation_Organisationone) 
j*prior("normal(0, 10)")*Organisation_Organisationtwo + k*prior("normal(0, 10)")*Organisation_Organisationthree + 

# ERI # ERI_est
g*prior("normal(0.2738432, 0.02124379)")*ERI_n_scaled + 

# Social support # SS_comp_mean
h*prior("normal(-0.1129176, 0.02163945)")*SS_comp_scaled + 

# Personal insecurity # PS_2_est
i*prior("normal(-0.005403912, 0.01666645)")*PS_2_num_scaled  # My work puts me in dangerous situations

### Mediated effects ###

# ERI and position (RL = Organisation_Organisationone)
ERI_n_scaled ~ a*prior("normal(0, 10)")*Organisation_Organisationtwo + d*prior("normal(0, 10)")*Organisation_Organisationthree 

# Social support and position 
SS_comp_scaled ~ b*prior("normal(0, 10)")*Organisation_Organisationtwo + e*prior("normal(0, 10)")*Organisation_Organisationthree 

# Dangerous situation and position 
PS_2_num_scaled ~ c*prior("normal(0, 10)")*Organisation_Organisationtwo + f*prior("normal(0, 10)")*Organisation_Organisationthree 

### Position (RL = position_simple_Researchacademia)
K6_est ~ prior("normal(0, 10)")*position_simple_Frontlinepractitioner  +prior("normal(0, 10)")*position_simple_Otherunknown + prior("normal(0, 10)")*position_simple_Office_basedpractitioner

###### Total effects 
#### Org 2 and distress
# Direct effects = j

# Indirect effects 
ag := a*g # Org 2 -> ERI -> distress 
bh := b*h # Org 2 -> Social support -> distress 
ci := c*i # Org 2 -> Danger -> distress 

# Total effects 
Total.Front := ag + bh + ci + j

#### Office-based and distress
# Direct effects = k

# Indirect effects 
dg := d*g # Office -> ERI -> distress 
eh := e*h # Office -> Social support -> distress 
fi := f*i # Office -> Danger -> distress 

# Total effects 
Total.Office := dg + eh + fi + k
'

# Model list
fit_sub_list_ORG <- list()

# Number of iterations 
iters <- 10 # (change back to run full model)

# Run the model for each imputed dataset
for (i in seq_along(1:iters)){
  
  # Bayesian SEM
  fit_sub_list_ORG[[i]] <- bsem(model_7, data=mice.imp.DEEP[[i]], fixed.x = FALSE, n.chains = chains, burnin = burnin, sample = sample, seed = seed, target = "stan",  bcontrol = list(cores = 4))
}

# Save the model
saveRDS(fit_sub_list_ORG, "fit_sub_list_ORG.rds")
fit_sub_list_ORG <- readRDS("fit_sub_list_ORG.rds")

# Extract the variable names 
summary_DF_names_sub_ORG <- paste0(fit_sub_list_ORG[[1]]@ParTable$lhs, fit_sub_list_ORG[[1]]@ParTable$op , fit_sub_list_ORG[[1]]@ParTable$rhs)

# Create a DF from the MCMC draws 
fit_sub_list_sample_ORG <-list()
fit_2_MCMCbinded_DF_ORG <- list()
for (i in seq_along(1:iters)) {
  fit_sub_list_sample_ORG[[i]] <- blavInspect(fit_sub_list_ORG[[i]], what="mcmc")
  fit_sub_list_sample_ORG[[i]] <- as.matrix(fit_sub_list_sample_ORG[[i]])
  fit_2_MCMCbinded_DF_ORG[[i]] <- data.frame(fit_sub_list_sample_ORG[[i]])
}

# Pool posterior 
fit_2_MCMCbinded_DF_ORG <- do.call("rbind", fit_2_MCMCbinded_DF_ORG)

# Rename the columns (blavInspect does not appear to return the operator defined parameters, so only returning the 1:n names from 'summary_DF_names'. If needed these can be calculated manually: https://groups.google.com/g/blavaan/c/69ukdYLpHXI). 
colnames(fit_2_MCMCbinded_DF_ORG) <- make.unique(summary_DF_names_sub_ORG)[1:length(colnames(fit_2_MCMCbinded_DF_ORG))]

# Even nicer names for key variables
colnames(fit_2_MCMCbinded_DF_ORG) <- revalue(colnames(fit_2_MCMCbinded_DF_ORG), c("K6_est~LOTR_total_scaled" = "Distress ~ Dispositional optimism", #
                                                                                    "K6_est~gender_simple_Male" = "Distress ~ Male", #
                                                                                    "K6_est~age_year_scaled" = "Distress ~ Age", #
                                                                                    "K6_est~years_cons_scaled" = "Distress ~ Years in conservation", #
                                                                                    "K6_est~health_num_scaled" = "Distress ~ Physical health", #
                                                                                    "K6_est~Organisation_Organisationtwo" = "Distress ~ Organisation two", #
                                                                                    "K6_est~Organisation_Organisationthree" = "Distress ~ Organisation three", #
                                                                                    "K6_est~ERI_n_scaled" = "Distress ~ Effort-reward imbalance", #
                                                                                    "K6_est~SS_comp_scaled" = "Distress ~ Social support", #
                                                                                    "K6_est~PS_2_num_scaled" = "Distress ~ Dangerous situations", #
                                                                                    "ERI_n_scaled~Organisation_Organisationtwo" = "Effort-reward imbalance ~ Organisation two", #
                                                                                    "ERI_n_scaled~Organisation_Organisationthree" = "Effort-reward imbalance ~ Organisation three", #
                                                                                    "SS_comp_scaled~Organisation_Organisationtwo" = "Social support ~ Organisation two", #
                                                                                    "SS_comp_scaled~Organisation_Organisationthree" = "Social support ~ Organisation three", #
                                                                                    "PS_2_num_scaled~Organisation_Organisationtwo" = "Dangerous situations ~ Organisation two", #
                                                                                    "PS_2_num_scaled~Organisation_Organisationthree" = "Dangerous situations ~ Organisation three", #
                                                                                    "K6_est~position_simple_Frontlinepractitioner" = "Distress ~ Field-based", #
                                                                                  "K6_est~position_simple_Otherunknown" = "Distress ~ Unknown position", #
                                                                                  "K6_est~position_simple_Office_basedpractitioner" = "Distress ~ Office-based" #
                                                                                    
                                                                                    
))

# Changing the variable order, so it matches Table 1 in the manuscript.
fit_2_MCMCbinded_DF_ORG <- fit_2_MCMCbinded_DF_ORG[c("Effort-reward imbalance ~ Organisation two", "Effort-reward imbalance ~ Organisation three",
                                                       "Social support ~ Organisation two" , "Social support ~ Organisation three", 
                                                       "Dangerous situations ~ Organisation two", "Dangerous situations ~ Organisation three",
                                                       "Distress ~ Organisation two", "Distress ~ Organisation three",
                                                       "Distress ~ Effort-reward imbalance", "Distress ~ Social support", "Distress ~ Dangerous situations",
                                                       "Distress ~ Dispositional optimism",
                                                       "Distress ~ Age",
                                                       "Distress ~ Years in conservation",
                                                       "Distress ~ Male",
                                                       "Distress ~ Physical health",
                                                       "Distress ~ Field-based",
                                                       "Distress ~ Office-based", 
                                                       colnames(fit_2_MCMCbinded_DF_ORG)[13:length(colnames(fit_2_MCMCbinded_DF_ORG))] )]

# Save pooled samples
saveRDS(fit_2_MCMCbinded_DF_ORG, "fit_2_MCMCbinded_DF_ORG.rds")




