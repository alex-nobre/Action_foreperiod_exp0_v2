
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(forcats)
library(simr)

# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

#=================== Simular usando o resultado do exp 1 ===================
# Usamos os betas do exp 1 pra construir um modelo no simr

betas <- readRDS("./Analysis/trimlogfplmm.RDS")

inter_beta <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "scaledNumForeperiod:condition1") |>
  pull(Estimate)

inter_beta_sd <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "scaledNumForeperiod:condition1") |>
  rename(std_error = "Std. Error") |>
  mutate(std_dev = std_error*sqrt(summary(betas)$ngrps)) |>
  pull(std_dev) 




# rodar n simulações e calcular a porcentagem que fica significativa.


# Com isso, pegamos o valor do beta que é detectável 80% das vezes.


# Read exp 1 data
goData2 <- read_csv("E:/Post-doc_data/Action_foreperiod/Action_fp_exp1/Analysis/goData2.csv")

#============= Build model without data ==================
ID_column <- seq(1:38)
fp_column <- c(1.0, 1.6, 2.2, 2.8)
condition_column <- c("external", "action")
oneBackFP_column <- c(1.0, 1.6, 2.2, 2.8)

beta_data <- expand.grid(ID = ID_column,
                         scaledNumForeperiod = fp_column,
                         condition = condition_column,
                         oneBackFP = oneBackFP_column) %>%
  mutate(ID = as_factor(ID), 
         oneBackFP = as_factor(oneBackFP))


beta_coefs <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  pull(Estimate)

beta_vcov <- betas |>
  vcov() |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname, Intercept = `(Intercept)`) |>
  as_tibble() |>
  filter(coefficient %in% c("(Intercept)", "scaledNumForeperiod", "condition1")) |>
  select(Intercept, scaledNumForeperiod, condition1) |>
  t() |>
  #as.vector() |>
  as.list()

beta_ressd <- attr(betas, "sigma")

# Get formula: attr(betas, "call")$formula

sim_model <- makeLmer(logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                        condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                        scaledNumForeperiod:oneBackFP:condition + (1 + condition + 
                                                                     scaledNumForeperiod | ID),
                      fixef = beta_coefs,
                      VarCorr = beta_vcov,
                      sigma = beta_ressd,
                      data = beta_data)

#==========================================================================#
#=================== Simular usando o resultado do exp 2 ===================
#==========================================================================#
# Read exp 2 data
load("./Analysis/data2.Rdata")

# Usamos os betas do exp 2 pra construir um modelo no simr

betas <- readRDS("./Analysis/fplmm.rds")

inter_beta <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "foreperiod.L:condition1") |>
  pull(Estimate)

inter_beta_sd <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "foreperiod.L:condition1") |>
  rename(std_error = "Std. Error") |>
  mutate(std_dev = std_error*sqrt(summary(betas)$ngrps)) |>
  pull(std_dev) 



# rodar n simulações e calcular a porcentagem que fica significativa.
# Com isso, pegamos o valor do beta que é detectável 80% das vezes.

# variamos o valor dos betas, e fazemos como eles recomendam no artigo
n_betas <- 20
#betas_sim <- rnorm(n_betas, inter_beta, 3*inter_beta_sd)
low_beta <- 0 #-0.01 #min(betas_sim)
high_beta <- 0.01 #max(betas_sim)

betas_sim <- seq(low_beta, high_beta, length.out = n_betas)

# Simulate for each beta
n_sim <- 1000


# Function to simulate
simulate_beta <- function(eff_size) {
  model <- betas
  fixef(model)["foreperiod.L:condition1"] <- eff_size
  sim <- powerSim(model, nsim = n_sim, test = fixed("foreperiod.L:condition1", "z"))
  sim_pvalues <- sim$pval
}

pvals_list <- lapply(betas_sim, simulate_beta)


load(file = "./Analysis/pvals_list_1000sim_20betas_positive.Rdata")
power <- sapply(pvals_list, function(x) sum(x < 0.05)/length(x))
names(power) <- round(betas_sim,4)

plot(round(betas_sim,4), power,
     xlab = "beta value",
     ylab = "power")
lines(round(betas_sim,4), power)
abline(v = inter_beta, lty = "dashed")
abline(h = 0.8)

# function to convert rts to ms
10^betas_sim


xtabs(logRT ~ scaledNumForeperiod, data2)
xtabs(RT ~ foreperiod, data2)

mean(data2$RT)

by(data2$RT, data2$foreperiod, mean)
by(data2$logRT, data2$foreperiod, mean)
by(data2$logRT, data2$scaledNumForeperiod, mean)


data3 <- data2 |>
  mutate(scaledNumForeperiod = round(scaledNumForeperiod,4))

mean(data2[data2$numForeperiod==1,]$logRT) - mean(data2[data2$numForeperiod==1.6,]$logRT)
mean(data3[data3$scaledNumForeperiod==-0.9091,]$logRT) - mean(data3[data3$scaledNumForeperiod==0.8909,]$logRT)

summary(betas)$coefficients # change in 1 in FP = change of -0.008 in logRT

10^(-0.395-0.0081797726)
10^-0.0081797726


#============= Build model without data ==================
ID_column <- seq(1:38)
fp_column <- c(1.0, 1.6, 2.2, 2.8)
condition_column <- c("external", "action")
oneBackFP_column <- c(1.0, 1.6, 2.2, 2.8)

beta_data <- expand.grid(ID = ID_column,
                         scaledNumForeperiod = fp_column,
                         condition = condition_column,
                         oneBackFP = oneBackFP_column) %>%
  mutate(ID = as_factor(ID), 
         oneBackFP = as_factor(oneBackFP))

# Saved model parameters
betas <- readRDS("./Analysis/trimlogfplmm.RDS")


inter_beta <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "scaledNumForeperiod:condition1") |>
  pull(Estimate)

inter_beta_sd <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  filter(coefficient == "scaledNumForeperiod:condition1") |>
  rename(std_error = "Std. Error") |>
  mutate(std_dev = std_error*sqrt(summary(betas)$ngrps)) |>
  pull(std_dev) 

# Build model from scratch  
beta_coefs <- summary(betas)$coefficients |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname) |>
  as_tibble() |>
  pull(Estimate)

beta_vcov <- betas |>
  vcov() |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column() |>
  rename(coefficient = rowname, Intercept = `(Intercept)`) |>
  as_tibble() |>
  filter(coefficient %in% c("(Intercept)", "scaledNumForeperiod", "condition1", "scaledNumForeperiod:condition1")) |>
  select(Intercept, scaledNumForeperiod, condition1, `scaledNumForeperiod:condition1`) |>
  as.matrix()
  #t() |>
  #as.list()

beta_ressd <- attr(betas, "sigma")

# Get formula: attr(betas, "call")$formula

sim_model <- makeLmer(logRT ~ 1 + scaledNumForeperiod + condition + scaledNumForeperiod:condition + 
                        oneBackFP + scaledNumForeperiod:oneBackFP + condition:oneBackFP + 
                        scaledNumForeperiod:condition:oneBackFP + (1 + condition + scaledNumForeperiod + scaledNumForeperiod:condition | ID),
                      fixef = beta_coefs,
                      VarCorr = beta_vcov,
                      sigma = beta_ressd,
                      data = beta_data)

# We used the save model parameters to choose effect sizes from which to simulate
n_betas <- 10
betas_sim <- rnorm(n_betas, inter_beta, inter_beta_sd)
low_beta <- min(betas_sim)
high_beta <- max(betas_sim)

betas_sim <- seq(low_beta, high_beta, length.out = n_betas)

# Simulate for each beta
n_sim <- 100

pvals_list <- vector(mode = "list", length = n_betas)

for(eff_size in 1:n_betas) {
  model <- betas
  fixef(model)["scaledNumForeperiod:condition1"] <- betas_sim[eff_size]
  
  sim <- powerSim(model, nsim = n_sim, test = fcompare(logRT ~ scaledNumForeperiod:condition))
  
  sim_pvalues <- sim$pval
  pvals_list[[eff_size]] <- sim_pvalues
}

power <- sapply(pvals_list, function(x) sum(x < 0.05)/length(x))


#=================================================================================================================================#
# sim_test <- powerSim(betas, nsim = 10, test = fcompare(logRT ~ scaledNumForeperiod:condition))
# 
# beta_test <- betas
# fixef(beta_test)["scaledNumForeperiod:condition1"] <- betas_sim[1]