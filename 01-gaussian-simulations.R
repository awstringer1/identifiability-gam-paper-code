### Simulations ###
# Gaussian response
# For paper: Identifiability Constraints in Generalized Additive Models
# Section 4.1
# Alex Stringer
# 2022/05

## Global parameters ----
n <- c(200,500,1000)
f1 <- function(x) sin(2*pi*x) # Smooth function 1
f2 <- function(x) cos(2*pi*x) # Smooth function 2
alpha <- 2 # Intercept
sigma <- 1 # Gaussian SD
B <- 1e03 # Number of simulated datasets to produce
# Reproducible random number generation with mclapply
RNGkind("L'Ecuyer-CMRG")
## Setup ----
# Packages
library(mgcv)
library(tidyverse)
library(parallel)
options(mc.cores = parallel::detectCores())

# Paths
resultspath <- tempdir()



## Simulation function ----
do_a_simulation <- function(lst) {
  # lst: list with elements:
  #      n: sample size
  #      id: unique identifier of the specific simulation, id = 1...B
  
  n <- lst$n
  id <- lst$id
  
  # Generate Data #
  x1 <- runif(n,0,1)
  x2 <- runif(n,0,1)
  mu <- alpha + f1(x1) + f2(x2)
  y <- rnorm(n,mu,sigma)
  dat <- data.frame(y=y,x1=x1,x2=x2)

  # Fit the models #
  pcx1_1 <- mean(x1)
  pcx1_2 <- 0.1
  pcx2_1 <- mean(x2)
  pcx2_2 <- 0.1
  
  mod_stz <- gam(y ~ s(x1)+s(x2),data = dat,method = "REML")
  mod_pc1 <- gam(y ~ s(x1,pc=pcx1_1)+s(x2,pc=pcx2_1),data = dat,method = "REML")
  mod_pc2 <- gam(y ~ s(x1,pc=pcx1_2)+s(x2,pc=pcx2_2),data = dat,method = "REML")
  
  # Get summary information #
  
  preddat <- data.frame(
    x1 = c(sort(x1),rep(0,n)),
    x2 = c(rep(0,n),sort(x2))
  )
  
  f1_stz_c <- mean(f1(preddat$x1[1:n]))
  f2_stz_c <- mean(f2(preddat$x2[(n+1):(2*n)]))
  f1_pc1_c <- f1(pcx1_1)
  f1_pc2_c <- f1(pcx1_2)
  f2_pc1_c <- f2(pcx2_1)
  f2_pc2_c <- f2(pcx2_2)
  
  pred_stz <- predict(mod_stz,type='terms',newdata = preddat,se.fit=TRUE,unconditional=FALSE)
  pred_pc1 <- predict(mod_pc1,type='terms',newdata = preddat,se.fit=TRUE,unconditional=FALSE)
  pred_pc2 <- predict(mod_pc2,type='terms',newdata = preddat,se.fit=TRUE,unconditional=FALSE)
  
  fit_stz_x1 <- pred_stz$fit[1:n,1]
  fit_stz_x2 <- pred_stz$fit[(n+1):(2*n),2]
  fit_pc1_x1 <- pred_pc1$fit[1:n,1]
  fit_pc1_x2 <- pred_pc1$fit[(n+1):(2*n),2]
  fit_pc2_x1 <- pred_pc2$fit[1:n,1]
  fit_pc2_x2 <- pred_pc2$fit[(n+1):(2*n),2]
  
  se_stz_x1 <- pred_stz$se.fit[1:n,1]
  se_stz_x2 <- pred_stz$se.fit[(n+1):(2*n),2]
  se_pc1_x1 <- pred_pc1$se.fit[1:n,1]
  se_pc1_x2 <- pred_pc1$se.fit[(n+1):(2*n),2]
  se_pc2_x1 <- pred_pc2$se.fit[1:n,1]
  se_pc2_x2 <- pred_pc2$se.fit[(n+1):(2*n),2]
  
  truth_stz_x1 <- f1(preddat$x1[1:n]) - f1_stz_c
  truth_stz_x2 <- f2(preddat$x2[(n+1):(2*n)]) - f2_stz_c
  truth_pc1_x1 <- f1(preddat$x1[1:n]) - f1_pc1_c
  truth_pc1_x2 <- f2(preddat$x2[(n+1):(2*n)]) - f2_pc1_c
  truth_pc2_x1 <- f1(preddat$x1[1:n]) - f1_pc2_c
  truth_pc2_x2 <- f2(preddat$x2[(n+1):(2*n)]) - f2_pc2_c
  
  # Compute summary statistics #
  RMSE_stz_f1 <- sqrt(mean((fit_stz_x1 - truth_stz_x1)^2))
  RMSE_pc1_f1 <- sqrt(mean((fit_pc1_x1 - truth_pc1_x1)^2))
  RMSE_pc2_f1 <- sqrt(mean((fit_pc2_x1 - truth_pc2_x1)^2))
  RMSE_stz_f2 <- sqrt(mean((fit_stz_x2 - truth_stz_x2)^2))
  RMSE_pc1_f2 <- sqrt(mean((fit_pc1_x2 - truth_pc1_x2)^2))
  RMSE_pc2_f2 <- sqrt(mean((fit_pc2_x2 - truth_pc2_x2)^2))
  
  SE_stz_f1 <- mean(se_stz_x1)
  SE_pc1_f1 <- mean(se_pc1_x1)
  SE_pc2_f1 <- mean(se_pc2_x1)
  SE_stz_f2 <- mean(se_stz_x2)
  SE_pc1_f2 <- mean(se_pc1_x2)
  SE_pc2_f2 <- mean(se_pc2_x2)
  
  covr_stz_f1 <- mean(truth_stz_x1 >= fit_stz_x1 - 2*se_stz_x1 & truth_stz_x1 <= fit_stz_x1 + 2*se_stz_x1)
  covr_pc1_f1 <- mean(truth_pc1_x1 >= fit_pc1_x1 - 2*se_pc1_x1 & truth_pc1_x1 <= fit_pc1_x1 + 2*se_pc1_x1)
  covr_pc2_f1 <- mean(truth_pc2_x1 >= fit_pc2_x1 - 2*se_pc2_x1 & truth_pc2_x1 <= fit_pc2_x1 + 2*se_pc2_x1)
  covr_stz_f2 <- mean(truth_stz_x2 >= fit_stz_x2 - 2*se_stz_x2 & truth_stz_x2 <= fit_stz_x2 + 2*se_stz_x2)
  covr_pc1_f2 <- mean(truth_pc1_x2 >= fit_pc1_x2 - 2*se_pc1_x2 & truth_pc1_x2 <= fit_pc1_x2 + 2*se_pc1_x2)
  covr_pc2_f2 <- mean(truth_pc2_x2 >= fit_pc2_x2 - 2*se_pc2_x2 & truth_pc2_x2 <= fit_pc2_x2 + 2*se_pc2_x2)
  
  # Return #
  df1 <- data.frame(
    id = id,
    n = n,
    var = 'x1',
    RMSE_stz = RMSE_stz_f1,RMSE_pc1 = RMSE_pc1_f1,RMSE_pc2 = RMSE_pc2_f1,
    SE_stz = SE_stz_f1,SE_pc1 = SE_pc1_f1,SE_pc2 = SE_pc2_f1,
    covr_stz = covr_stz_f1,covr_pc1 = covr_pc1_f1,covr_pc2 = covr_pc2_f1
  )
  df2 <- data.frame(
    id = id,
    n = n,
    var = 'x2',
    RMSE_stz = RMSE_stz_f2,RMSE_pc1 = RMSE_pc1_f2,RMSE_pc2 = RMSE_pc2_f2,
    SE_stz = SE_stz_f2,SE_pc1 = SE_pc1_f2,SE_pc2 = SE_pc2_f2,
    covr_stz = covr_stz_f2,covr_pc1 = covr_pc1_f2,covr_pc2 = covr_pc2_f2
  )
  bind_rows(df1,df2)
}

## Run the simulations ##
simstodo <- expand.grid(n=n,id=1:B)
simlist <- vector(mode='list',length=nrow(simstodo))
for (i in 1:nrow(simstodo)) simlist[[i]] <- simstodo[i, ]
cat("Doing",nrow(simstodo),"total simulations...\n")
tm <- Sys.time()
set.seed(476329)
mc.reset.stream()
sims <- mclapply(simlist,do_a_simulation)
simframe <- bind_rows(sims) %>% as_tibble()
dt <- round(as.numeric(difftime(Sys.time(),tm,units='secs')))
cat("Finished simulations, they took",dt,"seconds.\n")

## Summarize results ##
results <- simframe %>%
  group_by(var,n) %>%
  summarize(across(-contains(c("id","var")),list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) %>%
  mutate(across(contains("covr"),~.x*100))

# Print nicely for paper
knitr::kable(results,format = 'markdown',digits = 2)
# 
# 
# |var |    n| RMSE_stz_mean| RMSE_stz_sd| RMSE_pc1_mean| RMSE_pc1_sd| RMSE_pc2_mean| RMSE_pc2_sd| SE_stz_mean| SE_stz_sd| SE_pc1_mean| SE_pc1_sd| SE_pc2_mean| SE_pc2_sd| covr_stz_mean| covr_stz_sd| covr_pc1_mean| covr_pc1_sd| covr_pc2_mean| covr_pc2_sd|
# |:---|----:|-------------:|-----------:|-------------:|-----------:|-------------:|-----------:|-----------:|---------:|-----------:|---------:|-----------:|---------:|-------------:|-----------:|-------------:|-----------:|-------------:|-----------:|
# |x1  |  200|          0.14|        0.05|          0.18|        0.07|          0.19|        0.08|        0.15|      0.01|        0.20|      0.01|        0.21|      0.02|         96.30|        7.92|         97.65|        6.27|         96.34|        9.44|
# |x1  |  500|          0.10|        0.03|          0.13|        0.05|          0.13|        0.05|        0.11|      0.00|        0.14|      0.01|        0.14|      0.01|         96.75|        6.31|         97.24|        7.12|         96.99|        7.93|
# |x1  | 1000|          0.07|        0.02|          0.10|        0.04|          0.10|        0.04|        0.08|      0.00|        0.10|      0.00|        0.11|      0.00|         96.68|        5.92|         97.10|        7.45|         96.98|        7.38|
# |x2  |  200|          0.14|        0.04|          0.19|        0.08|          0.19|        0.08|        0.15|      0.01|        0.19|      0.01|        0.20|      0.01|         95.71|        8.22|         95.74|        9.66|         95.11|       10.77|
# |x2  |  500|          0.10|        0.03|          0.13|        0.05|          0.13|        0.05|        0.10|      0.00|        0.13|      0.01|        0.14|      0.01|         96.07|        6.86|         96.26|        8.60|         95.77|        9.37|
# |x2  | 1000|          0.07|        0.02|          0.09|        0.03|          0.10|        0.04|        0.08|      0.00|        0.10|      0.00|        0.10|      0.00|         96.39|        6.44|         96.77|        7.46|         96.27|        8.77|
# 
# 
# |var |    n| RMSE_stz_mean| RMSE_stz_sd| RMSE_pc1_mean| RMSE_pc1_sd| RMSE_pc2_mean| RMSE_pc2_sd| SE_stz_mean| SE_stz_sd| SE_pc1_mean| SE_pc1_sd| SE_pc2_mean| SE_pc2_sd| covr_stz_mean| covr_stz_sd| covr_pc1_mean| covr_pc1_sd| covr_pc2_mean| covr_pc2_sd|
# |:---|----:|-------------:|-----------:|-------------:|-----------:|-------------:|-----------:|-----------:|---------:|-----------:|---------:|-----------:|---------:|-------------:|-----------:|-------------:|-----------:|-------------:|-----------:|
# |x1  |  200|          0.14|        0.05|          0.18|        0.07|          0.19|        0.08|        0.15|      0.01|        0.20|      0.01|        0.21|      0.02|         96.30|        7.92|         97.65|        6.27|         96.34|        9.44|
# |x1  |  500|          0.10|        0.03|          0.13|        0.05|          0.13|        0.05|        0.11|      0.00|        0.14|      0.01|        0.14|      0.01|         96.75|        6.31|         97.24|        7.12|         96.99|        7.93|
# |x1  | 1000|          0.07|        0.02|          0.10|        0.04|          0.10|        0.04|        0.08|      0.00|        0.10|      0.00|        0.11|      0.00|         96.68|        5.92|         97.10|        7.45|         96.98|        7.38|
# |x2  |  200|          0.14|        0.04|          0.19|        0.08|          0.19|        0.08|        0.15|      0.01|        0.19|      0.01|        0.20|      0.01|         95.71|        8.22|         95.74|        9.66|         95.11|       10.77|
# |x2  |  500|          0.10|        0.03|          0.13|        0.05|          0.13|        0.05|        0.10|      0.00|        0.13|      0.01|        0.14|      0.01|         96.07|        6.86|         96.26|        8.60|         95.77|        9.37|
# |x2  | 1000|          0.07|        0.02|          0.09|        0.03|          0.10|        0.04|        0.08|      0.00|        0.10|      0.00|        0.10|      0.00|         96.39|        6.44|         96.77|        7.46|         96.27|        8.77|

# Save
write_csv(results,file = file.path(resultspath,"gaussian-simulation-results.csv"))

cat("Done. View the results at",resultspath,"\n")




