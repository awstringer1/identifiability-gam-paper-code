### Simulations: Bernoulli Data, Optimal Constraints ###
# Find the optimal constraints for simulated Poisson data, and asses
# how much lower the variance is.
# For paper: Identifiability Constraints in Generalized Additive Models
# Section 4.2
# Alex Stringer
# 2022/05

## Load Libraries ----
library(tidyverse)
library(here)

library(mgcv)
library(TMB)
precompile()
library(Matrix)
library(parallel)
options(mc.cores = parallel::detectCores())

# Reproducible random number generation
RNGkind("L'Ecuyer-CMRG")
## Global parameters ----
# Number of simulations to do
B <- 100
# Sample sizes to use
n <- c(50,100,200)
# True smooth function
truefunc <- function(x) 3*sin(2*pi*x)
# Intercept
alpha <- 0
# Number of points to predict at, for calculating variance
pn <- 100
# Smooth specification- thin plate regression spline with defaults
smooth <- list(s(x1))

## Set Paths ----

resultspath <- tempdir()
if (!dir.exists(resultspath)) dir.create(resultspath)
tmbpath <- file.path(here::here(),"tmb/bernoulli_gam")

## Compile TMB Code ----
compile(paste0(tmbpath,".cpp"))
dyn.load(dynlib(tmbpath))

## Function to fit a Bernoulli GAM ----
# Get constraint matrix
getZ <- function(cc) {
  if (is.matrix(cc)) {
    dd <- ncol(cc)
    out <- vector(mode='list',length=dd)
    for (j in 1:dd) out[[j]] <- getZ(cc[ ,j])
    return(out)
  }
  p <- length(cc)
  if (norm(cc,type='2') < sqrt(.Machine$double.eps)) return(diag(p)[ ,2:p])
  normc <- norm(cc,type='2')
  r <- c(normc,rep(0,p-1))
  Q <- diag(p) - 2*tcrossprod(cc-r)/norm(cc-r,type='2')^2
  Q[ ,2:p]
}
fit_bernoulli_gam <- function(smooth,dat,preddat=dat,constr,method=c('BFGS','trust'),...) {
  ## SETUP CONTROL ##
  method <- method[1]
  absorbcons <- FALSE # Do NOT change this.
  verbose <- FALSE # why not
  ## END SETUP CONTROL ##
  
  ## MODEL SETUP ##
  # 1: SMOOTHS #
  Xr <- Xf <- mam::newsparsemat(nrow(dat),0)
  Xrpred <- Xfpred <- mam::newsparsemat(nrow(preddat),0)
  numsmooth <- 0
  r <- 0
  if (!is.null(smooth)) {
    tm <- Sys.time()
    if (verbose) cat("Constructing smooths... ")
    if (!inherits(smooth,'list')) smooth <- list(smooth)
    # Conditional model
    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=dat,absorb.cons = absorbcons),'[[',1)
    XXpredlist <- lapply(SS,mgcv::PredictMat,data = preddat)
    numsmooth <- length(smooth) # Number of smooth terms
    
    # Design matrices
    Xlist <- lapply(SS,'[[','X')
    
    ## CONSTRAINTS ##
    
    # Form the constrained design matrix and penalty matrix
    if (is.character(constr)) {
      if (constr == 'stz') {
        # Sum to zero
        cc <- lapply(Xlist,colSums)
        Zlist <- lapply(cc,getZ)
      }
    } else if (is.numeric(constr)) {
      # User-provided numeric vector or matrix
      if (is.matrix(constr)) {
        Zlist <- getZ(constr) # Returns a list
      } else {
        Zlist <- list(getZ(constr)) # Force list, for downstream consistency
      }
    } else {
      stop("Unknown provision of constraints.")
    }
    # List of constrained design matrices
    XZlist <- mapply(function(x,y) x %*% y,Xlist,Zlist,SIMPLIFY = FALSE)
    XZpredlist <- mapply(function(x,y) x %*% y,XXpredlist,Zlist,SIMPLIFY = FALSE)
    
    # List of constrained penalty matrices
    Plist <- lapply(lapply(SS,'[[','S'),'[[',1)
    PZlist <- mapply(function(x,y) crossprod(x,crossprod(y,x)),Zlist,Plist,SIMPLIFY = FALSE)
    P <- bdiag(Plist)
    PZ <- bdiag(PZlist)
    
    X <- Reduce(cbind,Xlist)
    XZ <- Reduce(cbind,XZlist)
    Xpred <- Reduce(cbind,XXpredlist)
    XZpred <- Reduce(cbind,XZpredlist)
    
    ## END CONSTRAINTS ##
    
    EE <- lapply(PZlist,eigen)
    
    p <- sapply(lapply(EE,'[[','vectors'),ncol)
    r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps))
    m <- p-r
    URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
    UFlist <- mapply(
      function(x,y,z) {
        if (y<z) return(x[ ,(1+y):z])
        mam::newsparsemat(z,z)
      },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
    URlist <- lapply(URlist,cbind) # Ensure they stay matrices
    UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices
    
    UR <- Matrix::bdiag(URlist)
    UF <- Matrix::bdiag(UFlist)
    # if m=1 UF gets coerced to numeric
    if (!is.matrix(UF)) UF <- cbind(UF)
    
    Dpi <- Matrix::Diagonal(sum(r),1 / sqrt(Reduce(c,lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps]))))
    
    Xr <- as.matrix(XZ %*% UR %*% Dpi)
    Xf <- as.matrix(XZ %*% UF)
    # ADD the intercept (!)
    Xf <- cbind(1,Xf)
    
    
    # Prediction
    # Use the SAME Ur and Uf
    Xrpred <- as.matrix(XZpred %*% UR %*% Dpi)
    Xfpred <- cbind(1,as.matrix(XZpred %*% UF))
    Xpred <- cbind(Xfpred,Xrpred)
    
    dt <- difftime(Sys.time(),tm,units = 'secs')
    if (verbose) cat("finished, took",round(dt),"seconds.\n")
  }
  # END SMOOTHS #
  ## END MODEL SETUP ##
  
  ## MODEL ##
  
  # Fit with TMB
  tm <- Sys.time()
  if (verbose) cat("Fitting model... ")
  tmbdat <- list(
    XF = as.matrix(Xf),
    XR = as.matrix(Xr),
    y = dat$y, # Response
    p = as.integer(numsmooth), # Number of smooth terms
    r = r # Rank of each smooth
  )
  tmbparams <- with(tmbdat,list(
    betaF = rep(0,ncol(XF)),
    bR = rep(0,ncol(XR)),
    logsmoothing = rep(0,p)
  ))
  rand <- c('bR','betaF')
  
  template <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = rand,
    silent = TRUE,
    DLL = 'bernoulli_gam'
  )
  
  if (method == 'BFGS') {
    opt <- with(template,stats::optim(par,fn,gr,method='BFGS',hessian=TRUE))
    # Rename output to match trustOptim
    opt$solution <- opt$par
  } else if (method == 'trust') {
    utils::capture.output(opt <- with(template,trustOptim::trust.optim(
      x = par,
      fn = fn,
      gr = gr,
      hs = function(x) as(numDeriv::jacobian(gr,x),'dgCMatrix'),
      method = 'Sparse'
    )))
  } else {
    stop(paste0("Unrecognized optimization method:",method,"\n"))
  }
  
  # Point estimates
  tmbcoefs <- with(template$env,last.par[random])
  tmbbetaF <- tmbcoefs[names(tmbcoefs)=='betaF']
  tmbbR <- tmbcoefs[names(tmbcoefs)=='bR']
  
  condcoefs <- c(tmbbetaF,tmbbR)
  condest <- as.numeric(Xpred[ ,-1] %*% condcoefs[-1])
  
  # compute variances & SEs
  H <- TMB::sdreport(template,getJointPrecision = TRUE)$jointPrecision
  condidx <- which(rownames(H) %in% c('betaF','bR'))[-1] # No intercept- H is in the order of (betaF,bR)
  estvar <- diag(Xpred[ ,-1] %*% solve(H)[condidx,condidx] %*% t(Xpred[ ,-1]))
  
  dt <- difftime(Sys.time(),tm,units = 'secs')
  if (verbose) cat("finished, took",round(dt),"seconds.\n")
  ## END CONDITIONAL MODEL ##
  list(
    est = condest,
    sd = sqrt(estvar)
  )
}


## Optimize over constraints ##
# Function to do a single optimization
fit_model_constr <- function(cc,dat,preddat) {
  cc <- cc / norm(cc,type='2')
  # Fit model subject to t(cc)%*%beta = 0
  vv <- tryCatch(fit_bernoulli_gam(smooth,dat,preddat,constr = cc,method='BFGS')$sd^2,error=function(e) e)
  if (inherits(vv,'condition')) return(999)
  mean(sqrt(vv))
}

doopt <- function(lst) {
  # lst: list with elements:
  #      n: sample size
  #      id: unique identifier of the specific simulation, id = 1...B
  
  n <- lst$n
  id <- lst$id
  cat("Doing simulation, id =",id,"n =",n,"\n")
  
  # Generate a data set
  x1 <- runif(n)
  eta <- truefunc(x1) + alpha
  y <- rbinom(n,1,mam::ilogit(eta))
  dat <- data.frame(x1=x1,y=y)
  preddat <- data.frame(
    x1 = seq(min(x1),max(x1),length.out=pn)
  )
  
  # Get the sum-to-zero
  SS <- lapply(lapply(smooth,mgcv::smoothCon,data=dat,absorb.cons = FALSE),'[[',1)
  Xlist <- lapply(SS,'[[','X')
  X <- Reduce(cbind,Xlist)
  
  opt <- tryCatch(optim(
    colSums(X),
    fit_model_constr,
    method = 'Nelder-Mead',
    control = list(trace=0,maxit=1000,reltol=1e-08),
    dat = dat,
    preddat = preddat
  ),error = function(e) e)
  if (inherits(opt,'condition')) {
    stz <- optimal <- -1
  } else {
    stz <- fit_model_constr(colSums(X),dat,preddat)
    optimal <- opt$value
  }
  data.frame(
    id = id,
    n = n,
    stz = stz,
    optimal = optimal
  )
}

# Do the simulations
simstodo <- expand.grid(n=n,id=1:B)
simlist <- vector(mode='list',length=nrow(simstodo))
for (i in 1:nrow(simstodo)) simlist[[i]] <- simstodo[i, ]
set.seed(5206726)
mc.reset.stream()
cat("Doing",nrow(simstodo),"total simulations...\n")
tm <- Sys.time()
sims <- mclapply(simlist,doopt)
simframe <- bind_rows(sims) %>% as_tibble()
dt <- round(as.numeric(difftime(Sys.time(),tm,units='secs')))
cat("Finished simulations, they took",dt,"seconds.\n")

# Save results
write_csv(simframe,file.path(resultspath,"bernoulli-optconstraint-sims.csv"))
# simframe <- read_csv(file.path(resultspath,"bernoulli-optconstraint-sims.csv"))

# Summarize results
results <- simframe %>%
  filter(stz > -1) %>%
  mutate(sddiff = stz - optimal) %>%
  group_by(n) %>%
  summarize(mn = mean(sddiff),se = sd(sddiff)) %>%
  mutate(across(mn:se,~.x)) %>%
  knitr::kable(
    digits = 5,
    format = 'markdown'
  )

results
# 
# 
# |   n|      mn|      se|
# |---:|-------:|-------:|
# |  50| 0.65382| 5.66328|
# | 100| 0.01050| 0.03926|
# | 200| 0.00143| 0.00102|
