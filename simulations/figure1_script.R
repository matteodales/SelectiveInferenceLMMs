## script for running simulation of time dependency to number of variables on server


args <- commandArgs(trailingOnly = TRUE)  # Read command-line arguments
if (length(args) == 0) {
  stop("No arguments provided!")
}
itnum <- as.numeric(args[1])  # Convert the first argument to numeric


library(lme4)
library(lmerTest)
library(expm)
library(ggplot2)
library(dplyr)
library(reshape2)
library(mvtnorm)
library(MASS)
library(R.utils)
library(MultBiplotR)
library(foreach)
library(parallel)
library(mgcv)
library(cowplot)
library(glmmLasso)
library(cv.glmmLasso)
library(remotes)

library(grid)
library(glmnet)
library(gridExtra)

library(postcAIC)
library(ks)
library(tmg)


# simulation settings 

plot_results = TRUE

seed = 0

sim_num <- 100 # number of simulations
its <- 1 # number of iterations to compute fdr in each simulation
tot_its <- its*sim_num

fdr_level <- 0.05
set.seed(itnum)

selection_frac <- 0.5 # fraction of observations to use for selecting the variables

#####################################################
######## generating the random intercept data #######
#####################################################

n <- 200
nrsubj <- 40
nrobs <- 5

list_p <- list(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 35, 40, 50, 60, 75)

time_results_df <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(time_results_df) <- c('method', 'p', 'time')












for(p in list_p){
  
  
  naive_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(naive_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  naive_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(naive_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  datasplitting_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(datasplitting_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  datasplitting_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(datasplitting_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  UVILassoLMM_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(UVILassoLMM_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  UVILassoLMM_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(UVILassoLMM_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  selfmadelasso_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(selfmadelasso_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  selfmadelasso_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(selfmadelasso_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  selfmadestep_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(selfmadestep_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  selfmadestep_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(selfmadestep_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  postcAIC_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(postcAIC_results_df) <- c('method', 'tpr', 'fdr', 'fdr', 'coverage', 'avg_ci')
  
  postcAIC_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(postcAIC_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  
  
  n <- 200
  nrsubj <- 40
  nrobs <- 5
  p_rel = 3 # number of non-zero coefficients
  q = 0 # number of mixed effects, other than the intercept
  
  print(paste0("p ", p))
  
  p_rel = 3 # number of non-zero coefficients
  q = 0 # number of mixed effects, other than the intercept
  
  
  
  SNR = 1 #signal-to-noise ratio
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
  var_names <- c(paste("X", 1:p, sep=""))
  
  
  
  #####################################################
  ################## generating data  #################
  #####################################################
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
  var_names <- c(paste("X", 1:p, sep=""))
  
  
  Xs <- list()
  real_betas <- list()
  yss <- list()
  Zs <- list()
  
  
  for(i in 1:sim_num){
    
    real_beta = c(sample (c(1,2,-1,-2), size=p_rel, replace=T), rep(0, p-p_rel))
    X <- cbind(mvrnorm(n = n, mu = rep(0,p), Sigma=Sigmap))
    colnames(X) <- var_names
    
    sig_e = sd(X%*%real_beta)/SNR #sd of random noise
    sig_v = sig_e
    
    ys <- matrix(0, ncol = its, nrow = n)
    
    for(j in 1:its){
      
      eij <- rnorm(n, mean = 0, sd = sig_e)
      if(q==0){
        vi <- rnorm(nrsubj, mean = 0, sd = sig_v)
        y <- X%*%real_beta + eij + vi[subjind]
        ys[,j] <- y
        Z <- matrix(nrow=n,ncol=1)
      }else{
        b <- mvrnorm(n = nrsubj, mu = rep(0,q+1), Sigma=sig_v^2*varcovar)  # Simulate random effects
        Z <- cbind(rep(1,n), mvrnorm(n = n, mu = rep(0,q), Sigma=Sigmaq)) # Make the random design matrix
        colnames(Z) <- c("(Intercept)", paste("Z", 1:q, sep=""))
        Zb <- rowSums(b[subjind,]*Z)
        y <- X%*%real_beta + eij + Zb
        ys[,j] <- y
      }
      
    }
    
    Xs[[i]] <- X[,]
    real_betas[[i]] <- real_beta
    Zs[[i]] <- Z[,-1]
    yss[[i]] <- ys
    
  }
  
  
  
  
  form <- as.formula(
    paste(" ~ -1 + ", paste("X", 1:p, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
  ) 
  
  yform <- as.formula(
    paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
  )
  
  if(q>0){
    form <- as.formula(
      paste(" ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
    )
    
    yform <<- as.formula(
      paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
    )
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  for(jj in 1:sim_num){
    
    
    n <- 200
    nrsubj <- 40
    nrobs <- 5
    p_rel = 3 # number of non-zero coefficients
    q = 0 # number of mixed effects, other than the intercept
    
    
    
    SNR = 1 #signal-to-noise ratio
    
    Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
    Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
    
    subjind <- as.factor(rep(1:nrsubj, each=nrobs))
    varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
    var_names <- c(paste("X", 1:p, sep=""))
    
    
    
    ### generate the data
    
    Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
    Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
    
    subjind <- as.factor(rep(1:nrsubj, each=nrobs))
    varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
    var_names <- c(paste("X", 1:p, sep=""))
    
    
    
    
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)
    
    
    
    
    
    
    
    print('Naive inference')
    
    y <- ys[,1]
    
    dat <- data.frame(X, Z, y, subjind)
    
    t1 <- Sys.time()
    
    sx <- as.matrix(scale(X))
    sy <- as.vector(y)
    
    lambdaMax <- max(abs(colSums(sx*sy)))
    lambda_vec <- exp(seq(from = log(lambdaMax),
                          to = log(lambdaMax * 0.0001),
                          length.out = 50))
    
    PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose=FALSE)
    Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
    Q.start<-as.numeric(VarCorr(PQL)[1,1])
    
    BIC_vec <- rep(Inf, length(lambda_vec))
    sel_vec <- matrix(0, p, length(lambda_vec))
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = dat,
        lambda = lambda_vec[l],
        control = list(center = FALSE, start=Delta.start[l,], q_start = Q.start[l])
      ), silent=TRUE)
      
      if (!inherits(glm3, "try-error")) {
        BIC_vec[l] <- glm3$bic
        sel_vec[,l] <- as.numeric(glm3$coefficients!=0)
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      }
    }
    
    opt3 <- which.min(BIC_vec)
    
    selected <- sel_vec[,opt3]
    
    
    selected_tot_names <- var_names[selected==1]
    
    selected_tot_yform <- as.formula(
      paste(" y ~ -1 +", paste(selected_tot_names, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
    )
    
    if(q>0){
      selected_tot_yform <- as.formula(
        paste(" y ~ -1 +", paste(selected_tot_names, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
      )
    }
    
    suppressWarnings(suppressMessages(sel_mod <- lmer(formula = selected_tot_yform, data = dat)))
    
    
    if(sum(selected)==1){
      
      pvals_tot <- as.data.frame(t(coef(summary(sel_mod))[,c(1,5)]))
      suppressMessages(confint_tot <- confint(sel_mod, method='Wald', parm = 'beta_'))
      confint_tot <- as.data.frame(confint_tot)
      
    }else{
      
      pvals_tot <- as.data.frame(coef(summary(sel_mod))[,c(1,5)])
      suppressMessages(confint_tot <- confint(sel_mod, method='Wald', parm = 'beta_'))
      confint_tot <- as.data.frame(confint_tot)
      
    }
    
    k <- 1
    for(i in 1:(p)){
      
      if(selected[i] == 1){
        naive_pvals_df[nrow(naive_pvals_df) + 1,] <- c('Naive', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], pvals_tot[k,1], pvals_tot[k,2],
                                                       confint_tot[k,1], confint_tot[k,2], as.numeric(confint_tot[k,1]<=real_beta[i] & confint_tot[k,2]>=real_beta[i]) )
        k = k + 1
      }
    }
    
    naive_pvals_df[,3:9] <- lapply(naive_pvals_df[,3:9],as.numeric)
    
    coverage <- mean(naive_pvals_df[(nrow(naive_pvals_df)-sum(selected)):nrow(naive_pvals_df), 9])
    avg_CI <- mean(naive_pvals_df[(nrow(naive_pvals_df)-sum(selected)):nrow(naive_pvals_df), 8] - naive_pvals_df[(nrow(naive_pvals_df)-sum(selected)):nrow(naive_pvals_df), 7])
    
    
    corrected_pvals_tot <- p.adjust(pvals_tot[,2], method = 'BH')
    
    selected[selected==1] <- corrected_pvals_tot<=fdr_level
    tot_metrics <- metrics(selected, real_beta!=0)
    
    naive_results_df[nrow(naive_results_df) + 1,] <- c('Naive', tot_metrics$tpr, tot_metrics$fdr, tot_metrics$fdr, coverage, avg_CI)
    
    time_results_df[nrow(time_results_df) + 1,] <- c('Naive', p, difftime(Sys.time(),t1,units = "secs"))
    
    
    
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    y <- ys[,1]
    dat <- data.frame(X, Z, y, subjind)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    print('Data Splitting')
    
    #data splitting: the data is split in two, using the first half for selection and the second for inference
    
    t2 <- Sys.time()
    
    selection_sample <- sample(1:nrsubj, round(nrsubj*selection_frac))
    selection_dat <- dat[dat$subjind %in% selection_sample,]
    inference_dat <- dat[!(dat$subjind %in% selection_sample),]
    
    X_sel <- as.matrix(selection_dat[,1:p])
    y <- selection_dat$y
    
    sx <- as.matrix(scale(X_sel))
    sy <- as.vector(y)
    
    lambdaMax <- max(abs(colSums(sx*sy)))
    lambda_vec <- exp(seq(from = log(lambdaMax),
                          to = log(lambdaMax * 0.001),
                          length.out = 50))
    
    selection_dat$subjind <- droplevels(selection_dat$subjind)
    
    
    PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=selection_dat,verbose=FALSE)
    Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
    Q.start<-as.numeric(VarCorr(PQL)[1,1])
    
    BIC_vec <- rep(Inf, length(lambda_vec))
    sel_vec <- matrix(0, p, length(lambda_vec))
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = selection_dat,
        lambda = lambda_vec[l],
        control = list(start = Delta.start[l, ], q_start = Q.start[l], center = FALSE)
      ), silent=TRUE)
      
      if (!inherits(glm3, "try-error")) {
        BIC_vec[l] <- glm3$bic
        sel_vec[,l] <- as.numeric(glm3$coefficients!=0)
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      }
    }
    
    opt3 <- which.min(BIC_vec)
    
    selected <- sel_vec[,opt3]
    
    selected_split_names <- var_names[selected==1]
    
    selected_split_yform <- as.formula(
      paste(" y ~ -1 +", paste(selected_split_names, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
    )
    
    if(q>0){
      selected_split_yform <- as.formula(
        paste(" y ~ -1 +", paste(selected_split_names, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
      )
    }
    
    suppressWarnings(suppressMessages(sel_mod <- lmer(formula = selected_split_yform, data = inference_dat)))
    
    if(sum(selected)==1){
      
      pvals_split <- as.data.frame(t(coef(summary(sel_mod))[,c(1,5)]))
      suppressMessages(confint_split <- confint(sel_mod, method='Wald', parm='beta_'))
      confint_split  <- as.data.frame(confint_split )
      
    }else{
      
      pvals_split  <- as.data.frame(coef(summary(sel_mod))[,c(1,5)])
      suppressMessages(confint_split <- confint(sel_mod, method='Wald', parm='beta_'))
      confint_split  <- as.data.frame(confint_split )
      
    }
    
    
    k <- 1
    for(i in 1:(p)){
      
      if(selected[i] == 1){
        datasplitting_pvals_df[nrow(datasplitting_pvals_df) + 1,] <- c('Data Splitting', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], pvals_split[k,1], pvals_split[k,2],
                                                                       confint_split[k,1], confint_split[k,2], as.numeric(confint_split[k,1]<=real_beta[i] & confint_split[k,2]>=real_beta[i]) )
        k = k + 1
      }
    }
    
    datasplitting_pvals_df[,3:9] <- lapply(datasplitting_pvals_df[,3:9],as.numeric)
    
    
    
    coverage <- mean(datasplitting_pvals_df[(nrow(datasplitting_pvals_df)-sum(selected)):nrow(datasplitting_pvals_df), 9])
    avg_CI <- mean(datasplitting_pvals_df[(nrow(datasplitting_pvals_df)-sum(selected)):nrow(datasplitting_pvals_df), 8] - datasplitting_pvals_df[(nrow(datasplitting_pvals_df)-sum(selected)):nrow(datasplitting_pvals_df), 7])
    
    corrected_pvals_split <- p.adjust(pvals_split[,2], method = 'BH')
    
    selected[selected==1] <- corrected_pvals_split<=fdr_level
    split_metrics <- metrics(selected, real_beta!=0)
    
    datasplitting_results_df[nrow(datasplitting_results_df) + 1,] <- c('Data splitting', split_metrics$tpr, split_metrics$fdr, split_metrics$fdr, coverage, avg_CI)
    
    time_results_df[nrow(time_results_df) + 1,] <- c('Data splitting', p, difftime(Sys.time(),t2,units = "secs"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    

  ncp_fun <- function(p, Cinv, lambdas) {
    max_val <- -Inf  # Initialize the maximum value

    # Iterate over all possible sign combinations (2^p combinations)
    for (i in 0:(2^p - 1)) {
      # Generate the sign vector using binary representation
      sign_vector <- ifelse(as.integer(intToBits(i))[1:p] == 1, 1, -1)

      # Calculate the function for the current sign combination
      current_val <- t(sign_vector) %*% Cinv %*% sign_vector

      # Update the maximum value if needed
      if (current_val > max_val) {
        max_val <- current_val
      }
    }

    # Multiply by the lambda^2 factor
    ncp <- lambdas[1]^2 * max_val
    return(ncp)
  }


  print('UVILassoLMM')

    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    y <- ys[,1]

    dat <- data.frame(X, Z, y, subjind)

    t3 <- Sys.time()


      mod <- lmer(formula = yform, data = dat)
      V <- getME(mod,'Z')%*%(sigma(mod)^2*getME(mod,'Lambda')%*%getME(mod,'Lambdat'))%*%t(getME(mod,'Z'))+ diag(sigma(mod)^2,n)
      V_inv <- solve(V)

      V_menunmezz <- matrixsqrtinv(V)
      cv_mod <- cv.glmnet(x = V_menunmezz%*%X, y = V_menunmezz%*%y, intercept=FALSE, standardize=FALSE)
      beta <- coef(cv_mod, s = cv_mod$lambda.min)[-1]

      # what is C
      C = as.matrix(t(X)%*%V_inv%*%X)/n
      Cinv <- solve(C)

      kram_lambda <- cv_mod$lambda.min*n

      lambdas <- c(rep(kram_lambda/sqrt(n),p))

      ncp <- ncp_fun(p,Cinv,lambdas)


      pvals <- rep(0,p)
      k <- 1

      for(i in seq_len(p)){


        test_stat <- beta[i]**2 / Cinv[i,i] * n
        pvals[k] <- pchisq(test_stat, 1, ncp=ncp, lower.tail = FALSE)


        cil <- beta[i] - sqrt(qchisq(1 - fdr_level, 1,  ncp=ncp)*Cinv[i,i]/n)
        ciu <- beta[i] + sqrt(qchisq(1 - fdr_level, 1, ncp=ncp)*Cinv[i,i]/n)


        UVILassoLMM_pvals_df[nrow(UVILassoLMM_pvals_df)+1,] <- c('UVILassoLMM', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i],
                                                               beta[i], pvals[k], cil, ciu, as.numeric(cil<=real_beta[i] & ciu>=real_beta[i]))

        k <- k+1
      }

      UVILassoLMM_pvals_df[,3:9] <- lapply(UVILassoLMM_pvals_df[,3:9],as.numeric)

      coverage <- mean(UVILassoLMM_pvals_df[(nrow(UVILassoLMM_pvals_df)-p+1):nrow(UVILassoLMM_pvals_df), 9])
      avg_CI <- mean(UVILassoLMM_pvals_df[(nrow(UVILassoLMM_pvals_df)-p):nrow(UVILassoLMM_pvals_df), 8] - UVILassoLMM_pvals_df[(nrow(UVILassoLMM_pvals_df)-p):nrow(UVILassoLMM_pvals_df), 7])

      pvals <- p.adjust(pvals, method = 'BH')

      selected <- pvals<=fdr_level
      UVILassoLMM_metrics <- metrics(selected, real_beta!=0)

      UVILassoLMM_results_df[nrow(UVILassoLMM_results_df) + 1,] <- c('UVILassoLMM', UVILassoLMM_metrics$tpr,UVILassoLMM_metrics$fdr,
                                                                   UVILassoLMM_metrics$fdr, coverage, avg_CI)

      time_results_df[nrow(time_results_df) + 1,] <- c('UVILassoLMM', p, difftime(Sys.time(),t3,units = "secs"))




        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    fixed_form <<- as.formula(
      paste("y ~ - 1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")
    )
    
    modFun <- function(yy)
    {
      dat$y <- as.numeric(yy)
      lmer(yform, REML = FALSE, data = dat)
    }
    
    selFun <- function(yy){
      
      dat$y <- as.numeric(yy)
      
      
      
      PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose='FALSE')
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])
      BIC_vec <- c(rep(Inf,length(lambda_vec)))
      sel_vec <- matrix(0,length(lambda_vec),p)
      
      for (l in 1:length(lambda_vec)) {
        glm3 <- try(glmmLasso(
          fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
          rnd = list(subjind = ~1),
          data = dat,
          lambda = lambda_vec[l],
          control = list(center = FALSE, start=Delta.start[l,], q_start = Q.start[l])
        ), silent=TRUE)
        
        if (!inherits(glm3, "try-error")) {
          BIC_vec[l] <- glm3$bic
          sel_vec[l,] <- as.numeric(glm3$coefficients != 0)
          Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
          Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
        }
      }
      
      opt <- which.min(BIC_vec)
      
      names_vec = var_names[sel_vec[opt,]==1]
      
      return(c(names_vec))
    }
    
    
    
    print('selfmade-lasso')
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    y <- ys[,1]
    
    dat <- data.frame(X, Z, y, subjind)
    
    t4 <- Sys.time()
    
    
    sx <- as.matrix(scale(X))
    sy <- as.vector(y)
    
    lambdaMax <- max(abs(colSums(sx*sy)))
    lambda_vec <<- exp(seq(from = log(lambdaMax),
                           to = log(lambdaMax * 0.0001),
                           length.out = 50))
    
    
    names_vec = selFun(y)
    
    if(length(names_vec)==0){
      selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', 0, 0, 0, NA, NA)
    }else{
      
      
      sel_form <- as.formula(
        paste("y ~ -1 + ", paste(names_vec[1:length(names_vec)], collapse = "+"), "+ (1 |subjind)")
      )
      if(q>0){
        sel_form <- as.formula(
          paste("y ~ -1 + ", paste(names_vec[1:length(names_vec)], collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)")
        )
      }
      
      final_model <- suppressWarnings(suppressMessages(lmer(formula = sel_form, data = dat)))
      
      beta <- fixef(final_model)
      selection = c(names(fixef(final_model)))
      
      
      #define function which checks congruency
      checkFun <- function(yb){
        
        setequal( selFun(yy = yb),
                  selection )
        
      }
      
      r <- suppressWarnings(mocasin(mod = final_model,
                                    this_y = y,
                                    checkFun = checkFun,
                                    nrSamples = 100,
                                    which = 1:length(selection),
                                    bayesian = FALSE,
                                    conditional = FALSE,
                                    efficient = TRUE,
                                    varForSampling = 'est',
                                    VCOV_sampling = NULL,
                                    VCOV_vT = NULL,
                                    trace=FALSE))
      
      r <- do.call("rbind", r$selinf)
      r$variable <- selection
      
      selected <- rep(0,p)
      
      k <- 1
      for(i in seq_len(p)){
        
        if(var_names[i] %in% selection){
          
          selfmadelasso_pvals_df[nrow(selfmadelasso_pvals_df)+1,] <- c('selfmade-lasso', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], beta[k],
                                                                     r$pval[k], r$cil[k], r$ciu[k], as.numeric(r$cil[k]<=real_beta[i] & r$ciu[k]>=real_beta[i]))
          
          selected[i] <- 1
          k <- k +1
          
        }
      }
      
      selfmadelasso_pvals_df[,3:9] <- lapply(selfmadelasso_pvals_df[,3:9],as.numeric)
      selfmadelasso_pvals_df[is.infinite(selfmadelasso_pvals_df[,7]),7] <- NA
      selfmadelasso_pvals_df[is.infinite(selfmadelasso_pvals_df[,8]),8] <- NA
      
      coverage <- mean(selfmadelasso_pvals_df[(nrow(selfmadelasso_pvals_df)-sum(selected)):nrow(selfmadelasso_pvals_df), 9], na.rm=TRUE)
      avg_CI <- mean(selfmadelasso_pvals_df[(nrow(selfmadelasso_pvals_df)-sum(selected)):nrow(selfmadelasso_pvals_df), 8] - selfmadelasso_pvals_df[(nrow(selfmadelasso_pvals_df)-sum(selected)):nrow(selfmadelasso_pvals_df), 7], na.rm=TRUE)
      
      pvals <- p.adjust(r$pval, method = 'BH')
      
      selected[selected==1] <- pvals<=fdr_level
      selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
      rugamer_metrics <- metrics(selected, real_beta!=0)
      
      selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fdr, coverage, avg_CI)
    }
    
    time_results_df[nrow(time_results_df) + 1,] <- c('selfmade-lasso', p, difftime(Sys.time(),t4,units = "secs"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    modFun <- function(yy, dat = dat) 
    {
      dat$y <- as.numeric(yy)
      suppressWarnings(suppressMessages(lmer(as.formula(paste("y ~ -1 + ", paste("X", 1:p, sep = "", collapse = "+"), "+ (1 |subjind)")), REML = FALSE, data = dat)))
    }
    
    selFun <- function(mod)
    {
      
      suppressWarnings(suppressMessages(attr(lmerTest:::step.lmerModLmerTest(mod, reduce.random = FALSE), "model")))
      
    }
    
    extractSelFun <- function(this_mod){
      
      if(class(this_mod)=="lm") 
        return(attr(this_mod$coefficients, "names")) else
          return(c(names(fixef(this_mod))))
      
    }
    
    
    # define function which checks congruency
    checkFun <- function(yb){ 
      
      setequal( extractSelFun(selFun(modFun(yy = yb, dat))), 
                selection )
      
    }
    
    
    print('selfmade-step')
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    
    
    y = ys[,1]
    
    dat <<- data.frame(X, Z, y, subjind)
    t5 <- Sys.time()
    
    mod <- modFun(yy=y, dat=dat)
    
    
    # suppressWarnings(suppressMessages(final_model <- attr(
    #     # use lmerTest:::step.lmerModLmerTest directly
    #     # as the packages overloading of step
    #     # does not work for the latest version
    #     lmerTest:::step.lmerModLmerTest(mod, reduce.random = FALSE), "model")))
    
    final_model <- selFun(mod)
    beta <- fixef(final_model)
    
    
    selection = extractSelFun(selFun(mod))
    
    
    # define function which checks congruency
    checkFun <- function(yb){ 
      
      setequal( extractSelFun(selFun(modFun(yy = yb, dat=dat))), 
                selection )
      
    }
    
    
    
    r <- suppressWarnings(mocasin(mod = final_model, 
                                  this_y = y,
                                  checkFun = checkFun, 
                                  nrSamples = 100,
                                  which = 1:length(selection),
                                  bayesian = FALSE,
                                  conditional = FALSE,
                                  efficient = TRUE,
                                  varForSampling = 'est',
                                  VCOV_sampling = NULL,
                                  VCOV_vT = NULL,
                                  trace=FALSE))
    
    
    r <- do.call("rbind", r$selinf)
    r$variable <- selection
    
    selected <- rep(0,p)
    
    k <- 1
    for(i in seq_len(p)){
      
      if(var_names[i] %in% selection){
        
        selfmadestep_pvals_df[nrow(selfmadestep_pvals_df)+1,] <- c('selfmade-step', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], beta[k], 
                                                             r$pval[k], r$cil[k], r$ciu[k], as.numeric(r$cil[k]<=real_beta[i] & r$ciu[k]>=real_beta[i]))
        
        selected[i] <- 1
        k <- k +1
        
      }
    }
    
    selfmadestep_pvals_df[,3:9] <- lapply(selfmadestep_pvals_df[,3:9],as.numeric)
    selfmadestep_pvals_df[is.infinite(selfmadestep_pvals_df[,7]),7] <- NA
    selfmadestep_pvals_df[is.infinite(selfmadestep_pvals_df[,8]),8] <- NA
    
    coverage <- mean(selfmadestep_pvals_df[(nrow(selfmadestep_pvals_df)-sum(selected)):nrow(selfmadestep_pvals_df), 9], na.rm=TRUE)
    avg_CI <- mean(selfmadestep_pvals_df[(nrow(selfmadestep_pvals_df)-sum(selected)):nrow(selfmadestep_pvals_df), 8] - selfmadestep_pvals_df[(nrow(selfmadestep_pvals_df)-sum(selected)):nrow(selfmadestep_pvals_df), 7], na.rm=TRUE)
    
    pvals <- p.adjust(r$pval, method = 'BH')
    
    selected[selected==1] <- pvals<=fdr_level
    selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
    rugamer_metrics <- metrics(selected, real_beta!=0)
    
    selfmadestep_results_df[nrow(selfmadestep_results_df) + 1,] <- c('selfmade-step', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fdr, coverage, avg_CI)
    time_results_df[nrow(time_results_df) + 1,] <- c('selfmade-step', p, difftime(Sys.time(),t5,units = "secs"))
    


    
    
    
    
    
    
    
    
    
    
    
    


    print('postcAIC')


    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]


      timeout <<- 0
      t6 <- Sys.time()

      y = ys[,1]


      tryCatch({
        cAIC_model_set <- withTimeout({
          compute_cAIC_for_model_set(
            X,
            y,
            subjind,
            model = "NERM",
            covariate_selection_matrix = NULL,
            modelset  = "all_subsets",
            common = NULL,
            intercept = TRUE
          )

        }, timeout = 10000)
      }, TimeoutException = function(ex) {
        message("Timeout on modelset. Skipping.")
        timeout<<-1
      })


      if(timeout==0){

      cAIC_min = cAIC_model_set$cAIC_min
      degcAIC_models = cAIC_model_set$degcAIC_models
      X_full = cAIC_model_set$X_full
      X_cluster_full = cAIC_model_set$X_cluster_full

      sig_u_full = cAIC_model_set$sig_u_full
      sig_e_full = cAIC_model_set$sig_e_full

      beta_sel = cAIC_model_set$beta_sel
      mu_sel = cAIC_model_set$mu_sel

      modelset_matrix = cAIC_model_set$modelset_matrix
      x_beta_lin_com = cAIC_model_set$X_cluster_full



      p_full = ncol(X_full)
      model='NERM'
      clusterID = subjind
      Z = create_Z(model, clusterID)
      n = ncol(Z)
      C_cluster_full = cbind(X_cluster_full, diag(n))
      R_full = sig_e_full * diag(nrow(X_full))
      invR_full = 1 / sig_e_full * diag(nrow(X_full))
      G_full = sig_u_full * diag(n)
      n_cluster_units = as.data.frame(table(clusterID))$Freq

      V_full_list <- list()
      invV_full_list <- list()

      for (i in 1:n) {
        V_full_list[[i]] <- sig_e_full * diag(n_cluster_units[i]) +
          sig_u_full * matrix(1, nrow = n_cluster_units[i], ncol = n_cluster_units[i])
        invV_full_list[[i]] <- solve(V_full_list[[i]])
      }

      tryCatch({
        postcAIC_CI_results <- withTimeout({
          postcAIC_CI(
            cAIC_min,
            degcAIC_models,

            X_full,
            X_cluster_full,
            sig_u_full,
            sig_e_full,
            model = "NERM",
            subjind,

            beta_sel,
            mu_sel,

            n_samples = 1000,

            modelset_matrix,
            scale_mvrnorm = 10,
            n_starting_points = 5,
            x_beta_lin_com = NULL
          )

        }, timeout = 10000)
      }, TimeoutException = function(ex) {
        message("Timeout. Skipping.")
        results_df[nrow(results_df) + 1,] <<- c('postcAIC', NA, NA, NA, NA, NA)
        timeout<<-1
      })
      }


      if(timeout==0){

      beta_postcAIC_CI_up <- postcAIC_CI_results$beta_postcAIC_CI_up
      beta_postcAIC_CI_do <- postcAIC_CI_results$beta_postcAIC_CI_do
      beta_postcAIC_pval <- postcAIC_CI_results$beta_postcAIC_pval

      selected <- modelset_matrix[cAIC_min,]

      k <- 1
      for(i in seq_len(p)){

        if(selected[i]!=0){

          postcAIC_pvals_df[nrow(postcAIC_pvals_df)+1,] <- c('postcAIC', paste0('X',i), as.numeric(real_beta[i]!=0), real_beta[i], beta_sel[k],
                                           beta_postcAIC_pval[k], beta_postcAIC_CI_do[i], beta_postcAIC_CI_up[i],
                                           as.numeric(beta_postcAIC_CI_do[i]<=real_beta[i] & beta_postcAIC_CI_up[i]>=real_beta[i]))

          k <- k + 1

        }
      }

      postcAIC_pvals_df[,3:9] <- lapply(postcAIC_pvals_df[,3:9],as.numeric)

      coverage <- mean(postcAIC_pvals_df[(nrow(postcAIC_pvals_df)-sum(selected)):nrow(postcAIC_pvals_df), 9])
      avg_CI <- mean(postcAIC_pvals_df[(nrow(postcAIC_pvals_df)-sum(selected)):nrow(postcAIC_pvals_df), 8] - postcAIC_pvals_df[(nrow(postcAIC_pvals_df)-sum(selected)):nrow(postcAIC_pvals_df), 7])

      pvals <- p.adjust(beta_postcAIC_pval, method = 'BH')

      selected[selected==1] <- pvals<=fdr_level
      postcAIC_metrics <- metrics(selected, real_beta!=0)

      postcAIC_results_df[nrow(postcAIC_results_df) + 1,] <- c('postcAIC', postcAIC_metrics$tpr, postcAIC_metrics$fdr, postcAIC_metrics$fdr, coverage, avg_CI)
      time_results_df[nrow(time_results_df) + 1,] <- c('postcAIC', p, difftime(Sys.time(),t6,units = "secs"))
      }else{time_results_df[nrow(time_results_df) + 1,] <- c('postcAIC', p, NA)}


    
  }
  
  saveRDS(time_results_df, paste0("time_results_df_p",p,"_itnum", itnum,".RDS"))
  
  saveRDS(naive_pvals_df, file = paste0("naive_pvals_p",p,"_itnum", itnum,".RDS"))
  saveRDS(naive_results_df, file = paste0("naive_results_p",p,"_itnum", itnum,".RDS"))
  
  saveRDS(datasplitting_pvals_df, file = paste0("datasplitting_pvals_p",p,"_itnum", itnum,".RDS"))
  saveRDS(datasplitting_results_df, file = paste0("datasplitting_results_p",p,"_itnum", itnum,".RDS"))
  
  saveRDS(selfmadelasso_pvals_df, file = paste0("selfmadelasso_pvals_p",p,"_itnum", itnum,".RDS"))
  saveRDS(selfmadelasso_results_df, file = paste0("selfmadelasso_results_p",p,"_itnum", itnum,".RDS"))
  
  saveRDS(selfmadestep_pvals_df, file = paste0("selfmadestep_pvals_p",p,"_itnum", itnum,".RDS"))
  saveRDS(selfmadestep_results_df, file = paste0("selfmadestep_results_p",p,"_itnum", itnum,".RDS"))
  
}

