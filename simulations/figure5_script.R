setwd("C:/Users/matteda/OneDrive - Universitetet i Oslo/Skrivebord/phd/ConditioningApproach/src")

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
library(grid)
library(glmnet)
library(gridExtra)

library(postcAIC)
library(ks)
library(tmg)

# simulation settings

seed = 0
nlambda <- 10

its <- 2 # number of iterations to compute fdr in each simulation
fdr_level <- 0.05
set.seed(seed)

selection_frac <- 0.5 # fraction of observations to use for selecting the variables























#####################################################
################## generating data  #################
#####################################################

n = 200
nrsubj = 40
nrobs = 5
p = 20 # number of fixed effects, other than the intercept
q = 0
p_rel = 5 # number of non-zero coefficients, other than the intercept

SNR = 1 #signal-to-noise ratio

Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables

subjind <- as.factor(rep(1:nrsubj, each=nrobs))

var_names <- c(paste("X", 1:p, sep=""))

yform <- as.formula(
  paste("y ~ - 1 +", paste("X", 1:p, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
)



Xs <- list()
real_betas <- list()
yss <- list()
Zs <- list()


for(i in 1:sim_num){
  
  print(i)
  
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


saveRDS(Xs, paste0('Xs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs, '_SNR',SNR,'_changinglambda.RDS'))
saveRDS(real_betas, paste0('real_betas_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'_changinglambda.RDS'))
saveRDS(Zs, paste0('Zs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'_changinglambda.RDS'))
saveRDS(yss, paste0('yss_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'_changinglambda.RDS'))




















lambdaMax <- 600

lambda_vec <- exp(seq(from = log(lambdaMax),
                      to = log(lambdaMax * 0.0001),
                      length.out = 10))



#####################################################
###########       applying methods        ###########
#####################################################


results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(results_df) <- c('method', 'tpr', 'fdr', 'coverage', 'avg_ci', 'lambda', 'lambda_cv','nonzero','cvnonzero')

pvals_df <- data.frame(matrix(nrow = 0, ncol = 10))
colnames(pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered', 'lambda')


print('Naive inference')


for(jj in itnum*5+c(-4:0)){
  
  X <- Xs[[jj]]
  real_beta <- real_betas[[jj]]
  Z <- Zs[[jj]]
  ys <- yss[[jj]]
  print(jj)
  
  for(j in c(1)){
    
    
    y <- ys[,j]
    
    dat <- data.frame(X, Z, y, subjind)
    
    sx <- as.matrix(scale(X))
    sy <- as.vector(y)
    lambdaMax <- max(abs(colSums(sx*sy)))
    
    
    PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose=FALSE)
    Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
    Q.start<-as.numeric(VarCorr(PQL)[1,1])
    
    BIC_vec <- rep(Inf, length(lambda_vec))
    sel_vec <- matrix(NA,p,length(lambda_vec))
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = dat,
        lambda = lambda_vec[l],
        control = list(center = FALSE, start=Delta.start[l,], q_start = Q.start[l])
      ), silent=TRUE)
      
      if (!inherits(glm3, "try-error")) {
        sel_vec[,l] <- as.numeric(glm3$coefficients!=0)
        BIC_vec[l] <- glm3$bic
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      }
    }
    
    opt3 <- which.min(BIC_vec)
    lambda_cv <- lambda_vec[opt3]
    cvnonzero <- sum(sel_vec[,opt3])
    
    for (l in 1:length(lambda_vec)) {
      
      if(any(is.na(sel_vec[,l]))) next
      
      selected <- sel_vec[,l]
      nonzero=sum(selected)
      
      selected_tot_names <- var_names[selected==1]
      
      selected_tot_yform <- as.formula(
        paste(" y ~ -1 + ", paste(selected_tot_names, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
      )
      
      if(q>0){
        selected_tot_yform <- as.formula(
          paste(" y ~ -1 +", paste(selected_tot_names, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
        )
      }

      
      suppressWarnings(suppressMessages(sel_mod <- lmer(formula = selected_tot_yform, data = dat)))
      
      
      if(sum(selected)==1){
        
        pvals_tot <- as.data.frame(t(coef(summary(sel_mod))[,c(1,5)]))
        
      }else{
        
        pvals_tot <- as.data.frame(coef(summary(sel_mod))[,c(1,5)])
        
        
      }
      
      suppressMessages(confint_tot <- confint(sel_mod, method='Wald', parm = 'beta_'))
      confint_tot <- as.data.frame(confint_tot)
      
      k <- 1
      for(i in 1:(p)){
        
        if(selected[i] == 1){
          pvals_df[nrow(pvals_df) + 1,] <- c('Naive', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], pvals_tot[k,1], pvals_tot[k,2],
                                             confint_tot[k,1], confint_tot[k,2], as.numeric(confint_tot[k,1]<=real_beta[i] & confint_tot[k,2]>=real_beta[i]), lambda_vec[l])
          k = k + 1
        }
      }
      
      pvals_df[,3:10] <- lapply(pvals_df[,3:10],as.numeric)
      
      coverage <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 9])
      avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 8] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7])
      
      
      corrected_pvals_tot <- p.adjust(pvals_tot[,2], method = 'BH')
      
      selected[selected==1] <- corrected_pvals_tot<=fdr_level
      tot_metrics <- metrics(selected, real_beta!=0)
      
      results_df[nrow(results_df) + 1,] <- c('Naive', tot_metrics$tpr, tot_metrics$fdr, coverage, avg_CI, lambda_vec[l], lambda_cv, nonzero, cvnonzero)
    }
  }
}

results_df[,2:9] <- lapply(results_df[,2:9],as.numeric)

saveRDS(pvals_df, file = paste0('naive_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))
saveRDS(results_df, file = paste0('naive_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))



















results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(results_df) <- c('method', 'tpr', 'fdr', 'coverage', 'avg_ci', 'lambda', 'lambda_cv','nonzero','cvnonzero')

pvals_df <- data.frame(matrix(nrow = 0, ncol = 10))
colnames(pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered', 'lambda')


print('Data Splitting')

for(jj in itnum*5+c(-4:0)){

  X <- Xs[[jj]]
  real_beta <- real_betas[[jj]]
  Z <- Zs[[jj]]
  ys <- yss[[jj]]
  print(jj)

  for(j in c(1)){


    y <- ys[,j]

    dat <- data.frame(X, Z, y, subjind)

    # data splitting: the data is split in two, using the first half for selection and the second for inference

    selection_sample <- sample(1:nrsubj, round(nrsubj*selection_frac))
    selection_dat <- dat[dat$subjind %in% selection_sample,]
    inference_dat <- dat[!(dat$subjind %in% selection_sample),]

    X_sel <- as.matrix(selection_dat[,1:p])
    y <- selection_dat$y


    selection_dat$subjind <- droplevels(selection_dat$subjind)

    PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=selection_dat,verbose=FALSE)
    Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
    Q.start<-as.numeric(VarCorr(PQL)[1,1])

    BIC_vec <- rep(Inf, length(lambda_vec))
    sel_vec <- matrix(NA,p,length(lambda_vec))

    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = selection_dat,
        lambda = lambda_vec[l],
        control = list(center = FALSE, start=Delta.start[l,], q_start = Q.start[l])
      ), silent=TRUE)

      if (!inherits(glm3, "try-error")) {
        sel_vec[,l] <- as.numeric(glm3$coefficients!=0)
        BIC_vec[l] <- glm3$bic
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      }
    }

    opt3 <- which.min(BIC_vec)
    lambda_cv <- lambda_vec[opt3]
    cvnonzero <- sum(sel_vec[,opt3])

    for (l in 1:length(lambda_vec)) {

      if(any(is.na(sel_vec[,l]))) next

      selected <- sel_vec[,l]
      nonzero=sum(selected)


      selected_split_names <- var_names[selected==1]

      selected_split_yform <- as.formula(
        paste(" y ~ -1 + ", paste(selected_split_names, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
      )

      if(q>0){
        selected_split_yform <- as.formula(
          paste(" y ~ -1 +", paste(selected_split_names, sep = "", collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)", sep = "")
        )
      }

      suppressWarnings(suppressMessages(sel_mod <- lmer(formula = selected_split_yform, data = inference_dat)))

      if(sum(selected)==1){

        pvals_split <- as.data.frame(t(coef(summary(sel_mod))[,c(1,5)]))

      }else{

        pvals_split  <- as.data.frame(coef(summary(sel_mod))[,c(1,5)])

      }

      suppressMessages(confint_split <- confint(sel_mod, method='Wald', parm='beta_'))
      confint_split  <- as.data.frame(confint_split )

      k <- 1
      for(i in 1:(p)){

        if(selected[i] == 1){
          pvals_df[nrow(pvals_df) + 1,] <- c('Data Splitting', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], pvals_split[k,1], pvals_split[k,2],
                                             confint_split[k,1], confint_split[k,2], as.numeric(confint_split[k,1]<=real_beta[i] & confint_split[k,2]>=real_beta[i]), lambda_vec[l])
          k = k + 1
        }
      }

      pvals_df[,3:10] <- lapply(pvals_df[,3:10],as.numeric)



      coverage <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 9])
      avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 8] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7])

      corrected_pvals_split <- p.adjust(pvals_split[,2], method = 'BH')

      selected[selected==1] <- corrected_pvals_split<=fdr_level
      split_metrics <- metrics(selected, real_beta!=0)

      results_df[nrow(results_df) + 1,] <- c('Data splitting', split_metrics$tpr, split_metrics$fdr, coverage, avg_CI, lambda_vec[l], lambda_cv, nonzero, cvnonzero)

    }
  }
}

results_df[,2:9] <- lapply(results_df[,2:9],as.numeric)

saveRDS(pvals_df, file = paste0('datasplitting_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))
saveRDS(results_df, file = paste0('datasplitting_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))





















results_df <- data.frame(matrix(nrow = 0, ncol = 11))
colnames(results_df) <- c('method', 'tpr', 'fdr', 'coverage', 'avg_ci', 'lambda', 'lambda_cv','nonzero','cvnonzero','survivors_signal','survivors_noise')

pvals_df <- data.frame(matrix(nrow = 0, ncol = 12))
colnames(pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered', 'lambda','survivors_signal','survivors_noise')

nsteps <- 5

print('selfmade-lasso')

for (index in 1:length(lambda_vec)) {

  print(index)


  lambda <- lambda_vec[index]

  modFun <- function(yy)
  {
    dat$y <- as.numeric(yy)
    lmer(yform, REML = FALSE, data = dat)
  }

  selFun <- function(yy){

    lambda_vec_small <- exp(seq(from = log(lambdaMax),
                                to = log(lambda),
                                length.out = nsteps))

    dat$y <- as.numeric(yy)

    PQL<-try(glmmPQL(y~-1,random = ~1|subjind, family='gaussian',data=dat,verbose='FALSE'), silent = TRUE)

    if(inherits(PQL,"try-error")){
      Delta.start<-NULL
      Q.start<-NULL
    }else{
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])
    }


    for (l in 1:nsteps) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = dat,
        lambda = lambda_vec_small[l],
        control = list(center = FALSE, start = Delta.start[nrow(Delta.start), ], q_start = Q.start[nrow(Delta.start)])
      ), silent = TRUE)

      if (!inherits(glm3, "try-error")) {

        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      }
    }

    if(inherits(glm3, "try-error")){
      names_vec <- NA}else{
        names_vec = paste("X", 1:p, sep = "")[glm3$coefficients != 0]
      }



    return(c(names_vec))
  }


  for(jj in itnum*5+c(-4:0)){

    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)

    for(j in c(1)){

      #print(j)

      y = ys[,j]

      dat <- data.frame(X, Z, y, subjind)

      y = ys[,j]

      dat <- data.frame(X, Z, y, subjind)
      X <- as.matrix(dat[,1:p])

      names_vec = selFun(y)

      nonzero = length(names_vec)

      if(nonzero==0){
        results_df[nrow(results_df) + 1,] <- c('selfmade-lasso', 0, 0, NA, NA, lambda, NA, 0, NA, NA, NA)
        next
      }

      sel_form <- as.formula(
        paste("y ~ - 1 + ", paste(names_vec[1:length(names_vec)], collapse = "+"), "+ (1 |subjind)")
      )
      if(q>0){
        sel_form <- as.formula(
          paste("y ~ - 1 + ", paste(names_vec[1:length(names_vec)], collapse = "+"), "+ (1 + ",paste0("Z",1:q,collapse = '+'),"|subjind)")
        )
      }

      final_model <- suppressWarnings(suppressMessages(lmer(formula = sel_form, data = dat)))

      beta <- fixef(final_model)
      selection = c(names(fixef(final_model)))

      # define function which checks congruency
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
                                    trace = FALSE))

      r <- do.call("rbind", r$selinf)
      r$variable <- selection

      survivors_signal <- mean(r$nrsurv[r$variable %in% signal_var_names])
      survivors_noise <- mean(r$nrsurv[!(r$variable %in% signal_var_names)])

      selected <- rep(0,p)

      k <- 1
      for(i in seq_len(p)){

        if(var_names[i] %in% selection){

          pvals_df[nrow(pvals_df)+1,] <- c('selfmade-lasso', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], beta[k],
                                           r$pval[k], r$cil[k], r$ciu[k], as.numeric(r$cil[k]<=real_beta[i] & r$ciu[k]>=real_beta[i]), lambda, survivors_signal, survivors_noise)

          selected[i] <- 1
          k <- k +1

        }
      }

      pvals_df[,3:10] <- lapply(pvals_df[,3:10],as.numeric)
      pvals_df[is.infinite(pvals_df[,7]),7] <- NA
      pvals_df[is.infinite(pvals_df[,8]),8] <- NA

      coverage <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 9], na.rm=TRUE)
      avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 8] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7], na.rm=TRUE)

      pvals <- p.adjust(r$pval, method = 'BH')

      selected[selected==1] <- pvals<=fdr_level
      selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
      rugamer_metrics <- metrics(selected, real_beta!=0)

      results_df[nrow(results_df) + 1,] <- c('selfmade-lasso', rugamer_metrics$tpr, rugamer_metrics$fdr, coverage, avg_CI, lambda, NA, nonzero, NA, survivors_signal, survivors_noise)

    }

  }

  results_df[,2:11] <- lapply(results_df[,2:11],as.numeric)

  saveRDS(pvals_df, file = paste0('selfmadelasso_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))
  saveRDS(results_df, file = paste0('selfmadelasso_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_SNR',SNR,'.RDS'))

}
