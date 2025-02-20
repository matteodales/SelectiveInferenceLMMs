
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

plot_results = TRUE

seed = 0

sim_num <- 100 # number of simulations
its <- 1 # number of iterations to compute FWER in each simulation
tot_its <- its*sim_num

fdr_level <- 0.05
set.seed(seed)

selection_frac <- 0.5 # fraction of observations to use for selecting the variables

























#####################################################
################## generating data  #################
#####################################################

list_settings <- list(list(n = 200, nrsubj = 20, nrobs = 10, p = 25, p_rel=10),
                      list(n = 200, nrsubj = 20, nrobs = 10, p = 50, p_rel=10),
                      list(n = 200, nrsubj = 20, nrobs = 10, p = 150, p_rel=10),
                      list(n = 200, nrsubj = 20, nrobs = 10, p = 250, p_rel=10))


for(setting in list_settings){
  
  n <- setting$n
  nrsubj <- setting$nrsubj
  nrobs <- setting$nrobs
  p <- setting$p
  p_rel <- setting$p_rel
  
  
  SNR = 2 #signal-to-noise ratio
  
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
  
  
  saveRDS(Xs, paste0('Xs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs, '_SNR',SNR,'.RDS'))
  saveRDS(real_betas, paste0('real_betas_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(Zs, paste0('Zs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(yss, paste0('yss_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
}




















#####################################################
###########       applying methods        ###########
#####################################################



for(setting in list_settings){
  
  n <- setting$n
  nrsubj <- setting$nrsubj
  nrobs <- setting$nrobs
  p <- setting$p
  p_rel <- setting$p_rel
  ratio <- if(p<=n) 0.0001 else 0.01
  print(paste0("p ", p))
  
  
  SNR = 2 #signal-to-noise ratio
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
  var_names <- c(paste("X", 1:p, sep=""))
  
  
  Xs <- readRDS(paste0('Xs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs, '_SNR',SNR,'.RDS'))
  real_betas <- readRDS(paste0('real_betas_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  Zs <- readRDS(paste0('Zs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  yss <- readRDS(paste0('yss_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
  
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
  


  
  
  
  
  
  


  naive_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(naive_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')

  naive_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(naive_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')


  print('Naive inference')

  for(jj in itnum*5+c(-4:0)){

    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)

    for(j in 1:its){


      y <- ys[,j]

      dat <- data.frame(X, Z, y, subjind)

      sx <- as.matrix(scale(X))
      sy <- as.vector(y)

      lambdaMax <- max(abs(colSums(sx*sy)))

      lambda_vec <- exp(seq(from = log(lambdaMax),
                            to = log(lambdaMax * ratio),
                            length.out = 50))

      PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose=FALSE)
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])

      BIC_vec <- rep(Inf, length(lambda_vec))
      sel_vec <- matrix(0,length(lambda_vec),p)

      # Early stopping parameters
      patience_limit <- 3  # Max consecutive iterations without improvement
      patience_count <- 0  # Counter for consecutive non-improvements
      kk <- 0
      best_bic <- Inf      # Keep track of the best BIC value
      current_bic <- Inf

      for (l in 1:length(lambda_vec)) {
        glm3 <- try(glmmLasso(
          fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
          rnd = list(subjind = ~1),
          data = dat,
          lambda = lambda_vec[l],
          control = list(center = FALSE, start = Delta.start[nrow(Delta.start), ], q_start = Q.start[nrow(Delta.start)])
        ), silent = TRUE)

        if (!inherits(glm3, "try-error")) {
          current_bic <- glm3$bic
          BIC_vec[l] <- current_bic
          sel_vec[l, ] <- as.numeric(glm3$coefficients != 0)
          Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
          Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])

          # Check for improvement
          if (!is.numeric(current_bic) || is.nan(current_bic)) next
          if (current_bic < 1.1*best_bic) {
            if (current_bic < best_bic) best_bic <- current_bic
            patience_count <- 0  # Reset patience
          } else {
            patience_count <- patience_count + 1  # Increment patience
          }

          #Early stopping condition
          if (patience_count >= patience_limit) break

        }
      }


      opt <- which.min(BIC_vec)

      selected <- sel_vec[opt,]

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

      naive_results_df[nrow(naive_results_df) + 1,] <- c('Naive', tot_metrics$tpr, tot_metrics$fdr, tot_metrics$fwer, coverage, avg_CI)

    }
  }

  naive_results_df[,2:6] <- lapply(naive_results_df[,2:6],as.numeric)

  saveRDS(naive_pvals_df, file = paste0('naive_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(naive_results_df, file = paste0('naive_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



  datasplitting_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(datasplitting_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')

  datasplitting_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(datasplitting_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')


  print('Data Splitting')

  for(jj in itnum*5+c(-4:0)){

    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)

    for(j in 1:its){

      y <- ys[,j]

      dat <- data.frame(X, Z, y, subjind)

      #data splitting: the data is split in two, using the first half for selection and the second for inference

      selection_sample <- sample(1:nrsubj, round(nrsubj*selection_frac))
      selection_dat <- dat[dat$subjind %in% selection_sample,]
      inference_dat <- dat[!(dat$subjind %in% selection_sample),]

      X_sel <- as.matrix(selection_dat[,1:p])
      y <- selection_dat$y

      sx <- as.matrix(scale(X))
      sy <- as.vector(y)

      lambdaMax <- max(abs(colSums(sx*sy)))


      lambda_vec <- exp(seq(from = log(lambdaMax),
                            to = log(lambdaMax * ratio),
                            length.out = 50))

      selection_dat$subjind <- droplevels(selection_dat$subjind)

      PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=selection_dat,verbose=FALSE)
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])

      BIC_vec <- rep(Inf, length(lambda_vec))
      sel_vec <- matrix(0,length(lambda_vec),p)

      # Early stopping parameters
      patience_limit <- 3  # Max consecutive iterations without improvement
      patience_count <- 0  # Counter for consecutive non-improvements
      best_bic <- Inf      # Keep track of the best BIC value
      current_bic <- Inf

      for (l in 1:length(lambda_vec)) {
        glm3 <- try(glmmLasso(
          fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
          rnd = list(subjind = ~1),
          data = selection_dat,
          lambda = lambda_vec[l],
          control = list(center = FALSE, start = Delta.start[nrow(Delta.start), ], q_start = Q.start[nrow(Delta.start)])
        ), silent = TRUE)

        if (!inherits(glm3, "try-error")) {
          current_bic <- glm3$bic
          BIC_vec[l] <- current_bic
          sel_vec[l, ] <- as.numeric(glm3$coefficients != 0)
          Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
          Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])

          # Check for improvement
          if (!is.numeric(current_bic) || is.nan(current_bic)) next
          if (current_bic < 1.1*best_bic) {
            if (current_bic < best_bic) best_bic <- current_bic
            patience_count <- 0  # Reset patience
          } else {
            patience_count <- patience_count + 1  # Increment patience
          }

          #Early stopping condition
          if (patience_count >= patience_limit) break

        }
      }

      opt <- which.min(BIC_vec)

      selected <- sel_vec[opt,]

      if(sum(selected)>100) next

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

      datasplitting_results_df[nrow(datasplitting_results_df) + 1,] <- c('Data splitting', split_metrics$tpr, split_metrics$fdr, split_metrics$fwer, coverage, avg_CI)

    }
  }

  datasplitting_results_df[,2:6] <- lapply(datasplitting_results_df[,2:6],as.numeric)

  saveRDS(datasplitting_pvals_df, file = paste0('datasplitting_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(datasplitting_results_df, file = paste0('datasplitting_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #dataframes to save results
  
  selfmadelasso_results_df <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(selfmadelasso_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci','surv')
  
  selfmadelasso_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(selfmadelasso_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
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
    
    PQL<-try(glmmPQL(y~-1,random = ~1|subjind, family='gaussian',data=dat,verbose='FALSE'), silent = TRUE)
    
    if(inherits(PQL,"try-error")){
      Delta.start<-NULL
      Q.start<-NULL
    }else{
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])
    }
    
    BIC_vec <- c(rep(Inf,length(lambda_vec)))
    sel_vec <- matrix(0,length(lambda_vec),p)
    
    # Early stopping parameters
    patience_limit <- 3  # Max consecutive iterations without improvement
    patience_count <- 0  # Counter for consecutive non-improvements
    best_bic <- Inf      # Keep track of the best BIC value
    current_bic <- Inf
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ -1 +", paste("X", 1:p, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = dat,
        lambda = lambda_vec[l],
        control = list(center = FALSE, start = Delta.start[nrow(Delta.start), ], q_start = Q.start[nrow(Delta.start)])
      ), silent = TRUE)
      
      if (!inherits(glm3, "try-error")) {
        current_bic <- glm3$bic
        BIC_vec[l] <- current_bic
        sel_vec[l, ] <- as.numeric(glm3$coefficients != 0)
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
        
        # Check for improvement
        if (!is.numeric(current_bic) || is.nan(current_bic)) next
        if (current_bic < 1.1*best_bic) {
          if (current_bic < best_bic) best_bic <- current_bic
          patience_count <- 0  # Reset patience
        } else {
          patience_count <- patience_count + 1  # Increment patience
        }
        
        #Early stopping condition
        if (patience_count >= patience_limit) break
      }
    }
    
    opt <- which.min(BIC_vec)
    
    names_vec = paste("X", 1:p, sep = "")[sel_vec[opt,]==1]
    
    return(c(names_vec))
  }
  
  
  print('selfmade-lasso')
  
  for(jj in itnum*5+c(-4:0)){
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)
    
    
    for(j in c(1)){
      
      
      print(j)
      
      y = ys[,1]
      
      dat <- data.frame(X, Z, y, subjind)
      
      
      sx <- as.matrix(scale(X))
      sy <- as.vector(y)
      
      lambdaMax <- max(abs(colSums(sx*sy)))
      lambda_vec <<- exp(seq(from = log(lambdaMax),
                             to = log(lambdaMax * ratio),
                             length.out = 50))
      
      
      names_vec = selFun(y)
      
      if(length(names_vec)==0){
        selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', 0, 0, 0, NA, NA)
        next
      }
      
      
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
      surv <- mean(r$nrsurv, na.rm=TRUE)
      
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
      
      selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fwer, coverage, avg_CI, surv)
      
    }
  }
  
  selfmadelasso_results_df[,2:6] <- lapply(selfmadelasso_results_df[,2:6],as.numeric)
  
  saveRDS(selfmadelasso_pvals_df, file = paste0('selfmadelasso_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(selfmadelasso_results_df, file = paste0('selfmadelasso_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  selfmadestep_results_df <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(selfmadestep_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci', 'surv')
  
  selfmadestep_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(selfmadestep_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  
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
  for(jj in itnum*5+c(-4:0)){
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)
    
    
    for(j in c(1)){
      
      
      print(j)
      
      y = ys[,1]
      
      dat <- data.frame(X, Z, y, subjind)
      
      
      mod <- modFun(yy=y, dat=dat)
      
      final_model <- selFun(mod)
      beta <- fixef(final_model)
      print(final_model)
      
      selection = extractSelFun(selFun(mod))
      
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
                                    trace=TRUE))
      
      
      r <- do.call("rbind", r$selinf)
      r$variable <- selection
      surv <- mean(r$nrsurv, na.rm=TRUE)
      
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
      
      selfmadestep_results_df[nrow(selfmadestep_results_df) + 1,] <- c('selfmade-step', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fwer, coverage, avg_CI, surv)
    }
  }
  
  selfmadestep_results_df[,2:6] <- lapply(selfmadestep_results_df[,2:6],as.numeric)
  
  saveRDS(selfmadestep_pvals_df, file = paste0('selfmadestep_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(selfmadestep_results_df, file = paste0('selfmadestep_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
  
}
