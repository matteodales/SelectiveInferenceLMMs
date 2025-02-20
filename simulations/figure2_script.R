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

seed = 0

sim_num <- 25 # number of simulations
its <- 20 # number of iterations to compute FWER in each simulation
tot_its <- its*sim_num

fdr_level <- 0.05
set.seed(seed)

selection_frac <- 0.5 # fraction of observations to use for selecting the variables

myCluster <- makeCluster(25, outfile="" )
registerDoParallel(myCluster)


















#####################################################
################## generating data  #################
#####################################################

list_settings <- list(list(n = 100, nrsubj = 20, nrobs = 5),
                      list(n = 200, nrsubj = 40, nrobs = 5),
                      list(n = 250, nrsubj = 50, nrobs = 5),
                      list(n = 100, nrsubj = 10, nrobs = 10),
                      list(n = 200, nrsubj = 10, nrobs = 20),
                      list(n = 250, nrsubj = 10, nrobs = 25))



for(setting in list_settings){
  
  n <- setting$n
  nrsubj <- setting$nrsubj
  nrobs <- setting$nrobs
  
  # n = 250
  # nrsubj = 10
  # nrobs = 25
  p = 6 # number of fixed effects, other than the intercept
  p_rel = 3 # number of non-zero coefficients, other than the intercept
  q = 0 # number of mixed effects, other than the intercept
  
  
  SNR = 1 #signal-to-noise ratio
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
  var_names <- c(paste("X", 1:p, sep=""))
  
  
  Xs <- list()
  real_betas <- list()
  yss <- list()
  Zs <- list()
  
  real_beta <- c(1,-2,-1,0,0,0)
  
  
  for(i in 1:sim_num){
    
    #real_beta = c(sample (c(1,2,-1,-2), size=p_rel, replace=T), rep(0, p-p_rel))
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
########### computing modelsets postcAIC  ###########
#####################################################




for(setting in list_settings){
  
  n <- setting$n
  nrsubj <- setting$nrsubj
  nrobs <- setting$nrobs
  
  
  SNR = 1 #signal-to-noise ratio
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  
  
  Xs <- readRDS(paste0('Xs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs, '_SNR',SNR,'.RDS'))
  real_betas <- readRDS(paste0('real_betas_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  Zs <- readRDS(paste0('Zs_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  yss <- readRDS(paste0('yss_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
  
  cAIC_model_setss <- list()
  
  print('postcAIC')
  
  for(jj in 1:sim_num){
    
    print(jj)
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    
    cAIC_model_setss[[jj]] <- list()
    
    for(j in 1:its){
      
      print(j)
      
      y = ys[,j]
      
      timeout <<- 0
      
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
          
        }, timeout = 600)
      }, TimeoutException = function(ex) {
        message("Timeout on modelset. Skipping.")
        timeout<<-1
      })
      
      if(timeout==1){
        cAIC_model_setss[[jj]][[j]] <- NA
        next}
      
      cAIC_model_setss[[jj]][[j]] <- cAIC_model_set
      
    }
    
  }
  
  
  saveRDS(cAIC_model_setss, paste0('cAIC_model_setss_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
}













#####################################################
###########       applying methods        ###########
#####################################################


for(setting in list_settings){
  
  n <- setting$n
  nrsubj <- setting$nrsubj
  nrobs <- setting$nrobs
  
  # n = 250
  # nrsubj = 10
  # nrobs = 25
  p = 6 # number of fixed effects, other than the intercept
  p_rel = 3 # number of non-zero coefficients, other than the intercept
  q = 0 # number of mixed effects, other than the intercept
  
  
  SNR = 1 #signal-to-noise ratio
  
  Sigmap = matrix(0.3, p, p) + 0.7*diag(p) # covariance matrix for the covariates
  Sigmaq = matrix(0.3, q, q) + 0.7*diag(q) # covariance matrix for the random variables
  
  subjind <- as.factor(rep(1:nrsubj, each=nrobs))
  varcovar = matrix(0.3, q+1, q+1) + 0.7*diag(q+1)
  var_names <- c(paste("X", 1:p, sep=""))
  
  lambda_cv_list <- list()
  selected_cv_list <- list()
  glmmlasso_cvmod_list <- list()

  
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
  
  for(jj in 1:sim_num){
    
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
                            to = log(lambdaMax * 0.001),
                            length.out = 50))
      
      PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose=FALSE)
      Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
      Q.start<-as.numeric(VarCorr(PQL)[1,1])
      
      BIC_vec <- rep(Inf, length(lambda_vec))
      
      for (l in 1:length(lambda_vec)) {
        glm3 <- try(glmmLasso(
          fix = y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6,
          rnd = list(subjind = ~1),
          data = dat,
          lambda = lambda_vec[l],
          control = list(center = FALSE, start=Delta.start[l,], q_start = Q.start[l])
        ), silent=TRUE)
        
        if (!inherits(glm3, "try-error")) {
          BIC_vec[l] <- glm3$bic
          Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
          Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
        }
      }
      
      opt3 <- which.min(BIC_vec)
      
      glm3_final <- try(glmmLasso(fix = y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6, rnd = list(subjind=~1), data = dat, lambda=lambda_vec[opt3], switch.NR=FALSE,final.re=FALSE,
                                  control=list(center=FALSE, start=Delta.start[opt3,], q_start = Q.start[opt3])), silent=TRUE)
      
      
      selected <- as.numeric(glm3_final$coefficients!=0)

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
      
      
      corrected_pvals_tot <- p.adjust(pvals_tot[,2], method = 'holm')
      
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
  
  for(jj in 1:sim_num){
    
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
      
      for (l in 1:length(lambda_vec)) {
        glm3 <- try(glmmLasso(
          fix = y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6,
          rnd = list(subjind = ~1),
          data = selection_dat,
          lambda = lambda_vec[l],
          control = list(start = Delta.start[l, ], q_start = Q.start[l], center = FALSE)
        ), silent=TRUE)
        
        if (!inherits(glm3, "try-error")) {
          BIC_vec[l] <- glm3$bic
          Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
          Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
        }
      }
      
      opt3 <- which.min(BIC_vec)
      
      glm3_final <- try(glmmLasso(fix = y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6, rnd = list(subjind=~1), data = selection_dat, lambda=lambda_vec[opt3], switch.NR=FALSE,final.re=FALSE,
                                  control=list(start=Delta.start[opt3,],q_start=Q.start[opt3], center=FALSE)), silent=TRUE)
      
    
      
      selected <- as.numeric(glm3_final$coefficients!=0)
      
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
      
      corrected_pvals_split <- p.adjust(pvals_split[,2], method = 'holm')
      
      selected[selected==1] <- corrected_pvals_split<=fdr_level
      split_metrics <- metrics(selected, real_beta!=0)
      
      datasplitting_results_df[nrow(datasplitting_results_df) + 1,] <- c('Data splitting', split_metrics$tpr, split_metrics$fdr, split_metrics$fwer, coverage, avg_CI)
      
    }
  }
  
  datasplitting_results_df[,2:6] <- lapply(datasplitting_results_df[,2:6],as.numeric)
  
  saveRDS(datasplitting_pvals_df, file = paste0('datasplitting_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(datasplitting_results_df, file = paste0('datasplitting_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  

  
  
  
  
  
  
  
  
  





  print('postcAIC')

  output <- foreach(jj=1:sim_num, .packages = c("glmnet","lme4",'lmerTest','nlme','R.utils','postcAIC','ks','tmg','mgcv')) %dopar% {


      results_df <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')

      pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
      colnames(pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')

      print(jj)

      X <- Xs[[jj]]
      real_beta <- real_betas[[jj]]
      Z <- Zs[[jj]]
      ys <- yss[[jj]]

      for(j in 1:its){

          timeout <<- 0

          y = ys[,j]

          print(j)

          if(any(is.na(cAIC_model_setss[[jj]][[j]]))) next

          cAIC_model_set <- cAIC_model_setss[[jj]][[j]]

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

          t1 <- Sys.time()


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

              }, timeout = 600)
              }, TimeoutException = function(ex) {
              message("Timeout. Skipping.")
              results_df[nrow(results_df) + 1,] <<- c('postcAIC', NA, NA, NA, NA, NA)
              timeout<<-1
              })

          time_elapsed <- difftime(Sys.time(),t1,units='secs')



          if(timeout==1) next

          beta_postcAIC_CI_up <- postcAIC_CI_results$beta_postcAIC_CI_up
          beta_postcAIC_CI_do <- postcAIC_CI_results$beta_postcAIC_CI_do
          beta_postcAIC_pval <- postcAIC_CI_results$beta_postcAIC_pval

          selected <- modelset_matrix[cAIC_min,]

          k <- 1
          for(i in seq_len(p)){

              if(selected[i]!=0){

                  pvals_df[nrow(pvals_df)+1,] <- c('postcAIC', paste0('X',i), as.numeric(real_beta[i]!=0), real_beta[i], beta_sel[k],
                  beta_postcAIC_pval[k], beta_postcAIC_CI_do[i], beta_postcAIC_CI_up[i],
                  as.numeric(beta_postcAIC_CI_do[i]<=real_beta[i] & beta_postcAIC_CI_up[i]>=real_beta[i]))

                  k <- k + 1

              }
          }

          pvals_df[,3:9] <- lapply(pvals_df[,3:9],as.numeric)

          coverage <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 9])
          avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 8] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7])

          pvals <- p.adjust(beta_postcAIC_pval, method = 'holm')

          selected[selected==1] <- pvals<=fdr_level
          postcAIC_metrics <- metrics(selected, real_beta!=0)

          results_df[nrow(results_df) + 1,] <- c('postcAIC', postcAIC_metrics$tpr, postcAIC_metrics$fdr, postcAIC_metrics$fwer, coverage, avg_CI)

      }

      print(paste0('Done ',jj))

      results_df[,2:6] <- lapply(results_df[,2:6],as.numeric)

      return(list(results = results_df, pvals = pvals_df))

  }

  results_df <- do.call('rbind', lapply(output,function(x){x[[1]]}))
  pvals_df <- do.call('rbind', lapply(output,function(x){x[[2]]}))


  results_df[,2:6] <- lapply(results_df[,2:6],as.numeric)

  saveRDS(pvals_df, file = paste0('postcAIC_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR', SNR,'.RDS'))
  saveRDS(results_df, file = paste0('postcAIC_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR', SNR,'.RDS'))

  
  
  

  
  selfmadelasso_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(selfmadelasso_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')
  
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
    
    
    
    PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose='FALSE')
    Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
    Q.start<-as.numeric(VarCorr(PQL)[1,1])
    BIC_vec <- c(rep(Inf,length(lambda_vec)))
    sel_vec <- matrix(0,length(lambda_vec),p)
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = y ~ -1 + X1 + X2 + X3 + X4 + X5 + X6,
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
    
    names_vec = names(glm3$coefficients)[sel_vec[opt,]==1]
    
    return(c(names_vec))
  }

  
  
  print('selfmade-lasso')
  
  for(jj in 1:sim_num){
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    print(jj)
    
    for(j in 1:its){
      
      print(j)
      
      y = ys[,j]
      
      dat <- data.frame(X, Z, y, subjind)
      
      
      sx <- as.matrix(scale(X))
      sy <- as.vector(y)
      
      lambdaMax <- max(abs(colSums(sx*sy)))
      lambda_vec <<- exp(seq(from = log(lambdaMax),
                             to = log(lambdaMax * 0.001),
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
      
      pvals <- p.adjust(r$pval, method = 'holm')
      
      selected[selected==1] <- pvals<=fdr_level
      selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
      rugamer_metrics <- metrics(selected, real_beta!=0)
      
      selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fwer, coverage, avg_CI)
      
    }
    
  }
  
  selfmadelasso_results_df[,2:6] <- lapply(selfmadelasso_results_df[,2:6],as.numeric)
  
  saveRDS(selfmadelasso_pvals_df, file = paste0('selfmadelasso_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(selfmadelasso_results_df, file = paste0('selfmadelasso_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
  
  
  
  print('selfmade-step')
  
  output <- foreach(jj=1:sim_num, .packages = c("glmnet","lme4",'lmerTest')) %dopar% {
    
    modFun <- function(yy, dat = dat) 
    {
      dat$y <- as.numeric(yy)
      suppressWarnings(suppressMessages(lmer(y ~ - 1 + X1 + X2 + X3 + X4 + X5 + X6 + (1 |subjind), REML = FALSE, data = dat)))
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
    
    
    
    results_df <- data.frame(matrix(nrow = 0, ncol = 6))
    colnames(results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')
    
    pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
    colnames(pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
    
    print(jj)
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    
    for(j in 1:its){
      
      print(j)
      
      y = ys[,j]
      
      dat <<- data.frame(X, Z, y, subjind)
      
      mod <- modFun(yy=y, dat=dat)
  
      
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
          
          pvals_df[nrow(pvals_df)+1,] <- c('selfmade-step', var_names[i], as.numeric(real_beta[i]!=0), real_beta[i], beta[k], 
                                           r$pval[k], r$cil[k], r$ciu[k], as.numeric(r$cil[k]<=real_beta[i] & r$ciu[k]>=real_beta[i]))
          
          selected[i] <- 1
          k <- k +1
          
        }
      }
      
      pvals_df[,3:9] <- lapply(pvals_df[,3:9],as.numeric)
      pvals_df[is.infinite(pvals_df[,7]),7] <- NA
      pvals_df[is.infinite(pvals_df[,8]),8] <- NA
      
      coverage <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 9], na.rm=TRUE)
      avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 8] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7], na.rm=TRUE)
      
      pvals <- p.adjust(r$pval, method = 'holm')
      
      selected[selected==1] <- pvals<=fdr_level
      selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
      rugamer_metrics <- metrics(selected, real_beta!=0)
      
      results_df[nrow(results_df) + 1,] <- c('selfmade-step', rugamer_metrics$tpr, rugamer_metrics$fdr, rugamer_metrics$fwer, coverage, avg_CI)
      
    }
    
    print(paste0('Done ',jj))
    return(list(pvals_df, results_df))
  }
  
  results_df <- do.call('rbind', lapply(output,function(x){x[[2]]}))
  pvals_df <- do.call('rbind', lapply(output,function(x){x[[1]]}))
  
  results_df[,2:6] <- lapply(results_df[,2:6],as.numeric)
  
  saveRDS(pvals_df, file = paste0('selfmadestep_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR', SNR,'.RDS'))
  saveRDS(results_df, file = paste0('selfmadestep_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR', SNR,'.RDS'))
  
  
  
  
  #dataframes to save results
  
  UVILassoLMM_results_df <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(UVILassoLMM_results_df) <- c('method', 'tpr', 'fdr', 'fwer', 'coverage', 'avg_ci')
  
  UVILassoLMM_pvals_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(UVILassoLMM_pvals_df) <- c('method', 'variable', 'signal', 'real_beta', 'estimate', 'pval', 'cil', 'ciu', 'covered')
  
  print('UVILassoLMM')
  
  for(jj in 1:sim_num){
    
    X <- Xs[[jj]]
    real_beta <- real_betas[[jj]]
    Z <- Zs[[jj]]
    ys <- yss[[jj]]
    y <- ys[,1]
    print(jj)
    
    dat <- data.frame(X, Z, y, subjind)
    
    for(j in 1:its){
      
      y = ys[,j]
      dat$y <- y
      
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

      dmat <- expand.grid(rep(list(c(-1,1)),p))[1:2**(p-1),]
      ncp = lambdas[1]**2 * max(apply(dmat, 1, FUN = function(a){t(a) %*% Cinv %*% a}))
      
      
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
      
      pvals <- p.adjust(pvals, method = 'holm')
      
      selected <- pvals<=fdr_level
      UVILassoLMM_metrics <- metrics(selected, real_beta!=0)
      
      UVILassoLMM_results_df[nrow(UVILassoLMM_results_df) + 1,] <- c('UVILassoLMM', UVILassoLMM_metrics$tpr,UVILassoLMM_metrics$fdr,
                                                                   UVILassoLMM_metrics$fwer, coverage, avg_CI)
      
    }
  }
  
  
  UVILassoLMM_results_df[,2:6] <- lapply(UVILassoLMM_results_df[,2:6],as.numeric)
  
  saveRDS(UVILassoLMM_pvals_df, file = paste0('UVILassoLMM_pvals_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  saveRDS(UVILassoLMM_results_df, file = paste0('UVILassoLMM_results_n',n,'_p',p,'_p_rel',p_rel,'_q',q,'_nrsubj',nrsubj,'_nrobs',nrobs,'_SNR',SNR,'.RDS'))
  
}





stopCluster(myCluster)
