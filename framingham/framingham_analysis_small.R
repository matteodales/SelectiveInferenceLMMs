library(lme4)
library(lmerTest)
library(expm)
library(dplyr)
library(reshape2)
library(mvtnorm)
library(MultBiplotR)
library(MASS)
library(R.utils)
library(glmnet)
library(gridExtra)
library(scales)
library(postcAIC)
library(ks)
library(tmg)

# simulation settings 

seed = 0

sim_num <- 500 # number of simulations

fdr_level <- 0.05
set.seed(seed)

selection_frac <- 0.5 # fraction of observations to use for selecting the variables







#####################################################
########### preparing the framingham data ###########
#####################################################

data <- read.csv('framingham.csv')

n <- nrow(data)
nrsubj <- length(unique(data$newid))

data$ID <- NULL
data$X <- NULL

data$age <- scale(data$age)
data$year <- (data$year-5)/10


pnoise = 5 # number of fixed effects, other than the intercept
p <- 3 + pnoise

var_names <- c('sex','age','year', paste("X", 1:pnoise, sep=""))
y <- scale(data$cholst, scale=FALSE)/100

SNR <- 0.5

varnoise <- as.numeric(var(y))/SNR*5

Sigmap = matrix(0.3, pnoise, pnoise) + 0.7*diag(pnoise) # covariance matrix for the covariates

noisevars <- mvrnorm(n = n, mu = rep(0,pnoise), varnoise*Sigmap)

subjind <- data$newid

X <- cbind(as.matrix(data[,c('sex','age','year')]), noisevars)
colnames(X) <- c('sex','age','year', paste("X", 1:pnoise, sep=""))

yform <- as.formula(paste("y ~ -1 + sex + age + year + ", paste("X", 1:pnoise, sep = "", collapse = "+"), "+ (1|subjind)"))

Xs <- list()

for(j in 1:sim_num){
  
  noisevars <- mvrnorm(n = n, mu = rep(0,pnoise), varnoise*Sigmap)
  
  X <- cbind(as.matrix(data[,c('sex','age','year')]), noisevars)
  colnames(X) <- c('sex','age','year', paste("X", 1:pnoise, sep=""))
  
  Xs[[j]] <- X
  
}










#####################################################
########### computing modelsets postcAIC  ###########
#####################################################


cAIC_model_sets <- list()

print('postcAIC')

for(jj in 1:sim_num){
  
  print(jj)
  
  X <- Xs[[jj]]
  
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
    cAIC_model_sets[[jj]] <- NA
    next}
  
  cAIC_model_sets[[jj]] <- cAIC_model_set
  
}





#####################################################
###########       applying methods        ###########
#####################################################



print('Naive inference')

naive_results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(naive_results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')

naive_pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(naive_pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')

for(jj in 1:sim_num){
  
  X <- Xs[[jj]]
  print(jj)
  
  dat <- data.frame(X, y, subjind)
  
  sx <- as.matrix(scale(X))
  sy <- as.vector(y)
  
  lambdaMax <- max(abs(colSums(sx*sy)))
  lambda_vec <- exp(seq(from = log(lambdaMax),
                        to = log(lambdaMax * 0.0001),
                        length.out = 50))
  
  PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=dat,verbose=FALSE)
  Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
  Q.start<-as.numeric(VarCorr(PQL)[1,1])
  sel_vec <- matrix(0,p,length(lambda_vec))
  
  BIC_vec <- rep(Inf, length(lambda_vec))
  # Early stopping parameters
  patience_limit <- 10  # Max consecutive iterations without improvement
  patience_count <- 0  # Counter for consecutive non-improvements
  best_bic <- Inf      # Keep track of the best BIC value
  current_bic <- Inf
  
  for (l in 1:length(lambda_vec)) {
    glm3 <- try(glmmLasso(
      fix = y ~ -1 + sex + age + year + X1 + X2 + X3 + X4 + X5,
      rnd = list(subjind = ~1),
      data = dat,
      lambda = lambda_vec[l],
      control = list(start = Delta.start[l, ], q_start = Q.start[l], center = FALSE, standardize=FALSE)
    ), silent=TRUE)
    
    
    if (!inherits(glm3, "try-error")) {
      current_bic <- glm3$bic
      BIC_vec[l] <- current_bic
      sel_vec[,l] <- as.numeric(glm3$coefficients != 0)
      Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
      Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      
      # Check for improvement
      if (!is.numeric(current_bic) || is.nan(current_bic)) next
      if (current_bic < best_bic) {
        if (current_bic < best_bic) best_bic <- current_bic
        patience_count <- 0  # Reset patience
      } else {
        patience_count <- patience_count + 1  # Increment patience
      }
      
      #Early stopping condition
      if (patience_count >= patience_limit) break
      
    }
  }
  
  opt3 <- which.min(BIC_vec)
  selected <- sel_vec[,opt3]
  
  
  selected_tot_names <- var_names[selected==1]
  
  selected_tot_yform <- as.formula(
    paste(" y ~ -1 + ", paste(selected_tot_names, sep = "", collapse = "+"), "+ (1|subjind)", sep = "")
  )

  
  suppressWarnings(suppressMessages(sel_mod <- lmer(formula = selected_tot_yform, data = dat)))
  
  
  if(sum(selected)==1){
    
    pvals_tot <- as.data.frame(t(coef(summary(sel_mod))[,c(1,5)]))
    row.names(pvals_tot) <- selected_tot_names
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
      naive_pvals_df[nrow(naive_pvals_df) + 1,] <- c('Naive', var_names[i], as.numeric(real_beta[i]!=0), pvals_tot[k,1], pvals_tot[k,2],
                                                     confint_tot[k,1], confint_tot[k,2])
      k = k + 1
    }
  }
  
  naive_pvals_df[,3:7] <- lapply(naive_pvals_df[,3:7],as.numeric)
  
  avg_CI <- mean(naive_pvals_df[(nrow(naive_pvals_df)-sum(selected)):nrow(naive_pvals_df), 7] - naive_pvals_df[(nrow(naive_pvals_df)-sum(selected)):nrow(naive_pvals_df), 6])
  
  corrected_pvals_tot <- p.adjust(pvals_tot[,2], method = 'holm')
  
  selected[selected==1] <- corrected_pvals_tot<=fdr_level
  tot_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = pvals_tot['sex',1] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = pvals_tot['age',1] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = pvals_tot['year',1] else coef_year = 0
  
  naive_results_df[nrow(naive_results_df) + 1,] <- c('Naive', tot_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
}

naive_results_df[,2:9] <- lapply(naive_results_df[,2:9],as.numeric)

naive_results_df[,2:9] <- lapply(naive_results_df[,2:9],as.numeric)
colMeans(naive_results_df[,2:9])

saveRDS(naive_pvals_df, file = paste0('framingham_naive_pvals_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(naive_results_df, file = paste0('framingham_naive_results_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))







print('Data Splitting')

datasplitting_results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(datasplitting_results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')

datasplitting_pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(datasplitting_pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')


for(jj in 1:sim_num){
  
  print(jj)
  X <- Xs[[jj]]
  
  dat <- data.frame(X, y, subjind)
  
  # data splitting: the data is split in two, using the first half for selection and the second for inference
  
  selection_sample <- sample(1:nrsubj, round(nrsubj*selection_frac))
  selection_dat <- dat[dat$subjind %in% selection_sample,]
  inference_dat <- dat[!(dat$subjind %in% selection_sample),]
  
  X_sel <- as.matrix(selection_dat[,1:p])
  y_sel <- selection_dat$y
  
  sx <- as.matrix(scale(X_sel))
  sy <- as.vector(y_sel)
  
  lambdaMax <- max(abs(colSums(sx*sy)))
  lambda_vec <- exp(seq(from = log(lambdaMax),
                        to = log(lambdaMax * 0.0001),
                        length.out = 50))
  
  selection_dat$subjind <- droplevels(selection_dat$subjind)
  
  
  PQL<-glmmPQL(y~-1,random = ~1|subjind,family='gaussian',data=selection_dat,verbose=FALSE)
  Delta.start<-as.matrix(t(c(rep(0,p),as.numeric((PQL$coef$random$subjind)))))
  Q.start<-as.numeric(VarCorr(PQL)[1,1])
  
  BIC_vec <- rep(Inf, length(lambda_vec))
  sel_vec <- matrix(0, p, length(lambda_vec))
  
  # Early stopping parameters
  patience_limit <- 10  # Max consecutive iterations without improvement
  patience_count <- 0  # Counter for consecutive non-improvements
  best_bic <- Inf      # Keep track of the best BIC value
  current_bic <- Inf
  
  for (l in 1:length(lambda_vec)) {
    glm3 <- try(glmmLasso(
      fix = y ~ -1 + sex + age + year + X1 + X2 + X3 + X4 + X5,
      rnd = list(subjind = ~1),
      data = selection_dat,
      lambda = lambda_vec[l],
      control = list(start = Delta.start[l, ], q_start = Q.start[l], center = FALSE, standardize=FALSE)
    ), silent=TRUE)
    
    
    if (!inherits(glm3, "try-error")) {
      current_bic <- glm3$bic
      BIC_vec[l] <- current_bic
      sel_vec[, l] <- as.numeric(glm3$coefficients != 0)
      Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
      Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
      
      # Check for improvement
      if (!is.numeric(current_bic) || is.nan(current_bic)) next
      if (current_bic < best_bic) {
        if (current_bic < best_bic) best_bic <- current_bic
        patience_count <- 0  # Reset patience
      } else {
        patience_count <- patience_count + 1  # Increment patience
      }
      
      #Early stopping condition
      if (patience_count >= patience_limit) break
      
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
    row.names(pvals_split) <- selected_split_names
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
      datasplitting_pvals_df[nrow(datasplitting_pvals_df) + 1,] <- c('Data Splitting', var_names[i], as.numeric(real_beta[i]!=0), pvals_split[k,1], pvals_split[k,2],
                                                                     confint_split[k,1], confint_split[k,2] )
      k = k + 1
    }
  }
  
  datasplitting_pvals_df[,3:7] <- lapply(datasplitting_pvals_df[,3:7],as.numeric)
  
  
  
  avg_CI <- mean(datasplitting_pvals_df[(nrow(datasplitting_pvals_df)-sum(selected)):nrow(datasplitting_pvals_df), 7] - datasplitting_pvals_df[(nrow(datasplitting_pvals_df)-sum(selected)):nrow(datasplitting_pvals_df), 6])
  
  corrected_pvals_split <- p.adjust(pvals_split[,2], method = 'holm')
  
  selected[selected==1] <- corrected_pvals_split<=fdr_level
  split_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = pvals_split['sex',1] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = pvals_split['age',1] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = pvals_split['year',1] else coef_year = 0
  
  datasplitting_results_df[nrow(datasplitting_results_df) + 1,] <- c('Data Splitting', split_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
  
}

# 
datasplitting_results_df[,2:9] <- lapply(datasplitting_results_df[,2:9],as.numeric)
colMeans(datasplitting_results_df[,2:9])
# 
saveRDS(datasplitting_pvals_df, file = paste0('framingham_datasplitting_pvals_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(datasplitting_results_df, file = paste0('framingham_datasplitting_results_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))







print('postcAIC')

output <- foreach(j=1:sim_num, .packages = c("glmnet","lme4",'lmerTest','nlme','R.utils','postcAIC','ks','tmg','mgcv')) %dopar% {
  
  results_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')
  
  pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')
  
  
  print(j)
  
  X <- Xs[[j]]
  
  
  timeout <<- 0
  
  
  if(any(is.na(cAIC_model_setss[[j]]))) next
  
  cAIC_model_set <- cAIC_model_setss[[j]]
  
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
  
  # Post-cAIC CI for mixed and fixed parameters -------------------------------------
  
  t1 <- Sys.time()
  
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
      
    }, timeout = 6000)
  }, TimeoutException = function(ex) {
    message("Timeout. Skipping.")
    results_df[nrow(results_df) + 1,] <<- c('postcAIC', NA, NA, NA, NA, NA)
    timeout<<-1
  })
  
  
  if(timeout==1) next
  
  beta_postcAIC_CI_up <- postcAIC_CI_results$beta_postcAIC_CI_up
  beta_postcAIC_CI_do <- postcAIC_CI_results$beta_postcAIC_CI_do
  beta_postcAIC_pval <- postcAIC_CI_results$beta_postcAIC_pval
  
  selected <- modelset_matrix[cAIC_min,]
  
  k <- 1
  for(i in seq_len(p)){
    
    if(selected[i]!=0){
      
      pvals_df[nrow(pvals_df)+1,] <- c('postcAIC', var_names[i], as.numeric(real_beta[i]!=0), beta_sel[k],
                                       beta_postcAIC_pval[k], beta_postcAIC_CI_do[i], beta_postcAIC_CI_up[i])
      
      k <- k + 1
      
    }
  }
  
  pvals_df[,3:7] <- lapply(pvals_df[,3:7],as.numeric)
  
  avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 6])
  
  pvals <- p.adjust(beta_postcAIC_pval, method = 'holm')
  
  selected[selected==1] <- pvals<=fdr_level
  postcAIC_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = pvals_df[pvals_df$variable=='sex','estimate'] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = pvals_df[pvals_df$variable=='age','estimate'] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = pvals_df[pvals_df$variable=='year','estimate'] else coef_year = 0
  
  results_df[nrow(results_df) + 1,] <- c('postcAIC', postcAIC_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
  return(list(results = results_df, pvals = pvals_df))
  
}

results_df <- do.call('rbind', lapply(output,function(x){x[[1]]}))
pvals_df <- do.call('rbind', lapply(output,function(x){x[[2]]}))


results_df[,2:9] <- lapply(results_df[,2:9],as.numeric)

saveRDS(pvals_df, file = paste0('framingham_postcAIC_pvals_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(results_df, file = paste0('framingham_postcAIC_results_pnoise',pnoise,'_SNR',SNR,'.RDS'))

stopCluster(myCluster)




#dataframes to #save results

print('selfmade-lasso')

selfmadelasso_results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(selfmadelasso_results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')

selfmadelasso_pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(selfmadelasso_pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')

for(jj in 1:sim_num){
  
  print(jj)
  X <- Xs[[jj]]
  
  
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
    patience_limit <- 5  # Max consecutive iterations without improvement
    patience_count <- 0  # Counter for consecutive non-improvements
    best_bic <- Inf      # Keep track of the best BIC value
    current_bic <- Inf
    
    for (l in 1:length(lambda_vec)) {
      glm3 <- try(glmmLasso(
        fix = as.formula(paste("y ~ - 1 +", paste(var_names, sep = "", collapse = "+"), sep = "")),
        rnd = list(subjind = ~1),
        data = dat,
        lambda = lambda_vec[l],
        control = list(center = FALSE, start = Delta.start[nrow(Delta.start), ], q_start = Q.start[nrow(Delta.start)], standardize=FALSE)
      ), silent = TRUE)
      
      if (!inherits(glm3, "try-error")) {
        current_bic <- glm3$bic
        BIC_vec[l] <- current_bic
        sel_vec[l, ] <- as.numeric(glm3$coefficients != 0)
        Delta.start <- rbind(Delta.start, glm3$Deltamatrix[glm3$conv.step, ])
        Q.start <- c(Q.start, glm3$Q_long[[glm3$conv.step + 1]])
        
        # Check for improvement
        if (!is.numeric(current_bic) || is.nan(current_bic)) next
        if (current_bic < 1.5*best_bic) {
          if (current_bic < best_bic) best_bic <- current_bic
          patience_count <- 0  # Reset patience
        } else {
          patience_count <- patience_count + 1  # Increment patience
        }
        
        # Early stopping condition
        if (patience_count >= patience_limit) break
      }
    }
    
    opt <- which.min(BIC_vec)
    
    names_vec = var_names[sel_vec[opt,]==1]
    
    return(c(names_vec))
  }
  
  trace <- FALSE
  
  dat <- data.frame(X, y, subjind)
  
  sx <- as.matrix(scale(X))
  sy <- as.vector(y)
  
  lambdaMax <- max(abs(colSums(sx*sy)))
  
  lambda_vec <- exp(seq(from = log(lambdaMax),
                        to = log(lambdaMax * 0.0001),
                        length.out = 50))
  
  names_vec = selFun(y)
  
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
  
  
  
  selected <- rep(0,p)
  
  k <- 1
  for(i in seq_len(p)){
    
    if(var_names[i] %in% selection){
      
      selfmadelasso_pvals_df[nrow(selfmadelasso_pvals_df)+1,] <- c('selfmade-lasso', var_names[i], as.numeric(real_beta[i]!=0), beta[k],
                                                                 r$pval[k], r$cil[k], r$ciu[k])
      
      selected[i] <- 1
      k <- k +1
      
    }
  }
  
  selfmadelasso_pvals_df[,3:7] <- lapply(selfmadelasso_pvals_df[,3:7],as.numeric)
  selfmadelasso_pvals_df[is.infinite(selfmadelasso_pvals_df[,6]),6] <- NA
  selfmadelasso_pvals_df[is.infinite(selfmadelasso_pvals_df[,7]),7] <- NA
  
  avg_CI <- mean(selfmadelasso_pvals_df[(nrow(selfmadelasso_pvals_df)-sum(selected)):nrow(selfmadelasso_pvals_df), 7] - selfmadelasso_pvals_df[(nrow(selfmadelasso_pvals_df)-sum(selected)):nrow(selfmadelasso_pvals_df), 6], na.rm=TRUE)
  
  pvals <- p.adjust(r$pval, method = 'holm')
  
  selected[selected==1] <- pvals<=fdr_level
  selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
  rugamer_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = r[r$variable=='sex','tstat'] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = r[r$variable=='age','tstat'] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = r[r$variable=='year','tstat'] else coef_year = 0
  
  selfmadelasso_results_df[nrow(selfmadelasso_results_df) + 1,] <- c('selfmade-lasso', rugamer_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
}

selfmadelasso_results_df[,2:9] <- lapply(selfmadelasso_results_df[,2:9],as.numeric)

saveRDS(selfmadelasso_pvals_df, file = paste0('framingham_selfmadelasso_pvals_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(selfmadelasso_results_df, file = paste0('framingham_selfmadelasso_results_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))







print('selfmade-step')

output <- foreach(j=1:sim_num, .packages = c("glmnet","lme4",'lmerTest')) %dopar% {
  
  modFun <- function(yy, dat = dat)
  {
    dat$y <- as.numeric(yy)
    suppressWarnings(suppressMessages(lmer(y ~ -1 + sex + age + year + X1 + X2 + X3 + X4 + X5 + (1|subjind), REML = FALSE, data = dat)))
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
  

  
  results_df <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')
  
  pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')
  
  print(j)
  
  X <- Xs[[j]]
  
  
  dat <<- data.frame(X, y, subjind)
  
  mod <- modFun(yy=y, dat=dat)
  
  final_model <- selFun(mod)
  beta <- fixef(final_model)
  
  
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
                                trace=FALSE))
  
  
  r <- do.call("rbind", r$selinf)
  r$variable <- selection
  
  selected <- rep(0,p)
  
  k <- 1
  for(i in seq_len(p)){
    
    if(var_names[i] %in% selection){
      
      pvals_df[nrow(pvals_df)+1,] <- c('selfmade-step', var_names[i], as.numeric(real_beta[i]!=0), beta[k],
                                       r$pval[k], r$cil[k], r$ciu[k])
      
      selected[i] <- 1
      k <- k +1
      
    }
  }
  
  pvals_df[,3:7] <- lapply(pvals_df[,3:7],as.numeric)
  pvals_df[is.infinite(pvals_df[,6]),6] <- NA
  pvals_df[is.infinite(pvals_df[,7]),7] <- NA
  
  avg_CI <- mean(pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 7] - pvals_df[(nrow(pvals_df)-sum(selected)):nrow(pvals_df), 6], na.rm=TRUE)
  
  pvals <- p.adjust(r$pval, method = 'holm')
  
  selected[selected==1] <- pvals<=fdr_level
  selected[is.na(selected)] <- rep(0,sum(is.na(selected)))
  rugamer_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = r[r$variable=='sex','tstat'] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = r[r$variable=='age','tstat'] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = r[r$variable=='year','tstat'] else coef_year = 0
  
  results_df[nrow(results_df) + 1,] <- c('selfmade-step', rugamer_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
  return(list(pvals_df, results_df))
}

results_df <- do.call('rbind', lapply(output,function(x){x[[2]]}))
pvals_df <- do.call('rbind', lapply(output,function(x){x[[1]]}))

results_df[,2:9] <- lapply(results_df[,2:9],as.numeric)

saveRDS(pvals_df, file = paste0('framingham_selfmadestep_pvals_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(results_df, file = paste0('framingham_selfmadestep_results_pnoise',pnoise,'_SNR',SNR,'.RDS'))









print('UVILassoLMM')

UVILassoLMM_results_df <- data.frame(matrix(nrow = 0, ncol = 9))
colnames(UVILassoLMM_results_df) <- c('method', 'fwer', 'avg_ci', 'sel_sex', 'coef_sex', 'sel_age', 'coef_age', 'sel_year', 'coef_year')

UVILassoLMM_pvals_df <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(UVILassoLMM_pvals_df) <- c('method', 'variable', 'signal', 'estimate', 'pval', 'cil', 'ciu')


for(jj in 1:sim_num){
  
  print(jj)
  X <- Xs[[jj]]
  
  
  dat <- data.frame(X, y, subjind)
  
  suppressWarnings(suppressMessages(mod <- lmer(formula = yform, data = dat)))
  V <- getME(mod,'Z')%*%(sigma(mod)^2*getME(mod,'Lambda')%*%getME(mod,'Lambdat'))%*%t(getME(mod,'Z'))+ diag(sigma(mod)^2,n)
  V_inv <- solve(V)
  
  V_menunmezz <- matrixsqrtinv(V)
  cv_mod <- cv.glmnet(x = V_menunmezz%*%X, y = V_menunmezz%*%y, intercept=FALSE, standardize=FALSE)
  beta <- coef(cv_mod, s = cv_mod$lambda.min)[-1]
  
  selected <- as.numeric(beta!=0)
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
    
    UVILassoLMM_pvals_df[nrow(UVILassoLMM_pvals_df)+1,] <- c('UVILassoLMM', var_names[i], as.numeric(real_beta[i]!=0),
                                                           beta[i], pvals[k], cil, ciu)
    
    k <- k +1
    
  }
  
  UVILassoLMM_pvals_df[,3:7] <- lapply(UVILassoLMM_pvals_df[,3:7],as.numeric)
  
  avg_CI <- mean(UVILassoLMM_pvals_df[(nrow(UVILassoLMM_pvals_df)-sum(selected)):nrow(UVILassoLMM_pvals_df), 7] - UVILassoLMM_pvals_df[(nrow(UVILassoLMM_pvals_df)-sum(selected)):nrow(UVILassoLMM_pvals_df), 6])
  
  pvals <- p.adjust(pvals, method = 'holm')
  
  selected <- as.numeric(pvals<=fdr_level)
  UVILassoLMM_metrics <- metrics(selected, real_beta!=0)
  
  sel_sex = selected[1]
  if (sel_sex==1) coef_sex = beta[1] else coef_sex = 0
  sel_age = selected[2]
  if (sel_age==1) coef_age = beta[2] else coef_age = 0
  sel_year = selected[3]
  if (sel_year==1) coef_year = beta[3] else coef_year = 0
  
  UVILassoLMM_results_df[nrow(UVILassoLMM_results_df) + 1,] <- c('UVILassoLMM', UVILassoLMM_metrics$fwer, avg_CI, sel_sex, coef_sex, sel_age, coef_age, sel_year, coef_year)
  
}


UVILassoLMM_results_df[,2:9] <- lapply(UVILassoLMM_results_df[,2:9],as.numeric)

saveRDS(UVILassoLMM_pvals_df, file = paste0('framingham_UVILassoLMM_pvals_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))
saveRDS(UVILassoLMM_results_df, file = paste0('framingham_UVILassoLMM_results_framingham_pnoise',pnoise,'_SNR',SNR,'.RDS'))

