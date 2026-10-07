set.seed(2082026)
library(tidyverse)
# function to simulate data (replicates the null simulation in paper, where no variant has a direct effect on mortality, that is not through incidence, i.e logHR_m=0)
# to replicate simulation where 20% of mortality logHRs are 1.05, set the argument logHR_m=0.05 
# for null simulation, with increased heritability, set logHR_m=0, H2=0.5

sim_data <- function(N=10000, prev=0.1,   logHR_l=0.2,   logHR_m=0, H2=0.2552){
  
  ## 451 SNPs.  Uniform minor allele frequence between 1% and 50%
  MAF <- runif(451,0.01,0.5)
  
  # SNP effects (measured by log odds ratios) inversely proportional to MAF*(1-MAF) implying SNP heritabiltiy roughly independent of MAF # or bigger effects for smaller MAF
  beta <- rnorm(451,0,sd=sqrt(MAF*(1-MAF))) # I've made all the minor variants deleterious here.  This ensures that the 20\% of SNPs have signs agreeing with the PRS
  
  
  # scale logORs so SNPs represent a heritability of 16% 
  S_naive <- sum(beta^2*2*MAF*(1-MAF)/1.8^2)
  S_true <- H2/(1-H2)
  true_scale <- sqrt(S_true/S_naive)
  
   beta <- beta*true_scale
  
  ## Calculate beta0 so that overall prevalence is about 0.1 (This is assuming latent liability model)
  ## The 1.8 constant is for differneces between the standard normal and standard logistic distribuiton
  beta0 = -1.8*qnorm(1-prev)*sqrt(1+sum(2*MAF*(1-MAF)*beta^2/1.8^2))
  
  ## simulate genetic data for 10000 people.
  G <- matrix(nrow=N,ncol=451)
  for(j in 1:451) G[,j] <- rbinom(N,size=2, prob=MAF[j])-2*MAF[j] # centering the genotype variables
  eta <- beta0+G%*%beta #(this is the population PRS!)
  
  # use logistic variables to simulate case status on latent scale (this allows using the non-genetic component to specify selection bias, but is completely equivalent to simulating with a binomial).
  latent_logistic <- rlogis(N, location=0,scale=1)
  cases <- as.numeric(latent_logistic + eta > 0)
  
  # now restrict to cases. 
  data_case <- G[cases==1,]
  n_case <- nrow(data_case)
  
  logHR_l <- 0.2
  
  p <- 0.2

  mortality_snps <- abs(451*p)
  betasign <- 1 - 2*as.numeric(beta<0)
  eta <- eta/sd(eta) # scale per sd
  
  data_case <- data.frame(time=rexp(n_case,rate=exp(log(1/10)+logHR_l*latent_logistic[cases==1] + logHR_m*apply(as.matrix(data_case[,1:mortality_snps]),1,function(x){sum(x*betasign[1:mortality_snps])}))),status = rep(1, n_case),data_case, PRS=eta[cases==1],PRS_mortality_true=logHR_m*apply(as.matrix(data_case[,1:mortality_snps]),1,function(x){sum(x*betasign[1:mortality_snps])}))
  
  ### create a score of the genetic variables weighted by their effects (i.e. log HRs) on mortality
  ### We use case-only data to do this
  ### We estimate per-SNP HRs with half the data, and use the other half to estimate the HR, per sd unit of the mortality PRS (using the estimated logHRs from the first half of the data as coefficients) on survival.    In cases where the SNPs do have a population effect on progression, we also estimate the HR for the 'population' mortality PRS, per SD unit, where the coefficients are the population log HRs.
  split_n <- abs(nrow(data_case)/2)
  data_case1 <- data_case[1:split_n,]
  data_case2 <- data_case[(split_n+1):nrow(data_case),]
  log_HR_vec1 <- numeric(451)
  for(i in 1:451){
    
    form <- as.formula(paste0("Surv(time, status) ~ X", i))
    model <- coxph(form, data = data_case1)
    log_HR_vec1[i] <- summary(model)$coefficients[1]
    
  }
  log_HR_vec1[is.na(log_HR_vec1)] <- 0
  data_case2$PRS_mortality <- as.matrix(data_case2[,3:453])%*%log_HR_vec1
  data_case2$PRS_mortality <- data_case2$PRS_mortality/sd(data_case2$PRS_mortality)
  

 if(logHR_m != 0) data_case2$PRS_mortality_true <- data_case2$PRS_mortality_true/sd(data_case2$PRS_mortality_true)
  
  # estimated mortality association with incidence PRS.  Association should be negative, due to event index bias
  a1 <- summary(coxph(Surv(time, status) ~ PRS, data = data_case2))$coefficients
  
  ## estimated corrected mortality association with incidence PRS, corrected for the latent environmental factor causing event index bias.
  correction <- latent_logistic[cases==1][(split_n+1):nrow(data_case)]
  a2 <- summary(coxph(Surv(time, status) ~ PRS+correction, data = data_case2))$coefficients
  
  # estimated mortality association with estimated mortality PRS
  a3 <- summary(coxph(Surv(time, status) ~ PRS_mortality, data = data_case2))$coefficients
  
  # estimated mortality association with estimated mortality PRS, corrected for event index bias
  a4 <- summary(coxph(Surv(time, status) ~ PRS_mortality+correction, data = data_case2))$coefficients
  
  # estimated mortality association with population mortality PRS
  a5 <- summary(coxph(Surv(time, status) ~ PRS_mortality_true, data = data_case2))$coefficients
  # estimated mortality association with population mortality PRS, corrected for event index bias.  However, there should be no event index bias if the true mortality associations (log HRs) are independent of the true incidence associations
  a6 <- summary(coxph(Surv(time, status) ~ PRS_mortality_true+correction, data = data_case2))$coefficients
  #browser()
  return(list(p_values = c(a1[5],a2[1,5],a3[5],a4[1,5],a5[1,5],a6[1,5]),HRs=exp(c(a1[1],a2[1,1],a3[1],a4[1,1],a5[1],a6[1,1]))))
}

## parallelize the above function to run 1000 simulations quickly
library(parallel)
library(survival)

num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("sim_data"))
clusterEvalQ(cl, {
  library(survival)
})

results_list <- parLapply(cl, 1:1000, function(i) {
  sim_data()
})
stopCluster(cl)


df_hr <- as.data.frame(do.call(rbind, lapply(results_list, function(x) x$HRs)))

# Name columns for plotting
col_names <- c("1: PRS", "2: PRS + Latent", "3: PRS_Mort", "4: PRS_Mort + Latent", "5: PRS_Mort_true", "6: PRS_Mort_true + Latent")
colnames(df_hr) <- col_names

# Pivot data into "long" format for ggplot
df_hr_long <- pivot_longer(df_hr, cols = everything(), names_to = "Model", values_to = "Value")


# Theme for HRs: Removes horizontal (y) grids, keeps vertical (x) grids
clean_theme_hr <- theme(
  legend.position = "none",
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major.y = element_blank(), # Blanks only horizontal major lines
  panel.grid.minor.y = element_blank()  # Blanks only horizontal minor lines
)

# get rid of NAs when logHR_m=0
df_hr_long <- df_hr_long[!is.na(df_hr_long$Value),]

# Create the Hazard Ratio plot
hr_plot <- ggplot(df_hr_long, aes(x = Value, fill = Model)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Model, scales = "free", ncol = 2) + 
  theme_minimal() + # Set the base theme first
  theme(plot.title = element_text(hjust = 0.5)) + # Apply customization second
  labs(title = "Mortality Hazard Ratios", x = "Hazard Ratio", y = "Density") +
  clean_theme_hr







