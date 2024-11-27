##=================================================================##
##                                                                 ##
##                 Fit integrated population model                 ##
##                                                                 ##
##=================================================================##
pacman::p_load(here,tidyverse,rstan,gtools,salmonIPM)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

covar_effects<-TRUE ## TRUE/FALSE
biastest<-FALSE ## TRUE/FALSE

##============================================================## data
fish_dat<-read.csv(here("data","IPM_fish_dat_all.csv"))
if(biastest) fish_dat<-read.csv(here("data","IPM_fish_dat_biastest.csv"))
fish_dat<-fish_dat %>% dplyr::select(-F_rate_NA)  

if(covar_effects) {
   covar_dat<-read.csv(here("data","IPM_covar_dat_selected.csv"))
   fish_dat<-fish_dat %>%
      left_join(covar_dat,by='year') %>%
      na.omit() ## no NAs allowed in covariate data
   ## covariate models
   par_models<-list(s_SS~SST+pinks+av_CFS_min,R~NPGO_2+SST_4+pinks_4)
}else{
   par_models<-NULL
}

##=========================================================## fit IPM
IPM_fit<-salmonIPM(
   life_cycle="SSiter",
   pool_pops=TRUE,
   SR_fun="Ricker", 
   par_models=par_models, 
   fish_data=fish_dat,
   prior=list(
      mu_alpha~normal(1.5,0.5),
      mu_p~dirichlet(c(1,2,47,44,5,1)),
      mu_SS~beta(1.5,3)),
   pars=c(stan_pars("IPM_SSiter_pp")),
   chains=3,
   iter=2000,
   warmup=1000, 
   control=list(adapt_delta=0.98,max_treedepth=12))

##=====================================## save model data and results
if(covar_effects) { 
   write.csv(fish_dat,here("output","IPM_fish_dat_with_covars.csv"),
             row.names=F)
   if(biastest){
      saveRDS(IPM_fit,here("output","IPM_fit_with_covars_biastest.rds"))
   }else{
      saveRDS(IPM_fit,here("output","IPM_fit_with_covars.rds"))  
   }
}else{ 
   write.csv(fish_dat,here("output","IPM_fish_dat_without_covars.csv") ,
             row.names=F)
   saveRDS(IPM_fit,here("output","IPM_fit_without_covars.rds"))
}

##=================================================================##