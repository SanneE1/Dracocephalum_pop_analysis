
run_ipm_analyses <- function(VR_FLM, state_independent_variables, lag = 48){


# environmental params set to mean
mean_env_params <- list(
  lags = lag,
  yr_increase_temp = 0,
  yr_increase_prcp = 0,
  herb_shading = 2,
  shrub_shading = 0,
  slope = 0,
  soil = 3
)

# param/model list 
params <- list(
  surv_mod = VR_FLM$surv,
  s_int = coef(VR_FLM$surv)[1],
  s_stems = coef(VR_FLM$surv)[2],
  s_site_CR = 0,
  s_site_HK = coef(VR_FLM$surv)[3],
  s_site_KS = coef(VR_FLM$surv)[4],
  s_site_RU = coef(VR_FLM$surv)[5],
  s_herb = coef(VR_FLM$surv)[6],
  grow_mod = VR_FLM$growth,
  g_int = coef(VR_FLM$growth)[1],
  g_stems = coef(VR_FLM$growth)[2],
  g_herb = coef(VR_FLM$growth)[3],
  grow_sd = sd(resid(VR_FLM$growth)),
  pflower_mod = VR_FLM$flower_p,
  fp_int = coef(VR_FLM$flower_p)[1],
  fp_stems = coef(VR_FLM$flower_p)[2],
  fp_site_CR = 0,
  fp_site_HK = coef(VR_FLM$flower_p)[3],
  fp_site_KS = coef(VR_FLM$flower_p)[4],
  fp_site_RU = coef(VR_FLM$flower_p)[5],
  fp_shrub = coef(VR_FLM$flower_p)[6],
  fp_slope = coef(VR_FLM$flower_p)[7],
  fp_soil = coef(VR_FLM$flower_p)[8],
  pabort_mod = VR_FLM$abort_p,
  ab_int = coef(VR_FLM$abort_p)[1],
  ab_stems = coef(VR_FLM$abort_p)[2],
  ab_site_CR = 0,
  ab_site_HK = coef(VR_FLM$abort_p)[3],
  ab_site_KS = coef(VR_FLM$abort_p)[4],
  ab_site_RU = coef(VR_FLM$abort_p)[5],
  ab_shrub = coef(VR_FLM$abort_p)[6],
  ab_soil = coef(VR_FLM$abort_p)[7],
  nseed_mod = VR_FLM$n_seeds,
  ns_int = coef(VR_FLM$n_seeds)[1],
  ns_stems = coef(VR_FLM$n_seeds)[2],
  ns_site_CR = 0,
  ns_site_HK = coef(VR_FLM$n_seeds)[3],
  ns_site_KS = coef(VR_FLM$n_seeds)[4],
  ns_site_RU = coef(VR_FLM$n_seeds)[5],
  ns_slope = coef(VR_FLM$n_seeds)[6],
  seed_surv1 = 0.578,  ## get a few more digits 
  seed_surv2 = 0.154,  ## get a few more digits
  seed_surv3 = 0.667,  ## get a few more digits
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     ## Using mean for now
  # germ_alpha = mean(state_independent_variables$est_germination_rate$alpha), ## Switch to using beta dist. later
  # germ_kappa = mean(state_independent_variables$est_germination_rate$kappa),
  sdl_surv_mod = state_independent_variables$sdl_surv,
  sdl_s_int = lme4::fixef(state_independent_variables$sdl_surv)[1],
  sdl_s_site_CR = 0,
  sdl_s_site_HK = lme4::fixef(state_independent_variables$sdl_surv)[2],
  sdl_s_site_KS = lme4::fixef(state_independent_variables$sdl_surv)[3],
  sdl_s_site_RU = lme4::fixef(state_independent_variables$sdl_surv)[4],
  sdl_s_herb = lme4::fixef(state_independent_variables$sdl_surv)[5],
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_d_site_CR = 0,
  sdl_d_site_HK = lme4::fixef(state_independent_variables$sdl_size_d)[2],
  sdl_d_site_KS = lme4::fixef(state_independent_variables$sdl_size_d)[3],
  sdl_d_site_RU = lme4::fixef(state_independent_variables$sdl_size_d)[4],
  sdl_d_herb = lme4::fixef(state_independent_variables$sdl_size_d)[5],
  sdl_size_d_sd = sd(resid(state_independent_variables$sdl_size_d))
)



## Set integration params
L <- min(VR_FLM$growth$model$ln_stems_t0, na.rm = T)
U <- max(VR_FLM$growth$model$ln_stems_t1, na.rm = T) * 1.2
n = 100


##------------------------------------------------------------------------------
## Calculate lambda for mean env_variables
##------------------------------------------------------------------------------

mean_ipm_CR <- run_ipm(params = params,
                    env_params = mean_env_params,
                    locality = "CR", n_it = 1000,
                    U = U, L = L, n = n)

mean_ipm_HK <- run_ipm(params = params,
                       env_params = mean_env_params,
                       locality = "HK", n_it = 1000,
                       U = U, L = L, n = n)

mean_ipm_KS <- run_ipm(params = params,
                       env_params = mean_env_params,
                       locality = "KS", n_it = 1000,
                       U = U, L = L, n = n)

mean_ipm_RU <- run_ipm(params = params,
                       env_params = mean_env_params,
                       locality = "RU", n_it = 1000,
                       U = U, L = L, n = n)

saveRDS(list(mean_ipm_CR = mean_ipm_CR,
             mean_ipm_HK = mean_ipm_HK,
             mean_ipm_KS = mean_ipm_KS,
             mean_ipm_RU = mean_ipm_RU),
        file = "results/rds/mean_ipms_per_loc.rds")


mean_lambda_CR <- lambda(mean_ipm_CR)
mean_lambda_HK <- lambda(mean_ipm_HK)
mean_lambda_KS <- lambda(mean_ipm_KS)
mean_lambda_RU <- lambda(mean_ipm_RU)

print("mean lambda for CR")
print(mean_lambda_CR)

print("mean lambda for HK")
print(mean_lambda_HK)

print("mean lambda for KS")
print(mean_lambda_KS)

print("mean lambda for RU")
print(mean_lambda_RU)


##------------------------------------------------------------------------------
## Calculate lambda at different co-variate levels
##------------------------------------------------------------------------------

### Loop through different populations and env_param levels 
localities <- c("CR", "HK", "KS", "RU")
herb_shading <- seq(0,10, length.out = 11)
shrub_shading <- seq(0,10, length.out = 11)
slope <- seq(0, 80, length.out = 9)
soil <- seq(0, 15, length.out = 7)

df_env <- rbind(expand.grid(localities = localities, 
                      herb_shading = herb_shading, 
                      shrub_shading = shrub_shading,
                      slope = 0,
                      soil = 0),
                expand.grid(localities = localities,
                            herb_shading = 0,
                            shrub_shading = 0,
                            slope = slope,
                            soil = 0),
                expand.grid(localities = localities,
                            herb_shading = 0,
                            shrub_shading = 0,
                            slope = 0,
                            soil = soil)
) %>% mutate(localities = as.character(localities))

df_env$lambda <- NA


### Set up parallel
cl <- makeCluster(detectCores() - 3)
clusterExport(cl, c("df_env", "ipm_loop", "run_ipm", "params", 
                    "U", "L", "n", "lag",
                    "sampling_env", "FLM_clim_predict"))
clusterEvalQ(cl, c(library("ipmr"), library("dplyr")))

df <- parLapplyLB(cl,
                  as.list(c(1:nrow(df_env))),
                  function(x) ipm_loop(i = x, df_env = df_env,
                                       params = params, n_it = 1000, 
                                       U = U, L = L, n = n)) %>% bind_rows()

write.csv(df, file = "results/overview_lambda_env_levels.csv", 
          row.names = F)

##------------------------------------------------------------------------------
## Calculate parameter level sensitivity/elasticity
##------------------------------------------------------------------------------

pert_param <- list( 
  "s_int", "s_stems",
  "g_int", "g_stems", "grow_sd",
  "fp_ing", "fp_stems",
  "ab_int", "ab_stems",
  "ns_int", "ns_stems",
  "seed_surv1", "seed_surv2", "seed_surv3",
  "germ_mean", 
  "sdl_s_int", 
  "sdl_d_int", "sdl_size_d_sd"
  
)

clusterExport(cl, c("pert_param", "elast_ipm"))

elas <- parLapplyLB(cl,
                    pert_param,
                    function(x) elast_ipm(params = params,
                                          mean_env_params = mean_env_params,
                                          parameter_perturb = x,
                                          n_it = 1000, 
                                          U = U, L = L, n = n)) %>% bind_rows()

write.csv(elas, file = "results/parameter_sensitivity_elasticity.csv", 
          row.names = F)



return(list(
  mean_ipms = list(mean_ipm_CR = mean_ipm_CR,
                   mean_ipm_HK = mean_ipm_HK,
                   mean_ipm_KS = mean_ipm_KS,
                   mean_ipm_RU = mean_ipm_RU),
  lambda_overview_env_params = df,
  parameter_elasticity = elas
  )
)




}