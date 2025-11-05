#### PARAMETER SWEEP LISTS ####


######## EXPERIMENT 1 #########

set.seed(1234)

# NEUTRAL SNAPSHOT -------------------------------------------------------------

# Innovation rate
n_snap_mu_params <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.025, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.05, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.075, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.125, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.15, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.175, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.2, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000)
)

# Population size
n_snap_N_params <- list(
  list(N=10, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=50, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=150, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=200, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=250, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=300, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=350, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=400, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000)
)

# Time series
n_snap_time_params <- list(
  list(N=100, mu=0.01, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=200, timesteps=200, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=500, timesteps=500, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=750, timesteps=750, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=1500, timesteps=1500, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=2000, timesteps=2000, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=2500, timesteps=2500, p_value_lvl=0.05, n_runs=1000),
  list(N=100, mu=0.01, burnin=3000, timesteps=3000, p_value_lvl=0.05, n_runs=1000)
)


# MULTIVARIATE SWEEPS

# N x mu
n_snap_N_x_mu <- list(
  # low N x low mu
  list(N=25, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  
  # low N x high mu
  list(N=25, mu=0.4, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  
  # mid N x mid mu
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  
  # high N x low mu
  list(N=400, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000),
  
  # high N x high mu
  list(N=400, mu=0.4, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000)
)


# N x t
n_snap_N_x_time <- list(
  # low N x low t
  list(N=25, mu=0.05, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  
  # low N x high t
  list(N=25, mu=0.05, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000),
  
  # mid N x mid t
  list(N=100, mu=0.05, burnin=1000, timesteps=500, p_value_lvl=0.05, n_runs=1000),
  
  # high N x low t
  list(N=400, mu=0.05, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  
  # high N x high t
  list(N=400, mu=0.05, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000)
)


# mu x t
n_snap_mu_x_t <- list(
  # low mu x low t
  list(N=100, mu=0.01, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  
  # low mu x high t
  list(N=100, mu=0.01, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000),
  
  # mid mu x mid t
  list(N=100, mu=0.1, burnin=1000, timesteps=500, p_value_lvl=0.05, n_runs=1000),
  
  # high mu x low t
  list(N=100, mu=0.4, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000),
  
  # high mu x high t
  list(N=100, mu=0.4, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000)
)



# Run and store

library(purrr)
library(tibble)

# INNOVATION RATE (MU)
n_snap_mu_results <- map_dfr(n_snap_mu_params, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_mu_results, "tables/n_snap_output/n_snap_mu_params.xlsx") # mu


# POPULATION SIZE (N)
n_snap_N_results <- map_dfr(n_snap_N_params, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_N_results, "tables/n_snap_output/n_snap_N_params.xlsx") # N


# TIME SERIES
n_snap_time_results <- map_dfr(n_snap_time_params, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_time_results, "tables/n_snap_output/n_snap_time_params.xlsx") # time series

# N x mu
n_snap_N_x_mu_results <- map_dfr(n_snap_N_x_mu, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_N_x_mu_results, "tables/n_snap_output/int/n_snap_N_x_mu.xlsx")


# N x t
n_snap_N_x_time_results <- map_dfr(n_snap_N_x_time, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_N_x_time_results, "tables/n_snap_output/int/n_snap_N_x_t.xlsx")

# mu x t
n_snap_mu_x_t_results <- map_dfr(n_snap_mu_x_t, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_snap_mu_x_t_results, "tables/n_snap_output/int/n_snap_mu_x_t.xlsx")

# PRINT RESULTS
print(n_snap_mu_results)
print(n_snap_N_results)
print(n_snap_time_results)
print(n_snap_N_x_mu_results)
print(n_snap_N_x_time_results)
print(n_snap_mu_x_t_results)












# NEUTRAL TIME AVERAGING -------------------------------------------------------

# Innovation rate
n_ta_mu_params <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.025, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.05, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.075, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.125, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.15, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.175, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.2, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20)
)

# Population size
n_ta_N_params <- list(
  list(N=10, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=50, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=150, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=200, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=250, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=300, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=350, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=400, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20)
)

# Time series
n_ta_time_params <- list(
  list(N=100, mu=0.01, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=200, timesteps=200, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=500, timesteps=500, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=750, timesteps=750, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=1500, timesteps=1500, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=2000, timesteps=2000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=2500, timesteps=2500, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=3000, timesteps=3000, p_value_lvl=0.05, n_runs=1000, time_window = 20)
)

# Time window size
n_ta_tw_params <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 5),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 10),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 20),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 40),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 50),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 70),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 80),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 100),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window = 200)
)

# N x mu
n_ta_N_x_mu <- list(
  # low N x low mu
  list(N=25, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # low N x high mu
  list(N=25, mu=0.4, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # mid N x mid mu
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high N x low mu
  list(N=400, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high N x high mu
  list(N=400, mu=0.4, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=1000, time_window=20)
)


# N x t
n_ta_N_x_time <- list(
  # low N x low t
  list(N=25, mu=0.05, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # low N x high t
  list(N=25, mu=0.05, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # mid N x mid t
  list(N=100, mu=0.05, burnin=1000, timesteps=500, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high N x low t
  list(N=400, mu=0.05, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high N x high t
  list(N=400, mu=0.05, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000, time_window=20)
)


# mu x t
n_ta_mu_x_t <- list(
  # low mu x low t
  list(N=100, mu=0.01, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # low mu x high t
  list(N=100, mu=0.01, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # mid mu x mid t
  list(N=100, mu=0.1, burnin=1000, timesteps=500, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high mu x low t
  list(N=100, mu=0.4, burnin=1000, timesteps=100, p_value_lvl=0.05, n_runs=1000, time_window=20),
  
  # high mu x high t
  list(N=100, mu=0.4, burnin=1000, timesteps=2000, p_value_lvl=0.05, n_runs=1000, time_window=20)
)

# Run and store

# INNOVATION RATE (MU)
n_ta_mu_results <- map_dfr(n_ta_mu_params, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = sim$mean_accuracy,
         "FPR" = sim$mean_FPR,
         "%NA" = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_mu_results, "tables/n_ta_output/n_ta_mu_params2.xlsx") # mu


# POPULATION SIZE (N)
n_ta_N_results <- map_dfr(n_ta_N_params, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = sim$mean_accuracy,
         "FPR" = sim$mean_FPR,
         "%NA" = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_N_results, "tables/n_ta_output/n_ta_N_params2.xlsx") # N


# TIME SERIES
n_ta_time_results <- map_dfr(n_ta_time_params, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = sim$mean_accuracy,
         "FPR" = sim$mean_FPR,
         "%NA" = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_time_results, "tables/n_ta_output/n_ta_time_params2.xlsx") # time series


# TIME WINDOW (W)
n_ta_tw_results <- map_dfr(n_ta_tw_params, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = sim$mean_accuracy,
         "FPR" = sim$mean_FPR,
         "%NA" = sim$proportionNA,
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_tw_results, "tables/n_ta_output/n_ta_tw_params2.xlsx") # time window size

# N x mu
n_ta_N_x_mu_results <- map_dfr(n_ta_N_x_mu, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_N_x_mu_results, "tables/n_ta_output/int/n_ta_N_x_mu.xlsx")


# N x t
n_ta_N_x_time_results <- map_dfr(n_ta_N_x_time, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_N_x_time_results, "tables/n_ta_output/int/n_ta_N_x_time.xlsx")

# mu x t
n_ta_mu_x_t_results <- map_dfr(n_ta_mu_x_t, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "sd NDR" = round(sim$sd_NDR, 3),
         "FPR" = round(sim$mean_FPR, 3),
         "sd FPR" = round(sim$sd_FPR, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})
write_xlsx(n_ta_mu_x_t_results, "tables/n_ta_output/int/n_ta_mu_x_t.xlsx")


# PRINT RESULTS
print(n_ta_mu_results)
print(n_ta_N_results)
print(n_ta_time_results)
print(n_ta_tw_results)






