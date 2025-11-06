
# Beyond the Signals: An evaluation of the Signal Selection Test using simulated data

## Outline

- [Overview](#Overview)
- [Quick start](#Quick-start)
- [Methods](#Methods)
- [Reproduce experiments](#Reproduce-experiments)
- [Further remarks](#Further-remarks)
- [Final remark](#Final-remark)
- [Contributing](#Contributing)
- [Acknowledgments](#Acknowledgments)

## Overview

This repo includes all scripts, figures and tables created to carry out
the analysis of my dissertation, submitted for the MPhil in Human
Evolutionary Studies at the University of Cambridge. The aim of this
project is to critically assess the performance of the core statistical
component of the
[`signatselect`](https://github.com/benmarwick/signatselect), the
**Frequency Increase Test (FIT)**. To have access to the dissertation,
refer to the following email (creator’s email):
`albertok.cooper@gmail.com`.

Neutrality tests, such as the Ewens-Watterson and [Slatkin’s exact
test](https://github.com/mmadsen/slatkin-exact-tools), have yielded
insightful results. Nonetheless, the nature of the archaeological record
and cultural evolution raises critical questions about the inferences
from these tests: Does departure from neutrality always mean selection?
Can cultural transmission biases equate to selection? How reliable are
these tests given the poor time resolution of archaeological
assemblages?

We address this gap by evaluating the `signatselect` (SST), a
statistical method designed to detect selection in genetic time series
data [Feder et al. (2014)](https://doi.org/10.1534/genetics.113.158220).
It is a computationally efficient technique, employing a one-sample
*t*-test (the Frequency Increase Test, FIT) to compare the rescaled
allele frequency differences between two successive intervals against
theoretical expectations under neutrality. It requires at least three
time points (the **three-time-point** rule), and it flags selection at
$\alpha$ = 0.05.

Its efficiency for archaeological time series, shaped by evolutionary
and depositional processes, remains untested. Thus, our research has
three objectives:

- To evaluate neutral detection and Type I error (false positives) under
  unbiased transmission.
- To appraise its statistical power (SSR) when variants undergo
  selection via two cultural transmission biases: content and conformist
  bias.
- To critically gauge the confounding effect of temporal structure of
  assemblages (i.e., time averaging) on the test’s performance.

## Quick start

### Installation

``` r
# Install required packages
install.packages(c("dplyr","ggplot2","tidyr","gridExtra",
                   "purrr","tibble","writexl","pak","tidyverse",
                   "viridis","stringr","viridisLite",
                   "readxl","patchwork","scales"))

# Install signatselect
pak::pkg_install("benmarwick/signatselect")
```

``` r
# Read packages
pkgs <- c(
  "dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl", "tidyverse",
  "viridis", "stringr", "viridisLite",
  "readxl", "patchwork", "scales"
)
invisible(lapply(pkgs, library, character.only = TRUE))

# Read signatselect
library(signatselect)
```

## Methods

We conduct two simulation-based experiments modelling cultural
transmission. Both use base parameters: population size (`N`),
innovation rate (`μ`), a burn-in period (`B`), discrete time units
(`t`), and time window size (`w`) for averaging. Thus, we simulate three
transmission models:

1.  **Unbiased transmission**: naive individuals randomly adopt a
    variant from the population.
2.  **Content-biased transmission**: learners are more likely to adopt a
    variant with a higher “attractiveness” product (1 + `b`).
3.  **Frequency-biased transmission**: learners are more likely to adopt
    a variant the more common it becomes, modeled with a conformist
    exponent (1 + `c`).

Each experiment is paired with two temporal deposition schemes: a
**snapshot** regime (fine-grained ideal deposition)

<img src="man/figures/README-snapshot-1.png" width="100%" style="display: block; margin: auto;" />

and **time-averaging** regime (coarse-grained, palimpsest deposition).

<img src="man/figures/README-time-averaging-1.png" width="100%" style="display: block; margin: auto;" />

## Reproduce experiments

In order to reproduce the results reported in the dissertation the user
should refer to two types of
[scripts](https://github.com/Cooperr18/Dissertation-Code/tree/master/scripts)
and the following workflow:

1.  Script storing the functions defining the models and baseline
    parameters. Firstly, define the function of each model, and run the
    baseline levels at the end of the same script. E.g.:

``` r
# PARAMETERS -----------------------------------------------------------------
# N = Population size
# mu = innovation rate (between 0 and 1, inclusive)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs

set.seed(1234)

# PIPELINE -------------------------------------------------------------------

neutral_snapshot <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs) {
  
  accuracy_snapshot <- numeric(n_runs) # empty vector for accuracy tracking each run
  FPR_snapshot <- numeric(n_runs)
  fit_p_count <- vector("list", n_runs) # store p-values
  
  for (run in 1:n_runs) {
    ini <- 1:N # Initial cultural variants
    traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
    pop <- ini # initial population of variants to pop
    maxtrait <- N 
    
    # Burn-in stage
    for(i in 1:burnin) {
      # neutral transmission:
      pop <- sample(pop,replace=T)
      
      # Add innovations
      innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
    }
    
    # Observation period after equilibrium
    for (i in 1:timesteps) {
      pop <- sample(pop,replace=T)
      innovate <- which(runif(N)<mu) 
      if(length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      
      traitmatrix[i,] <- pop #record the variants of each individual
    }
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to frequency matrix
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- unique_variants # give names
    
    # Prepare FIT input
    freq_long <- as.data.frame(freq_mat) %>%
      mutate(time = 1:timesteps) %>%
      pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
      filter(freq > 0) %>% # remove zeros
      mutate(variant = as.integer(variant))
    
    # Filter data to 3 time points
    freq_long_filtered <- freq_long %>%
      group_by(variant) %>%
      filter(n_distinct(time) >= 3) %>%
      ungroup()
    
    # FIT application and storing
    fit_results <- freq_long_filtered %>%
      group_split(variant) %>%
      map_dfr(~ {
        df <- as.data.frame(.x)  # .x is a dataframe
        
        # check if we have enough data points
        if (nrow(df) < 3) {
          return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
        }
        
        # safely apply FIT
        res <- tryCatch(
          fit(time = df$time, v = df$freq),
          error = function(e) list(fit_p = NA)  # handle errors
        )
        
        data.frame(
          variant = df$variant[1],
          time_points = nrow(df),
          fit_p = res$fit_p,
          stringsAsFactors = FALSE
        )
      }) %>%
      mutate(
        sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"), # neutral if >0.05
        sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
      )
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    # Store metrics for this run
    total_variants <- sum(fit_results$sig != "NA")
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_snapshot[run] <- NDR  # track for each run
    FPR_snapshot[run] <- FPR
  }
  
  # Store results
  # MEAN
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_snapshot >= 0.95) / n_runs * 100
  mean_FPR <- mean(FPR_snapshot, na.rm = TRUE)
  
  # SD
  sd_NDR <- sd(accuracy_snapshot, na.rm = T) # sd across runs
  sd_FPR <- sd(FPR_snapshot, na.rm = T)
  
  all_pvals <- unlist(fit_p_count)
  
  # Proportion of NA from the simulation
  sumNA <- sum(is.na(all_pvals))
  proportionNA <- sumNA / length(all_pvals) * 100
  if(sumNA == 0) {
    proportionNA <- 0
  }
  
  # Store output
  return(list(
    
    # PARAMETERS
    N = N,
    mu = mu,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    
    # OUTPUTS
    accuracy_snapshot = accuracy_snapshot,
    mean_accuracy = mean_accuracy,
    mean_FPR = mean_FPR,
    sd_NDR = sd_NDR,
    sd_FPR = sd_FPR,
    high_accuracy_runs = high_accuracy_runs,
    sumNA = sumNA,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    all_pvals = all_pvals
    )
  )
}

# BASELINE SIMULATION
n_snap_sim <- neutral_snapshot(N = 100, mu = 0.01, burnin = 1000,
                              timesteps = 1000, p_value_lvl = 0.05, n_runs = 100)

# Store output across runs in a table
results_table_neutral_snapshot <- tibble(
  N  = n_snap_sim$N,
  "µ" = n_snap_sim$mu,
  "B" = n_snap_sim$burnin,
  "Time steps" = n_snap_sim$timesteps,
  "α" = n_snap_sim$p_value_lvl,
  NDR = n_snap_sim$mean_accuracy,
  FPR = n_snap_sim$mean_FPR,
  ">95% runs" = n_snap_sim$high_accuracy_runs,
  "%NA" = n_snap_sim$proportionNA,
  mean_p_value = mean(n_snap_sim$all_pvals, na.rm = TRUE),
  "Runs" = n_snap_sim$n_runs
)
results_table_neutral_snapshot
```

2.  Script including the lists of parameter combinations, and the
    specific function to run each set of parameters. E.g.:

``` r
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

# Run and store

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
```

And then plot the results (for a more detailed description of figure
reproducibility, see [Further remarks](#Further-remarks)):

``` r
# import sweep
n_snap_mu_df <- read_excel("tables/n_snap_output/n_snap_mu_params2.xlsx")
df <- n_snap_mu_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = `mu`)) +
  geom_col(aes(y = `%NA`), fill = "grey80", width = 0.0125) +
  geom_line(aes(y = NDR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = expression(mu), breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )
```

<img src="man/figures/README-neutral snap mu-1.png" width="100%" style="display: block; margin: auto;" />

### Experiment 1

In Experiment 1, we simulate neutral conditions under both deposition
schemes,
[`Unbiased_Snapshot.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Unbiased_Snapshot.R)
and
[`Unbiased_Time_Averaging.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Unbiased_Time_Averaged.R).
At equilibrium, we record variants and run the FIT over multiple
iterations. We then calculate the Neutral Detection Rate or `NDR` (True
Negatives) and the Type I Error Rate (or False Positive Rate, `FPR`),
and run through different parameter combinations, defined in
[`Experiment_1.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Experiment_1.R).

### Experiment 2

In Experiment 2, we introduce one focal variant with a content or
conformist bias at equilibrium, as put by the models in
[`Content_Bias_Snapshot.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Content_Bias_Snapshot.R),
[`Content_Bias_Time_Averaging.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Content_Bias_Time_Averaging.R),
[`Conformist_Bias_Snapshot.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Conformist_Bias_Snapshot.R),
[`Conformist_Bias_Time_Averaging.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Conformist_Bias_Time_Averaging.R).
We then follow the same procedure as before
([`Experiment_2.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Experiment_2.R)),
and compute the Signal of Selection Rate (statistical power, `SSR`) and
the Type II Error Rate (or False Negative Rate, `FNR`).

Both experiments are run with a base parameter set, followed by
parameter sweeps to observe the confounding effect of each parameter to
the test’s ability to flag selection. This approach allows us to capture
the interaction between cultural evolution and assemblage formation,
helping us pin down potential confounding agents in the archaeological
record.

## Further remarks

To provide full transparency of the decision-making process along the
realization of this work, refer to the
[`Report of daily progress.docx`](https://github.com/Cooperr18/Dissertation-Code/blob/master/Report%20of%20daily%20progress.docx).
It is a diary that shows full accountability for both included and not
included decisions in the final project. For instance,
[`Alt_Content_Unused.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/Alt_Content_Unused.R)
stores alternative approaches to modeling content bias.

Furthermore, all listed figures can be found at
[`model_sims_plots.R`](https://github.com/Cooperr18/Dissertation-Code/blob/master/scripts/model_sims_plots.R).
It is a script that breaks down the code for each numbered figure in the
dissertation. Again, to fully reproduce them the workflow should follow
[Reproduce experiments](#Reproduce-experiments): store the model
functions, the parameter combinations and run the functions for each
set; store the output, and import using `read_excel()`. In some cases,
referral to other scripts is not mandatory, as the functions are often
included within this script. E.g.:

``` r
neutral_snapshot <- function(N, mu, timesteps, p_value_lvl, n_runs) {
  runs <- vector("list", length = n_runs)
  
  for (run in 1:n_runs) {
    # Initialize population: 33% start with variant 1, others have unique variants
    initial_variant_count <- round(N * 0.33)
    pop <- c(rep(1, initial_variant_count), 2:(N - initial_variant_count + 1))
    
    traitmatrix <- matrix(NA, nrow = timesteps, ncol = N)
    maxtrait <- max(pop)
    
    for (i in 1:timesteps) {
      pop <- sample(pop, replace = TRUE)  # Neutral drift
      innovate <- which(runif(N) < mu)
      if (length(innovate) > 0) {
        new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
        pop[innovate] <- new_variants
        maxtrait <- max(pop)
      }
      traitmatrix[i, ] <- pop
    }
    
    # Track ONLY variant 1 (initial frequency = 0.33)
    freq_variant1 <- rowMeans(traitmatrix == 1, na.rm = TRUE)
    
    runs[[run]] <- data.frame(
      time = 1:timesteps,
      freq = freq_variant1,
      run = run,
      N = N,
      mu = mu,
      p_value_lvl = p_value_lvl,
      variant = "Variant 1 (Initial freq = 0.33)"
    )
  }
  bind_rows(runs)
}

n_snap_N_params <- list(
  list(N=10, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=50, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=100, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=250, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=500, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4),
  list(N=1000, mu=0,timesteps=200, p_value_lvl=0.05, n_runs=4)
)

# Run all simulations and bind into one big tibble
all_results <- map_dfr(n_snap_N_params, ~
                         neutral_snapshot(.x$N, .x$mu, .x$timesteps, 
                                          .x$p_value_lvl, .x$n_runs))

ggplot(all_results, aes(x = time, y = freq, color = factor(run))) +
  geom_line(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0.33, linetype = "dashed", color = "gray30", size = 0.8) +
  facet_wrap(~ N, scales = "fixed", ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Run") +
  labs(
    x = "Time Step",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 20) +  # Base font size increased to 20
  theme(
    axis.title.x = element_text(size = 22, margin = margin(t = 10)),  # Larger X-axis label
    axis.title.y = element_text(size = 22, margin = margin(r = 10)),  # Larger Y-axis label
    axis.text = element_text(size = 18),  # Larger tick labels
    strip.text = element_text(size = 20, face = "bold"),  # Larger facet labels
    legend.title = element_text(size = 20),  # Larger legend title
    legend.text = element_text(size = 18),  # Larger legend labels
    legend.position = "bottom",
    panel.spacing = unit(1.5, "lines")  # Extra space between facets
  )
```

<img src="man/figures/README-unbiased transmission-1.png" width="100%" style="display: block; margin: auto;" />

But in other figures the user will need to import the data from
spreadsheets available at
[tables](https://github.com/Cooperr18/Dissertation-Code/tree/master/tables),
with sub-folders for each one of the models:

``` r
# FIGURE 15 --------------------------------------------------------------------
n_snap_N_df <- read_excel("tables/n_snap_output/n_snap_N_params2.xlsx")
df <- n_snap_N_df
df$`%NA` <- df$`%NA`/100

# NDR
ggplot(df, aes(x = N)) +
  geom_col(aes(y = `%NA`), fill = "grey80") +
  geom_line(aes(y = NDR),   color = plasma(5)[4], linewidth = 1.2) +
  geom_point(aes(y = NDR),  color = plasma(5)[4], size = 3) +
  scale_x_continuous(name = "N", breaks = pretty_breaks(8)) +
  scale_y_continuous(
    name     = "NA",
    sec.axis = sec_axis(~ ., name = "NDR")) +
  labs(
    caption = paste0(
      "Runs = ", 1000, "    |    ",
      "Grey = NA, Orange = NDR")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y.left  = element_text(size = 20),
    axis.title.y.right = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 20, hjust = 0)
  )
```

<img src="man/figures/README-fig 15-1.png" width="100%" style="display: block; margin: auto;" />

## Final remark

This investigation is driven by the lack of systematic research of the
`signatselect` in archaeology. There are few sources which explore the
scope of this technique: the original publication, [Feder et
al. (2014)](https://doi.org/10.1534/genetics.113.158220); an application
to language time series data [Newberry et
al. (2017)](https://doi.org/10.1038/nature24455); a brief overview of
softwares to detect selection [Vlachos et
al. (2020)](https://doi.org/10.1186/s13059-019-1770-8), and an
evaluation of Newberry et al.’s results plus the effect of time binning
on the test [Karjus et al. (2019)](https://doi.org/10.5334/gjgl.909).
However, the only source which explicitly applies and evaluates the
`signatselect` with archaeological data is a public repo posted by Ben
Marwick, Hezekiah A. Bacovcin and Sergey Kryazhimskiy, at
<https://github.com/benmarwick/signatselect>. This source has served of
great inspiration to this study, thus I am truly thankful to the
authors.

## Contributing

We welcome contributions. Please:

1.  Fork the repository
2.  Create a feature branch (`git checkout -b amazing-feature`)
3.  Commit your changes (`git commit -m "Add amazing feature"`)
4.  Push to the branch (`git push origin amazing-feature`)
5.  Open a Pull Request

## Acknowledgments

I am grateful to my supervisor, [Dr. Enrico
Crema](https://ercrema.github.io/), for his unwavering encouragement,
combined with his superb coding expertise and keen insights into the
literature, all of which have profoundly benefited my work. In addition,
I appreciate his patient guidance in helping me improve my coding and
analytical skills, which have grown significantly under his mentorship.
Moreover, I am thankful for the guidance and inspiration provided by
[Dr. Daniel García
Rivero](https://investigacion.us.es/sisius/sis_showpub.php?idpers=8645),
for I owe him much of my interest in evolutionary archaeology.

I also thank Ben Marwick, Hezekiah Akiva Bacovcin and Sergey
Kryazhimskiy for their contributions in the application of the
[`signatselect`](https://github.com/benmarwick/signatselect) to
archaeological data,which served as a great initial push for this work.
Finally, to the HES and BAS cohort, and to Justin and Jose, for easing
this sometimes rather steep climb and putting up with my complaints.
