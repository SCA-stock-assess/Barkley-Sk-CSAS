# Packages ----------------------------------------------------------------

pkgs <- c("here", "tidyverse", "gsl", "rstan", "ggridges")
#install.packages(pkgs)


# Load packages
library(here)
library(tidyverse); theme_set(theme_bw())
library(ggridges)
library(gsl)
library(rstan)


set.seed(19)



# Load fitted AR1 models  -------------------------------------------------


# Load fitted models for the three CUs
AR1_fits <- here(
  "3. outputs",
  "Stock-recruit modelling"
) |> 
  list.files(
    pattern = "_AR1.rds",
    full.names = TRUE
  ) |> 
  purrr::set_names(~str_extract(.x, ".{3}(?=_AR1.rds)")) |> 
  map(readRDS)


# Functions to calculate benchmarks ---------------------------------------


# Both taken from Dylan Glaser's Fraser Pink analysis:
# https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/functions.R
get_Smsy <- function(a, b){
  Smsy <- (1 - lambert_W0(exp(1 - a))) / b
  if(Smsy <0){Smsy <- 0.001} #dumb hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy) {Sgen * a * exp(-b * Sgen) - Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}


# Wrap benchmark bootstrapping into a function for iteration
get_benchmarks <- function(AR1_model) {
  
  model_pars <- rstan::extract(AR1_model)

  bench <- matrix(
    NA,
    1000,
    4,
    dimnames = list(
      seq(1:1000), 
      c("Sgen","Smsy","Umsy", "Seq")
    )
  )
  
  for(i in 1:1000){
    r <- sample(seq(1,1000),1,replace=TRUE)
    ln_a <- model_pars$lnalpha[r]
    b <- model_pars$beta[r]
    
    bench[i,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[i,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[i,2]) #S_gen
    bench[i,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[i,4] <- ln_a/b #S_eq
  }
  
  
  return(bench)
}


# Bootstrap reference points from AR1 parameters --------------------------


# Iterate over CUs
ref_pts <- AR1_fits |> 
  map(get_benchmarks) |> 
  map(as_tibble) |> 
  list_rbind(names_to = "stock") |> 
  mutate(Smsy_0.8 = 0.8*Smsy) |> 
  pivot_longer(
    !stock,
    names_to = "ref_pt",
    values_to = "value"
  ) |> 
  mutate(
    stock = factor(
      stock, 
      levels = c("GCL", "SPR", "HED"),
      labels = c("GCL", "SPR", "HUC")
    )
  )


# Median values for each stock and reference point
ref_pts_median <- ref_pts |> 
  summarize(
    .by = c(stock, ref_pt),
    median = median(value)
  ) |> 
  mutate(median = if_else(ref_pt == "Umsy", round(median, 2), round(median, 0)))


# Plot distributions of bootstrapped reference points
(ref_pts_plot <- ref_pts |> 
  filter(!ref_pt == "Smsy") |> 
  ggplot(aes(x = value, y = fct_rev(stock))) +
  facet_wrap(
    ~ ref_pt, 
    scales = "free_x"
  ) +
  geom_density_ridges(
    fill = "grey",
    quantile_lines = TRUE,
    quantiles = 2,
    vline_colour = "red",
    scale = 1.2,
    colour = alpha("black", 0.5),
    alpha = 0.75
  ) +
  geom_text(
    data = filter(ref_pts_median, !ref_pt == "Smsy"),
    aes(
      label = median,
      x = median
    ),
    hjust = -0.2, 
    vjust = -1.5
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 3),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .7))) + 
  labs(y = "Conservation Unit")
)


# Save plot
ggsave(
  ref_pts_plot,
  filename = here(
    "3. outputs",
    "Plots",
    "Barkley-Sk_bootstrapped_ref_pts.png"
  ),
  width = 6.5,
  units = "in",
  dpi = "print"
)


# Build summary dataframe of reference points to export -------------------


# Summary stats for each reference point
ref_pts_summary <- ref_pts |> 
  summarize(
    .by = c(stock, ref_pt),
    quant = list(quantile(value, c(0.05, 0.25, 0.5, 0.75, 0.95))),
    mean = mean(value),
    sd = sd(value)
  ) |> 
  unnest_wider(quant) |> 
  rename("median" = `50%`) |> 
  rename_with(
    .cols = matches("\\d"),
    \(x) str_remove_all(paste0("q", x), "[:punct:]")
  )


# Export reference points
write.csv(
  ref_pts_summary,
  file = here(
    "3. outputs",
    "Stock-recruit modelling",
    "Barkley-Sk_CU_ref_pts_summary.csv"
  ),
  row.names = FALSE
)

