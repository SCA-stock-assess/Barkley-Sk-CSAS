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
  map(readRDS)


# Use model_name attribute to assign names to the list of models
names(AR1_fits) <- map(AR1_fits, \(x) attr(x, "model_name")) |> list_c()


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


# Bootstrap reference points from AR1 parameters --------------------------


# Add AR1 model alpha and beta posterior draws to a dataframe
AR1_frame <- AR1_fits |> 
  enframe(name = "spec", value = "model") |> 
  rowwise() |> 
  mutate(parameter = list(str_subset(names(model), "lnalpha|beta"))) |> 
  unnest_longer(parameter) |> 
  rowwise() |> 
  mutate(
    posterior = list(as.data.frame(rstan::extract(model, pars = parameter))),
    stock = str_extract(spec, "GCL|SPR|HUC"),
    fert = case_when(
      str_detect(parameter, "_fert") ~ 1, 
      str_detect(parameter, "_unfert") ~ 0,
      .default = NA
    ),
    data_scope = if_else(str_detect(spec, "trim"), "trim", "full"),
    parameter = str_remove_all(parameter, "_.*")
  ) |>
  # Ensure '_fert' and '_unfert' suffixes are excluded from dataframe column names
  mutate(posterior = list(rename_with(.data = posterior, ~str_remove_all(.x, "_.*")))) |> 
  select(-model) |> 
  pivot_wider(names_from = parameter, values_from = posterior) |> 
  rowwise() |> 
  mutate(
    parameters = list(cbind(lnalpha, beta)),
    .keep = "unused"
  ) |> 
  ungroup()


# Use posterior draws to estimate reference points
ref_pts <- AR1_frame |> 
  rowwise() |> 
  # Take 5000 random samples from the posterior draws
  mutate(parameters = list(slice_sample(parameters, n = 5000, replace = TRUE))) |> 
  unnest(parameters) |> 
  rowwise() |> 
  # Calculate reference points from the 5000 random samples
  mutate(
    Smsy = get_Smsy(lnalpha, beta),
    Smsy0.8 = 0.8*Smsy,
    Sgen = get_Sgen(exp(lnalpha), beta, -1, 1/beta*2, Smsy),
    Umsy = (1 - lambert_W0(exp(1 - lnalpha))),
    Seq = lnalpha/beta
  ) |> 
  ungroup()


# Summary values for each model and reference point
ref_pts_summary <- ref_pts |> 
  summarize(
    .by = !c(lnalpha, beta, Smsy, Smsy0.8, Sgen, Umsy, Seq),
    across(
      c(lnalpha, beta, Smsy, Smsy0.8, Sgen, Umsy, Seq),
      .fns = c(
        `5` = ~quantile(.x, 0.05),
        `50` = ~quantile(.x, 0.5),
        `95` = ~quantile(.x, 0.95)
      ),
      .names = "{.col}_{.fn}"
    )
  ) |> 
  pivot_longer(
    cols = matches("_\\d+"),
    names_sep = "_",
    names_to = c("ref_pt", "quantile"),
    values_to = "value"
  ) |> 
  pivot_wider(names_from = quantile, names_prefix = "q") |>
  mutate(
    long_name = factor(
      spec,
      levels = c(
        "GCL", 
        "SPR", 
        "HUC_full_nofert", 
        "HUC_trim_nofert", 
        "HUC_full_fert", 
        "HUC_trim_fert"
      ),
      labels = c(
        "Great Central Lake",
        "Sproat Lake",
        "Hucuktlis Lake (all data)",
        "Hucuktlis Lake (excl. 1993)",
        "Hucuktlis Lake (w/fertilization)",
        "Hucuktlis Lake (w/fertilization; excl. 1993)"
      ) 
    ),
    stock = factor(
      stock,
      levels = c("GCL", "SPR", "HUC"),
      labels = c("Great Central", "Sproat", "Hucuktlis")
    )
  ) 


# Plot distributions of bootstrapped reference points
(ref_pts_plot <- ref_pts_summary |> 
    mutate(
      label = case_when(
        str_detect(ref_pt, "^S") ~ as.character(round(q50, 0)),
        ref_pt == "lnalpha" ~ as.character(round(q50, 2)),
        ref_pt == "Umsy" ~ paste0(round(q50*100, 0), "%"),
        ref_pt == "beta" ~ as.character(format(q50, digits = 2, scientific = TRUE))
      )
    ) |>
    ggplot(aes(x = q50, y = fct_rev(long_name))) +
    facet_wrap(
      ~ ref_pt, 
      strip.position = "bottom",
      scales = "free_x"
    ) +
    geom_pointrange(
      aes(
        colour = factor(fert),
        xmin = q5, 
        xmax = q95
      ),
      position = position_dodge(width = 0.7),
      orientation = "y"
    ) +
    geom_text(
      aes(
        label = label,
        colour = factor(fert)
      ),
      hjust = -0.2, 
      vjust = -0.3,
      position = position_dodge(width = 0.7),
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 3),
      expand = c(0, 0)
    ) +
    scale_colour_manual(values = c("red", "blue")) +
    guides(
      colour = guide_legend(
        title = "Fertilization state",
        theme = theme(
          legend.direction = "horizontal",
          legend.title.position = "top",
          legend.text.position = "bottom"
        )
      )
    ) +
    theme(
      axis.title = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing.y = unit(1, "lines"),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0),
      legend.justification.inside = c(1, 0),
      legend.background = element_rect(fill = alpha("grey", 0.25))
    )
)
# Note that Umsy values look low but Umsy is a function only of 
# alpha, which reflects productivity at low spawner abundances but not
# overall carrying capacity


# Save plot
ggsave(
  ref_pts_plot,
  filename = here(
    "3. outputs",
    "Plots",
    "Barkley-Sk_bootstrapped_ref_pts.png"
  ),
  width = 10,
  height = 10,
  units = "in",
  dpi = "print"
)


# Export summary dataframe of reference points -------------------


# Export reference points
ref_pts_summary |> 
  arrange(long_name, desc(ref_pt), data_scope, fert) |> 
  write.csv(
    file = here(
      "3. outputs",
      "Stock-recruit modelling",
      "Barkley-Sk_CU_ref_pts_summary.csv"
    ),
    row.names = FALSE
  )

