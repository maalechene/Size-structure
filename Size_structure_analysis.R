library(tidyr)
library(dplyr)
library(ggplot2)
library(effects)   #
library(ggeffects) #
library(emmeans)   #
library(tidyverse) #
library(rstanarm)   #
library(brms)       #
library(ggmcmc)     #
library(DHARMa)     #
library(rstan)      # version 2.26.22
library(tidybayes)  #
library(HDInterval) #
library(readr)
library(posterior)  #
library(patchwork)  #
library(reshape2)
library(coda)
library(RColorBrewer)
library(cowplot)


## DHARMa like function to check simulated residuals

make_brms_dharma_res <- function(brms_model, seed = 10, ...) {
  # equivalent to `simulateResiduals(lme4_model, use.u = FALSE)`
  # cores are set to 1 just to ensure reproducibility
  options(mc.cores = 1)
  on.exit(options(mc.cores = parallel::detectCores()))
  response <- brms::standata(brms_model)$Y
  ndraws <- nrow(as_draws_df(brms_model))
  manual_preds_brms <- matrix(0, ndraws, nrow(brms_model$data))
  random_terms <- insight::find_random(
    brms_model, split_nested = TRUE, flatten = TRUE
  )
  # for this to have a similar output to `glmmTMB`'s default, we need to
  #   create new levels in the hierarchical variables, so then we can
  #   use `allow_new_levels = TRUE` and `sample_new_levels = "gaussian"` in
  #   `brms::posterior_epred`. This is equivalent to
  #   `simulateResiduals(lme4_model, use.u = FALSE)`. See details in
  #   `lme4:::simulate.merMod` and `glmmTMB:::simulate.glmmTMB`
  new_data <- brms_model$data |>
    dplyr::mutate(across(
      all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
    ))
  set.seed(seed)
  brms_sims <- brms::posterior_predict(
    brms_model, re_formula = NULL, newdata = new_data,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  ) |>
    t()
  fitted_median_brms <- apply(brms_sims, 1, median)
  ## fitted_median_brms <- apply(
  ##     t(brms::posterior_epred(brms_model, ndraws = ndraws, re.form = NA)),
  ##     1,
  ##     mean)
  DHARMa::createDHARMa(
    simulatedResponse = brms_sims,
    observedResponse = response,
    fittedPredictedResponse = fitted_median_brms,
    ...
  )
}


plot0.1 <- theme(axis.title.x = 
                   element_text(
                     margin = margin(t = 10),
                     color = "Black", 
                     size = 15), 
                 axis.title.y = 
                   element_text(
                     margin = margin(r = 20),
                     color = "Black", 
                     size = 15), 
                 axis.text = 
                   element_text(
                     color = "Black",
                     size = 12), 
                 axis.line = 
                   element_line(
                     color = "Black",
                     linewidth = 0.5), 
                 axis.ticks = element_line(color = "Black"),
                 panel.background = element_rect(fill = "grey98"),
                 panel.grid.minor = element_line(color = "White"),
                 plot.title = 
                   element_text(
                     size = 15,
                     face = "bold",
                     hjust = 0.5),
                 plot.subtitle = 
                   element_text(
                     size = 15,
                     hjust = 0.5),
                 strip.text.x = 
                   element_text(
                     size = 15,
                     color = "black", 
                     face = "bold"),
                 strip.background = 
                   element_rect(
                     color = "white",
                     fill = "white", 
                     linewidth = 1, 
                     linetype = "solid"),
                 legend.background = element_blank(), 
                 legend.title = 
                   element_text(
                     colour = "black",
                     size = 15, 
                     face = "bold"), 
                 legend.text = 
                   element_text(
                     colour = "black",
                     size = 12,
                     face = "bold.italic"),
                 legend.key.width = unit(1.5, "lines"), 
                 legend.key = element_rect(fill = NA))

# Same theme without a legend:
plot0 <- plot0.1 + theme(legend.position = "none")



###############################################################################
###################           CLUSTER ANALYSIS           ######################
###############################################################################

# Import data
clu_data <- read_csv("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Data/CSV/cluster_data.csv")

clu_data <- clu_data |> 
  mutate(Cluster = factor(Cluster),
         Class = factor(Class))

# Model fitting
clu_m2 <- brm(bf(log10(Area_cm2) ~ Class*Cluster + (1 | Site/Plot),
                 quantile = 0.95),
              family = asym_laplace(),
              data = clu_data,
              cores = 3,
              chains = 3,
              warmup = 5000,
              thin = 12,
              iter = 20000,
              control = list(max_treedepth = 12))

clu_m3 <- brm(bf(log10(Area_cm2) ~ Class*Cluster + (1 | Site/Plot),
                 quantile = 0.5),
              family = asym_laplace(),
              data = clu_data,
              cores = 3,
              chains = 3,
              warmup = 5000,
              thin = 10,
              iter = 20000,
              control = list(max_treedepth = 12))

clu_m4 <- brm(bf(log10(Area_cm2) ~ Class*Cluster + (1 | Site/Plot),
                 quantile = 0.1),
              family = asym_laplace(),
              data = clu_data,
              cores = 3,
              chains = 3,
              warmup = 5000,
              thin = 12,
              iter = 20000,
              control = list(max_treedepth = 12))

# Models diagnostics and validation

clu_m2$fit |> stan_trace(inc_warmup = TRUE)
clu_m2$fit |> stan_ac()
clu_m2$fit |> stan_rhat()
clu_m2$fit |> stan_ess()
clu_m2 |> pp_check(type = "dens_overlay", ndraws = 100)
clu_m2 |> conditional_effects()
resids <- make_brms_dharma_res(clu_m2, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(clu_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

clu_m3$fit |> stan_trace(inc_warmup = TRUE)
clu_m3$fit |> stan_ac()
clu_m3$fit |> stan_rhat()
clu_m3$fit |> stan_ess()
clu_m3 |> pp_check(type = "dens_overlay", ndraws = 100)
clu_m3 |> conditional_effects()
resids <- make_brms_dharma_res(clu_m3, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(clu_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

clu_m4$fit |> stan_trace(inc_warmup = TRUE)
clu_m4$fit |> stan_ac()
clu_m4$fit |> stan_rhat()
clu_m4$fit |> stan_ess()
clu_m4 |> pp_check(type = "dens_overlay", ndraws = 100)
clu_m4 |> conditional_effects()
resids <- make_brms_dharma_res(clu_m4, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(clu_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

# Extract posterior distributions

load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/clu_m2.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/clu_m3.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/clu_m4.RData")

post10_c <- clu_m4 |> emmeans(~ Cluster|Class) |> gather_emmeans_draws()
post50_c <- clu_m3 |> emmeans(~ Cluster|Class) |> gather_emmeans_draws()
post95_c <- clu_m2 |> emmeans(~ Cluster|Class) |> gather_emmeans_draws()
post10_c$Size = "smallest"
post50_c$Size = "median"
post95_c$Size = "largest"

d_all = rbind(post10_c, post50_c)
d_all2 = rbind(d_all, post95_c)

#Re-order levels based on expectation that coral size increases with latitude and is larger inshore
d_all2 <- d_all2 |> 
  mutate(Cluster = factor(Cluster, levels = c("Torres Strait", "Offshore North", "Offshore Central", "Inshore Central", "Inshore South", "Offshore South")),
         Class = factor(Class),
         Size = factor(Size, levels = c("smallest", "median", "largest")))

library(RColorBrewer)

# Manually select the more saturated colors from the "Blues" palette
blues_palette <- brewer.pal(6, "Blues")[4:9]  # Select the darker, more saturated colors

clu_A <- ggplot(d_all2, aes(x = .value, y = Cluster, fill = Size)) +
  stat_slab(aes(fill = Size), alpha = 0.6) +
  stat_pointinterval(
    aes(color = Size), 
    point_interval = "median_hdi", # Use HPDI for the intervals
    .width = c(0.66, 0.95),
    position = position_dodge(width = .4, preserve = "single")
  ) +
  scale_fill_manual(values = blues_palette) +
  scale_color_manual(values = blues_palette) +
  scale_x_continuous(
    breaks = log10(c(seq(10, 90, 10), seq(100, 900, 100), seq(1000, 9000, 1000), seq(10000, 30000, 10000))),
    labels = c(10, "", "", "", "", "", "", "", "",
               100, "", "", "", "", "", "", "", "",
               1000, "", "", "", "", "", "", "", "",
               10000, "", 30000),
    limits = c(log10(10), log10(30000))
  ) +
  labs(title = "A",
       x = expression(Colony~size~(cm^2)),
       y = "",
       fill = "Confidence interval") +
  facet_grid(~Class) +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16)
  ) +
  background_grid(major = "x", color.major = "grey90") +
  scale_y_discrete(limits = rev(levels(d_all2$Cluster)))

# Planned contrasts between clusters

## Pairwise contrasts based on expectation that larger corals will be found in 
## higher latitude order of Cluster 
## Default order of cluster in clu_data is alphabetical

cluster_mat <- list('Inshore South vs. Offshore South' = c(0, 1, 0, 0, -1, 0),
                    'Inshore South vs. Inshore Central' = c(-1, 1, 0, 0, 0, 0),
                    'Inshore South vs. Offshore Central' = c(0, 1, -1, 0, 0, 0),
                    'Inshore South vs. Offshore North' = c(0, 1, 0, -1, 0, 0),
                    'Inshore South vs. Torres Strait' = c(0, 1, 0, 0, 0, -1),
                    'Offshore South vs. Inshore Central' = c(-1, 0, 0, 0, 1, 0),
                    'Offshore South vs. Offshore Central' = c(0, 0, -1, 0, 1, 0),
                    'Offshore South vs. Offshore North' = c(0, 0, 0, -1, 1, 0),
                    'Offshore South vs. Torres Strait' = c(0, 0, 0, 0, 1, -1),
                    'Inshore Central vs. Offshore Central' = c(1, 0, -1, 0, 0, 0),
                    'Inshore Central vs. Offshore North' = c(1, 0, 0, 0, -1, 0),
                    'Inshore Central vs. Torres Strait' = c(1, 0, 0, 0, 0, -1),
                    'Offshore Central vs. Offshore North' = c(0, 0, 1, -1, 0, 0),
                    'Offshore Central vs. Torres Strait' = c(0, 0, 1, 0, 0, -1),
                    'Offshore North vs. Torres Strait' = c(0, 0, 0, 1, 0, -1))

contrast_50_c <- emmeans(clu_m3, ~ Cluster|Class) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_50_df <- as.data.frame(contrast_50_c)

contrast_50_df <- contrast_50_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_50_df <- contrast_50_df |> 
  mutate(Size = "median")

contrast_10_c <- emmeans(clu_m4, ~ Cluster|Class) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_10_df <- as.data.frame(contrast_10_c)

contrast_10_df <- contrast_10_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_df <- contrast_10_df |> 
  mutate(Size = "smallest")


contrast_95_c <- emmeans(clu_m2, ~ Cluster|Class) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_95_df <- as.data.frame(contrast_95_c)

contrast_95_df <- contrast_95_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_df <- contrast_95_df |> 
  mutate(Size = "largest")

contrast_Clu_df <- rbind(contrast_10_df, contrast_50_df, contrast_95_df)
contrast_Clu_df <- contrast_Clu_df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Clu_df <- contrast_Clu_df |> 
  mutate(Sig = factor(Sig),
         Class = factor(Class, levels = c("Atab", "Amil", "Pdam", "Pmas")),
         contrast = factor(contrast, levels = c(
           'Inshore South vs. Offshore South',
           'Inshore South vs. Inshore Central',
           'Inshore South vs. Offshore Central',
           'Inshore South vs. Offshore North',
           'Inshore South vs. Torres Strait',
           'Offshore South vs. Inshore Central',
           'Offshore South vs. Offshore Central',
           'Offshore South vs. Offshore North',
           'Offshore South vs. Torres Strait',
           'Inshore Central vs. Offshore Central',
           'Inshore Central vs. Offshore North',
           'Inshore Central vs. Torres Strait',
           'Offshore Central vs. Offshore North',
           'Offshore Central vs. Torres Strait',
           'Offshore North vs. Torres Strait'
         )))

##Aggregate the data by contrast and Class
agg_contrast_df <- contrast_Clu_df %>%
  group_by(contrast, Class) %>%
  summarise(
    positive_count = sum(Sig == "positive difference"),
    negative_count = sum(Sig == "negative difference"),
    total_sig = positive_count + negative_count,
    .groups = 'drop'
  ) %>%
  mutate(
    dominant_sig = case_when(
      positive_count > 0 & negative_count == 0 ~ "positive difference",
      negative_count > 0 & positive_count == 0 ~ "negative difference",
      positive_count > 0 & negative_count > 0 ~ "mixed difference",
      TRUE ~ "no difference"
    ),
    opacity_level = factor(total_sig, levels = 0:3, labels = c("0/3", "1/3", "2/3", "3/3"))
  )

##Adjust the plot with fill color and discrete opacity levels
Clu_B <- ggplot(agg_contrast_df, aes(x = Class, y = contrast)) +
  geom_tile(aes(fill = dominant_sig, alpha = opacity_level), color = "white", size = 0.2) +
  scale_fill_manual(values = c(
    "positive difference" = "deepskyblue3",
    "negative difference" = "orange",
    "mixed difference" = "purple",
    "no difference" = "gray90"
  )) +
  scale_alpha_manual(
    values = c("0/3" = 0, "1/3" = 1/3, "2/3" = 2/3, "3/3" = 1),
    guide = guide_legend(title = "Number of Significant Differences")
  ) +
  labs(y = "", x = "Class") +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  scale_y_discrete(limits = rev(levels(agg_contrast_df$contrast)))

Fig3 <- clu_A + clu_B + plot_layout(nrow = 2, ncol = 1)
ggsave("Fig3_clusters.png", plot = Fig3, width = 20, height = 18, dpi = 300, bg = "transparent")

#Supplementary figure

contrast_50_90 <- emmeans(clu_m3, ~ Cluster|Class, level = 0.9) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_50_90_df <- as.data.frame(contrast_50_90)

contrast_50_90_df <- contrast_50_90_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_50_90_df <- contrast_50_90_df |> 
  mutate(Size = "median")

contrast_10_90 <- emmeans(clu_m4, ~ Cluster|Class, level = 0.9) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_10_90_df <- as.data.frame(contrast_10_90)

contrast_10_90_df <- contrast_10_90_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_90_df <- contrast_10_90_df |> 
  mutate(Size = "smallest")


contrast_95_90 <- emmeans(clu_m2, ~ Cluster|Class, level = 0.9) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_95_90_df <- as.data.frame(contrast_95_90)

contrast_95_90_df <- contrast_95_90_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_90_df <- contrast_95_90_df |> 
  mutate(Size = "largest")

contrast_Clu90_df <- rbind(contrast_10_90_df, contrast_50_90_df, contrast_95_90_df)
contrast_Clu90_df <- contrast_Clu90_df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Clu90_df <- contrast_Clu90_df |> 
  mutate(Sig = factor(Sig),
         Class = factor(Class, levels = c("Atab", "Amil", "Pdam", "Pmas")),
         contrast = factor(contrast, levels = c(
           'Inshore South vs. Offshore South',
           'Inshore South vs. Inshore Central',
           'Inshore South vs. Offshore Central',
           'Inshore South vs. Offshore North',
           'Inshore South vs. Torres Strait',
           'Offshore South vs. Inshore Central',
           'Offshore South vs. Offshore Central',
           'Offshore South vs. Offshore North',
           'Offshore South vs. Torres Strait',
           'Inshore Central vs. Offshore Central',
           'Inshore Central vs. Offshore North',
           'Inshore Central vs. Torres Strait',
           'Offshore Central vs. Offshore North',
           'Offshore Central vs. Torres Strait',
           'Offshore North vs. Torres Strait'
         )))

contrast_50_66 <- emmeans(clu_m3, ~ Cluster|Class, level = 0.66) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_50_66_df <- as.data.frame(contrast_50_66)

contrast_50_66_df <- contrast_50_66_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_50_66_df <- contrast_50_66_df |> 
  mutate(Size = "median")

contrast_10_66 <- emmeans(clu_m4, ~ Cluster|Class, level = 0.66) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_10_66_df <- as.data.frame(contrast_10_66)

contrast_10_66_df <- contrast_10_66_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_66_df <- contrast_10_66_df |> 
  mutate(Size = "smallest")


contrast_95_66 <- emmeans(clu_m2, ~ Cluster|Class, level = 0.66) |>
  contrast(method = cluster_mat) |>
  summary()
contrast_95_66_df <- as.data.frame(contrast_95_66)

contrast_95_66_df <- contrast_95_66_df |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_66_df <- contrast_95_66_df |> 
  mutate(Size = "largest")

contrast_Clu66_df <- rbind(contrast_10_66_df, contrast_50_66_df, contrast_95_66_df)
contrast_Clu66_df <- contrast_Clu66_df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Clu66_df <- contrast_Clu66_df |> 
  mutate(Sig = factor(Sig),
         Class = factor(Class, levels = c("Atab", "Amil", "Pdam", "Pmas")),
         contrast = factor(contrast, levels = c(
           'Inshore South vs. Offshore South',
           'Inshore South vs. Inshore Central',
           'Inshore South vs. Offshore Central',
           'Inshore South vs. Offshore North',
           'Inshore South vs. Torres Strait',
           'Offshore South vs. Inshore Central',
           'Offshore South vs. Offshore Central',
           'Offshore South vs. Offshore North',
           'Offshore South vs. Torres Strait',
           'Inshore Central vs. Offshore Central',
           'Inshore Central vs. Offshore North',
           'Inshore Central vs. Torres Strait',
           'Offshore Central vs. Offshore North',
           'Offshore Central vs. Torres Strait',
           'Offshore North vs. Torres Strait'
         )))

clu_95supp <- ggplot(contrast_Clu_df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "A",
       y = "", 
       x = "Estimated Marginal Means (HPD 95)") +
  facet_grid(~Class) +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  ) +
  background_grid(major = "x", color.major = "grey90") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  scale_y_discrete(limits = rev(levels(contrast_Clu_df$contrast)))

clu_90supp <- ggplot(contrast_Clu90_df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "B",
       y = "", 
       x = "Estimated Marginal Means (HPD 90)") +
  facet_grid(~Class) +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  ) +
  background_grid(major = "x", color.major = "grey90") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  scale_y_discrete(limits = rev(levels(contrast_Clu_df$contrast)))

clu_66supp <- ggplot(contrast_Clu66_df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "C",
       y = "", 
       x = "Estimated Marginal Means (HPD 66)") +
  facet_grid(~Class) +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  ) +
  background_grid(major = "x", color.major = "grey90") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  scale_y_discrete(limits = rev(levels(contrast_Clu_df$contrast)))



Fig.S1 <- clu_95supp + clu_90supp + clu_66supp + plot_layout(ncol = 1)
ggsave("FigS1_clusters.png", plot = Fig.S2, width = 20, height = 18, dpi = 300, bg = "transparent")

###############################################################################
###############           INSHORE/OFFSHORE ANALYSIS           #################
###############################################################################

# Import data
io_data <- read_csv("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Data/CSV/io_data.csv")

##Set inshore as reference

io_data <- io_data |> 
  mutate(
    CSP = factor(CSP),
    CSP = forcats::fct_relevel(CSP, 'Inshore')) 

# Model fitting

io_m2 <- brm(bf(log10(Area_cm2) ~ CSP*Class*Lat_cluster + (1 | Site/Plot),
                quantile = 0.95),
             family = asym_laplace(),
             data = io_data,
             cores = 3,
             chains = 3,
             warmup = 5000,
             thin = 12,
             iter = 20000,
             control = list(adapt_delta = 0.9, max_treedepth = 15))

io_m3 <- brm(bf(log10(Area_cm2) ~ CSP*Class*Lat_cluster + (1 | Site/Plot),
                quantile = 0.5),
             family = asym_laplace(),
             data = io_data,
             cores = 3,
             chains = 3,
             warmup = 5000,
             thin = 12,
             iter = 20000,
             control = list(adapt_delta = 0.9, max_treedepth = 15))

io_m4 <- brm(bf(log10(Area_cm2) ~ CSP*Class*Lat_cluster + (1 | Site/Plot),
                quantile = 0.1),
             family = asym_laplace(),
             data = io_data,
             cores = 3,
             chains = 3,
             warmup = 5000,
             thin = 20,
             iter = 50000,
             control = list(adapt_delta = 0.9, max_treedepth = 15))

# Models diagnostics and validation

load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/InshoreOffshore_analysis/io_m2.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/InshoreOffshore_analysis/io_m3.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/InshoreOffshore_analysis/io_m4.RData")

io_m2$fit |> stan_trace(inc_warmup = TRUE)
io_m2$fit |> stan_ac()
io_m2$fit |> stan_rhat()
io_m2$fit |> stan_ess()
io_m2 |> pp_check(type = "dens_overlay", ndraws = 100)
io_m2 |> conditional_effects()
resids <- make_brms_dharma_res(io_m2, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(io_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

io_m3$fit |> stan_trace(inc_warmup = TRUE)
io_m3$fit |> stan_ac()
io_m3$fit |> stan_rhat()
io_m3$fit |> stan_ess()
io_m3 |> pp_check(type = "dens_overlay", ndraws = 100)
io_m3 |> conditional_effects()
resids <- make_brms_dharma_res(io_m3, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(io_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

io_m4$fit |> stan_trace(inc_warmup = TRUE)
io_m4$fit |> stan_ac()
io_m4$fit |> stan_rhat()
io_m4$fit |> stan_ess()
io_m4 |> pp_check(type = "dens_overlay", ndraws = 100)
io_m4 |> conditional_effects()
resids <- make_brms_dharma_res(io_m4, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(io_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

# Extract posterior distributions

post95 <- io_m2 |> emmeans(~ CSP * Class * Lat_cluster) |> gather_emmeans_draws()
post50 <- io_m3 |> emmeans(~ CSP * Class * Lat_cluster) |> gather_emmeans_draws()
post10 <- io_m4 |> emmeans(~ CSP * Class * Lat_cluster) |> gather_emmeans_draws()
post95$Size = "largest"
post50$Size = "median"
post10$Size = "smallest"

io_all = rbind(post10, post50)
io_all2 = rbind(io_all, post95)
io_all2 <- io_all2 |> mutate(CSP = factor(CSP, levels = c("Inshore", "Offshore")),
                           Size = factor(Size, levels = c("smallest", "median", "largest")))


blues_palette <- brewer.pal(6, "Blues")[4:9]

io_A <- ggplot(io_all2, aes(x = .value, y = CSP, fill = Size)) +
  stat_slab(aes(fill = Size), alpha = 0.6) +
  stat_pointinterval(
    aes(color = Size), 
    point_interval = "median_hdi", # Use HPDI for the intervals
    .width = c(0.66, 0.95),
    position = position_dodge(width = .4, preserve = "single")
  ) +
  scale_fill_manual(values = blues_palette) +
  scale_color_manual(values = blues_palette) +
  scale_x_continuous(
    breaks = c(log10(seq(10, 90, 10)), log10(seq(100, 900, 100)), log10(seq(1000, 9000, 1000)), log10(seq(10000, 50000, 10000))),
    labels = c(10, "", "", "", "", "", "", "", "",
               100, "", "", "", "", "", "", "", "",
               1000, "", "", "", "", "", "", "", "",
               10000, "", "", "", 50000),
    limits = c(log10(10), log10(100000))
  ) +
  facet_grid(~ Class + Lat_cluster) +
  labs(title = "A",
       x = expression(Colony~size~(cm^2)),
       y = "",
       fill = "Confidence interval",
       color = "Size") +
  theme_light() +
  background_grid(major = "x", color.major = "grey90") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.spacing.y = unit(-0.5, "lines")
  )+
  scale_y_discrete(limits = rev(levels(io_all2$CSP)))

# Pairwise contrasts between inshore offshore

contrast_95 <- emmeans(io_m2, ~ CSP | Class * Lat_cluster) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_95 <- as.data.frame(contrast_95)

contrast_df_95 <- contrast_df_95 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_95 <- contrast_df_95 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "largest")

contrast_50 <- emmeans(io_m3, ~ CSP | Class * Lat_cluster) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_50 <- as.data.frame(contrast_50)

contrast_df_50 <- contrast_df_50 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_50 <- contrast_df_50 |> 
  mutate(Sig = factor(Sig)) |> 
  mutate(Size = "median")

contrast_10 <- emmeans(io_m4, ~ CSP | Class * Lat_cluster) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_10 <- as.data.frame(contrast_10)

contrast_df_10 <- contrast_df_10 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_10 <- contrast_df_10 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "smallest")

contrast_df <- rbind(contrast_df_10, contrast_df_50, contrast_df_95)
contrast_df <- contrast_df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

io_B <- ggplot(contrast_df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "B",
       y = "",
       x = "Estimated Marginal Means") +
  theme_light() +
  background_grid(major="none", minor="none") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  facet_grid(~ Class + Lat_cluster)


Fig4 <- io_A + io_B + plot_layout(nrow = 2, ncol = 1)
ggsave("Fig4_inshoreoffshore.png", plot = Fig4, width = 14, height = 10, dpi = 300, bg = "transparent")

#Supplementary figure
contrast_95_90 <- emmeans(io_m2, ~ CSP | Class * Lat_cluster, level = 0.9) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_95_90 <- as.data.frame(contrast_95_90)

contrast_df_95_90 <- contrast_df_95_90 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_95_90 <- contrast_df_95_90 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "largest")

contrast_50_90 <- emmeans(io_m3, ~ CSP | Class * Lat_cluster, level = 0.9) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_50_90 <- as.data.frame(contrast_50_90)

contrast_df_50_90 <- contrast_df_50_90 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_50_90 <- contrast_df_50_90 |> 
  mutate(Sig = factor(Sig)) |> 
  mutate(Size = "median")

contrast_10_90 <- emmeans(io_m4, ~ CSP | Class * Lat_cluster, level = 0.9) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_10_90 <- as.data.frame(contrast_10_90)

contrast_df_10_90 <- contrast_df_10_90 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_10_90 <- contrast_df_10_90 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "smallest")

contrast_90df <- rbind(contrast_df_10_90, contrast_df_50_90, contrast_df_95_90)
contrast_90df <- contrast_90df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_95_66 <- emmeans(io_m2, ~ CSP | Class * Lat_cluster, level = 0.66) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_95_66 <- as.data.frame(contrast_95_66)

contrast_df_95_66 <- contrast_df_95_66 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_95_66 <- contrast_df_95_66 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "largest")

contrast_50_66 <- emmeans(io_m3, ~ CSP | Class * Lat_cluster, level = 0.66) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_50_66 <- as.data.frame(contrast_50_66)

contrast_df_50_66 <- contrast_df_50_66 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_50_66 <- contrast_df_50_66 |> 
  mutate(Sig = factor(Sig)) |> 
  mutate(Size = "median")

contrast_10_66 <- emmeans(io_m4, ~ CSP | Class * Lat_cluster, level = 0.66) |> 
  contrast(method = "pairwise") |> 
  summary()
contrast_df_10_66 <- as.data.frame(contrast_10_66)

contrast_df_10_66 <- contrast_df_10_66 |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_df_10_66 <- contrast_df_10_66 |> 
  mutate(Sig = factor(Sig))|> 
  mutate(Size = "smallest")

contrast_66df <- rbind(contrast_df_10_66, contrast_df_50_66, contrast_df_95_66)
contrast_66df <- contrast_66df |> 
  mutate(Size = factor(Size, levels = c("smallest", "median", "largest")))

io_90 <- ggplot(contrast_90df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "A",
       y = "",
       x = "Estimated Marginal Means (HPD 90)") +
  theme_light() +
  background_grid(major="none", minor="none") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  facet_grid(~ Class + Lat_cluster)

io_66 <- ggplot(contrast_66df, aes(x = estimate, y = contrast, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig, shape = Size, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), size = 3) +
  geom_errorbar(aes(col = Sig, group = interaction(contrast, Size)), position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "B",
       y = "",
       x = "Estimated Marginal Means (HPD 66)") +
  theme_light() +
  background_grid(major="none", minor="none") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_shape_manual(values = c("smallest" = 15, "median" = 17, "largest" = 19)) +
  facet_grid(~ Class + Lat_cluster)


Fig.S2 <- io_90 + io_66 + plot_layout(nrow = 2, ncol = 1)
ggsave("FigS2_inshoreoffshore.png", plot = Fig.S3, width = 14, height = 10, dpi = 300, bg = "transparent")

###############################################################################
###################           HABITAT ANALYSIS           ######################
###############################################################################

# Import data
hab_data <- read_csv("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Data/CSV/hab_data.csv")

hab_data <- hab_data |> 
  mutate(Habitat = factor(Habitat),
         Class = factor(Class),
         Cluster = factor(Cluster))

# Model fitting
m2_habitat <- brm(bf(log10(Area_cm2) ~ Habitat*Class + Cluster + (1 | Reef/Site/Plot),
             quantile = 0.95),
          family = asym_laplace(),
          data = hab_data,
          cores = 3,
          chains = 3,
          warmup = 5000,
          thin = 10,
          iter = 20000,
          control = list(max_treedepth = 12)
)

hab_m3 <- brm(bf(log10(Area_cm2) ~ Habitat*Class + Cluster + (1 | Reef/Site/Plot),
             quantile = 0.5),
          family = asym_laplace(),
          data = hab_data,
          cores = 3,
          chains = 3,
          warmup = 5000,
          thin = 12,
          iter = 20000,
          control = list(max_treedepth = 12)
)

hab_m4 <- brm(bf(log10(Area_cm2) ~ Habitat*Class + Cluster + (1 | Reef/Site/Plot),
             quantile = 0.1),
          family = asym_laplace(),
          data = hab_data,
          cores = 3,
          chains = 3,
          warmup = 5000,
          thin = 12,
          iter = 20000,
          control = list(max_treedepth = 12)
)

# Model diagnostics and validation


m2_habitat$fit |> stan_trace(inc_warmup = TRUE)
m2_habitat$fit |> stan_ac()
m2_habitat$fit |> stan_rhat()
m2_habitat$fit |> stan_ess()
m2_habitat |> pp_check(type = "dens_overlay", ndraws = 100)
m2_habitat |> conditional_effects()
resids <- make_brms_dharma_res(m2_habitat, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(hab_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

hab_m3$fit |> stan_trace(inc_warmup = TRUE)
hab_m3$fit |> stan_ac()
hab_m3$fit |> stan_rhat()
hab_m3$fit |> stan_ess()
hab_m3 |> pp_check(type = "dens_overlay", ndraws = 100)
hab_m3 |> conditional_effects()
resids <- make_brms_dharma_res(hab_m3, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(hab_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

hab_m4$fit |> stan_trace(inc_warmup = TRUE)
hab_m4$fit |> stan_ac()
hab_m4$fit |> stan_rhat()
hab_m4$fit |> stan_ess()
hab_m4 |> pp_check(type = "dens_overlay", ndraws = 100)
hab_m4 |> conditional_effects()
resids <- make_brms_dharma_res(hab_m4, integerResponse = FALSE)
wrap_elements(~testUniformity(resids))+
  wrap_elements(~plotResiduals(resids, form = factor(rep(1, nrow(hab_data)))))+
  wrap_elements(~plotResiduals(resids))+
  wrap_elements(~testDispersion(resids))

# Extract posterior distributions


load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/m2_habitat.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/hab_m3.RData")
load("C:/Users/jc815475/OneDrive - James Cook University/PhD data/2_SizeFreqDistribution/GitHub/Model_output/Spatial_analysis/hab_m4.RData")

post95_h <- m2_habitat |> emmeans(~ Habitat|Class) |> gather_emmeans_draws()
post50_h <- hab_m3 |> emmeans(~ Habitat|Class) |> gather_emmeans_draws()
post10_h <- hab_m4 |> emmeans(~ Habitat|Class) |> gather_emmeans_draws()
post10_h$Size = "smallest"
post50_h$Size = "median"
post95_h$Size = "largest"

post_Hab_df <- rbind(post10_h, post50_h, post95_h)

post_Hab_df <- post_Hab_df |> 
  mutate(Habitat = factor(Habitat, levels = c("BA_D", "FL_D", "FR_D", "LA_S", "BA_S", "FL_S", "FR_S")),
         Size = factor(Size, levels = c("smallest", "median", "largest")))

post_Hab_df <- post_Hab_df |> 
  mutate(Habitat2 = fct_recode(Habitat,
                               `Low exposure (D)` = "BA_D", 
                               `Mid exposure (D)` = "FL_D", 
                               `High exposure (D)` ="FR_D", 
                               `Lowest exposure` ="LA_S", 
                               `Low exposure (S)` ="BA_S", 
                               `Mid exposure (S)` ="FL_S", 
                               `High exposure (S)` ="FR_S")) |> 
  mutate(Habitat2 = factor(Habitat2, levels = c("Lowest exposure",
                                                "Low exposure (S)",
                                                "Mid exposure (S)",
                                                "High exposure (S)",
                                                "Low exposure (D)",
                                                "Mid exposure (D)",
                                                "High exposure (D)")))

# Plot posterior distributions of habitat models

blues_palette <- brewer.pal(6, "Blues")[4:9]

hab1 <- ggplot(post_Hab_df, aes(x = .value, y = Habitat2, fill = Size)) +
  stat_slab(aes(fill = Size), alpha = 0.6) +
  stat_pointinterval(
    aes(color = Size), 
    point_interval = "median_hdi", # Use HPDI for the intervals
    .width = c(0.66, 0.95),
    position = position_dodge(width = .4, preserve = "single")
  ) +
  scale_fill_manual(values = blues_palette) +
  scale_color_manual(values = blues_palette) +
  scale_x_continuous(
    breaks = log10(c(seq(10, 90, 10), seq(100, 900, 100), seq(1000, 9000, 1000), seq(10000, 20000, 10000))),
    labels = c(10, "", "", "", "", "", "", "", "",
               100, "", "", "", "", "", "", "", "",
               1000, "", "", "", "", "", "", "", "",
               10000, 20000),
    limits = c(log10(10), log10(20000))
  ) +
  facet_wrap(~Class, ncol = 3) +
  labs(x = expression(Colony~size~(cm^2)),
       y = "",
       fill = "Confidence interval") +
  theme_light() +
  background_grid(major = "x", color.major = "grey90") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.spacing.y = unit(-0.5, "lines")
  )+
  scale_y_discrete(limits = rev(levels(post_Hab_df$Habitat2)))


# Planned contrasts between habitats

#Create contrast matrix based on expectation that larger corals will be found in least/less exposed habitats (so put habitat where larger corals are expected in first position)

habitat_mat <- list('BA vs FR (S)' = c(0, 1, 0, 0, 0, -1, 0),
                    'BA vs FR (D)' = c(1, 0, 0, 0, -1, 0, 0),
                    'BA vs FL (S)' = c(0, 1, 0, -1, 0, 0, 0),
                    'BA vs FL (D)' = c(1, 0, -1, 0, 0, 0, 0),
                    'LA vs BA (S)' = c(0, -1, 0, 0, 0, 0, 1),
                    'BA (D) vs BA (S)' = c(1, -1, 0, 0, 0, 0, 0),
                    'FL vs FR (S)' = c(0 ,0 , 0, 1, 0, -1, 0),
                    'FL vs FR (D)' = c(0, 0, 1, 0, -1, 0, 0),
                    'LA vs FL (S)' = c(0, 0, 0, -1, 0, 0, 1),
                    'LA vs FR (S)' = c(0, 0, 0, 0, 0, -1, 1),
                    'FR (D) vs FR (S)' = c(0, 0, 0, 0, 1, -1, 0),
                    'FL (D) vs FL (S)' = c(0, 0, 1, -1, 0, 0, 0))

contrast_50_h <- emmeans(hab_m3, ~ Habitat|Class) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_50_h <- as.data.frame(contrast_50_h)

contrast_50_h <- contrast_50_h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_h <- emmeans(hab_m4, ~ Habitat|Class) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_10_h <- as.data.frame(contrast_10_h)

contrast_10_h <- contrast_10_h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_h <- emmeans(m2_habitat, ~ Habitat|Class) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_95_h <- as.data.frame(contrast_95_h)

contrast_95_h <- contrast_95_h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_h <- contrast_10_h |> 
  mutate(Size = "smallest")
contrast_50_h <- contrast_50_h |> 
  mutate(Size = "median")
contrast_95_h <- contrast_95_h |> 
  mutate(Size = "largest")

contrast_Hab_df <- rbind(contrast_10_h, contrast_50_h, contrast_95_h)
contrast_Hab_df <- contrast_Hab_df |> 
  mutate(Sig = factor(Sig),
         Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Hab_df <- contrast_Hab_df |> 
  mutate(contrast2 = fct_recode(contrast,
                                `Low vs High exposure (S)` = "BA vs FR (S)", 
                                `Low vs High exposure (D)` = "BA vs FR (D)", 
                                `Low vs Mid exposure (S)` = "BA vs FL (S)", 
                                `Low vs Mid exposure (D)` = "BA vs FL (D)", 
                                `Mid vs High exposure (S)` = "FL vs FR (S)", 
                                `Mid vs High exposure (D)` = "FL vs FR (D)",
                                `Lowest vs Low exposure (S)` = "LA vs BA (S)", 
                                `Lowest vs High exposure (S)` = "LA vs FR (S)",
                                `Lowest vs Mid exposure (S)` = "LA vs FL (S)",
                                `Deep vs Shallow (Low exposure)` = "BA (D) vs BA (S)",
                                `Deep vs Shallow (Mid exposure)` = "FL (D) vs FL (S)",
                                `Deep vs Shallow (High exposure)` = "FR (D) vs FR (S)"
  ))|> 
  mutate(contrast2 = factor(contrast2, levels = c("Lowest vs Low exposure (S)", 
                                                  "Lowest vs Mid exposure (S)",
                                                  "Lowest vs High exposure (S)",
                                                  "Low vs Mid exposure (D)", 
                                                  "Low vs Mid exposure (S)", 
                                                  "Low vs High exposure (D)", 
                                                  "Low vs High exposure (S)", 
                                                  "Mid vs High exposure (D)", 
                                                  "Mid vs High exposure (S)",
                                                  "Deep vs Shallow (Low exposure)",
                                                  "Deep vs Shallow (Mid exposure)",
                                                  "Deep vs Shallow (High exposure)")))

# Custom function to set y-axis labels for each facet
set_y_axis_labels <- function(data) {
  data$y_axis_labels <- ifelse(data$Class == unique(data$Class)[1] | data$Class == unique(data$Class)[4], as.character(data$contrast2), "")
  return(data)
}

# Apply the custom function to the dataframe
contrast_Hab_df <- set_y_axis_labels(contrast_Hab_df)

# Aggregate the data by contrast and Class
agg_contrast_df <- contrast_Hab_df |> 
  group_by(contrast2, Class) |> 
  summarise(
    positive_count = sum(Sig == "positive difference"),
    negative_count = sum(Sig == "negative difference"),
    total_sig = positive_count + negative_count,
    .groups = 'drop'
  ) |> 
  mutate(
    dominant_sig = case_when(
      positive_count > 0 & negative_count == 0 ~ "positive difference",
      negative_count > 0 & positive_count == 0 ~ "negative difference",
      positive_count > 0 & negative_count > 0 ~ "mixed difference",
      TRUE ~ "no difference"
    ),
    opacity_level = factor(total_sig, levels = 0:3, labels = c("0/3", "1/3", "2/3", "3/3"))
  )

hab2 <- ggplot(agg_contrast_df, aes(x = Class, y = contrast2)) +
  geom_tile(aes(fill = dominant_sig, alpha = opacity_level), color = "white", size = 0.2) +
  scale_fill_manual(values = c(
    "positive difference" = "deepskyblue3",
    "negative difference" = "orange",
    "mixed difference" = "purple",
    "no difference" = "gray90"
  )) +
  scale_alpha_manual(
    values = c("0/3" = 0, "1/3" = 1/3, "2/3" = 2/3, "3/3" = 1),
    guide = guide_legend(title = "Number of Significant Differences")
  ) +
  labs(y = "", x = "Class") +
  theme_light() +
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  scale_y_discrete(limits = rev(levels(agg_contrast_df$contrast2)))

Fig5 <- hab1 + hab2 + plot_layout(ncol = 1)
ggsave("Fig5_habitats.png", plot = Fig5, width = 20, height = 10, dpi = 300, bg = "transparent")

#Supplementary figure

## HPD 90%
contrast_50_90h <- emmeans(hab_m3, ~ Habitat|Class, level = 0.9) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_50_90h <- as.data.frame(contrast_50_90h)

contrast_50_90h <- contrast_50_90h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_90h <- emmeans(hab_m4, ~ Habitat|Class, level = 0.9) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_10_90h <- as.data.frame(contrast_10_90h)

contrast_10_90h <- contrast_10_90h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_90h <- emmeans(m2_habitat, ~ Habitat|Class, level = 0.9) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_95_90h <- as.data.frame(contrast_95_90h)

contrast_95_90h <- contrast_95_90h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_90h <- contrast_10_90h |> 
  mutate(Size = "smallest")
contrast_50_90h <- contrast_50_90h |> 
  mutate(Size = "median")
contrast_95_90h <- contrast_95_90h |> 
  mutate(Size = "largest")

contrast_Hab_90df <- rbind(contrast_10_90h, contrast_50_90h, contrast_95_90h)
contrast_Hab_90df <- contrast_Hab_90df |> 
  mutate(Sig = factor(Sig),
         Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Hab_90df <- contrast_Hab_90df |> 
  mutate(contrast2 = fct_recode(contrast,
                                `Low vs High exposure (S)` = "BA vs FR (S)", 
                                `Low vs High exposure (D)` = "BA vs FR (D)", 
                                `Low vs Mid exposure (S)` = "BA vs FL (S)", 
                                `Low vs Mid exposure (D)` = "BA vs FL (D)", 
                                `Mid vs High exposure (S)` = "FL vs FR (S)", 
                                `Mid vs High exposure (D)` = "FL vs FR (D)",
                                `Lowest vs Low exposure (S)` = "LA vs BA (S)", 
                                `Lowest vs High exposure (S)` = "LA vs FR (S)",
                                `Lowest vs Mid exposure (S)` = "LA vs FL (S)",
                                `Deep vs Shallow (Low exposure)` = "BA (D) vs BA (S)",
                                `Deep vs Shallow (Mid exposure)` = "FL (D) vs FL (S)",
                                `Deep vs Shallow (High exposure)` = "FR (D) vs FR (S)"
  ))|> 
  mutate(contrast2 = factor(contrast2, levels = c("Lowest vs Low exposure (S)", 
                                                  "Lowest vs Mid exposure (S)",
                                                  "Lowest vs High exposure (S)",
                                                  "Low vs Mid exposure (D)", 
                                                  "Low vs Mid exposure (S)", 
                                                  "Low vs High exposure (D)", 
                                                  "Low vs High exposure (S)", 
                                                  "Mid vs High exposure (D)", 
                                                  "Mid vs High exposure (S)",
                                                  "Deep vs Shallow (Low exposure)",
                                                  "Deep vs Shallow (Mid exposure)",
                                                  "Deep vs Shallow (High exposure)")))

set_y_axis_labels <- function(data) {
  data$y_axis_labels <- ifelse(data$Class == unique(data$Class)[1] | data$Class == unique(data$Class)[4], as.character(data$contrast2), "")
  return(data)
}

contrast_Hab_90df <- set_y_axis_labels(contrast_Hab_90df)

## HPD 66%
contrast_50_66h <- emmeans(hab_m3, ~ Habitat|Class, level = 0.66) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_50_66h <- as.data.frame(contrast_50_66h)

contrast_50_66h <- contrast_50_66h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_66h <- emmeans(hab_m4, ~ Habitat|Class, level = 0.66) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_10_66h <- as.data.frame(contrast_10_66h)

contrast_10_66h <- contrast_10_66h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_95_66h <- emmeans(m2_habitat, ~ Habitat|Class, level = 0.66) |> 
  contrast(method = habitat_mat) |> 
  summary()
contrast_95_66h <- as.data.frame(contrast_95_66h)

contrast_95_66h <- contrast_95_66h |> 
  mutate(Sig = case_when(
    lower.HPD < 0 & upper.HPD > 0 ~ 'no difference',
    upper.HPD >= 0 & lower.HPD >= 0 ~ 'positive difference',
    upper.HPD <= 0 & lower.HPD <= 0 ~ 'negative difference'
  ))

contrast_10_66h <- contrast_10_66h |> 
  mutate(Size = "smallest")
contrast_50_66h <- contrast_50_66h |> 
  mutate(Size = "median")
contrast_95_66h <- contrast_95_66h |> 
  mutate(Size = "largest")

contrast_Hab_66df <- rbind(contrast_10_66h, contrast_50_66h, contrast_95_66h)
contrast_Hab_66df <- contrast_Hab_66df |> 
  mutate(Sig = factor(Sig),
         Size = factor(Size, levels = c("smallest", "median", "largest")))

contrast_Hab_66df <- contrast_Hab_66df |> 
  mutate(contrast2 = fct_recode(contrast,
                                `Low vs High exposure (S)` = "BA vs FR (S)", 
                                `Low vs High exposure (D)` = "BA vs FR (D)", 
                                `Low vs Mid exposure (S)` = "BA vs FL (S)", 
                                `Low vs Mid exposure (D)` = "BA vs FL (D)", 
                                `Mid vs High exposure (S)` = "FL vs FR (S)", 
                                `Mid vs High exposure (D)` = "FL vs FR (D)",
                                `Lowest vs Low exposure (S)` = "LA vs BA (S)", 
                                `Lowest vs High exposure (S)` = "LA vs FR (S)",
                                `Lowest vs Mid exposure (S)` = "LA vs FL (S)",
                                `Deep vs Shallow (Low exposure)` = "BA (D) vs BA (S)",
                                `Deep vs Shallow (Mid exposure)` = "FL (D) vs FL (S)",
                                `Deep vs Shallow (High exposure)` = "FR (D) vs FR (S)"
  ))|> 
  mutate(contrast2 = factor(contrast2, levels = c("Lowest vs Low exposure (S)", 
                                                  "Lowest vs Mid exposure (S)",
                                                  "Lowest vs High exposure (S)",
                                                  "Low vs Mid exposure (D)", 
                                                  "Low vs Mid exposure (S)", 
                                                  "Low vs High exposure (D)", 
                                                  "Low vs High exposure (S)", 
                                                  "Mid vs High exposure (D)", 
                                                  "Mid vs High exposure (S)",
                                                  "Deep vs Shallow (Low exposure)",
                                                  "Deep vs Shallow (Mid exposure)",
                                                  "Deep vs Shallow (High exposure)")))

set_y_axis_labels <- function(data) {
  data$y_axis_labels <- ifelse(data$Class == unique(data$Class)[1] | data$Class == unique(data$Class)[4], as.character(data$contrast2), "")
  return(data)
}

contrast_Hab_66df <- set_y_axis_labels(contrast_Hab_66df)


## Supplementary figure
hab_A95 <- ggplot(contrast_Hab_df[contrast_Hab_df$Size == "smallest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Smallest colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_B95 <- ggplot(contrast_Hab_df[contrast_Hab_df$Size == "median", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Median colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_C95 <- ggplot(contrast_Hab_df[contrast_Hab_df$Size == "largest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Largest colonies",
       y = "",
       x = "Estimated Marginal Means (HPD 95)") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

HPD95 <- hab_A95 + hab_B95 + hab_C95 + plot_layout(nrow = 3, ncol = 1)

hab_A90 <- ggplot(contrast_Hab_90df[contrast_Hab_90df$Size == "smallest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Smallest colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_B90 <- ggplot(contrast_Hab_90df[contrast_Hab_90df$Size == "median", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Median colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_C90 <- ggplot(contrast_Hab_90df[contrast_Hab_90df$Size == "largest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Largest colonies",
       y = "",
       x = "Estimated Marginal Means (HPD 90)") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

HPD90 <- hab_A90 + hab_B90 + hab_C90 + plot_layout(nrow = 3, ncol = 1)


hab_A66 <- ggplot(contrast_Hab_66df[contrast_Hab_66df$Size == "smallest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Smallest colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_B66 <- ggplot(contrast_Hab_66df[contrast_Hab_66df$Size == "median", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Median colonies",
       y = "",
       x = "") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

hab_C66 <- ggplot(contrast_Hab_66df[contrast_Hab_66df$Size == "largest", ], aes(x = estimate, y = contrast2, xmin = lower.HPD, xmax = upper.HPD, fill = Sig)) +
  geom_point(aes(col = Sig)) +
  geom_errorbar(aes(col = Sig, width = 0.2)) +
  labs(title = "Largest colonies",
       y = "",
       x = "Estimated Marginal Means (HPD 66)") +
  theme_light() +
  background_grid(major="none",minor="none")+
  theme(
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "gray28", colour = "white", size = 1),
    legend.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  facet_wrap(~ Class, nrow = 1)+
  scale_fill_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange")) +
  scale_color_manual(values = c("no difference" = "gray", "positive difference" = "deepskyblue3", "negative difference" = "orange"))+
  scale_y_discrete(limits = rev(levels(contrast_Hab_df$contrast2)))

HPD66 <- hab_A66 + hab_B66 + hab_C66 + plot_layout(nrow = 3, ncol = 1)


plot_HPD95 <- ggdraw() + draw_plot(HPD95)
plot_HPD90 <- ggdraw() + draw_plot(HPD90)
plot_HPD66 <- ggdraw() + draw_plot(HPD66)

Fig.S3 <- plot_grid(
  plot_HPD95, plot_HPD90, plot_HPD66,
  nrow = 3, align = "v", labels = c("A", "B", "C")
)

ggsave("FigS3_habitats.png", plot = Fig.S3, width = 30, height = 70, units = "cm", dpi = 300, bg = "transparent")
