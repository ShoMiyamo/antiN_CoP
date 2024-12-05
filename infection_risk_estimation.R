library('brms')
library('ggplot2')
library('gratia')
library('dplyr')
library('tidyr')
library('viridis')
library('scales')
library('cmdstanr')

# Step 1: Predict baseline infection probability using a Poisson regression model
# The model adjusts for demographic and regional factors to reduce bias and estimate baseline infection risk.

model1 <- glm(new_infection ~ aged + male_sex + comorbidity 
              + vaccination_count + prior_infect + bivalent +
                sitea1 + sitea2 + sitea3 + sitem1 + sitem2 + sitem3 + 
                siteo1 + siteo2 + siteo3 + sitet1 + sitet2 + sitet3 + 
                sitef1 + sitef2 + sitef3,
              family=poisson, data = data1)

# Predict infection probability based on the baseline characteristics
data1$p_s <- predict(model1, newdata = data1, type = 'response')

# Calculate inverse probability weights (IPW) to adjust for baseline characteristics
data1 <- data1 %>% mutate(ipw = 1 / p_s)

# Step 2: Bayesian generalized additive model (GAM) for immune correlates analysis
# The GAM is used to explore the nonlinear relationship between log-transformed antibody titers and infection risk.
# Cubic spline smoothing (k = 3) is applied to allow flexibility in the functional form of antibody effects.
gam1 <- brm(
  bf(new_infection | weights(ipw) ~ s(anti_s, bs = "cr", k = 3) + s(anti_n, bs = "cr", k = 3)),
  data = data1, 
  family = poisson(), 
  cores = 4, 
  seed = 1234,  
  iter = 5000, 
  warmup = 1000, 
  thin = 5,  
  refresh = 500, 
  control = list(adapt_delta = 0.99, max_treedepth = 15),  
  backend = 'cmdstanr'  
)

# Step 3: Function to visualize infection risk predictions and antibody effects
# This function generates three types of plots:
# 1. Predicted risk vs. anti-S antibody titers (pgams)
# 2. Predicted risk vs. anti-N antibody levels (pgamn)
# 3. Absolute risk distribution based on antibody combinations (AR)

predict_and_plot <- function(model, data, colp = c("#D98825", "#347aa5"), base_size = 12) {
  # Calculate conditional effects to explore the model's smooth terms
  ce <- brms::conditional_effects(model, surface = TRUE, ordinary = TRUE, resolution = 1000, ask = FALSE, re_formula = NULL)
  
  # Extract conditional effects for anti-S and anti-N
  ce_s <- ce$`anti_s`
  ce_n <- ce$`anti_n`
  
  # Predict infection risk for the given data
  predicted_risk <- predict(model, newdata = data)
  if (nrow(predicted_risk) != nrow(data)) {
    stop("The number of rows in the predicted risk does not match the number of rows in the input data.")
  }
  data$.predicted_risk <- predicted_risk
  
  # Ensure minimum prediction threshold for plotting
  unique_values <- sort(unique(data$.predicted_risk[, "Estimate"]))
  MIN <- ifelse(length(unique_values) > 1, unique_values[2], unique_values[1])
  data <- data %>% mutate(Estimate = ifelse(.predicted_risk[, "Estimate"] < MIN, MIN / 2, .predicted_risk[, "Estimate"]))
  
  # Define overall risk reduction thresholds
  rr_overall <- 0.184  # Overall infection risk
  rr90 <- rr_overall * 0.1  # Threshold for 90% risk reduction
  
  # Plot: Anti-S antibody effects on predicted risk
  pgams <- ggplot() + theme_bw(base_size = base_size) +
    geom_hline(yintercept = rr_overall, color = "#A2A637", linetype = 1, alpha = 1, size = .8) +
    geom_hline(yintercept = rr90, color = "#A2A637", linetype = 3, alpha = 1, size = 1) +
    geom_ribbon(data = ce_s, aes(x = effect1__, ymin = lower__, ymax = upper__), alpha = 0.4, fill = "#D98825") +
    geom_line(data = ce_s, aes(x = effect1__, y = estimate__), size = 1, alpha = .9, color = "#D98825") +
    geom_point(data = data, aes(x = anti_s, y = Estimate), alpha = .1, size = 1, color = "#024873") +
    labs(x = bquote(paste(Log[10], " anti-S antibodies (BAU/mL)")), y = bquote(paste(Log[10], " absolute risk"))) +
    scale_y_log10(labels = trans_format("log10", math_format(.x))) +
    coord_cartesian(ylim = c(10^(-4.5), 10^(-0.5)))
  
  # Plot: Anti-N antibody effects on predicted risk
  pgamn <- ggplot() + theme_bw(base_size = base_size) +
    geom_hline(yintercept = rr_overall, color = "#A2A637", linetype = 1, alpha = 1, size = .8) +
    geom_hline(yintercept = rr90, color = "#A2A637", linetype = 3, alpha = 1, size = 1) +
    geom_ribbon(data = ce_n, aes(x = effect1__, ymin = lower__, ymax = upper__, fill = Hybrid), alpha = 0.4, fill = "#D98825") +
    geom_line(data = ce_n, aes(x = effect1__, y = estimate__, color = Hybrid), size = 1, alpha = .9, color = "#D98825") +
    geom_point(data = data, aes(x = anti_n, y = Estimate), alpha = .1, size = 1, color = "#024873") +
    labs(x = bquote(paste(Log[10], " anti-N antibodies (COI)")), y = bquote(paste(Log[10], " absolute risk"))) +
    scale_y_log10(labels = trans_format("log10", math_format(.x))) +
    coord_cartesian(ylim = c(10^(-4.5), 10^(-0.5)))
  
  # Plot: Risk distribution based on anti-S and anti-N combinations
  len <- 100
  anti_s_new <- seq(min(data$anti_s), max(data$anti_s), length.out = len)
  anti_n_new <- seq(min(data$anti_n), max(data$anti_n), length.out = len)
  newdata <- expand.grid(anti_s = anti_s_new, anti_n = anti_n_new)
  newdata$fit <- fitted(model, newdata = newdata)
  newdata <- newdata %>% mutate(logfit = log10(fit[, "Estimate"]),
                                fit2 = fit[, "Estimate"],
                                rr = fit2 / max(fit2),
                                protection = 1 - rr)
  
  AR <- ggplot(newdata, aes(x = anti_n, y = anti_s)) + theme_bw(base_size = 12) +
    geom_raster(aes(fill = logfit), alpha = 0.9) +
    geom_contour(aes(z = logfit), colour = "white", binwidth = 1, alpha = .9, linetype = "solid") +
    geom_contour(aes(z = logfit), colour = "white", binwidth = .5, alpha = .5, linetype = "dashed") +
    labs(fill = bquote(paste(Log[10], " absolute risk")),
         x = bquote(paste(Log[10], " anti-N antibodies (COI)")),
         y = bquote(paste(Log[10], " anti-S antibodies (BAU/mL)")))
  
  return(list(pgams = pgams, pgamn = pgamn, AR = AR))
}

# Example usage of the function to generate and view plots
plots <- predict_and_plot(gam1, data1)
plots$pgams  # Plot: Anti-S antibody effects
plots$pgamn  # Plot: Anti-N antibody effects
plots$AR     # Plot: Risk distribution based on antibody combinations
