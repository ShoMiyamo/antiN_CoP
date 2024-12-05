## Modeling the anti-N antibody response

library(rstan) 

# Model for primary infection
## Load data for individuals with primary infection
d <- d_primaryinfect 
N <- nrow(d)  # Number of observations

## Prepare data for the Stan model
# `t` is the time of measurement (days post-diagnosis)
# `y` is the log-transformed anti-N antibody levels with an offset added for numerical stability
data <- list(N = N, 
             t = d$days_post_diagnosis, 
             y = d$anti_n + 1.5)

## Fit the Bayesian model using Stan
# The model file `nmodel01.stan` implements the cumulative gamma function with decay rate
fit1 <- stan(file = 'nmodel01.stan', data = data, iter = 4000, thin = 5, seed = 1234)

# Model for reinfection
## Load data for individuals with reinfection
d <- d_reinfect #%>% filter(DPE <= 365) # Optionally filter for data within 365 days post-diagnosis
N <- nrow(d)  # Number of observations

## Prepare data for the Stan model
# `t` is the time of measurement (days post-diagnosis)
# `y` is the log-transformed anti-N antibody levels with an offset added for numerical stability
data <- list(N = N, 
             t = d$days_post_diagnosis, 
             y = d$anti_n + 1)

## Fit the Bayesian model using Stan
# The model file `nmodel01.stan` implements the cumulative gamma function with decay rate
fit2 <- stan(file = 'nmodel01.stan', data = data, iter = 4000, thin = 5, seed = 1234)

# Extract posterior samples from the primary infection model
posterior_samples1 <- rstan::extract(fit1)

## Define a function to compute predicted antibody responses
compute_y_pred1 <- function(t) {
  # Prediction based on the cumulative gamma function and decay
  posterior_samples1$h * pgamma(t, posterior_samples1$alpha, posterior_samples1$beta) * exp(-posterior_samples1$lambda * t)
}

# Generate predictions for a sequence of time points
t_values <- seq(0, 1100)  # Time range for predictions
predictions1 <- sapply(t_values, compute_y_pred1)

## Compute the 95% credible interval and median of the predictions
ci_lower1 <- apply(predictions1, 2, function(y_pred) quantile(y_pred, 0.025))
med1 <- apply(predictions1, 2, function(y_pred) quantile(y_pred, 0.50))
ci_upper1 <- apply(predictions1, 2, function(y_pred) quantile(y_pred, 0.975))

## Prepare data for plotting (adjusting for the initial offset)
data_to_plot <- data.frame(t = t_values, 
                           ci_lower = ci_lower1 - 1.5, 
                           ci_upper = ci_upper1 - 1.5, 
                           med = med1 - 1.5)

# Extract posterior samples from the reinfection model
posterior_samples2 <- rstan::extract(fit2)

## Define a function to compute predicted antibody responses
compute_y_pred2 <- function(t) {
  # Prediction based on the cumulative gamma function and decay
  posterior_samples2$h * pgamma(t, posterior_samples2$alpha, posterior_samples2$beta) * exp(-posterior_samples2$lambda * t)
}

# Generate predictions for a sequence of time points
predictions2 <- sapply(t_values, compute_y_pred2)

## Compute the 95% credible interval and median of the predictions
ci_lower2 <- apply(predictions2, 2, function(y_pred) quantile(y_pred, 0.025))
med2 <- apply(predictions2, 2, function(y_pred) quantile(y_pred, 0.50))
ci_upper2 <- apply(predictions2, 2, function(y_pred) quantile(y_pred, 0.975))

## Prepare data for plotting (adjusting for the initial offset)
data_to_plot2 <- data.frame(t = t_values, 
                            ci_lower = ci_lower2 - 1, 
                            ci_upper = ci_upper2 - 1, 
                            med = med2 - 1)

# Plotting antibody response
library(ggplot2)  # For visualization

p <- ggplot() + theme_classic(base_size = 12) +
  # Reinfection: 95% credible interval and median
  geom_ribbon(data = data_to_plot2, aes(x = t, ymin = ci_lower, ymax = ci_upper), fill = "#024873", alpha = 0.2) +
  geom_line(data = data_to_plot2, aes(x = t, y = med), color = "#024873", size = 1) +
  
  # Primary infection: 95% credible interval and median
  geom_ribbon(data = data_to_plot, aes(x = t, ymin = ci_lower, ymax = ci_upper), fill = "#D98825", alpha = 0.4) +
  geom_line(data = data_to_plot, aes(x = t, y = med), color = "#D98825", size = 1) +
  
  # Observed data points for primary infection and reinfection
  geom_point(data = d_primaryinfect, aes(x = days_post_diagnosis, y = anti_n), color = "#D98825", alpha = 0.2, size = 1.2) +
  geom_point(data = d_reinfect, aes(x = days_post_diagnosis, y = anti_n), color = "#024873", alpha = 0.3, size = 1.2) +
  geom_hline(yintercept = 0, color = "black", linetype = 3, alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = -1, color = "gray40", linetype = 3, alpha = 0.8, size = 0.8) +
  scale_y_continuous(breaks = seq(-2, 4, 1), limits = c(-1.5, 3)) +
  scale_x_continuous(breaks = seq(0, 700, 100), limits = c(0, 750)) +
  labs(y = bquote(paste(Log[10], " anti-N antibodies (COI)")), x = "Days post diagnosis") +
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.position = "none")
p
