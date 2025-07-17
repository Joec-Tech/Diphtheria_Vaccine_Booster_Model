library(deSolve)
library(ggplot2)

# Define the model
SIRBVVW_model <- function(t, state, parameters) {
  S <- state["S"]; I <- state["I"]; R <- state["R"]
  B <- state["B"]; V1 <- state["V1"]; V2 <- state["V2"]; W <- state["W"]
  
  with(as.list(parameters), {
    N <- S + I + R + V1 + V2
    
    dS  <- Lambda - (beta_h * I / N + (beta0 * B / (kappa + B))) * S - nu_s * S - mu * S
    dI  <- (beta_h * I / N + (beta0 * B / (kappa + B))) * S - (1 - xi_v_1) * V1 + (1 - xi_w) * W - mu * I - delta * I - gamma * I
    dR  <- gamma * I - mu * R - zeta_R * R
    dB  <- xi_I * I - mu_B * B + alpha
    dV1 <- nu_s * S - mu * V1 - psi_v_1 * V1 - (1 - xi_v_1) * V1
    dV2 <- psi_v_1 * V1 - rho_v_2 * V2 - mu * V2 + rho_W * W
    dW  <- zeta_R * R - mu * W - rho_W * W - (1 - xi_w) * W + rho_v_2 * V2
    
    list(c(dS, dI, dR, dB, dV1, dV2, dW))
  })
}

# Initial states
initial_state <- c(S = 6846408, I = 10661, R = 10292, B = 500, V1 = 272917, V2 = 0, W = 0)

# Parameters (excluding rho_W which we will vary)
base_params <- c(
  gamma = 0.0061884, kappa = 500000, Lambda = 38, beta_h = 0.00003,
  beta0 = 0.234, psi_v_1 = 0, mu = 0.0003465, rho_v_2 = 0.0083,
  nu_s = 1 / 385, zeta_R = 1 / 15, xi_I = 0.8, alpha = 0.014,
  mu_B = 0.0345, delta = 0.00022187, xi_v_1 = 0.9917, xi_w = 0.8925
)

# Time
times <- seq(0, 155, by = 1)



# Booster rates to test
booster_rates <- c(0, 0.005, 0.01, 0.02)
colors <- c("black", "blue", "green", "red")

# Store output
results <- list()

for (i in seq_along(booster_rates)) {
  param_set <- c(base_params, rho_W = booster_rates[i])
  out <- ode(y = initial_state, times = times, func = SIRBVVW_model,
             parms = param_set, atol = 1e-6, rtol = 1e-6)
  results[[i]] <- as.data.frame(out)
  results[[i]]$Scenario <- paste0("rho_W = ", booster_rates[i])
}

# Combine all outputs
all_results <- do.call(rbind, results)

# Plot Infectious Population (I)
ggplot(all_results, aes(x = time, y = I, color = Scenario)) +
  geom_line(size = 1) +
  labs(title = "Effect of Booster Rate for Waning sub population on Infectious Population",
       x = "Time (days)", y = "Number of Infected Individuals") +
  theme_minimal()

# Plot Boosted Individuals (W)
ggplot(all_results, aes(x = time, y = W, color = Scenario)) +
  geom_line(size = 1) +
  labs(title = "Effect of Booster Rate for waning sub population on Boosted Individuals (W)",
       x = "Time (days)", y = "Number of Boosted Individuals") +
  theme_minimal()

