# Load libraries
library(deSolve)
library(ggplot2)

# Define the diphtheria model
SIRBVVW_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + A + R + V1 + V2 + W
    lambda <- (beta_h * (I + eta * A) / N) + (beta0 * B / (kappa + B))
    lambda_v1 <- (1 - xi_v_1) * lambda
    lambda_w  <- (1 - xi_w) * lambda
    
    dS  <- Lambda - lambda * S - mu * S - nu_s * S
    dE  <- lambda * S + lambda_v1 * V1 + lambda_w * W - sigma * E - mu * E
    dI  <- sigma * rho * E - gamma * I - delta * I - mu * I
    dA  <- sigma * (1 - rho) * E - tau_A * A - mu * A
    dR  <- gamma * I + tau_A * A - mu * R - zeta_R * R
    dB  <- xi_I * I + pi_A * A + pi_I * I - mu_B * B + alpha
    dV1 <- nu_s * S - mu * V1 - psi_v_1 * V1 - lambda_v1 * V1 - rho_v_1 * V1
    dV2 <- psi_v_1 * V1 - rho_v_2 * V2 - mu * V2 + rho_W * W
    dW  <- zeta_R * R - mu * W - rho_W * W - lambda_w * W + rho_v_2 * V2
    
    list(c(dS, dE, dI, dA, dR, dB, dV1, dV2, dW))
  })
}

# Initial state
initial_state <- c(
  S = 92846408, E = 0, I = 21, A = 500, R = 3,
  B = 500, V1 = 57040000, V2 = 0, W = 0
)

# Fixed parameters
base_params <- c(
  gamma = 0.0061884, sigma = 0.285714, rho = 0.2, tau_A = 0.071429,
  kappa = 500000, Lambda = 38, beta_h = 0.00003, beta0 = 0.234,
  psi_v_1 = 0, mu = 0.0003465, rho_v_2 = 0.0083, nu_s = 1 / 385,
  zeta_R = 1 / 15, xi_I = 0.8, alpha = 0.014, mu_B = 0.0345,
  delta = 0.00022187, xi_v_1 = 0.9917, xi_w = 0.8925,
  eta = 0.5, pi_I = 0.8, pi_A = 0.9
)

# Booster scenarios
booster_ages <- c(1, 5, 10)
coverage_levels <- c(0.3, 0.6, 0.8)
times <- seq(0, 155, by = 1)

# Storage
results_df <- data.frame()

# Run scenarios
for (age in booster_ages) {
  for (cov in coverage_levels) {
    rho_v_1 <- cov / (age * 52)
    scenario_params <- c(base_params, rho_v_1 = rho_v_1, rho_W = 0.05)
    
    out <- ode(y = initial_state, times = times,
               func = SIRBVVW_model, parms = scenario_params,
               atol = 1e-6, rtol = 1e-6)
    
    total_I <- sum(out[, "I"])
    
    results_df <- rbind(results_df, data.frame(
      Booster_Age = age,
      Coverage = cov,
      Total_I = total_I
    ))
  }
}

# Baseline: Age 10, 30% coverage
baseline <- results_df[results_df$Coverage == 0.3 & results_df$Booster_Age == 10, "Total_I"]

# % Cases Averted with 4-digit precision
results_df$Cases_Averted_Pct <- round((baseline - results_df$Total_I) / baseline * 100, 4)

# Scenario label
results_df$Scenario <- paste0("Age ", results_df$Booster_Age, " (", results_df$Coverage * 100, "%)")

# Reorder levels: Age 1, 5, 10
results_df$Scenario <- factor(results_df$Scenario, 
                              levels = c(
                                "Age 1 (30%)", "Age 1 (60%)", "Age 1 (80%)",
                                "Age 5 (30%)", "Age 5 (60%)", "Age 5 (80%)",
                                "Age 10 (30%)", "Age 10 (60%)", "Age 10 (80%)"
                              )
)

# Plot
ggplot(results_df, aes(x = Scenario, y = Cases_Averted_Pct, fill = factor(Booster_Age))) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(Cases_Averted_Pct, "%")), vjust = -0.5, size = 3.5) +
  scale_fill_brewer(palette = "Set1", name = "Booster Age") +
  theme_minimal() +
  labs(
    title = "Impact of Booster Timing and Coverage on Diphtheria Cases Averted",
    x = "Booster Scenario (Age & Coverage)",
    y = "Percentage of Cases Averted"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

