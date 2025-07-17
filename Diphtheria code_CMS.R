library(deSolve)

# Extended SEIR-VVW model with Exposed (E) and Asymptomatic (A)
SIRBVVW_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + A + R + V1 + V2 + W
    
    lambda <- (beta_h * (I + eta * A) / N) + (beta0 * B / (kappa + B))
    lambda_v1 <- (1 - xi_v_1) * lambda
    lambda_w  <- (1 - xi_w) * lambda
    
    dS <- Lambda - lambda * S - mu * S - nu_s * S
    dE <- lambda * S + lambda_v1 * V1 + lambda_w * W - sigma * E - mu * E
    dI <- sigma * rho * E - gamma * I - delta * I - mu * I
    dA <- sigma * (1 - rho) * E - tau_A * A - mu * A
    dR <- gamma * I + tau_A * A - mu * R - zeta_R * R
    dB <- xi_I * I + pi_A * A + pi_I * I - mu_B * B + alpha
    dV1 <- nu_s * S - mu * V1 - psi_v_1 * V1 - lambda_v1 * V1 - rho_v_1 * V1
    dV2 <- psi_v_1 * V1 - rho_v_2 * V2 - mu * V2 + rho_W * W
    dW  <- zeta_R * R - mu * W - rho_W * W - lambda_w * W + rho_v_2 * V2
    
    list(c(dS, dE, dI, dA, dR, dB, dV1, dV2, dW))
  })
}

initial_state <- c(
  S = 92846408, E = 0, I = 21, A = 500, R = 3,
  B = 500, V1 = 57040000, V2 = 0, W = 0
)

parameters <- c(
  gamma = 0.0061884, sigma = 0.285714, rho = 0.2, tau_A = 0.071429,
  kappa = 500000, Lambda = 38, beta_h = 0.00003, beta0 = 0.234,
  psi_v_1 = 0, mu = 0.0003465, rho_v_2 = 0.0083, nu_s = 1 / 385,
  zeta_R = 1 / 15, xi_I = 0.8, alpha = 0.014, mu_B = 0.0345,
  delta = 0.00022187, xi_v_1 = 0.9917, xi_w = 0.8925,
  rho_v_1 = 0.000548, rho_W = 0.05, eta = 0.5, pi_I = 0.8, pi_A = 0.9
)

times <- seq(0,155, by = 1)

solution <- ode(y = initial_state, times = times,
                func = SIRBVVW_model, parms = parameters,
                atol = 1e-6, rtol = 1e-6)
solution_df <- as.data.frame(solution)

plot(solution_df$time, solution_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Time (days)", ylab = "Population",
     main = "Infected and Asymptomatic Cases Over Time")
lines(solution_df$time, solution_df$A, col = "blue", lwd = 2)
legend("topright", legend = c("Infected (I)", "Asymptomatic (A)"),
       col = c("red", "blue"), lty = 1, lwd = 2)
plot(solution)

# Solve the ODE as usual

# Rename columns to full names for clarity
colnames(solution_df) <- c("Time", 
                           "Susceptible", 
                           "Exposed", 
                           "Infected", 
                           "Asymptomatic", 
                           "Recovered", 
                           "Bacteria Load", 
                           "Vaccinated (1st Dose)", 
                           "Vaccinated (2nd Dose)", 
                           "Boosted (Waning Immunity)")

# Plot all compartments
matplot(solution_df$Time, solution_df[, -1], type = "l", lty = 1, lwd = 2,
        col = rainbow(ncol(solution_df) - 1),
        xlab = "Time (days)", ylab = "Population", 
        main = "Diphtheria Model with Full Compartment Labels")

legend("topright", legend = colnames(solution_df)[-1], 
       col = rainbow(ncol(solution_df) - 1), lty = 1, lwd = 2, cex = 0.6)

# Rename columns to full names for clarity
colnames(solution) <- c("Time", 
                           "Susceptible", 
                           "Exposed", 
                           "Infected", 
                           "Asymptomatic", 
                           "Recovered", 
                           "Bacteria Load", 
                           "Vaccinated (1st Dose)", 
                           "Vaccinated (2nd Dose)", 
                           "Boosted (Waning Immunity)")
plot(solution)
base_params <- c(
  gamma = 0.0061884, sigma = 0.285714, rho = 0.2, tau_A = 0.071429,
  kappa = 500000, Lambda = 38, beta_h = 0.00003, beta0 = 0.234,
  psi_v_1 = 0, mu = 0.0003465, rho_v_2 = 0.0083, nu_s = 1 / 385,
  zeta_R = 1 / 15, xi_I = 0.8, alpha = 0.014, mu_B = 0.0345,
  delta = 0.00022187, xi_v_1 = 0.9917, xi_w = 0.8925,
  rho_v_1 = 0.000548, eta = 0.5, pi_I = 0.8, pi_A = 0.9
)
# Parameters (excluding rho_W which we will vary)
#base_params <- c(
  #gamma = 0.0061884, kappa = 500000, Lambda = 38, beta_h = 0.00003,
 # beta0 = 0.234, psi_v_1 = 0, mu = 0.0003465, rho_v_2 = 0.0083,
 # nu_s = 1 / 385, zeta_R = 1 / 15, xi_I = 0.8, alpha = 0.014,
 # mu_B = 0.0345, delta = 0.00022187, xi_v_1 = 0.9917, xi_w = 0.8925
#)

# Time
times <- seq(0, 1095, by = 1)



# Booster rates to test
booster_rates <- c(0, 0.05, 0.2, 0.6)

colors <- c("black", "blue", "green", "red")

# Store output
results <- list()

for (i in seq_along(booster_rates)) {
  param_set <- c(base_params, rho_W = booster_rates[i])
  out <- ode(y = initial_state, times = times, func = SIRBVVW_model,
             parms = param_set, atol = 1e-6, rtol = 1e-6)
  results[[i]] <- as.data.frame(out)
  results[[i]]$Scenario <- paste0("Booster rate = ", booster_rates[i])
}
Total0 = sum(results[[1]]$I)
Total1 = sum(results[[2]]$I)
Total2 = sum(results[[3]]$I)
Total3 = sum(results[[4]]$I)

Cases_A0 <- Total0 - Total0
Cases_A1 <- Total0 - Total1
Cases_A2 <- Total0 - Total2
Cases_A3 <- Total0 - Total3

# Convert to percentage
Cases_A0_pct <- (Cases_A0 / Total0) * 100
Cases_A1_pct <- (Cases_A1 / Total0) * 100
Cases_A2_pct <- (Cases_A2 / Total0) * 100
Cases_A3_pct <- (Cases_A3 / Total0) * 100

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
       x = "Time (days)", y = "Number of Boosted Individuals") 
