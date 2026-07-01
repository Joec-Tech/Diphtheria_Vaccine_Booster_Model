#SEIR model
library(deSolve)
SIRBVVW_model<- function(t,state,parameters) {
  S<- state["S"]
  I<- state["I"]
  R<- state["R"]
  B<- state ["B"]
  V1<- state ["V1"]
  V2<- state ["V2"]
  W<- state ["W"]
  
  beta0 <- parameters  ["beta0"]
  gamma <- parameters ["gamma"]
  Lambda <- parameters ["Lambda"]
  rho_v_2 <- parameters["rho_v_2"]
  zeta_R <- parameters["zeta_R"]
  xi_I <- parameters["xi_I"]
  alpha <- parameters ["alpha"]
  rho_W <- parameters["rho_W"]
  nu_s <- parameters["nu_s"]
  psi_v_1 <- parameters["psi_v_1"]
  mu <- parameters ["mu"]
  beta_h <- parameters ["beta_h"]
  xi_w <- parameters ["xi_w"]
  xi_v_1 <- parameters ["xi_v_1"]
  delta <- parameters ["delta"]
  mu_B <- parameters ["mu_B"]
  kappa <- parameters ["kappa"]
  N <- parameters["N"]
  N <- S + I + R + V1 + V2
  dS <- Lambda-(beta_h * I/N + (beta0 *B /(kappa + B))) *S - nu_s*S - mu*S
  dI <- (beta_h*I/N + (beta0 * B/(kappa + B)))*S-(1-xi_v_1)*V1 + (1-xi_w)*W - mu* I - delta * I - gamma * I
  dR <- gamma * I - mu*R - zeta_R * R
  dB <- xi_I * I -mu_B * B + alpha
  dV1 <- nu_s * S - mu*V1 - psi_v_1 * V1- (1-xi_v_1)*V1
  dV2 <- psi_v_1 * V1 - rho_v_2 * V2 - mu * V2 + rho_W * W
  dW <- zeta_R * R - mu * W - rho_W * W - (1-xi_w)*W + rho_v_2 * V2
  
  list(c(dS,dI,dR,dB,dV1,dV2,dW))
}

initial_state <- (c(S=6846408,
                    I=10661, 
                    R=10292, 
                    B=500, 
                    V1=272917, 
                    V2=0, 
                    W =0)) 

parameters <- c(gamma=0.0061884,
                kappa=500000 , 
                Lambda = 38,
                beta_h= 0.00003,
                beta0= 0.234,
                psi_v_1 = 0 ,
                mu= 0.0003465,
                rho_v_2 = 0.0083,
                nu_s = 1/385,
                zeta_R = 1/15,
                xi_I= 0.8,
                alpha = 0.014 ,
                rho_W = 0,
                mu_B= 0.0345,
                delta =0.00022187,
                xi_v_1 =0.9917,
                xi_w=0.8925)
times<- seq(0, 153, by = 0.1)

solution <- ode (y = initial_state, times= times, func = SIRBVVW_model, parms = parameters ,atol = 1e-6, rtol = 1e-6)
# Convert to data frame and plot
solution_df <- as.data.frame(solution)
print(solution_df)
