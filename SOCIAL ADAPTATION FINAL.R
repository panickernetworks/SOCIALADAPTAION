library(igraph)
library(spatstat)
library(spatgraphs)
###############################################################################
#N = NUMBER OF AGNETS
#b = CHARACTERISTIC TRANSMISSION RADIUS
#tmax = Epidemic Run Time
# SIR_Maxmobility: Simulate Infectious Disease Spread with Mobility

SIR_Maxmobility <- function(N, b, tmax) {
  
  # Model Parameters
  beta <- 0.5     # Infection rate
  gamma <- 0.05   # Recovery rate
  
  # Initialize Agent Locations
  xa <- runif(N, 0, 1)  # x-coordinates
  ya <- runif(N, 0, 1)  # y-coordinates
  dfa <- data.frame(xa, ya)  # Data frame with agent locations
  
  # Create Initial Geometric Graph
  g1 <- spatgraph(dfa, "geometric", par = b)
  p1 <- g1[]
  gr1 <- graph_from_adj_list(p1$edges)  # Initial graph without mobility
  
  # Initialize States
  I <- matrix(rep(0, N), nrow = N, ncol = 1)  # Infected
  S <- matrix(rep(1, N), nrow = N, ncol = 1)  # Susceptible
  R <- matrix(rep(0, N), nrow = N, ncol = 1)  # Recovered
  
  # Initial Infection
  I1 <- sample(1:N, size = 1)
  I[I1, 1] <- 1
  S[I1, 1] <- 0
  
  t <- 1  # Time step
  sus <- c()  # Susceptible count over time
  inf <- c()  # Infected count over time
  rec <- c()  # Recovered count over time
  
  gr <- gr1  # Initial graph for mobility
  S_list <- list(S)  # List to store susceptible matrices
  I_list <- list(I)  # List to store infected matrices
  R_list <- list(R)  # List to store recovered matrices
  while (t <= tmax) {
    t <- t + 1
    
    # Calculate Infections
    infneigh <- gr[] %*% I[, t - 1]
    pinf <- 1 - (1 - beta) ** infneigh
    
    newI <- c()
    for (i in 1:N) {
      new <- rbinom(1, S[, t - 1][i], pinf[i])
      newI <- c(newI, new)
    }
    
    # Calculate Recoveries
    newR <- rbinom(N, I[, t - 1], gamma)
    
    # Update States
    nextS <- S[, t - 1] - newI
    nextI <- I[, t - 1] + newI - newR
    nextR <- R[, t - 1] + newR
    
    I <- cbind(I, nextI)
    S <- cbind(S, nextS)
    R <- cbind(R, nextR)
    
    # Update Counts
    inf <- c(inf, sum(nextI))
    sus <- c(sus, sum(nextS))
    rec <- c(rec, sum(nextR))
    
    S_list[[t]] <- nextS
    I_list[[t]] <- nextI
    R_list[[t]] <- nextR
    
    # Update Agent Locations and Graph
    xa <- runif(N, 0, 1)  # Random mobility of all agents
    ya <- runif(N, 0, 1)
    dfa <- data.frame(xa, ya)
    g1 <- spatgraph(dfa, "geometric", par = b)
    p1 <- g1[]
    gr <- graph_from_adj_list(p1$edges)
  }
  
  return(list(S = S_list, I = I_list, R = R_list, prevalence = inf))
}

SIR_data1 <- SIR_Maxmobility(N = 100, b = 0.05, tmax = 50)
# Access matrices for susceptible, infected, and recovered
#individuals at a specific time step
susceptible_at_t10 <- SIR_data1$S[[10]]
infected_at_t10 <- SIR_data1$I[[10]]
recovered_at_t10 <- SIR_data1$R[[10]]
####################################################################
# SIR_zeromobility: Simulate Infectious Disease Spread with Zero Mobility

SIR_zeromobility <- function(N, b, tmax) {
  
  # Model Parameters
  beta <- 0.5     # Infection rate
  gamma <- 0.05   # Recovery rate
  
  # Initialize Agent Locations
  xa <- runif(N, 0, 1)  # x-coordinates
  ya <- runif(N, 0, 1)  # y-coordinates
  dfa <- data.frame(xa, ya)  # Data frame with agent locations
  
  # Create Initial Geometric Graph
  g1 <- spatgraph(dfa, "geometric", par = b)
  p1 <- g1[]
  gr1 <- graph_from_adj_list(p1$edges)  # Initial graph without mobility
  
  # Initialize States
  I <- matrix(rep(0, N), nrow = N, ncol = 1)  # Infected
  S <- matrix(rep(1, N), nrow = N, ncol = 1)  # Susceptible
  R <- matrix(rep(0, N), nrow = N, ncol = 1)  # Recovered
  
  # Initial Infection
  I1 <- sample(1:N, size = 1)
  I[I1, 1] <- 1
  S[I1, 1] <- 0
  
  t <- 1  # Time step
  sus <- c()  # Susceptible count over time
  inf <- c()  # Infected count over time
  rec <- c()  # Recovered count over time
  
  gr <- gr1  # Initial graph (static geometric graph)
  S_list <- list(S)  # List to store susceptible matrices
  I_list <- list(I)  # List to store infected matrices
  R_list <- list(R)  # List to store recovered matrices
  while (t <= tmax) {
    t <- t + 1
    
    # Calculate Infections
    infneigh <- gr[] %*% I[, t - 1]
    pinf <- 1 - (1 - beta) ** infneigh
    
    newI <- c()
    for (i in 1:N) {
      new <- rbinom(1, S[, t - 1][i], pinf[i])
      newI <- c(newI, new)
    }
    
    # Calculate Recoveries
    newR <- rbinom(N, I[, t - 1], gamma)
    
    # Update states
    nextS <- S[, t - 1] - newI
    nextI <- I[, t - 1] + newI - newR
    nextR <- R[, t - 1] + newR
    
    I <- cbind(I, nextI)
    S <- cbind(S, nextS)
    R <- cbind(R, nextR)
    
    # Update counts
    inf <- c(inf, sum(nextI))
    sus <- c(sus, sum(nextS))
    rec <- c(rec, sum(nextR))
    
    S_list[[t]] <- nextS
    I_list[[t]] <- nextI
    R_list[[t]] <- nextR
  }
  return(list(S = S_list, I = I_list, R = R_list, prevalence = inf))
}

# Get prevalence data for a single trial
SIR_data2 <- SIR_zeromobility(N = 500, b = 0.05, tmax = 150)
# Access matrices for susceptible, infected, and recovered
#individuals at a specific time step
susceptible_at_t10 <- SIR_data2$S[[10]]
infected_at_t10 <- SIR_data2$I[[10]]
recovered_at_t10 <- SIR_data2$R[[10]]
####################################################################
# Define a function named SIR_StaticAdaptation with input parameters N, b, tmax, f, and th
#f  = social adaptation factor
#th = prevalence threshold for adaptation
SIR_StaticAdaptation = function(N, b, tmax, f, th) {
  # Set values for parameters beta and gamma
  beta = 0.5
  gamma = 0.05
  
  # Generate random initial coordinates for nodes
  xa = runif(N, 0, 1)
  ya = runif(N, 0, 1)
  
  # Create a data frame with the generated coordinates
  dfa = data.frame(xa, ya)
  
  # Create a static network graph (g1) without social adaptation
  g1 = spatgraph(dfa, "geometric", par = b)
  
  # Calculate a modified value for parameter b (b2) based on SA factor f
  b2 = b / f
  
  # Create a static network graph (g2) with social adaptation
  g2 = spatgraph(dfa, "geometric", par = b2)
  
  # Extract adjacency lists from the generated graphs
  p1 = g1[]
  p2 = g2[]
  
  # Convert adjacency lists into graph objects
  gr1 = graph_from_adj_list(p1$edges)
  gr2 = graph_from_adj_list(p2$edges)
  
  # Initialize matrices for infected (I), susceptibles (S), and removed (R) individuals
  I = matrix(rep(0, N), nrow = N, ncol = 1)
  S = matrix(rep(1, N), nrow = N, ncol = 1)
  R = matrix(rep(0, N), nrow = N, ncol = 1)
  
  # Select a random individual to be the first infected
  I1 = sample(1:N, size = 1)
  I[I1, 1] = 1
  S[I1, 1] = 0
  
  # Initialize time step
  t = 1
  
  # Initialize vectors to store counts of susceptible, infected, and recovered individuals
  sus = c()
  inf = c()
  rec = c()
  
  # Start with the first graph (no social adaptation)
  gr = gr1
  
  S_list <- list(S)  # List to store susceptible matrices
  I_list <- list(I)  # List to store infected matrices
  R_list <- list(R)  # List to store recovered matrices
  # Main simulation loop
  while (t <= tmax) {
    t = t + 1
    
    # Calculate the number of infected neighbors for each individual
    infneigh = gr[] %*% I[, t - 1]
    
    # Calculate the probability of infection for each individual
    pinf = 1 - (1 - beta)^infneigh
    
    # Initialize a vector to store the new infected individuals
    newI = c()
    
    # Generate new infections for each individual
    for (i in 1:N) {
      new = rbinom(1, S[, t - 1][i], pinf[i])
      newI = c(newI, new)
    }
    
    # Generate new recovered individuals
    newR = rbinom(N, I[, t - 1], gamma)
    
    # Calculate the next values for susceptible, infected, and recovered individuals
    nextS = S[, t - 1] - newI
    nextI = I[, t - 1] + newI - newR
    nextR = R[, t - 1] + newR
    
    # Update the vectors and matrices with new counts
    inf = c(inf, sum(nextI))
    sus = c(sus, sum(nextS))
    rec = c(rec, sum(nextR))
    I = cbind(I, nextI)
    S = cbind(S, nextS)
    R = cbind(R, nextR)
    
    # Check if prevalence exceeds the threshold
    if ((sum(nextI) / N) > th) {
      gr = gr2  # Switch to the graph with social adaptation
    } else {
      gr = gr1  # Use the graph without social adaptation
    }
    S_list[[t]] <- nextS
    I_list[[t]] <- nextI
    R_list[[t]] <- nextR
  }
  return(list(S = S_list, I = I_list, R = R_list, prevalence = inf))
}

SIR_data3=SIR_StaticAdaptation(N=500,b=0.1,tmax=50,f = 4,th = 0.1)
# Access matrices for susceptible, infected, and recovered
#individuals at a specific time step
susceptible_at_t10 <- SIR_data3$S[[10]]
infected_at_t10 <- SIR_data3$I[[10]]
recovered_at_t10 <- SIR_data3$R[[10]]
#########################################################################
# SIR_FullymixedAdaptation with input parameters N, b, tmax, f, and th
#adaptation in fully mixed network
#f  = social adaptation factor
#th = prevalence threshold for adaptation

SIR_FullyMixedAdaptation = function(N, b, tmax, f, th) {
  # Set values for parameters beta and gamma
  beta = 0.5
  gamma = 0.05
  
  # Generate random initial coordinates for nodes
  xa = runif(N, 0, 1)
  ya = runif(N, 0, 1)
  # Create a data frame with the generated coordinates
  dfa = data.frame(xa, ya)
  
  # Create a static network graph (g1) without social adaptation
  g1 = spatgraph(dfa, "geometric", par = b)
  
  # Extract adjacency lists from the generated graph
  p1 = g1[]
  
  # Convert adjacency lists into graph object
  gr1 = graph_from_adj_list(p1$edges)
  
  # Initialize matrices for infected (I), susceptibles (S), and removed (R) individuals
  I = matrix(rep(0, N), nrow = N, ncol = 1)
  S = matrix(rep(1, N), nrow = N, ncol = 1)
  R = matrix(rep(0, N), nrow = N, ncol = 1)
  
  # Select a random individual to be the first infected
  I1 = sample(1:N, size = 1)
  I[I1, 1] = 1
  S[I1, 1] = 0
  
  # Initialize time step
  t = 1
  
  # Initialize vectors to store counts of susceptible, infected, and recovered individuals
  sus = c()
  inf = c()
  rec = c()
  
  # Start with the first graph (no social adaptation)
  gr = gr1
  
  S_list <- list(S)  # List to store susceptible matrices
  I_list <- list(I)  # List to store infected matrices
  R_list <- list(R)  # List to store recovered matrices
  # Main simulation loop
  while (t <= tmax) {
    t = t + 1
    
    # Calculate the number of infected neighbors for each individual
    infneigh = gr[] %*% I[, t - 1]
    
    # Calculate the probability of infection for each individual
    pinf = 1 - (1 - beta)^infneigh
    
    # Initialize a vector to store the new infected individuals
    newI = c()
    
    # Generate new infections for each individual
    for (i in 1:N) {
      new = rbinom(1, S[, t - 1][i], pinf[i])
      newI = c(newI, new)
    }
    
    # Generate new recovered individuals
    newR = rbinom(N, I[, t - 1], gamma)
    
    # Calculate the next values for susceptible, infected, and recovered individuals
    nextS = S[, t - 1] - newI
    nextI = I[, t - 1] + newI - newR
    nextR = R[, t - 1] + newR
    
    # Update the vectors and matrices with new counts
    inf = c(inf, sum(nextI))
    sus = c(sus, sum(nextS))
    rec = c(rec, sum(nextR))
    I = cbind(I, nextI)
    S = cbind(S, nextS)
    R = cbind(R, nextR)
    # Check if prevalence exceeds the threshold
    if ((sum(nextI) / N) > th) {
      # Update Agent Locations 
      # and form graph with social adaptation
      xa <- runif(N, 0, 1)
      ya <- runif(N, 0, 1)
      dfa <- data.frame(xa, ya)
      g1 <- spatgraph(dfa, "geometric", par = b/f)
      p1 <- g1[]
      gr <- graph_from_adj_list(p1$edges)
    } else {
      # Update Agent Locations 
      # and form graph without social adaptation
      xa <- runif(N, 0, 1)  
      ya <- runif(N, 0, 1)
      dfa <- data.frame(xa, ya)
      g1 <- spatgraph(dfa, "geometric", par = b)
      p1 <- g1[]
      gr <- graph_from_adj_list(p1$edges)  # Use the graph without social adaptation
    }
    S_list[[t]] <- nextS
    I_list[[t]] <- nextI
    R_list[[t]] <- nextR
  }
  return(list(S = S_list, I = I_list, R = R_list, prevalence = inf))
}

SIR_data4=SIR_FullyMixedAdaptation(N=500,b=0.05,tmax=50,f = 4,th = 0.1)
# Access matrices for susceptible, infected, and recovered
#individuals at a specific time step
susceptible_at_t10 <- SIR_data4$S[[10]]
infected_at_t10 <- SIR_data4$I[[10]]
recovered_at_t10 <- SIR_data4$R[[10]]
#########################################################################
#function named SIRnumercial with input parameters N, tmax, and b
#meanfield
SIRnumercial_No_SA = function(N, tmax, b) {
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1) # Time vector
  
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  
  rho1 = N      # Density without social distancing
  inf = c()     # List to store daily infections
  
  beta = 0.5       # Transmission rate
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta)**(pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    # Break the loop if the number of Infected individuals becomes negative
    if (I[h + 1] < 0) {
      break
    }
    
    inf = c(inf, I[h + 1]) # Append daily infections to the list
  }
  
  # Return lists for Susceptible, Infected, and Recovered individuals
  return(list(S = S, I = I, R = R, inf))
}

# Example usage:
SIR_data5 <- SIRnumercial_No_SA(N = 500, tmax = 100, b = 0.05)
##########################################################################
SIRnumercial_SA = function(N, tmax, b, f, th) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 1  # Initial susceptible number
  I[1] <- 1      # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  #considering the equation of susceptibles 
  #changing density from rho into  rho/(f)**2 is numerically same as #reducing b into b/f 
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  
  beta = 0.5       # Transmission rate
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta)**(pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    sumI = sum(I[h + 1]) / N
    
    # Break the loop if the number of Infected individuals becomes negative
    if (I[h + 1] < 0) {
      break
    }
    
    if (sumI > 0.1) { #threshold = 0.1
      rho = rho2  # Apply social adaptation
    } else{
      rho = rho1  #else no SA
    }
    
    inf = c(inf, I[h + 1]) # Append daily infections to the list
  }
  
  # Return lists for Susceptible, Infected, and Recovered individuals
  return(list(S = S, I = I, R = R, inf))
}

# Example usage:
SIR_data6 <- SIRnumercial_SA(N = 500, tmax = 100, b = 0.05,f = 4,th =0.1)
##########################################################################
# Define a function named SIRDelay with input parameters N, t, r, l, f, th, and d
#f  = social adaptation factor
#th = prevalence threshold for adaptation
#d = delay
#l = length of patch
SIRnumerical_Delay = function(N, t, b, f, th, d, beta) {
  l = 1
  tmax = t
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, t, 1)  # Time vector
  
  S[1] <- N - 1 # Initial susceptible number
  I[1] <- 1     # Initial Infected number
  R[1] <- 0     # Initial Recovered number
  #considering the equation of susceptibles 
  #changing density from rho into  rho/(f)**2 is numerically same as #reducing b into b/f 
  rho1 = N / (l ** 2)           # Density without social adaptation
  rho2 = N / ((l * f) ** 2)     # Density with social adaptation
  sus = c()                       # List to store susceptible values
  inf = c()                      # List to store infected values
  rec = c()                       # List to store recovered values
  gamma = 0.05      # Recovery rate (gamma)
  rho = rho1    # Initialize density
  # Iterate over time steps
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta)**(pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - I[h] * gamma
    R[h + 1]  <- R[h] + I[h] * gamma
    
    sus = c(sus, S[h + 1])    # Append susceptible count
    inf = c(inf, I[h + 1])  # Append infected count
    rec = c(rec, R[h + 1])    # Append recovered count
    
    if (d == 0) {
      if ((inf[length(inf)] / N) > th) {  # no delay, prevalence>th
        rho = rho2
      } else {
        rho = rho1   #no delay, prevalence<th 
      }
    } else {  # delay
      if (h - 1 < (d)) {
        rho = rho1
      } else {
        v = (inf[h - d] / N)
        if (v > th) {
          rho = rho2
        } else {
          rho = rho1
        }
      }
    }
  }
  
  # Return the list of daily infected counts (prevalence)
  return(inf)
}

# Example usage:
SIR_data7 <- SIRnumerical_Delay(N = 500, t = 100, b = 0.05,
                                f = 4, th = 0.1, d = 3, beta = 0.1)
#####################################################################
# Define a function to get variation of peak infection with delay
#f  = social adaptation factor
#th = prevalence threshold for adaptation
#d = delay
#l = length of patch
getPeakInfectionWithDelay = function(N, t, b, f, th,beta) {
  peak_infections = c()  # List to store peak infection counts for each delay
  for (d in 1:30) {
    infections = SIRnumerical_Delay(N, t, b, f, th, d, beta)  # Get infection counts with the current delay
    peak_infections = c(peak_infections, max(infections))  # Store the peak infection count
  }
  return(peak_infections)
}

# Example usage:
SIR_data8 <- getPeakInfectionWithDelay(N = 500, t = 100, b = 0.03,
                                        f = 4, th = 0.1, beta=0.5)
##############################################################################
# Define a function named SIR1 with input parameters N, tmax, r, and f
SIRnumerical_double_threshold = function(N, tmax, b, f) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 1  # Initial susceptible number
  I[1] <- 1      # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  #considering the equation of susceptibles 
  #changing density from rho into  rho/(f)**2 is numerically same as.   #reducing b into b/f 
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  
  beta = 0.5       # Transmission rate
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > 0.2) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.05) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(list(S = S, I = I, R = R))
}

# Example usage:
SIR_data9 = SIRnumerical_double_threshold(N = 500, tmax = 150, b = 0.05, f = 6)
#####################################################################
