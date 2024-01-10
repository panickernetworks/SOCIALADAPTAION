par(mfrow = c(4, 2), mar = c(4, 5, 4, 2))
r=read.csv("countries-aggregated.csv")

Ge=data.frame(r$Confirmed[54673:55233],r$Recovered[54673:55233])
Is=data.frame(r$Confirmed[69361:69921],r$Recovered[69361:69921])
It=data.frame(r$Confirmed[70177:70737],r$Recovered[70177:70737])
Jp=data.frame(r$Confirmed[71809:72369],r$Recovered[71809:72369])
Po=data.frame(r$Confirmed[115057:115618],r$Recovered[115057:115618])
AR=data.frame(r$Confirmed[6529:7089],r$Recovered[6529:7089])
CD=data.frame(r$Confirmed[26929:27489],r$Recovered[26929:27489])


#canada
ICD=as.numeric(as.vector(CD$r.Confirmed.26929.27489.))
RCD=as.numeric(as.vector(CD$r.Recovered.26929.27489.))
PCD=ICD-RCD#prevalence

#armenia
IAR=as.numeric(as.vector(AR$r.Confirmed.6529.7089.))
RAR=as.numeric(as.vector(AR$r.Recovered.6529.7089.))
PAR=IAR-RAR #prevalence

#germany
IG=as.numeric(as.vector(Ge$r.Confirmed.54673.55233.))
RG=as.numeric(as.vector(Ge$r.Recovered.54673.55233.))
PG=IG-RG

#israel
IIs=as.numeric(as.vector(Is$r.Confirmed.69361.69921.))
RIs=as.numeric(as.vector(Is$r.Recovered.69361.69921.))
PIs=IIs-RIs

#italy
IIt=as.numeric(as.vector(It$r.Confirmed.70177.70737.))
RIt=as.numeric(as.vector(It$r.Recovered.70177.70737.))
PIt=IIt-RIt

#japan
IJp=as.numeric(as.vector(Jp$r.Confirmed.71809.72369.))
RJp=as.numeric(as.vector(Jp$r.Recovered.71809.72369.))
PJp=IJp-RJp

#poland
IPo=as.numeric(as.vector(Po$r.Confirmed.115057.115618.))
RPo=as.numeric(as.vector(Po$r.Recovered.115057.115618.))
PPo=IPo-RPo
SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N -5   # Initial susceptible number
  I[1] <- 5  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}

#####################################################################

plot(c(125,500),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(125:500,PAR[125:500]/(max(PAR[125:500])),col="red",pch=16,cex=0.5)
lines(125:500,PAR[125:500]/(max(PAR[125:500])),col="red",lwd=0.7)



SIR1 = SIRDT(N = 300, tmax = 376, b = 0.02, f = 2,0.4,0.15)
SIR2 = SIRDT(N = 300, tmax = 376, b = 0.02, f = 3,0.4,0.15)
SIR3 = SIRDT(N = 300, tmax = 376, b = 0.02, f = 4,0.4,0.15)
SIR4 = SIRDT(N = 300, tmax = 376, b = 0.02, f = 5,0.4,0.15)
SIR5 = SIRDT(N = 300, tmax = 376, b = 0.02, f = 6,0.4,0.15)


lines(125:500,SIR1/max(SIR1),col="blue",lwd=1)
lines(125:500,SIR2/max(SIR2),col="green4",lwd=1)
lines(125:500,SIR3/max(SIR3),col="black",lwd=1)
lines(125:500,SIR4/max(SIR4),col="brown",lwd=1)
lines(125:500,SIR5/max(SIR5),col="violet",lwd=1)
text(465,0.88,"ARMENIA",cex=1.1)
# text(360, 0.92, expression(paste("i", c[1],"= 0.15")),cex=1.1)  # for c1
# text(360, 0.8, expression(paste("i", c[2],"= 0.005")),cex=1.1)  # for c2


plot(c(270,550),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(270:550,PAU[270:550]/(max(PAU[270:550])),col="red",pch=16,cex=0.5)
lines(270:550,PAU[270:550]/(max(PAU[270:550])),col="red",lwd=0.7)



SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N -10   # Initial susceptible number
  I[1] <- 10  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}

legend("topright", legend=c("2", "3","4","5","6"),
       col=c("blue", "green4","black","brown","violet")
       , lty=1, lwd=1.9, cex=0.9, 
       title="SA Factor(f)")

SIR1 = SIRDT(N = 300, tmax = 281, b = 0.02, f = 2,0.4,0.2)
SIR2 = SIRDT(N = 300, tmax = 281, b = 0.02, f = 3,0.4,0.2)
SIR3 = SIRDT(N = 300, tmax = 281, b = 0.02, f = 4,0.4,0.2)
SIR4 = SIRDT(N = 300, tmax = 281, b = 0.02, f = 5,0.4,0.2)
SIR5 = SIRDT(N = 300, tmax = 281, b = 0.02, f = 6,0.4,0.2)


lines(270:550,SIR1/max(SIR1),col="blue",lwd=1)
lines(270:550,SIR2/max(SIR2),col="green4",lwd=1)
lines(270:550,SIR3/max(SIR3),col="black",lwd=1)
lines(270:550,SIR4/max(SIR4),col="brown",lwd=1)
lines(270:550,SIR5/max(SIR5),col="violet",lwd=1)
text(470,0.95,"AUSTRIA",cex=1.1)


plot(c(50,550),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(50:550,PCD[50:550]/(max(PCD[50:550])),col="red",pch=16,cex=0.5)
lines(50:550,PCD[50:550]/(max(PCD[50:550])),col="red",lwd=0.7)

SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N    # Initial susceptible number
  I[1] <- 1  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}


SIR1 = SIRDT(N = 300, tmax = 501, b = 0.02, f = 2,0.4,0.1)
SIR2 = SIRDT(N = 300, tmax = 501, b = 0.02, f = 3,0.4,0.1)
SIR3 = SIRDT(N = 300, tmax = 501, b = 0.02, f = 4,0.4,0.1)
SIR4 = SIRDT(N = 300, tmax = 501, b = 0.02, f = 5,0.4,0.1)
SIR5 = SIRDT(N = 300, tmax = 501, b = 0.02, f = 6,0.4,0.1)


lines(50:550,SIR1/max(SIR1),col="blue",lwd=1)
lines(50:550,SIR2/max(SIR2),col="green4",lwd=1)
lines(50:550,SIR3/max(SIR3),col="black",lwd=1)
lines(50:550,SIR4/max(SIR4),col="brown",lwd=1)
lines(50:550,SIR5/max(SIR5),col="violet",lwd=1)
text(520,0.96,"CANADA",cex=1.05)
# text(200, 1.2, expression(paste("i", c[1],"= 0.1")),cex=1.2)  # for c1
# text(350, 1.2, expression(paste("i", c[2],"= 0.005")),cex=1.2)  # for c2
# 







# 
# 
# plot(c(50,550),c(0,1.3),xlab="Day from Jan 22, 2020",
#      ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
# points(50:550,PCD[50:550]/(max(PCD[50:550])),col="red",pch=16,cex=0.5)
# lines(50:550,PCD[50:550]/(max(PCD[50:550])),col="red",lwd=0.7)
# 



plot(c(250,550),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(250:550,PG[250:550]/max(PG[280:550]),col="red",pch=16,cex=0.5)
lines(250:550,PG[250:550]/max(PG[280:550]),col="red",lwd=0.7)
SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 5  # Initial susceptible number
  I[1] <- 5  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}


# Example usage:
SIR1 = SIRDT(N = 300, tmax = 301, b = 0.02, f = 2,0.4,0.2)
SIR2 = SIRDT(N = 300, tmax = 301, b = 0.02, f = 3,0.4,0.2)
SIR3 = SIRDT(N = 300, tmax = 301, b = 0.02, f = 4,0.4,0.2)
SIR4 = SIRDT(N = 300, tmax = 301, b = 0.02, f = 5,0.4,0.2)
SIR5 = SIRDT(N = 300, tmax = 301, b = 0.02, f = 6,0.4,0.2)



lines(250:550,SIR1/max(SIR1),col="blue",lwd=1)
lines(250:550,SIR2/max(SIR2),col="green4",lwd=1)
lines(250:550,SIR3/max(SIR3),col="black",lwd=1)
lines(250:550,SIR4/max(SIR4),col="brown",lwd=1)
lines(250:550,SIR5/max(SIR5),col="violet",lwd=1)
text(518,0.93,"GERMANY",cex=1)

plot(c(220,475),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(220:475,PIs[220:475]/(max(PIs[220:475])),col="red",pch=16,cex=0.5)
lines(220:475,PIs[220:475]/(max(PIs[220:475])),col="red",lwd=0.7)



SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 15  # Initial susceptible number
  I[1] <- 15  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}

SIR1 = SIRDT(N = 300, tmax = 256, b = 0.02, f = 2,0.4,0.2)
SIR2 = SIRDT(N = 300, tmax = 256, b = 0.02, f = 3,0.4,0.2)
SIR3 = SIRDT(N = 300, tmax = 256, b = 0.02, f = 4,0.4,0.2)
SIR4 = SIRDT(N = 300, tmax = 256, b = 0.02, f = 5,0.4,0.2)
SIR5 = SIRDT(N = 300, tmax = 256, b = 0.02, f = 6,0.4,0.2)



lines(220:475,SIR1/max(SIR1),col="blue",lwd=1)
lines(220:475,SIR2/max(SIR2),col="green4",lwd=1)
lines(220:475,SIR3/max(SIR3),col="black",lwd=1)
lines(220:475,SIR4/max(SIR4),col="brown",lwd=1)
lines(220:475,SIR5/max(SIR5),col="violet",lwd=1)
text(450,0.92,"ISRAEL",cex=1.1)



SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 10  # Initial susceptible number
  I[1] <- 10  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}

SIR41 = SIRDT(N = 300, tmax = 231, b = 0.02, f = 2,0.4,0.2)
SIR42 = SIRDT(N = 300, tmax = 231, b = 0.02, f = 3,0.4,0.2)
SIR43 = SIRDT(N = 300, tmax = 231, b = 0.02, f = 4,0.4,0.2)
SIR44 = SIRDT(N = 300, tmax = 231, b = 0.02, f = 5,0.4,0.2)
SIR45 = SIRDT(N = 300, tmax = 231, b = 0.02, f = 6,0.4,0.2)



plot(c(270,500),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(270:500,PIt[270:500]/max(PIt[270:500]),col="red",pch=16,cex=0.3)
lines(270:500,PIt[270:500]/max(PIt[270:500]),col="red",lwd=0.7)

lines(270:500,SIR41/max(SIR41),col="blue",lwd=1)
lines(270:500,SIR42/max(SIR41),col="green4",lwd=1)
lines(270:500,SIR43/max(SIR41),col="black",lwd=1)
lines(270:500,SIR44/max(SIR41),col="brown",lwd=1)
lines(270:500,SIR45/max(SIR41),col="violet",lwd=1)
text(470,0.92,"ITALY",cex=1.2)


Jp=data.frame(r$Confirmed[71809:72369],r$Recovered[71809:72369])
IJp=as.numeric(as.vector(Jp$r.Confirmed.71809.72369.))
RJp=as.numeric(as.vector(Jp$r.Recovered.71809.72369.))
PJp=IJp-RJp

SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N-20   # Initial susceptible number
  I[1] <- 20  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}


plot(c(330,525),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(330:525,PJp[330:525]/(max(PJp[330:525])),col="red",pch=16,cex=0.5)
lines(330:525,PJp[330:525]/(max(PJp[330:525])),col="red",lwd=0.7)


SIR1 = SIRDT(N = 300, tmax = 196, b = 0.02, f = 2,0.4,0.2)
SIR2 = SIRDT(N = 300, tmax = 196, b = 0.02, f = 3,0.4,0.2)
SIR3 = SIRDT(N = 300, tmax = 196, b = 0.02, f = 4,0.4,0.2)
SIR4 = SIRDT(N = 300, tmax = 196, b = 0.02, f = 5,0.4,0.2)
SIR5 = SIRDT(N = 300, tmax = 196, b = 0.02, f = 6,0.4,0.2)



lines(330:525,SIR1/max(SIR1),col="blue",lwd=1)
lines(330:525,SIR2/max(SIR2),col="green4",lwd=1)
lines(330:525,SIR3/max(SIR3),col="black",lwd=1)
lines(330:525,SIR4/max(SIR4),col="brown",lwd=1)
lines(330:525,SIR5/max(SIR5),col="violet",lwd=1)
text(510,0.92,"JAPAN",cex=1.1)

SIRDT= function(N, tmax, b, f,beta,UT) {
  l = 1
  S = c()  # List to store Susceptible values
  I = c()  # List to store Infected values
  R = c()  # List to store Recovered values
  t <- seq(0, tmax, 1)  # Time vector
  
  S[1] <- N - 2  # Initial susceptible number
  I[1] <- 2  # Initial Infected number
  R[1] <- 0      # Initial Recovered number
  
  rho1 = N / (l ** 2)          # Density without social distancing
  rho2 = N / ((l * f) ** 2)    # Density with social distancing
  inf = c()                     # List to store daily infections
  gamma = 0.05      # Recovery rate
  rho = rho1    # Initialize density
  
  for (h in seq_len(length(t) - 1)) {
    dS <- S[h] * (1 - (1 - beta) ** (pi * b * b * rho * I[h] / N))
    S[h + 1]  <- S[h] - dS
    I[h + 1]  <- I[h] + dS - gamma * I[h]
    R[h + 1]  <- R[h] + gamma * I[h]
    
    sumI = sum(I[h + 1]) / N
    
    # Check upper and lower thresholds for adapting density
    if (sumI > UT) { #upper threshold = 0.2
      rho = rho2  # Apply social distancing
    } else if (sumI < 0.005) { #lower threshold = 0.05
      rho = rho1  # No social distancing
    }
    
    inf = c(inf, I[h + 1])  # Append daily infections to the list
  }
  
  # Return a list of S, I, and R time-series
  return(inf)
}

# par(mfrow = c(3, 2), mar = c(5, 5, 2, 2))

# par(mfrow = c(4, 2), mar = c(4, 5, 4, 2))


plot(c(255,500),c(0,1),xlab="Day from Jan 22, 2020",
     ylab="Prevalence",cex.lab=1.5,cex.axis=1.5,type="n",cex.title=1.5)
points(255:500,PPo[255:500]/max(PPo[225:500]),col="red",pch=16,cex=0.5)
lines(255:500,PPo[255:500]/max(PPo[225:500]),col="red",lwd=0.7)

# legend("topright", legend=c("2", "3","4","5","6"),
#        col=c("blue", "green4","black","brown","violet")
#        , lty=1, lwd=3, cex=0.9, horiz = TRUE,
#        title="SA Factor(f)")

SIR41 = SIRDT(N = 300, tmax = 246, b = 0.02, f = 2,0.4,0.2)
lines(255:500,SIR41/max(SIR41),col="blue",lwd=1)



SIR42 = SIRDT(N = 300, tmax = 246, b = 0.02, f = 3,0.4,0.2)
SIR43 = SIRDT(N = 300, tmax = 246, b = 0.02, f = 4,0.4,0.2)
SIR44 = SIRDT(N = 300, tmax = 246, b = 0.02, f = 5,0.4,0.2)
SIR45 = SIRDT(N = 300, tmax = 246, b = 0.02, f = 6,0.4,0.2)

# lines(255:500,pPOL[255:500],col="red",lwd=0.7)

# plot(c(250,500),c(0,300),type="n")
# max(SIR41)
max(SIR41)

lines(255:500,SIR41/max(SIR41),col="blue",lwd=1)
lines(255:500,SIR42/max(SIR42),col="green4",lwd=1)
lines(255:500,SIR43/max(SIR43),col="black",lwd=1)
lines(255:500,SIR44/max(SIR44),col="brown",lwd=1)
lines(255:500,SIR45/max(SIR45),col="violet",lwd=1)
text(480,0.9,"POLAND",cex=1.1)
