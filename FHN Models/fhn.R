rm(list=ls())

# Parameters
dt <- 0.1
Nt <- 5000
r  <- 0.1
K  <- 1.0
a  <- 0.1
I0 <- 0.0
b  <- 1
ep <- 0.005
v0 <- 0.05

# Define functions
f_v <- function(v, w) {
  r * v * (1 - v/K) * (v/a - 1) - b * w + I0
}
f_w <- function(v, w) {
  ep * (v - w)
}

# Initialize vectors
v <- numeric(Nt + 1)
w <- numeric(Nt + 1)
t <- seq(0, Nt * dt, by = dt)
v[1] <- v0
w[1] <- 0  # initial w

# Time evolution loop
for (n in 1:Nt) {
  v[n+1] <- v[n] + dt * f_v(v[n], w[n])
  w[n+1] <- w[n] + dt * f_w(v[n], w[n])
}

# Plot the time series
plot(t, v, type = 'l', col = 'blue', lwd = 2,
     xlab = "Time", ylab = "State", main = "FHN Model (Time Series)")
lines(t, w, col = 'red', lwd = 2)
legend("topright", legend = c("v", "w"), col = c("blue", "red"), lty = 1, lwd = 2)
