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
v0 <- 0.2

# Define model functions
f_v <- function(v, w) {
  r * v * (1 - v/K) * (v/a - 1) - b * w + I0
}
f_w <- function(v, w) {
  ep * (v - w)
}

# Initialize time series
v <- numeric(Nt + 1)
w <- numeric(Nt + 1)
t <- seq(0, Nt*dt, by = dt)
v[1] <- v0
w[1] <- 0  # assume initial w = 0

# Time evolution
for (n in 1:Nt) {
  v[n+1] <- v[n] + dt * f_v(v[n], w[n])
  w[n+1] <- w[n] + dt * f_w(v[n], w[n])
}

# Plot time series
par(mfrow = c(1,2))
plot(t, v, type = 'l', col = 'blue', lwd = 2,
     xlab = "Time", ylab = "v", main = "v over Time")
lines(t, w, col = 'red', lwd = 2)
legend("topright", legend = c("v", "w"), col = c("blue", "red"), lty = 1)

# Phase plane with nullclines
vv <- seq(-0.4, 1.2, by = 0.01)
# v-nullcline: dv/dt = 0 => w = [r*v*(1-v/K)*(v/a-1) + I0] / b
w_nullcline <- r * vv * (1 - vv/K) * (vv/a - 1) / b + I0/b
w_identity <- vv  # from dw/dt = ep*(v-w)=0 => v = w

plot(vv, w_nullcline, type = 'l', lty = 5, lwd = 3, col = 'cyan3',
     xlab = "v", ylab = "w", xlim = c(-0.4, 1.2), ylim = c(-0.05, 0.2),
     main = "Phase Plane & Nullclines")
lines(vv, w_identity, lty = 5, lwd = 3, col = 'coral3')
lines(v, w, col = 'black', lwd = 2)
legend("topright", legend = c("v-nullcline", "w-nullcline", "Trajectory"),
       col = c("cyan3", "coral3", "black"), lty = c(5,5,1), lwd = c(3,3,2))
