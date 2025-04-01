rm(list=ls())

# Parameters
dt <- 0.1
Nt <- 4000
Nx <- 100
dx <- 1.0
D <- 1.0
p <- D * dt / dx^2

# FHN parameters
r  <- 0.1
K  <- 1.0
a  <- 0.1
I0 <- 0
b  <- 1
ep <- 0.005
v0 <- 0.5

tplot <- 50
framepause <- 0.1

# Define model functions
f_v <- function(v, w) {
  r * v * (1 - v/K) * (v/a - 1) - b * w + I0
}
f_w <- function(v, w) {
  ep * (v - w)
}

# Initialize state vectors
v <- rep(0, Nx)
w <- rep(0, Nx)
v[1] <- v0

x <- ((1:Nx) - 0.5) * dx

# Time evolution
for (n in 1:Nt) {
  v_new <- numeric(Nx)
  w_new <- numeric(Nx)
  
  # Left boundary (no flux)
  v_new[1] <- v[1] + dt * f_v(v[1], w[1]) + p * (v[2] - v[1])
  w_new[1] <- w[1] + dt * f_w(v[1], w[1])
  
  # Interior points
  for (j in 2:(Nx-1)) {
    v_new[j] <- v[j] + dt * f_v(v[j], w[j]) + p * (v[j-1] - 2*v[j] + v[j+1])
    w_new[j] <- w[j] + dt * f_w(v[j], w[j])
  }
  
  # Right boundary (no flux)
  v_new[Nx] <- v[Nx] + dt * f_v(v[Nx], w[Nx]) + p * (v[Nx-1] - v[Nx])
  w_new[Nx] <- w[Nx] + dt * f_w(v[Nx], w[Nx])
  
  v <- v_new
  w <- w_new
  
  if (n %% tplot == 0) {
    plot(x, v, type = 'b', col = 'blue', ylim = c(-0.5, 1),
         xlab = "Space", ylab = "State",
         main = sprintf("FHN 1D: time = %.2f", n * dt))
    lines(x, w, type = 'b', col = 'red')
    Sys.sleep(framepause)
  }
}
