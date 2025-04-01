rm(list=ls())

# Parameters
dt <- 0.05       # time step
Nt <- 1500       # number of time steps
Nx <- 100        # grid size (Nx x Nx)
dx <- 1.0        # spatial step
D <- 1.0         # diffusion coefficient for v
p <- D * dt / dx^2  # diffusion factor

# FHN model parameters
r  <- 0.1        # rate parameter for v
K  <- 1.0        # carrying capacity
a  <- 0.1        # threshold parameter
I0 <- 0          # external input
b  <- 1          # coupling from w to v
ep <- 0.005      # recovery rate for w

tplot <- 20      # plot interval (time steps)
framepause <- 0.1

# Define model functions
f_v <- function(v, w) {
  r * v * (1 - v/K) * (v/a - 1) - b * w + I0
}
f_w <- function(v, w) {
  ep * (v - w)
}

# Initialize state matrices
v <- matrix(0, nrow = Nx, ncol = Nx)
w <- matrix(0, nrow = Nx, ncol = Nx)
v[1:10, 1:10] <- 0.5  # initial condition patch for v

# Function to create a 1D diffusion matrix (for no-flux boundaries)
create_diffusion_matrix <- function(N, p) {
  A <- matrix(0, nrow = N, ncol = N)
  for (i in 2:(N-1)) {
    A[i, i-1] <- p
    A[i, i]   <- -2 * p
    A[i, i+1] <- p
  }
  A[1,1] <- -p; A[1,2] <- p
  A[N,N-1] <- p; A[N,N] <- -p
  A
}
A_diff <- create_diffusion_matrix(Nx, p)

# Time evolution
for (n in 1:Nt) {
  # Compute diffusion in both dimensions: use separable Laplacian (A*v + v*A)
  v_diff <- A_diff %*% v + v %*% A_diff
  v_new <- v + dt * f_v(v, w) + v_diff
  w_new <- w + dt * f_w(v, w)
  
  v <- v_new
  w <- w_new
  
  if (n %% tplot == 0) {
    image(v, col = topo.colors(100), zlim = c(-0.5, 1.2),
          main = sprintf("FHN 2D: time = %.2f", n * dt),
          xlab = "x", ylab = "y")
    Sys.sleep(framepause)
  }
}
