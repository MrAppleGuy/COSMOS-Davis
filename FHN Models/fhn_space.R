#
# numerically solve FitzHugh-Nagumo equations on a 1D line with no flux boundaries
#    dv/dt = r*v*(1-v/K)*(v/a-1) - b*w + I0 + D v_{xx}
#    dw/dt = ep*(v-w)
#

# set parameters
# --------------
dt     = 0.1        # time step
Nt     = 4000       # number of time steps

dx     = 1.0        # space space
Nx     = 100        # number of space steps
D      = 1.0        # diffusion coefficient
p      = D*dt/dx^2  # fraction jumping left and right one cell per time step

r     = 0.1         # low v decay rate
K     = 1.0         # carrying capacity
a     = 0.1         # threshold
v0    = 0.5         # initial condition at left edge

ep    = 0.005       # recovery rate (inverse time scale)
b     = 1           # rate w lowers v
I0    = 0           # input current

tplot  = 50         # output time steps 
framepause = 0.1  # time to pause between plots

#  dv/dt = fv(v,w)
#  dw/dt = fw(v,w)
#
fv=function(v,w) r*v*(1-v/K)*(v/a-1) - b*w + I0
fw=function(v,w) ep*(v-w)

# initialize storage for solution 
#--------------------------------
v   = numeric(Nx)
w   = numeric(Nx)

vnew = numeric(Nx)
wnew = numeric(Nx)

# time and space variables
#
t   = (0:Nt)*dt
x  = ((1:Nx)-0.5)*dx

# set initial conditions
# ----------------------
v[1] = v0
v[1:2] = v0

# run simulation
# --------------

for (n in 1:Nt){
  
  # update left boundary points
  #
  vnew[1] = v[1] + p*(v[2]-v[1]) + dt*fv(v[1],w[1])
  wnew[1] = w[1]                 + dt*fw(v[1],w[1])

  # loop over interior points
  #
  for (j in 2:(Nx-1)){
    vnew[j]=v[j] + dt*fv(v[j],w[j]) + p*(v[j-1]-2*v[j]+v[j+1])
    wnew[j]=w[j] + dt*fw(v[j],w[j])
  }

  # update the right boundary point
  #
  vnew[Nx]=v[Nx]+p*(v[Nx-1]-v[Nx]) + dt*fv(v[Nx],w[Nx])
  wnew[Nx]=w[Nx]                   + dt*fw(v[Nx],w[Nx])

  # replace with new values
  #
  v = vnew
  w = wnew

  # make a plot of the solution
  #
  if(n %% tplot == 0) {
    plot(x,v,type='b',col='blue',xlim=c(0,Nx),ylim=c(-0.5,1),xlab='space')
    lines(x,w, type='b',col='red')
    text(0.85*Nx,1.0,"time = ")
    text(0.93*Nx,1.0,n*dt)

    # pause
    #
    Sys.sleep(framepause)

  }
  
}  # end loop in time

