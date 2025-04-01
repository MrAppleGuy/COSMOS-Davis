#
# numerically solve FitzHugh-Nagumo equations on a 1D line with no flux boundaries
#    dv/dt = r*v*(1-v/K)*(v/a-1) - b*w + I0 + D v_{xx}
#    dw/dt = ep*(v-w)
#

# set parameters
# --------------
dt     = 0.05        # time step
Nt     = 1500       # number of time steps

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

tplot  = 20         # output time steps 
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
v   =matrix(0,Nx,Nx)
vnew=matrix(0,Nx,Nx)
v[1:10,1:10]=v0
#v[1,1:Nx]=v0

w    = matrix(0,Nx,Nx)
wnew = matrix(0,Nx,Nx)

# initilize matrix corresponding to the diffusion operator -- no flux conditions
#
Adiff = matrix(0,Nx,Nx)
Adiff[1,1]   = -p
Adiff[1,2]   =  p
Adiff[Nx,Nx-1] =  p
Adiff[Nx,Nx]   = -p

for (j in 2:(Nx-1) ){

  Adiff[j,j-1] =    p
  Adiff[j,j  ] = -2*p
  Adiff[j,j+1] =    p
}




# run simulation
# --------------

for (n in 1:Nt){
  
  # update
  #
  vnew = v + fv(v,w) + (Adiff %*% v) + (v %*% Adiff)
  wnew = w + fw(v,w)

  # replace (v,w) with the new values
  #
  v=vnew
  w=wnew
  
  
  # make a plot of the solution
  #
  if(n %% tplot == 0) {
  
    image(v,zlim=c(-0.5,1.2),col=topo.colors(100),main=sprintf('time=%f',n*dt))
    #persp(v,zlim=c(-0.5,1.2))

    # pause
    #
    Sys.sleep(framepause)

  }
  
}  # end loop in time

