#

#
# Reaction--Diffusion system of Gierer and Meinhardt taken from
#  Edelstien-Keshet book 
#
#             2              (    2         2 )
#  da        a               ( d a        da  )
#  -- =  -------- - mu*a + Da( ---    +   --- )
#  dt    h(1+k*a^2)          (    2         2 )
#                            ( dx         dy  )
#
#
#                       2
#  dh    2             d h
#  -- = a  - nu*h + Dh --- 
#  dt                    2
#                      dx
#
#
#  homogeneous equilibrium
#        nu           nu
#   a = -----,   h = -----
#                       2
#        mu           mu
#
#  is stable when  nu > mu, without diffusion
#
#  
##########################################################################


# parameters
#
mu = 0.7     # self activation rate
nu = 1       # self inhibition rate
ks = 0.2     # saturation constant for activator
ks = 0.0

Da = 0.5    # diffusion coefficent of activator
Dh = 25     # diffusion coefficent of inhibitor

Nt = 50000   # number of time steps to take
dt = 0.005   # time step

tplot = 100    # make a plot every this many steps
tpause = 0.001 # system pause in animation

Nx = 100     # number of space steps 
dx = 1.0     # length of space step

# initial conditions -- equilibrium for k=0 or solve for equlibrium, then perturb
#
#a0 = nu/mu    
#h0 = nu/mu^2

# solve for equlibrium
#
sseqn=function(a)(nu/(1+ks*a*a) - mu*a)
a0 = uniroot(sseqn,lower=0,upper=nu/mu)
a0 = a0$root
h0 = a0^2/nu

A = matrix(a0,Nx,Nx) 
H = matrix(h0,Nx,Nx)
Anew = A
Hnew = H

# perturbation of uniform state
#
A = A+matrix(runif((Nx*Nx),-0.01,0.01),Nx,Nx) # perturb


# make reaction functions 
#
Ra = function(x,y) (x  /(y*(1+ks*x*x)) - mu)*x
Rh = function(x,y) (x^2/y - nu)*y 

# per time step hoping rates from diffusion
#
pa = dt*Da/dx^2
ph = dt*Dh/dx^2


# initilize matrix corresponding to the diffusion operator -- no flux conditions
#
Adiff = matrix(0,Nx,Nx)
Adiff[1,1]     = -1.0
Adiff[1,2]     =  1.0
Adiff[Nx,Nx-1] =  1.0
Adiff[Nx,Nx]   = -1.0
for (j in 2:(Nx-1) ){
  Adiff[j,j-1] =  1.0
  Adiff[j,j  ] = -2.0
  Adiff[j,j+1] =  1.0
}


# begin main loop in time
#
for( n in 1:Nt ){
  
  Anew = A + dt*Ra(A,H) + pa*(Adiff%*% A + A %*%Adiff)
  Hnew = H + dt*Rh(A,H) + ph*(Adiff%*% H + H %*%Adiff)
  
  A = Anew
  H = Hnew
  
  if(n %% tplot == 0) {
    
    # animate
    #
    zmax = max(A)
    zmin = min(A)    
    filled.contour(A,zlim=c(zmin,zmax),color.palette=heat.colors)
    title(sprintf("time=%g",n*dt))
    
    # pause
    #
    Sys.sleep(tpause)
  }
  
}  # end time loop

