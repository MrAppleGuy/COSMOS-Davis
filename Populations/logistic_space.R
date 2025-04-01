#
# simulate logistic equation with spatial diffusion
#


# set parameters
# --------------

dt     = 0.1        # time step
dx     = 1.0        # space space
Nx      = 100        # number of space steps
Nt      = 1000       # number of time steps
tplot  = 25         # output time steps 

D      = 1.0        # diffusion coefficient
p      = D*dt/dx^2  # fraction jumping left and right one cell per time step

r     = 0.5         # low v growth rate 
K     = 1.0         # carrying capacity
v0    = 0.01        # initial amount at the left edge
threshold = 0.2
count=1
framepause = 0.01   # time to pause between plots

# logistic growth function
#
fv=function(v) r*v*(1-v/K)


# need p<=1/2 for stability stop if violated
#
if( p > 0.5 ) stop( sprintf("p=%g, p must be less than 0.5, reduce the time step",p ) )

# initialize storage for solution
#--------------------------------
v   = numeric(Nx)
vnew= numeric(Nx)

# time and space variables
#
t   = (0:Nt)*dt
x  = ((1:Nx)-0.5)*dx

# set initial conditions
# ----------------------
v[1] = v0

# plot the initial condition
#
t=0
plot(dx*(1:Nx),v,type='b',col='blue',xlim=c(0,dx*Nx),ylim=c(0,1),xlab='space')

# long pause before main loop to show initial data
#
Sys.sleep(3)



# run simulation
# --------------

for (nt in 1:Nt){
  continue=TRUE
  
  # update left boundary point
  #
  vnew[1]=v[1]+p*(v[2]-v[1]) + dt*fv(v[1])

  # loop over interior points
  #
  for (j in 2:(Nx-1)){
	  vnew[j]=v[j]+p*(v[j-1]-2*v[j]+v[j+1]) + dt*fv(v[j])
  }
  
  # update the right boundary point
  #
  vnew[Nx]=v[Nx]+p*(v[Nx-1]-v[Nx]) + dt*fv(v[Nx])
  
  # replace u with the new values
  #
  v=vnew
  
  
  # make a plot of the solution
  #
  if(nt %% tplot == 0){
    plot(x,v,type='b',col='blue',xlim=c(0,dx*Nx),ylim=c(0,1),xlab='space')
    text(0.85*Nx,1.0,"time = ")
    text(0.93*Nx,1.0,t)
  }
  vcopy=v
  while(continue)
  {
    if(vcopy[count] >= 0.9)
    {
      print(nt)
      vcopy = vcopy[-count]
      count = count + 1
      continue = FALSE
    }
    continue = FALSE
  }
  # pause
  #
  Sys.sleep(framepause)
                                      
  
}  # end loop in time





