#
# simulate diffusion with a source term in the middle and constant decary rate
#


# set parameters
# --------------

dt     = 0.1        # time step
dx     = 1.0        # space space
Nx      = 51        # number of space steps
Js      = 24:26        # source location
Nt      = 1000      # number of time steps
tplot  = 25         # output time steps 

D      = 0.0        # diffusion coefficient
p      = D*dt/dx^2  # fraction jumping left and right one cell per time step

k     = 0.1         # decay rate 
                
framepause = 0.01   # time to pause between plots

# initialize storage for solution
#--------------------------------
v   = numeric(Nx)
vnew= numeric(Nx)

# make a source function
#  with no difussion steady state is just 1
#
S = numeric(Nx)
S[Js] = k

# time and space variables
#
t   = (0:Nt)*dt
x  = ((1:Nx)-0.5)*dx


# run simulation
# --------------

for (nt in 1:Nt){
  
  # update left boundary point
  #
  vnew[1]=v[1]+ dt*(S[1] - k*v[1]) + p*(v[2]-v[1])
  # loop over interior points
  #
  for (j in 2:(Nx-1)){
	  vnew[j]=v[j]+ dt*(S[j] - k*v[j])+ p*(v[j-1]-2*v[j]+v[j+1]) 
  }
  
  # update the right boundary point
  #
  vnew[Nx]=v[Nx]+dt*(S[Nx] - k*v[Nx]) + p*(v[Nx-1]-v[Nx]) 
  
  # replace v with the new values
  #
  v=vnew
  
  # make a plot of the solution
  #
  if(nt %% tplot == 0){
    plot(x,v,type='b',col='blue',xlim=c(0,dx*Nx),ylim=c(0,1),xlab='space')
  }

  # pause
  #
  Sys.sleep(framepause)
                                      
    
  
}  # end loop in time





