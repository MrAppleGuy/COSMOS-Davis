# program to simulate diffusion on a 1d grid
# with no flux (no leak) boundary conditions


# set parameters
# --------------

dt     = 0.1        # time step
dx     = 1.0        # space space
Nx     = 100        # number of space steps
Nt     = 1000       # number of time steps
tplot  = 10         # output every this many time steps 
D      = 5        # diffusion coefficient
p      = D*dt/dx^2  # fraction jumping left and right one cell per time step
Xj     = ((1:Nx)-1/2)/Nx # position of the box
thresholdLine = a40
threshold = 0.2
flag = TRUE

ictype = 4          #  used to pick the initial conditions
k = 4               #    1 = step in the middle
                    #    2 = random
                    #    3 = elevated cosine, k periods of oscillation

                    
framepause = 0.01   # time to pause between plots

# need p<=1/2 for stability stop if violated
#
if( p > 0.5 ) stop( sprintf("p=%g, p must be less than 0.5, reduce the time step",p ) )


# initialize storage of solution at current and next time step 
#
u   = numeric(Nx)
unew= numeric(Nx)

# set initial conditions
# ----------------------
# ictype is used to pick the initial conditions
#
# (1) step in middle
# (2) random
# (3) sinusoidal (k=frequency control)
#
if (ictype==1) {           # step
  halfwidth = Nx*0.1
  j1=(Nx/2)-halfwidth
  j2=(Nx/2)+halfwidth
  u[j1:j2]=1 
} else if (ictype==2){     # random
  u=runif(Nx,min=0,max=1)
} else if(ictype==3) {     # sinusoidal
  u=(1+cos(2*k*pi*Xj))/2 
} else if(ictype==4) # start all on one side
{
  u[0:20]=1
}

# plot the initial condition, and pause
#
t=0
plot(u,type='b',col='blue',xlim=c(0,Nx),ylim=c(0,1),xlab='space',ylab='concentration')
text(0.85*Nx,1.0,"time = ")
text(0.93*Nx,1.0,t)
Sys.sleep(3)



# run simulation
# --------------

for (n in 1:Nt){
  
  # update left boundary point
  #
  unew[1]=u[1]+p*(u[2]-u[1])

  # loop over interior points
  #
  for (j in 2:(Nx-1)){
	  unew[j]=u[j]+p*(u[j-1]-2.0*u[j]+u[j+1]) 
  }
  
  # update the right boundary point
  #
  unew[Nx]=u[Nx]+p*(u[Nx-1]-u[Nx])
  
  # replace u with the new values
  #
  u=unew
  
  # make a plot of the solution
  #
  if(n %% tplot == 0) {
    t=n*dt
    plot(u,type='b',col='blue',xlim=c(0,Nx),ylim=c(0,1),xlab='space',ylab='concentration')
    text(0.85*Nx,1.0,"time = ")
    text(0.93*Nx,1.0,t)
  }
  
  if(flag && u[thresholdLine] >= threshold)
  {
    print("passed!")
    print(n)
    flag = FALSE
  }
  # pause
  #
  Sys.sleep(framepause)
                                 
}  # end loop in time





