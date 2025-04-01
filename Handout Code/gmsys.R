#
# Reaction--Diffusion system of Gierer and Meinhardt taken from
#  Edelstien-Keshet book page 531.
#
#          2              2
#  da     a              d a
#  -- =  --- - mu*a + Da --- 
#  dt     h                2
#                        dx
#
#
#                        2
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
mu = 1.0     # decay rate rate activator
nu = 1.2     # decay rate rate inhibitor

Da = 1    # diffusion coefficent of activator
Dh = 20    # diffusion coefficent of inhibitor


Nt = 50000   # number of time steps to take
dt = 0.001   # time step

tplot = 100    # make a plot every this many steps
tpause = 0.1 # system pause in animation

Nx = 100     # number of space steps 
dx = 1.0     # length of space step


# initial conditions -- equilibrium, then perturb
#
a0 = nu/mu    
h0 = nu/mu^2
A = array(a0 , Nx) 
H = array(h0 , Nx)
Anew = A
Hnew = H

# perturbation of equilibrium
#
A = A + runif(Nx,-0.01,0.01) # perturb
#k = 1; A = A + 0.01*cos(k*pi*((1:Nx)-1/2)/Nx)

# make reaction functions 
#
Ra = function(x,y) (x  /y - mu)*x
Rh = function(x,y) (x^2/y - nu)*y

# per time step hopping rates from diffusion
#
pa = dt*Da/dx^2
ph = dt*Dh/dx^2

# plot the initial conditions and pause
#
ymax = max( max(A),max(H))
ymax = max(ymax, 1.005*max(a0,h0))
ymin = min( min(A),min(H))
ymin = min(ymin, 0.995*min(a0,h0))
plot(A,type='b',col='red',ylim=c(ymin,ymax),xlab="space",ylab="A in red,   H in black")
lines(H,type='b',col='black')
Sys.sleep(2.0)


# begin main loop in time
#
for( n in 1:Nt ){

    # no flux boundaries
    #
    Anew[1]  = A[1]  + dt*Ra(A[1], H[1] ) + pa * (A[2]   -A[1] )
    Hnew[1]  = H[1]  + dt*Rh(A[1], H[1] ) + ph * (H[2]   -H[1] )
    Anew[Nx] = A[Nx] + dt*Ra(A[Nx],H[Nx]) + pa * (A[Nx-1]-A[Nx])
    Hnew[Nx] = H[Nx] + dt*Rh(A[Nx],H[Nx]) + ph * (H[Nx-1]-H[Nx])

    
    # interior points
    #
    for( j in 2:(Nx-1) ){
        Anew[j] = A[j] + dt*Ra(A[j],H[j]) + pa * (A[j-1]-2.0*A[j]+A[j+1])
        Hnew[j] = H[j] + dt*Rh(A[j],H[j]) + ph * (H[j-1]-2.0*H[j]+H[j+1])
    }

    # update
    #
    A = Anew
    H = Hnew

    # make a plot
    #
    if(n %% tplot == 0) {

      # animate
      #
      ymax = max( max(A),max(H))
      ymax = max(ymax, 1.005*max(a0,h0))
      ymin = min( min(A),min(H))
      ymin = min(ymin, 0.995*min(a0,h0))
      
      plot(A,type='b',col='red',main=sprintf("time=%g",n*dt),xlab="space",ylab="A in red,   H in black",ylim=c(ymin,ymax))
      lines(H,type='b',col='black')
     
      # pause
      #
      Sys.sleep(tpause)
    }
    
}  # end time loop

