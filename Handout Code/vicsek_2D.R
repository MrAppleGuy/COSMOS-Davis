# Simulate the 2D Vicsek Model
#install.packages("rdist",)
library(rdist)

# System parameters ------
N = 20 # Number of agents
TF = 200 # Final Time
L = 10 # Length of domain
eta = 1  # Magnitude of noise
R = 2.5 # Interaction radius
v0 = 10 # Magnitude of Velocity

video_flag = 1 # 1 = make video, 0 = no video, just plot final time step

dt = 0.1  # Time step
K = TF/dt # Number of time steps

# Matrix set up
theta = matrix(0,N,K+1) # Orientations  
Xpos = matrix(0,N,K+1)  # X positions
Ypos = matrix(0,N,K+1)  # Y positions

# Set the initial conditions
Xpos[,1] = runif(N,min=0,max=L)
Ypos[,1] = runif(N,min=0,max=L)
theta[,1] = runif(N,min=0,max=2*pi)

# Time steps ------
for (k in 1:K){
  
  # Find pairwise distances between all particles
   Rdist = pdist(cbind(Xpos[,k],Ypos[,k]), metric = "euclidean", p = 2)
  
   Rdist[Rdist>R] = 0 # Set values of pairwise distances > R = 0
   Rdist[Rdist>0] = 1 # Set values of pairwise distances <= R = 1
   diag(Rdist) = 1 # Add the diagonal (yourself)
   
  # Update orientation angle to be average of those within neighborhood + noise
  theta[,k+1] = Rdist%*%theta[,k]/rowSums(Rdist) + runif(N,min=-eta/2,max=eta/2)
  
  # Update the positions
  Xpos[,k+1] = Xpos[,k] + v0*dt*cos(theta[,k+1])
  Ypos[,k+1] = Ypos[,k] + v0*dt*sin(theta[,k+1])
  
  # Periodic Boundary Conditions
  Xpos[,k+1] = Xpos[,k+1]%%L
  Ypos[,k+1] = Ypos[,k+1]%%L
}


# Compute the order parameter:
phi = abs(sum(v0*exp(1i*theta[,K+1])))/(N*v0)

print(paste('Order parameter: ', phi))

#################### Make Video -------
if (video_flag){
# install.packages("gganimate")
# install.packages("tidyr")
# install.packages("gifski")
library(gganimate)
library(tidyr)
library(gifski)

Xvals = as.data.frame(Xpos) #change to data frame
Yvals = as.data.frame(Ypos)
Xvals = cbind(seq(1:N),Xvals)
long_X = gather(Xvals, time, measurementX, V1:paste("V",sep="",(K+1))) #convert to long format
long_Y = gather(Yvals, time, measurementY, V1:paste("V",sep="",(K+1))) 
long_X = cbind(long_X, long_Y$measurementY)
long_X$time = rep(1:(K+1), each=N)
colnames(long_X) = c("agent", "time", "X", "Y")


g = ggplot(long_X, aes(x=X, y=Y)) + 
 geom_point() + 
  theme_bw() +
  # gganimate specific bits:
  transition_time(time) +
  ease_aes('linear') 

# Save at gif:
animate(g, fps=5, width = 800, height = 800, renderer = gifski_renderer())
anim_save("agentbased.gif")

file.remove(list.files(pattern=".png"))
}else{
  
  plot(Xpos[,K+1],Ypos[,K+1],xlab='x',ylab='y',xlim=c(0,L),ylim=c(0,L))
  alength =L/50
  arrows(Xpos[,K+1],Ypos[,K+1], x1=Xpos[,K+1]+alength*cos(theta[,K+1]), y1=Ypos[,K+1]+alength*sin(theta[,K+1]), 
         length=0.05, col="Black")
  
}
