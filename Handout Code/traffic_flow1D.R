# 1D agent-based model of traffic flow on a circle
# Optimal velocity model
library(gifski)
N = 50 # Number of cars
K = 500 # Number of time steps: Make < 500 if making gif
dt = 0.5 # Length of time step
L = 25 # Circumference of track
a = 0.1  # Sensitivity: lower corresponds to a slow reaction time

# Plotting commands
plot_number = 10 # Plot every plot_number time steps
video_flag = 1 # 1 makes a gif (longer), 0 plots final image (shorter)

# Setup for time evolutionx``
X = matrix(0,N,K+1)
V = matrix(0,N,K+1) # Velocity

X0 = seq(from=0,to=L,length.out=N+1)  # Set initial positions
X[,1] = X0[1:N] 
X[1,1] = X[1,1] + 0.1

for (j in 1:K){
  
  distCars = (c(X[2:N,j],X[1,j])-X[,j])%%L # Distance between each car
  optV = tanh(distCars-2) +tanh(2)
  V[,j+1] =V[,j] + a*dt*(optV - V[,j])

  X[,j+1] = X[,j] + dt*V[,j+1] 
  X[,j+1] = X[,j+1]%%L  # Keep cars on a circular track
  
}
  

if (video_flag){
  png(file="example%02d.png", width=900, height=900)
  par(bg="white")
  
  for (j in seq(from=1,to=(K+1), by = plot_number)){
    xx = L*cos(2*pi*X[,j]/L)/(2*pi)  # Put on a ring of circumference L
    yy = L*sin(2*pi*X[,j]/L)/(2*pi)
    rbPal = colorRampPalette(c('red','blue'))
    myCols = rbPal(10)[as.numeric(cut(V[,j],breaks = seq(0,2,0.2)))]
    plot(xx,yy,pch=20,xlab='',ylab='',col=myCols,cex=3)
    points(xx[1],yy[1],col='green',cex=3,pch=20)
    
  }
  dev.off()
  png_files <- list.files( pattern = ".*png$", full.names = TRUE)
  gifski(png_files, gif_file = "animation.gif", width = 800, height = 600, delay = 0.1)

 file.remove(list.files(pattern=".png"))
}else{
  xx = L*cos(2*pi*X[,K]/L)/(2*pi)  # Plot final time step on a ring of circumference L
  yy = L*sin(2*pi*X[,K]/L)/(2*pi)
  rbPal = colorRampPalette(c('red','blue'))
  myCols = rbPal(10)[as.numeric(cut(V[,j],breaks = seq(0,2,0.2)))]
  plot(xx,yy,pch=20,xlab='',ylab='',col=myCols,cex=3)
  points(xx[1],yy[1],col='green',cex=3,pch=20)
}

