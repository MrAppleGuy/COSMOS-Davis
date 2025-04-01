timesteps = 10000 # total timesteps
dt = 0.0001 #delta t, the size of timestep
a = 3 #constant water influx
m = 1 #constant biomass loss
k = 1000
L=50
dx=L/k

a1=-1*((-a+sqrt(a^2-4*m^2))/(2))
a2=((a+sqrt(a^2-4*m^2))/(2*m))
W = matrix(NA, timesteps, k)
N = matrix(NA, timesteps, k)
W[1,1:k] = 3
N[1,]=seq(0,5, length=k)
v = 182.5
waterProb = 0.1 #water diffusion 
d = 0.45 #plant diffusion

dwdt = function(w, n){
  return(a-w-w*n^2)
}
dndt = function(w, n){
  return(w*n^2-m*n)
}

A.veg=function(){
  A=matrix(0,k,k)
  for(i in 1:(k-1)){
    A[i,i]=-d*2
    A[i,i+1]=d
    A[i+1,i]=d
  }
  A[1,1]=-d*2
  A[k,k]=-d*2
  A[k,1]=d
  A[1,k]=d
  return(A)
}
AVeg=A.veg()

A.water=function(){
  A=matrix(0,k,k)
  for(i in 1:(k-1)){
    A[i,i]=-waterProb
    A[i+1,i]=waterProb
  }
  A[1,1]=waterProb
  A[k,k]=-waterProb
  A[1,1]=-waterProb
  A[1,k]=waterProb
    
  return(A)
}
AWater=A.water()
for (i in 1:(timesteps-1)){
  W[i+1,] = W[i,] + dt*(dwdt(W[i,], N[i,])) + (v*dt/dx)*AWater%*%W[i,]
  N[i+1,] = N[i,] + dt*(dndt(W[i,], N[i,])) + (dt/dx^2)*AVeg%*%N[i,]
}
print(c(a1,a2))
par(mfrow=c(1,2))
matplot(t(W),type="l",lty=1,col = topo.colors(timesteps),xlab="Space",ylab="Water Concentration (W)")
# abline(h=a1,lty="dashed",lwd=3)
matplot(t(N),type="l",lty=1,col = topo.colors(timesteps),xlab="Space",ylab="Biomass Concentration (N)")
# abline(h=a2,lty="dashed",lwd=3)
