  timesteps = 75000
  dt = 0.0001
  a = 4
  m = 2
  
  W = numeric(timesteps)
  N = numeric(timesteps)
  
  #Initial amounts for water and plant biomass
  W[1] = 7
  N[1] = 8
  a1=-1*((-a+sqrt(a^2-4*m^2))/(2))
  a2=((a+sqrt(a^2-4*m^2))/(2*m))
  
  #Function for plant biomass
  dNdt = function(w,n){
    return(w*n^2-m*n)
  }
  
  #Function for water
  dWdt = function(w,n){
    return(a-w-w*n^2)
  }
  
  #main loop
  for(i in 1:(timesteps-1)){
    W[i+1] = W[i]+dt*(dWdt(W[i], N[i])) 
    N[i+1] = N[i]+dt*(dNdt(W[i], N[i]))
  }
  amt_time = 1:length(W) #Number of time steps, used for plotting W and N
  
  par(mfrow=c(1,2)) #Makes 2 graphs instead of just 1
  plot(amt_time,W,type="l",col="red",ylim=c(min(W),max(W)),xlab="Time", ylab="Water", main="a = 2.1, m = 1")
  abline(h=a1,lty="dashed",lwd=3)
  plot(amt_time,N,col="blue",type="l",ylim=c(min(N),max(N)),xlab="Time", ylab="Plant Biomass", main="a = 2.1, m = 1")
  abline(h=a2,lty="dashed",lwd=3)
  # matplot(t(W),type="l",lty=1,col = topo.colors(timesteps))
  # matplot(t(N),type="l",lty=1,col = topo.colors(timesteps))
  
  
