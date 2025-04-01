# a simulator for the stochastic difference equation
# x[n+1]=x[n]*r[n+1]*exp(-a*x[n])
# where a is non-negative (strength of competition)
# x[n] density at time n
# r[n] is randomly draw from some list on numbers

# function to simulate
SPGFunction=function(x0=5,rs=c(0.9,1.5),a=0.001,finalN=100,howMany=10){
# set up and simulate! 
xMat=matrix(NA,nrow=finalN,ncol=howMany) # matrix to hold sims
xMat[1,]=x0 # setting first row to initial condition
for(i in 2:finalN){ # loop for the simulation
  r=sample(rs,size=howMany,replace = TRUE) # randomly choose r multipliers
  xMat[i,]=xMat[i-1,]*r*exp(-a*xMat[i-1,]) # update population densities. 
}
return(xMat)
}
## plot
xMat=SPGFunction()
matplot(xMat,type="l",lty=1)