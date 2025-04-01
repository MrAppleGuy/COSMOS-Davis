# a simulator for the stochastic difference equation
# x[n+1]=x[n]*r[n+1]*exp(-a*x[n])
# where a is non-negative (strength of competition)
# x[n] density at time n
# r[n] is randomly draw from some list on numbers - "fitness"
  
# definitions

x0=5 # initial density
rs=c(4,0.2) # list of r values
a=0.01 # strength of competition
final.n=1000 # length of run
how.many=50 # number of runs

# let the good times roll 
# set up and simulate! 
x.mat=matrix(NA,nrow=final.n,ncol=how.many) # matrix to hold sims
x.mat[1,]=x0 # setting first row to initial condition
for(i in 2:final.n){ # loop for the simulation
  r=sample(rs,size=how.many,replace = T) # randomly choose r multipliers
  x.mat[i,]=x.mat[i-1,]*r*exp(-a*x.mat[i-1,]) # update population densities. 
}

## plot
matplot(x.mat,type="l",lty=1)

