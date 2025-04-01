# a simulator for the stochastic difference equation
# x[n+1] = x[n]*r[n+1]*exp(-a*x[n])
# where a is non-negative (strength of competition)
# r[n] is randomly drawn from some list of numbers
# r[n] is randomly drawn from some list on numbers - "fitness"

# definitions
x0=5 # initial density
rs=c(4,0.2) # list of r values
a=0.01 # strength of competition
finalN=1000 # length of run
howMany = 500 # number of runs

# set up & simulate
xMat=matrix(NA,finalN,howMany)
xMat[1,]=x0
for(i in 2:finalN)
{
  r=sample(rs,howMany,TRUE)
  xMat[i,]=xMat[i-1,]*r*exp(-a*xMat[i-1,])
}
matplot(xMat,type="l")