# parameters
r=0.99
timeSteps = 30

# vector where population values will be stored
N=numeric(length=timeSteps)

# initial condition
N[1]=1

# use formula to calculate population at each time step
for(i in 2:timeSteps)
{
  N[i]=r*(1-N[i-1])*N[i-1]
}

# plot population N vs generation k
timeVector=seq(1,timeSteps,1)
plot(timeVector,N,type='b',col='blue',xlim=c(0,timeSteps),c(0,1))