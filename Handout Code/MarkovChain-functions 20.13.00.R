# simulator for discrete-time Markov chains in discrete time

MC=function(A,X0,t.end=100,reps=1,output=FALSE){
  k=dim(A)[1]
  X=matrix(NA,t.end,reps)
  X[1,]=X0
  for(i in 1:reps){
    for(t in 2:t.end){
      X[t,i]=sample.int(n=k,size=1,prob=A[,X[t-1,i]])
    }
  }
  matplot(X,type="l",col=topo.colors(reps),lty=1,lwd=2,bty="n")
  if(output)return(X)
}

# A for the SIS model with gamma=1

A.SIS=function(beta=2,N=15){
  A=matrix(NA,N+1,N+1)
  for(i in 1:(N+1)){
    A[,i]=dbinom(0:N,size=N+1-i,prob=1-exp(-beta*(i-1)/N))
  }
  return(A)
}

# A for the wright fisher model
# N is population size
A.wf=function(N=5){
  A=matrix(NA,N+1,N+1)
  for(i in 1:(N+1)){
    A[,i]=dbinom(0:N,size=N,prob=(i-1)/N)
  }
  return(A)
}

# the following creates random walk matrices of different types. inputs are 
# n number of locations
# d fraction of individuals moving 
# right fraction of moving individuals going to the right
A.walk=function(n=10,d=0.5,right=0.5){
  A=matrix(0,n,n)
  for(i in 1:(n-1)){
    A[i,i]=1-d
    A[i,i+1]=d*(1-right)
    A[i+1,i]=d*right
  }
  A[1,1]=1-d*right
  A[n,n]=1-d*(1-right)
  return(A)
}


# best response with errors A
# for simulating the number of 
# individuals playing strategy 1
A.best.response=function(population.size=10,payoffs=cbind(c(4,0),c(0,1)),error.rate=0.01){
  A=matrix(0,population.size+1,population.size+1)
  for(i in 1:(population.size-1)){
    one.payoff=(i/population.size)*payoffs[1,1]+(1-i/population.size)*payoffs[1,2]
    two.payoff=(i/population.size)*payoffs[2,1]+(1-i/population.size)*payoffs[2,2]
    best.response=sign(one.payoff-two.payoff)
    if(best.response!=0){
    A[i+1+best.response,i+1]=1-error.rate
    A[i+1-best.response,i+1]=error.rate}else{
      A[i+2,i+1]=0.5
      A[i,i+1]=0.5
    }
  }
  i=0
  one.payoff=(i/population.size)*payoffs[1,1]+(1-i/population.size)*payoffs[1,2]
  two.payoff=(i/population.size)*payoffs[2,1]+(1-i/population.size)*payoffs[2,2]
  if(one.payoff<two.payoff){A[1,1]=1-error.rate;A[2,1]=error.rate}else{A[1,1]=error.rate;A[2,1]=1-error.rate}
  i=population.size
  one.payoff=(i/population.size)*payoffs[1,1]+(1-i/population.size)*payoffs[1,2]
  two.payoff=(i/population.size)*payoffs[2,1]+(1-i/population.size)*payoffs[2,2]
  if(one.payoff>two.payoff){A[i+1,i+1]=1-error.rate;A[i,i+1]=error.rate}else{A[i+1,i+1]=error.rate;A[i,i+1]=1-error.rate}
  return(A)
}

# gambler's ruin
A.gamblers.ruin=function(N=5,p=0.5){
P=matrix(0,N+1,N+1)
p=0.5
for(i in 2:N)P[i,i+1]=p
for(i in 2:(N))P[i,i-1]=1-p
P[1,1]=1
P[N+1,N+1]=1
return(t(P))}