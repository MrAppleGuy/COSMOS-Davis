# simulator for the SIS model 
# I(n+1)=I(n)+p*C*I(n)*(1-I(n))-gamma*I(n)
# choose parameters
C=1 # contact rate per time step 
p=0.5 # transmission probability
gamma=0.1 # removal probability
Tf=15 # length of run 
I0=seq(0.1,0.9,length=10) # vector of initial infected frequencies


###################
# simulation code #
###################
cols=topo.colors(length(I0))
f=function(I)I+p*C*I*(1-I)-gamma*I # the update function for the recursion I(t)=f(I(t-1))
n=length(I0) # number of initial conditions
I.mat=matrix(NA,Tf,n) # matriI to hold the output of the simulations. Columns are different initial conditions and rows are time
I.mat[1,]=I0 # first row has the initial frequencies
for(n in 2:Tf){ # loop to simulate the model
 I.mat[n,]=f(I.mat[n-1,]) # vectorized updates i.e. simulate for all initial frequencies at once
  }

matplot(1:Tf,I.mat,type="l",lwd=3,lty=1,col=cols,bty="n",xlab="time",ylab="fraction infected",ylim=c(0,1)) # plot the matriI of simulations using the "matriI plot" command
