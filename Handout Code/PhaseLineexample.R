# Start by initializing the 1D system with r=1,a=3,K=10.
#install.packages("deSolve", dependencies = TRUE)
library(deSolve)
Fx=function(x,h,a,k,r){r*x*(1-x/k)*(x-a)-h*x}

#prepare to plot the phase line curves to study them. 
x=seq(0,12,.1)
#h is the harvesting/stocking parameter to play with
h=0
r=1
k=5
a=3
y=Fx(x,h,a,k,r)
plot(x,y,type="l",ylim=c(-5,5))
abline(h=0,col="black",)

#Now For the same values of h, lets see what the real dynamics are like! 

dx=function(t,x,h){list(x*(1-x/8)*(x-3)-h*x)}

t=seq(0,50)

#pick an initial condition
x0=6
#run the simulation
out = ode(x0,t,dx,h)
plot(out)
