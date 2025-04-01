# randomwalk_many.R -- simulate many random walks on a 1D lattice at
# the same time
#

K = 1000000      # number of particles in the system 
Nt = 100    # number of steps to take
p = 0.40    # probability of hoping left or right


X = array(0,c(Nt,K)) # array of positions at all time levels
T = 0:(Nt-1)         # array of time values


# loop in time
#
for (n in 1:(Nt-1)){
  r = runif(K,0,1)
  X[n+1,] = X[n,] + (r<p)*1.0 + (r>(1-p))*(-1.0) 
}


# plot all particle trajectories
#
#matplot(T,X,ylim=max(abs(X))*c(-1,1))
#matplot(X,T,xlim=max(abs(X))*c(-1,1))

# plot the particle distribution
#
#hist(X[Nt,],breaks=20,freq=FALSE)

# animate the histogram
#
# maxX =max(abs(X))
# nbins = 20
# bins = seq(-maxX,maxX,2*maxX/nbins)
# for (n in 1:(Nt-1)){
#  hist(X[n,],breaks=bins,freq=FALSE,
#       xlim=c(-maxX,maxX),ylim=c(0,0.2),
#       main=sprintf('Step=%i',n))
#  Sys.sleep(0.1)
# }


# plot the mean displacement as a function of time
#
Xmean = rowMeans(X)
XAbsMean = rowMeans(abs(X))
XSqr = rowMeans(X^2)
plot(T,XAbsMean,type='b',ylim=max(abs(XAbsMean))*c(0,1))
#plot(Xmean,T,type='b',xlim=max(abs(Xmean))*c(-1,1))

