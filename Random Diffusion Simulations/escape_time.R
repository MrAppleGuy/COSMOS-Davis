#
# escape_time --- compute the escape time for different length race courses
#

NumTrials = 10000      # number of trials to perform
MaxTime = 100000      # max length of time
p = 0.99             # probability of jumping

xmin = 0        # location of left most cell
Xstart = xmin   # cell location to start
Lmax = 20       # max length of domain to use


EscTime = numeric(Lmax)  # initialize array of escape times
CourseLengths = 1:Lmax   # array of course lengths to run

# loop over different length of domains
#
for(xmax in CourseLengths){

  cat(sprintf('Begin running length=%i for %i trials\n ',xmax,NumTrials))
  
  # initialize particles and individual escape times
  #
  X = array(Xstart,NumTrials)
  T = array(0,NumTrials)

  # loop in time
  #
  for (n in 1:MaxTime){
  
    # perform single step of random walk
    #
    r = runif(NumTrials,0,1)
    X = X + (r<p)*1.0 + (r>(1-p))*(-1.0) 

    # keep X>=xmin
    #
    X = X + (X<xmin)*1

    # check if this is the first exit and record
    #
    T = T + (X>=xmax & T ==0)*n 

    # exit loop if everyone escaped
    #
    if(all(T>0)) break
  
  } #

  EscTime[xmax] = mean(T)
 

}
print(p)
plot(CourseLengths,EscTime,'b',xlab='length',ylab='escape time')
title(main=sprintf('escape time, Ntrials = %i, p=%5.3f',NumTrials,p))
