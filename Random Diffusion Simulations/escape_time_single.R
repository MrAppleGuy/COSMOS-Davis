#
# escape_time_single -- perform random walk and compute the escape
#                       for a single particle by performing many trials
#
NumTrials = 5000    # number of trials to perform
MaxTime   = 100000  # max length of time
p = 0.5             # probability of jumping

xmin = 0            # location of left most cell
xmax = 1            # location of right most cell
Xstart = xmin       # cell location to start

T = numeric(NumTrials) # vector of exit times


# loop over trials
#
for (k in 1:NumTrials){

    # initialize
    #    
    X = Xstart  

    # begin loop in time
    #
    for (n in 1:MaxTime){

        # perform single step of random walk
        #
        r = runif(1,0,1)
 	      X = X + (r<p)*1.0 + (r>(1-p))*(-1.0) 
 
        # keep X>=xmin
        #
        X = X + (X<xmin)*1

	      # check for exit
        #
	      if( X >=xmax ){
	        T[k]=n
	        cat(sprintf("trial = %5i   T=%9.3f    mean(T)= %9.3f\n",k,T[k],mean(T[1:k])))
	        break	
        }	
        
    } # end time loop
    
} # end loop over trials


