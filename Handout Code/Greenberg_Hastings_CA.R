# Simulates a 2D cellular automaton in time
#install.packages('plot.matrix')
library(plot.matrix)

N = 100   # Number of time steps
numColumns = 100
numRows = 100
plotlive=0
defib=0
nfib=40
# Set up the state matrices
Z_old = matrix(0,nrow=numRows,ncol=numColumns) 
Z_new = matrix(0,nrow=numRows,ncol=numColumns) 

# Initial condition
Z_old[numColumns/2,(numColumns/2):numColumns] = 1
Z_old[(numColumns/2-1),(numColumns/2):numColumns] = -1

# Function to compute the sum of the neighborhood: Moore Neighborhood
nbhdSum = function(i,j,M) {  # Input the (i,j) coordinates of the cell, matrix M
  n = c(0,0,0,0)
  n[1] = M[j,i-1]
  n[2] = M[j+1,i]
  n[3] = M[j,i+1]
  n[4] = M[j-1,i]
  
  ns = 0
  for (m in 1:length(n)){
    if (n[m] >0){
      ns = ns+1
    }
  }
  return(ns)
}

for (n in 2:(N+1)){  # Evolve in time
  
  for (i in 2:(numColumns-1)){ # Loop over the columns
    
    for (j in 2:(numRows-1)){  # Loop over the rows
      
      ns_x = nbhdSum(i,j,Z_old) # Find the sum of the neighboring cells
      
      if (ns_x > 0 && Z_old[j,i]==0){  # Equilibrium -> Excited
        Z_new[j,i] = 1
      }
    
      if (defib==1 && n==nfib && Z_old[j,i]==0){Z_new[j,i] = 1}
      
      else if ( Z_old[j,i]==1){  # Recovery -> Equilibrium
        Z_new[j,i] = -1
      }
    } # End loop over rows
  }  # End loop over columns
  
 if (n ==(N+1)){  # Plot the last time step
    plot(Z_new,col=c('black','white','red'),border=NA)
    
 }
 if (plotlive==1){ 
  flush.console()
  plot(Z_new,col=c('black','white','red'),border=NA)
  Sys.sleep(.05)}
  
  # Prepare for the next time step by resetting the matrices
  Z_old = Z_new
  Z_new = matrix(0,nrow=numRows,ncol=numColumns) 
  
  
}  # End Time loop





