# Load libraries
library(ggplot2)
library(reshape2)


#INITIAL CONDITIONS
X = 1000000 # total timesteps
dt = 0.001 #delta t, the size of timestep
a = 0.21 #constant water influx
m = 0.09 #constant biomass loss
k = 100 # number of spatial locations
L = 50 # length of space

#total number of x values
timejump = 1000

N_all = matrix(NA,X/timejump,k)

dx = L/k #don't change this


#calculate a1
a1 = c(-1 * ((-1 * a + sqrt((a^2) -4*(m^2))) / 2), (a + sqrt((a^2) -4*(m^2))) / (2*m))

v = 182.5
WaterProb = 0.4 #water diffusion 
d = 0.8 #plant diffusion

#set up matrices
W = matrix(NA, X, k)
N = matrix(NA, X, k)

#starting water values
W[1,1:k] = 1

#STARTING biomass
N[1,]=seq(0,5, length=k)

#function for N
g = function(w, n){
  return(w*(n^2) - m * n)
}

#function for W
f = function(w, n){
  return(a - w - w*(n^2))
}



#################WATER ADVECTION##############################
# array of rotating values that will make matrix (used to make seq2)

#Probabilities

S = -1*WaterProb
R = WaterProb
wat1 = c(R, S, 0)
for (i in 4:(k)) {
  wat1[i] = 0
}

# initializes the array that will form the matrix
wat2 = numeric(k^2)

#prepares seq1 for use
for (i in 1:(k-1)){
  wat1 = c(tail(wat1, 1), head(wat1, -1))
}

# Fill seq2 with shifted versions of seq1
for (i in 1:k) {
  for (n in 1:k) {
    wat2[(i - 1) * k + n] = wat1[n]
  }
  wat1 = c(tail(wat1, 1), head(wat1, -1)) #wat1 is shifted every value of i
}


#Create matrix from seq2
mat.water = matrix(wat2, nrow = k, byrow = TRUE)


#################BIOMASS DIFFUSION############################

#plant matrix variables
p = -2*d #controls probability (don't worry about this)
l = d #controls probability (don't worry about this either)


# array of rotating values that will make matrix (used to make seq2)
seq1 = c(l, p, l)
for (i in 4:(k)) {
  seq1[i] = 0
}

# initializes the array that will form the matrix
seq2 = numeric(k^2)

#prepares seq1 for use
for (i in 1:(k-1)){
  seq1 = c(tail(seq1, 1), head(seq1, -1))
}

# Fill seq2 with shifted versions of seq1
for (i in 1:k) {
  for (n in 1:k) {
    seq2[(i - 1) * k + n] = seq1[n]
  }
  seq1 = c(tail(seq1, 1), head(seq1, -1)) #seq1 is shifted every value of i
}

#Create matrix from seq2
mat.veg = matrix(seq2, nrow = k, byrow = TRUE)

###################MAIN LOOP####################################
#water and plant biomass recursion
for (i in 1:(X-1)){
  W[i+1,] = W[i,] + dt*(f(W[i,], N[i,])) + (v*dt/dx)*mat.water%*%W[i,]
  N[i+1,] = N[i,] + dt*(g(W[i,], N[i,])) + (dt/dx^2)*mat.veg%*%N[i,]
  
    if (i%%timejump == 0){
      N_all[i/timejump,] = N[i,]
    }
}

#calculate a1
a1 = c(-1 * ((-1 * a + sqrt((a^2) -4*(m^2))) / 2), (a + sqrt((a^2) -4*(m^2))) / (2*m))
print(a1)

#Print plots
par(mfrow = c(1,2))

# matplot(N,type="l",lty=1,main=paste0('a=', a, ' m=',m))
# abline(h = a1[2],lty=2)
# matplot(W,type="l",lty=1,main=paste0('X=',X))
# abline(h=a1[1],lty=2)

# matplot(t(N),type="l",lty=1,col = topo.colors(X),xlab="spatial steps")
print("FOR LOOP complete")

# Transform the matrix into a data frame
data = melt(N_all)
colnames(data) = c("X1", "X2", "value")

# Create the heatmap
ggplot(data, aes(x = X1, y = X2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "beige", high = "darkgreen") +
  theme_minimal() +
  labs(x = "T", y = "X", fill = "Biomass")
