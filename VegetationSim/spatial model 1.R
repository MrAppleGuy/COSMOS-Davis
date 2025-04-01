  X = 2000 # total timesteps
  dt = 0.001 #delta t, the size of timestep
  a = 0.1 #constant water influx
  m = 0.45 #constant biomass loss
  k = 10 # number of spatial locations
  
  v = 182.5
  varray = rep(v, length=k)
  WaterProb = 0.2 #water diffusion 
  d = 0.45 #plant diffusion
  
  #set up matrices
  W = matrix(NA, X, k)
  N = matrix(NA, X, k)
  
  #starting water values
  W[1,1:k] = 1
  
  h = 0.1
  b = 10
  N[1,]=seq(0,5, length=k)
  
  #function for N
  g = function(w, n){
    return(w*(n^2) - m * n)
  }
  
  #function for W
  f = function(w, n){
    return(a - w - w*(n^2))
  }
  
  
  #################WATER DIFFUSION##############################
  # array of rotating values that will make matrix (used to make seq2)
  
  #Probabilities
  
  S = -1*WaterProb
  R = WaterProb
  wat1 = c(0, S, R)
  for (i in 4:(k+1)) {
    wat1[i] = 0
  }
  
  # initializes the array that will form the matrix 
  wat2 = numeric(k^2)
  
  #prepares seq1 for use
  for (i in 1:(k)){
    wat1 = c(tail(wat1, 1), head(wat1, -1))
  }
  
  # Fill seq2 with shifted versions of seq1
  for (i in 1:k) {
    for (n in 1:k) {
      wat2[(i - 1) * k + n] = wat1[n]
    }
    wat1 = c(tail(wat1, 1), head(wat1, -1)) #wat1 is shifted every value of i
  }
  
  #changes first and last values to work
  wat2[1] = -1*WaterProb
  wat2[k^2] = -1*WaterProb
  
  #Create matrix from seq2
  mat.water = matrix(wat2, nrow = k, byrow = TRUE)
  
  
  #################BIOMASS DIFFUSION############################
  
  #plant matrix variables
  d = 0.45 #probability of plant movement
  p = -2*d #controls probability (don't worry about this)
  l = d #controls probability (don't worry about this either)
  
  
  # array of rotating values that will make matrix (used to make seq2)
  seq1 = c(l, p, l)
  for (i in 4:(k+1)) {
    seq1[i] = 0
  }
  
  # initializes the array that will form the matrix 
  seq2 = numeric(k^2)
  
  #prepares seq1 for use
  for (i in 1:(k)){
    seq1 = c(tail(seq1, 1), head(seq1, -1))
  }
  
  # Fill seq2 with shifted versions of seq1
  for (i in 1:k) {
    for (n in 1:k) {
      seq2[(i - 1) * k + n] = seq1[n]
    }
    seq1 = c(tail(seq1, 1), head(seq1, -1)) #seq1 is shifted every value of i
  }
  
  #changes first and last values to be correct
  seq2[1] = -1*l
  seq2[k^2] = -1*l
  
  #Create matrix from seq2
  mat.veg = matrix(seq2, nrow = k, byrow = TRUE)
  ###################MAIN LOOP####################################
  
  #water and plant biomass recursion
  for (i in 1:(X-1)){
    W[i+1,] = W[i,] + dt*(f(W[i,], N[i,])) + varray%*%mat.water
    N[i+1,] = N[i,] + dt*(g(W[i,], N[i,])) + N[i,]%*%mat.veg
  }
  
  matplot(W,type="l",lty=1)
  matplot(t(N),type="l",lty=1,col = topo.colors(X)) #transposes modeling, time is COLOR
  
  
  #calculate a1
  a1 = c(-1 * ((-1 * a + sqrt((a^2) -4*(m^2))) / 2), (a + sqrt((a^2) -4*(m^2))) / (2*m))
  print(a1)
  
  varray
  mat.water
  mat.veg
  
  
