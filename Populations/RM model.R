#install.packages("phaseR", dependencies = TRUE)
library(deSolve)
library(phaseR)

predprey = function(t, y, parameters) {
  x <- y[1]
  y <- y[2]
  r=0.5
  a=4
  h=3
  e=0.5
  d=0.1
 #r <- parameters[1]
  K <- parameters#[2]
 # h <- parameters[3]
 # d <- parameters[4]
 # e <- parameters[5]
  dy <- numeric(2)
  dy[1] <- r*x*(1-x/K) - a*x*y/(1+a*h*x)
  
  dy[2] <- e*a*x*y/(1+a*h*x) - d*y
  list(dy)
}


#change this parameter K below
K=0.35

predprey.nullclines <- nullclines(predprey, xlim = c(0, .3), ylim = c(0, .3),
                                  parameters = K, points = 500,add=FALSE)


predprey.flowField <-
   flowField(predprey, xlim = c(0, .3), ylim = c(0, .3),
             parameters = K, points = 19, add = TRUE)



#y0 below represents the possible initial conditions given in pairs

y0 <- matrix(c(0.2, 0.2,0.1,0.1), ncol = 2, nrow = 2, byrow = TRUE)

predprey.trajectory = trajectory(predprey, y0 = y0, tlim=c(0,500),
               parameters = K, col = rep("black", 3))

 predprey.numericalSolution <-
   numericalSolution(predprey, y0 = c(.2, .2),tlim=c(0,2000), type = "one",
                       parameters = K, col = c("green", "orange"),
                       ylab = "x, y", ylim = c(0, .3))
 legend("bottomright", col = c("green", "orange"), legend = c("x", "y"), lty = 1)

EQ=findEquilibrium(
  predprey,
  y0=c(0.1,.1),
  parameters=k,
  summary=FALSE
)
EQ[5]
EQ[6]

