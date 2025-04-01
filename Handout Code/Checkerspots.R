#load our SPGF file
source(file="StochasticPopulationGrowthFunction.R")
# look at the linearized model and see if we can predict things about it using information about the rs

a = 0
rs = c(1.2,0.9)
x0 = 2
finalN = 100
howMany = 10000

# running function and storing into a matrix
xMat=SPGFunction(x0, rs, a, finalN, howMany)
#plot
matplot(0:(finalN-1),log(xMat),type="l",lty=1)
# what is going on at the end of the run?
y=log(xMat[finalN,])
hist(y,n=25,freq=FALSE)
M=mean(y)
SD=sd(y)
xs=seq(min(y),max(y),length=100)
ys=dnorm(xs,mean=M,sd=SD)
lines(xs,ys,lwd=3,col="darkred")

# what is the mean trajectory?
# geometric mean 
heart=prod(rs)^(1/length(rs))
# plot the line of mean trajectory
# abline(a=log(x0),b=log(heart),lwd=4,col="purple")