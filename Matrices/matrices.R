#the industrial revolution and its consequences have been a disaster for the human race
#define matrix
A=rbind(c(0,3),c(0.5,0.5))
#intitial condition
x0 = c(10000,0)
# how long
nEnd = 200
#store the output
xMat = matrix(NA,nEnd+1,2)
xMat[1,]=x0
# run the model
for(n in 1:nEnd) xMat[(n+1),] = A%*%xMat[n,]
#plot
matplot(xMat,type="b",lty=1)
# plot fraction of adults
plot(xMat[,2]/(xMat[,1]+xMat[,2]),type="l")
# change in pop from one time step to another
totalX=rowSums(xMat)
plot(totalX[-1] / totalX[-(nEnd+1)],type="l")