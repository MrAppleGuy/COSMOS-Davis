# the following function simulates a matrix model
# and produces one of the three types of plot
# type 1 plots X(t) versus t
# type 2 plots total abundance and fractions against t
# type 3 animates X over time
# inputs are the matrix A, the initial condition x0,
# the length of the run t.end, and the plot type
matrix.sim=function(A,x0,t.end=25,type=1,output=FALSE){
  # store all of the simulation
  # in a matrix with t.end+1 rows
  # and 2 columns
  X=matrix(NA,nrow=t.end+1,ncol=length(x0))
  # put in the initial row
  X[1,]=x0
  # run the simulation
  for(t in 1:t.end)X[t+1,]=A%*%X[t,]
  cols=topo.colors(length(x0))
  if(type==1){
    par(mfrow=c(1,1),mar=c(4,4,0,0))
    matplot(X,type="b",lty=1,bty="n",col=cols,bg=cols,pch=21)
  }
  if(type==2){
  # plot totals and fractions 
  totals=rowSums(X)
  par(mfrow=c(1,2),mar=c(4,4,0,0))
  plot(totals,type="b",lty=1,bty="n",pch=21,cex=0.75,xlab="years",bg="gray")
  fractions=diag(1/totals)%*%X
  matplot(fractions,type="b",lty=1,bty="n",col=cols,bg=cols,pch=21,xlab="years",cex=0.75)}
  if(type==3){
    par(mfrow=c(1,1),mar=c(2,2,0,0))
    barplot(X[1,],ylim=c(0,max(X)),xlim=c(1,length(x0)),space=0)
    Sys.sleep(1)
    for(t in 1:(t.end+1)){barplot(X[t,],ylim=c(0,max(X)),xlim=c(1,length(x0)),space=0);
      Sys.sleep(0.1)}
  }
  if(output)return(X)
}

multiply = 1.148535747738697

# the following is the loggerhead matrix of Crouse et al. 1987 Ecology
A.loggerhead=matrix(c(0, 0, 0,0, 127, 4, 80,
           0.6747, 0.737, 0, 0, 0, 0, 0,
           0, 0.0486, 0.6610, 0, 0, 0, 0,
           0, 0, 0.0147, 0.6907, 0, 0, 0,
           0, 0, 0, 0.0518, 0, 0, 0,
           0, 0, 0, 0, 0.8091, 0, 0,
           0, 0, 0, 0, 0, 0.8091, 0.8089),7,byrow=TRUE)

jhxt=matrix(c(0, 0, 0,0, 127, 4, 80,
                      0.6747, 0.737, 0, 0, 0, 0, 0,
                      0, 0.0486, 0.6610, 0, 0, 0, 0,
                      0, 0, 0.0147, 0.6907, 0, 0, 0,
                      0, 0, 0, 0.0518, 0, 0, 0,
                      0, 0, 0, 0, multiply*0.8091, 0, 0,
                      0, 0, 0, 0, 0, multiply*0.8091, multiply*0.8089),7,byrow=TRUE)

# the following is a matrix model of Northern Seal Data from Yorke 
b=c(0,0,0,0.105,0.267,0.286,0.315,0.315,0.315,0.315,0.315,0.315,0.315)
l=c(1,0.782,0.612,0.478,0.445,0.404,0.362,0.32,0.28,0.242,0.208,0.178,0.15)
A.seal=matrix(0,13,13)
A.seal[1,]=b*0.782
for(i in 2:13)A.seal[i,i-1]=l[i]/l[i-1]

# the following creates random walk matrices of different types. inputs are 
# n number of locations
# d fraction of individuals moving 
# right fraction of moving individuals going to the right
A.walk=function(n=10,d=0.5,right=0.5){
  A=matrix(0,n,n)
  for(i in 1:(n-1)){
  A[i,i]=1-d
  A[i,i+1]=d*(1-right)
  A[i+1,i]=d*right
  }
  A[1,1]=1-d*right
  A[n,n]=1-d*(1-right)
  return(A)
}

# age structured loggerhead
A.loggerhead.age=matrix(0,55,55)
A.loggerhead.age[2,1]=0.6747
# small juvs
for (i in 2:8){
  A.loggerhead.age[i+1,i]=0.7857
}
# large juvs
for (i in 9:16){
  A.loggerhead.age[i+1,i]=0.6758}
# subadults
for (i in 17:22){
  A.loggerhead.age[i+1,i]=0.7425}
# novice breeders
A.loggerhead.age[1,23]=127
A.loggerhead.age[24,23]=0.8091
# 1st remigrants
A.loggerhead.age[25,24]=0.8091
A.loggerhead.age[1,24]=4
# matures
for (i in 25:54){
  A.loggerhead.age[i+1,i]=0.8091
  A.loggerhead.age[1,i]=80}
A.loggerhead.age[1,55]=80

eigen(jhxt)$values

