SIR=function(beta=1,gamma=0.1,S0=0.99,I0=0.01,Tf=45){
  S=numeric(Tf)
  S[1]=S0
  I=numeric(Tf)
  I[1]=I0
  for(t in 1:(Tf-1)){
    S[t+1]=S[t]-beta*S[t]*I[t]
    I[t+1]=S[t]*beta*I[t]+(1-gamma)*I[t] 
  }
  R=1-S-I
  cols=c("black","red","gray")
  par(cex.axis=1,cex.lab=1)
  matplot(cbind(S,I,R),type="l",col=cols,xlab="days",bty="n",ylab="abundances",lwd=2,lty=1,ylim=c(0,max(S,I,R)))
  legend("right",c("Susceptible","Infected","Removed"),lty=1,col=cols,bty="n",cex=0.5,lwd=2)
  }

SIR.phase=function(beta=1,gamma=0.1,S0=seq(0,1,length=15),I0=rep(0.01,15),Tf=45){
  S=matrix(NA,Tf,length(I0))
  S[1,]=S0
  I=matrix(NA,Tf,length(I0))
  I[1,]=I0
  for(t in 1:(Tf-1)){
    S[t+1,]=S[t,]-beta*S[t,]*I[t,]
    I[t+1,]=S[t,]*beta*I[t,]+(1-gamma)*I[t,] 
  }
  R=1-S-I
  cols=colorRampPalette(colors = c("red","blue"))(length(I0))
  cols=topo.colors(length(I0))
  par(cex.axis=1,cex.lab=1)
  matplot(x=S,y=I,type="l",col=cols,xlab="Susceptible",bty="n",ylab="Infected",lwd=2,lty=1,ylim=c(0,1),xlim=c(0,1))
abline(v=gamma/beta,lty=2)
}


