# Modified from the original script pplane.r written by Daniel Kaplan,
# Dept. of Mathematics, Macalester College, kaplan@macalester.edu 

# Modifications by S. Ellner for use in connection with the textbook
# Dynamic Models in Biology by S.P. Ellner and J. Guckenheimer,
# Princeton University Press (2006)  

require(deSolve); 


phasearrows <- function(fun,xlims,ylims,resol=25, col='gray', add=F,parms=NULL,jitter=FALSE) {
  if (add==F) {par(mar=c(0.1,0.1,0.1,0.1))
    plot(1,xlim=xlims, ylim=ylims, type='n',xlab="",ylab="",xaxt="n",yaxt="n");
  }
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=T, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=F, resol, resol);
  npts <- resol*resol;
  if(jitter) {
    xspace <- abs(diff(xlims))/(resol*10);
    yspace <- abs(diff(ylims))/(resol*10);
    x <- x + matrix(runif(npts, -xspace, xspace),resol,resol);
    y <- y + matrix(runif(npts, -yspace, yspace),resol,resol);
  }
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  maxx <- max(abs(z1));
  maxy <- max(abs(z2));
  dt <- min( abs(diff(xlims))/maxx, abs(diff(ylims))/maxy)/resol;
  lens <- sqrt(z1^2 + z2^2);
  lens2 <- lens/max(lens); 
  arrows(c(x), c(y), c(x+dt*z1/((lens2)+.1)), c(y+dt*z2/((lens2)+.1)),length=.04, col=col);
}

showcontours <- function(fun,xlims, ylims,resol=250,add=F, colors=c('red', 'blue'),parms=NULL) {
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1, add=add, col=colors[1]);
  contour(x[,1],y[1,],z2, add=T, col=colors[2]); 
}

nullclines <- function(fun,xlims, ylims, resol=250, add=F,parms=NULL) {
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1,levels=c(0), drawlabels=F,add=add, col="blue",lwd=3,bty="n");
  contour(x[,1],y[1,],z2,levels=c(0), drawlabels=F,add=T, col="red",lwd=3); 
  #title(main="Blue=x nullcline, Green=y nullcline",cex=0.35); 
}

grid=function(fun,xlim,ylim,parms,ngrid,maxtime=50,add=F,color="black") {
  if (add==F) {
    plot(1,xlim=xlim, ylim=ylim, type='n',xlab="x",ylab="y");
  }
  xvals=seq(xlim[1],xlim[2],length=ngrid); 
  yvals=seq(ylim[1],ylim[2],length=ngrid); 
  for(i in 1:ngrid) {
    for(j in 1:ngrid) {
      out=lsoda(times=seq(0,maxtime,.2),y=c(xvals[i],yvals[j]),func=fun,parms=parms); 
      points(out[,2],out[,3],type="l",lwd=2,col=color);
      out=lsoda(times=-seq(0,maxtime,.2),y=c(xvals[i],yvals[j]),func=fun,parms=parms); 
      points(out[,2],out[,3],type="l",lwd=2,col=color);
      
    }}
}

# Newton's method to find equilibria of vector field.
# func() must have the same input arguments and returns as for lsoda/rk4.  
# Inputs: 
#   x0 = intial guess at equilibrium. If x0 is not supplied in the call, 
#        the user chooses it from the current graphics device via locator()
#         and the equilibrium is plotted to the same device. Plotting
#         symbol is closed/open=stable/unstable, circle/triangle=eigenvalues imaginary/real.   
#   tol= Convergence tolerance 
#   niter = Maximum number of iterations
#   inc = finite-difference increment for derivative estimates 
# Coded 5/25/06 by SPE based on Matlab toggle.m by JG 

newton=function(func,x0=NULL,parms=NULL,tol=1e-16,niter=40,inc=1e-6,plotit=TRUE) {
  x=x0; #initial x  
  if (is.null(x0)){x = locator(n=1); x=c(x$x,x$y)};
  nx = length(x); # length of state vector
  ######### Newton iteration loop: start  
  for(i in 1:niter){  
    y = func(0,x,parms)[[1]] 
    df = matrix(0,nx,nx); # Compute df
    for(j in 1:nx) {
      #Increment vector for estimating derivative wrt jth coordinate
      v=rep(0,nx); 
      v[j] = inc; 
      df[,j]=  (func(t,x+v,parms)[[1]] - func(t,x-v,parms)[[1]])/(2*inc) 
    }
    if (sum(y^2) < tol){  #check for convergence 
      if(plotit){
        ev=eigen(df)$values; pch1=1+as.numeric(Im(ev[1])!=0); pch2=1+as.numeric(max(Re(ev))<0);
        pchs=matrix( c(21,21,21,21),2,2,byrow=T); 
        bg.cols=c("white","red")
        cols=c("black","red")	
        points(x[1],x[2],type="p",pch=pchs[pch1,pch2],cex=3,col=cols[pch2],bg=bg.cols[pch2])
      }
      cat("Fixed point (x,y) = ",x,"\n"); 
      cat("Jacobian Df=","\n"); print(df);cat("Eigenvalues","\n"); print(eigen(df)$values)
      return(list(x=x,df=df))   
    }	
    x = x - solve(df,y) # one more step if needed 
    cat(i, x, "\n") #print out the next iterate 
  }
  ######### Newton iteration loop: end  
  cat("Convergence failed"); 
}

DrawManifolds=function(fun.lsoda,parms,x0=NULL,maxtime=100) {
  xbar=newton(fun.lsoda,x0=x0,parms=parms,plotit=FALSE);
  x=xbar$x; df=xbar$df; V=eigen(df)$vectors; ev=eigen(df)$values;  
  if (ev[1]*ev[2] > 0) {
    cat("Fixed point is not a saddle \n");
  }else{
    i1=which(ev>0); i2=which(ev<0); 
    v1=V[,i1]; v2=V[,i2]; eps=1e-3;  
    out1=lsoda(times=seq(0,maxtime,.1),y=x+eps*v1,func=fun.lsoda,parms=parms); points(out1[,2],out1[,3],type="l",lwd=2,col="red");
    out2=lsoda(times=seq(0,maxtime,.1),y=x-eps*v1,func=fun.lsoda,parms=parms); points(out2[,2],out2[,3],type="l",lwd=2,col="red");
    out3=lsoda(times=-seq(0,maxtime,.1),y=x+eps*v2,func=fun.lsoda,parms=parms); points(out3[,2],out3[,3],type="l",lwd=2,col="black");
    out4=lsoda(times=-seq(0,maxtime,.1),y=x-eps*v2,func=fun.lsoda,parms=parms); points(out4[,2],out4[,3],type="l",lwd=2,col="black");
    title(sub="Black=stable manifold, Red=unstable manifold"); 
  }
}

# Compute Jacobian of a planar vector field at a point (x,y),
# either input or chosen with locator().
jacobianAtXY <- function(fun,x=NULL, y=NULL,inc=1e-7){
  if (is.null(x)|is.null(y)) {
    x0 <- locator(n=1); x <- x0$x; y <- x0$y;  
  }
  foo <- fun(x,y); h = inc; 
  foox <- fun(x+h,y); fooy <- fun(x,y+h);
  A <- (foox[1] - foo[1])/h;
  B <- (fooy[1] - foo[1])/h;
  C <- (foox[2] - foo[2])/h;
  D <- (fooy[2] - foo[2])/h;
  return(matrix( c(A,B,C,D ),2,2,byrow=T))
}

# Cheap man's version of Polking's pplane.m for planar vector fields.
# The vector field fun() must be in the format of toggle() and ml() above: input
# arguments are (x,y,parms), the return is the vector field *as a vector*.   
Rpplane=function(fun,xlim,ylim,parms=NULL,add=F,ngrid=4,maxtime=100) {
  fun.lsoda=function(t,y,p) {dx=fun(y[1],y[2],parms=p); return(list(dx))}
  menu.go=1; while(menu.go>0) {
    j= menu(choices=c("Draw phase arrows", 
                      "Draw nullclines",
                      "Find fixed point (requires mouse click on plot)",
                      "Start Forward trajectory (requires mouse click on plot)",
                      "Start Backward trajectory (requires mouse click on plot)",
                      "Extend current trajectory (requires mouse click on plot)", 
                      "Local S/U manifolds for saddle (requires mouse click on plot)",
                      "Draw Grid of trajectories",
                      "Exit"),
            title = "R-pplane: select action & click OK",graphics=TRUE)
    if(j==1) {phasearrows(fun=fun,xlims=xlim,ylims=ylim,parms=parms,add=add); add=TRUE}
    if(j==2) {nullclines(fun=fun,xlims=xlim,ylims=ylim,parms=parms,add=add); add=TRUE};
    if(j==3) {xbar=newton(fun.lsoda,parms=parms); }
    if(j==4) {x=locator(n=1); out=lsoda(times=seq(0,maxtime,.1),y=c(x$x,x$y),func=fun.lsoda,parms=parms); 
    points(out[,2],out[,3],type="l",lwd=2,col="black"); } 
    if(j==5) {x=locator(n=1);out=lsoda(times=seq(0,-maxtime,-.1),y=c(x$x,x$y),func=fun.lsoda,parms=parms); 
    points(out[,2],out[,3],type="l",lwd=2,col="black");} 
    if(j==6) {nt=dim(out)[1]; x=out[nt,2:3]; times=out[,1]; 
    out=lsoda(times=out[,1],y=out[nt,2:3],func=fun.lsoda,parms=parms); 
    points(out[,2],out[,3],type="l",lwd=2,col="black");} 
    if(j==7) {DrawManifolds(fun.lsoda,parms=parms,maxtime=maxtime) }
    if(j==8) {grid(fun.lsoda,xlim=xlim,ylim=ylim,parms=parms,ngrid=ngrid,add=add); add=TRUE} 
    if(j==9) menu.go=0
  }
}



# NOTE: USE phase arrows first to get rid of axes etc.  
lotkaVolde=function(t,x,parms){
  a = parms[1]
  b = parms[2]
  c = parms[3]
  d = parms[4]
  dxdt=x[1]*(a-x[1] - b*x[2])
  dydt=x[2]*(c-x[2] - d*x[1])
  list(c(dxdt,dydt))
}
  atest=3
  mtest=1
lotkaVol=function(x,y,parms){
  a = parms[1]
  b = parms[2]
  c = parms[3]
  d = parms[4]
  dxdt=atest-x-x*y^2
  dydt=x*y^2-mtest*y
  c(dxdt,dydt)
}






## Start HERE
p=c(3,2,2,1)
t=seq(0,10,0.1)
# play with initial conditions here!
ini=c(x1=1,x2=5)
out = ode(ini,t,lotkaVolde,p)
# plot(out)



#Now here is the phase plane outputs

Rpplane(lotkaVol,c(0,5),c(0,5),c(3,2,2,1))