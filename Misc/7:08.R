z=8
z2=log10(z)

f = function(x) 2*sin(x)+1

g = function(x,y) 2*sin(x)+cos(y)

population = function(time) 1.5 * exp(0.126 * time)

k = 1:10
m=seq(1,10,by=0.1)
temp=rep(1,91)

population(m)

n=seq(1,10,by=0.5)

v=c(2,3,5,1,7,3.1415)
col=c("red", "blue", "green", "red")
threes=rep(3,101)

time=seq(0,10,by=0.1)
p=population(time)

plot(time,p,type='l',xlab='time (generations)',ylab='population size (individuals)')
lines(time, threes)
plot(time,threes)

##############################################################################

t=seq(0,10,by=0.5)
population1=function(time) 10*exp(0.1*time)
population2=function(time) 10*exp(-0.1*time)
p1=population1(t)
p2=population2(t)
plot(t,p1,xlab='time',ylab='pop size', col='blue', type='b',xlim=range(t),ylim=c(0,30))
lines(t,p2,col='red',type='b')

install.packages('plot3D')
library(plot3D)
M=mesh(seq(0,6*pi,length.out=50), seq(pi/3,pi,length.out=50))
u=M$x
v=M$y
x=v*cos(u)
y=v*sin(u)
z=10*u
surf3D(x,y,z,colvar=z,colkey=T,box=T,byu='b',phi=40,theta=0)

hist(runif(100, min=0, max=100))

x1=runif(1,0,100)
x1sum=0
for(i in 1:100)
{
  sum = sum+i
}
if(x1>20)
{
  yt=100  
} else
{
  yt=0
  yt
  
  initial = 2
  a = 1.5
  for(i in 1:10)
  {
    initial = a*initial
    print(initial)
  }
  
  ###########################################################################
  a=2
  iterations = 4
  x = rep(0, iterations)
  x[1]=2
  for(i in 1:(iterations - 1))
  {
    x[i+1] = x[i] * a
  }
  t=1:iterations
  plot(t,x,type='l')
  
  
  ###########################################################################
  
  my_data = read_excel("/cloud/project/sharks.xlsx")
  plot(my_data$Longitude, my_data$Latitude')