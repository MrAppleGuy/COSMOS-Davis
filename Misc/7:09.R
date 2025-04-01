a=2
b=1.0
iterations = 4
x = rep(0, iterations)
x[1]=0.5
for(i in 1:(iterations-1))
{
  x[i+1]=x[i] * a + b
}
t=1:iterations
plot(t,x,type='l')

#############################################################################

a = 0.5
influx = 2
iterations = 20
x = rep(0, iterations)
x[1] = 10000
for(i in 1:(iterations - 1))
{
  x[i+1] = x[i] * a
  x[i+1] = x[i+1] + influx
}
t=1:iterations
plot(t,x,type='l')