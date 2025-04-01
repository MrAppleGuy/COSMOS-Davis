jhxt=observations[which(observations$class=="Mammalia"),]
x=log(jhxt$body.mass)
y=log(jhxt$metabolic.rate)
plot(x,y,pch=21,bg="yellow")
slope=lm(y~x)
summary(slope)
abline(slope,lwd=3,col="red")