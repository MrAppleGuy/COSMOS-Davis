#extracting CA rows
CA = which(us_states$state=="California")
#getting number of cases in CA (cumulative)
cases = us_states$cases[CA]
#plot the cases
plot(cases,type='l')
#extract new cases by subtracting previous days
n = length(cases)
newCases = cases[-1]-cases[-n]
plot(newCases,type='l')
#get rid of negatives and adjust the zeros
newCasesCleaned = pmax(newCases,0)+0.5
#take logs
y=log(newCasesCleaned)
plot(y,type='l')
# first find a line that fits days 35-65
days=671:720
plot(days,y[days],type='p')
# fit a line
firstSegApp=lm(y[days]~days)
summary(firstSegApp)
# plot the line against data
a=firstSegApp$coefficients[1]
b=firstSegApp$coefficients[2]
lines(days,a+b*days,col="red",lwd=6)
# find the doubling/halving time
log(2)/b