# load SPGF file
source(file="StochasticPopulationGrowthFunction.R")
# load data
load("checkerspots.Rdata")
# plot the rs values
plot(yrs.data,rs.data,type="b") 

yrs.index=which(yrs.data>1971)
x.mat=SPGFunction(x0=1000,rs=rs.data[yrs.index],a=a,finalN=200,howMany=100)
matplot(x.mat,type="l",lty=1,log="y")
# geometric mean
prod(rs.data[yrs.index])^(1/(length(yrs.index)))

# >1990 0.9646526
# >=1971 0.9524209``
# <1971 1.041043