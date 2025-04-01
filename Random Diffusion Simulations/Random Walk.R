x = 0
r = runif(1)
walk = 0
count = 0
iterations = 1000
print(r)
maxDist = 20
minDist=1
currDist = 0
array=rep(1:iterations)
avgArray=rep(1:maxDist-minDist+1)
for(y in minDist:maxDist)
{
  for(x in 1:iterations)
  {
    while(walk < y)
    {
      r = runif(1)
      if(r < 0.5)
      {
        #print("h")
        walk = walk + 1
      }
      if(r>0.5)
      {
        #print("t")
        if(walk > 0)
        {
          walk = walk - 1
          
        }
      }
      count = count + 1
    }
    array[x] = count
    avgArray[y] = mean(array)
    count = 0
    walk=0
  }
  array=rep(1:iterations)
}
print(avgArray)
plot(minDist:maxDist,avgArray,xlab="Squared Traveled", ylab="Average Flips",type='b')