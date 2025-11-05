
rm(list=ls())

setwd("/mnt/DDrive/BONI-ROOT/Bio/Source/gt-cpp-ode-malariadrugpk-2019")
#setwd("D:/BONI-ROOT/Bio/Source/gt-cpp-ode-malariadrugpk-2019")

A = read.csv("out.dha.allpatients.csv")



pids = unique(A$PID)
allhours = unique(A$HOUR)



plot(-1,-1,xlim=c(0,75), ylim=c(-3, 2), xlab="hour", ylab="log mg DHA in central compt" )
for ( id in pids ) 
{
    B = subset( A, PID==id ) 
    #plot( B$HOUR, B$COMP2CONC, pch=19, xlim=c(1,9), ylim=c(-0.1, 0.7) )
    lines( B$HOUR, log10( B$COMP2CONC ), pch=19 )
}

plot(-1,-1,xlim=c(0,75), ylim=c(0, 27), xlab="hour", ylab="mg DHA in central compt" )
for ( id in pids ) 
{
  B = subset( A, PID==id ) 
  #plot( B$HOUR, B$COMP2CONC, pch=19, xlim=c(1,9), ylim=c(-0.1, 0.7) )
  lines( B$HOUR,  B$COMP2CONC, pch=19 )
}




num_cured = 0
plot(-1,-1,xlim=c(0,168), ylim=c(-0.1, 6), xlab="hour", ylab="log parasite density per ul" )
for ( id in pids ) 
{
  B = subset( A, PID==id ) 
  #plot( B$HOUR, B$COMP2CONC, pch=19, xlim=c(1,9), ylim=c(-0.1, 0.7) )
  lines( B$HOUR, log10( B$PARASITEDENSITY ), pch=19 )
  
  if( B$PARASITEDENSITY[168] < 1.0 )
  {
      num_cured=num_cured+1
  }
}

text(150,5,num_cured/ ( max(A$PID)+1 )  )

# plot(-1,-1,xlim=c(0,168), ylim=c(-0.1, 12000), xlab="hour", ylab="parasite density per ul" )
# for ( id in pids ) 
# {
#   B = subset( A, PID==id ) 
#   #plot( B$HOUR, B$COMP2CONC, pch=19, xlim=c(1,9), ylim=c(-0.1, 0.7) )
#   lines( B$HOUR, B$PARASITEDENSITY, pch=19 )
# }



# all the quantiles will go into this matrix; the first column is the hour
X = matrix(0,length(allhours),6)

r=1
for( hr in allhours )
{
    C = subset( A, HOUR==hr ) 
    qs=quantile( C$COMP2CONC , probs=c(.025,.25,.5,.75, .975) )  
    X[r,1] = hr
    X[r,2] = qs[1]
    X[r,3] = qs[2]
    X[r,4] = qs[3]
    X[r,5] = qs[4]
    X[r,6] = qs[5]
    
    r = r + 1
}


plot( -1, -1, xlim=c(0,96), ylim=c(-3.1, 1.8), xlab="hour", ylab="mg DHA in central compt", yaxt='n' )
grid(nx=NULL,ny=NULL)
#axis(2, at=c(0,0.4771213,1,1.4771213,2,2.4771213,3), labels = c("1","3","10","30","100","300","1000"), las=2) #,  las=2)  # 1 means bottom,2 is left
axis(2, at=c(-3,-2,-1,0,1), labels = c(".001",".01",".1","1","10"), las=2) 


lines( X[,1], log10( X[,2] ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,3] ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,4] ), pch=19, col="#AA4371", lwd=5 )
lines( X[,1], log10( X[,5] ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,6] ), pch=19, col="#BBBBBB", lwd=3 )


plot( -1, -1, xlim=c(0,96), ylim=c(-3.1, 1.8), xlab="hour", ylab="ng/ml PPQ in central compt", yaxt='n' )
grid(nx=NULL,ny=NULL)
axis(2, at=c(-3,-2,-1,0,1), labels = c(".001",".01",".1","1","10"), las=2) 

lines( X[,1], log10( X[,2]/5.5 ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,3]/5.5 ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,4]/5.5 ), pch=19, col="#AA4371", lwd=5 )
lines( X[,1], log10( X[,5]/5.5 ), pch=19, col="#BBBBBB", lwd=3 )
lines( X[,1], log10( X[,6]/5.5 ), pch=19, col="#BBBBBB", lwd=3 )




