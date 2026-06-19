### Calculate Life Expectancy for mortality rates by year
### Code by: Rafael Meza rmeza@umich.edu
### assumes mortality rates matrix as input
### See formulas: https://mathworld.wolfram.com/LifeExpectancy.html

LEcalc=function(deathrates){
  ages=dim(deathrates)[1] ### should be 100 if using CISNET's ages 0-99
  years=dim(deathrates)[2] ## number of years
  LE=matrix(0,nrow=ages,ncol=years) 
  
  for (j in 1:years){
    qx=deathrates[,j]
    LT=matrix(0,ages,6)
    LT[1,1]=10000
    LT[1,2]=LT[1,1]*(1-exp(-qx[1]))
    for (i in 2:ages){
      LT[i,1]=LT[i-1,1]-LT[i-1,2]
      LT[i,2]=LT[i,1]*(1-exp(-qx[i]))
    }
    LT[,3]=LT[,1]/LT[1,1]
    LT[1:(ages-1),4]=(LT[1:(ages-1),3]+LT[2:ages,3])/2
    LT[ages,4]= (LT[ages,3]+LT[ages,3]*(1-exp(-qx[ages])))/2#LT[100,3]/2#(LT[100,3]+LT[100,3]*(1-exp(-qx[100])))/2 ## I fixed this to add an extra year
    
    LT[1,5]=sum(LT[,4])
    for (i in 1:(ages-1)){
      LT[i+1,5]=LT[i,5]-LT[i,4]
    }
    LT[,6]=LT[,5]/LT[,4]
    
    LE[,j]=LT[,6]
  }

  return(LE)  
}

