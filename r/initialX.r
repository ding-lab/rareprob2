print("initialX1....")
if(omega=="")
	omega<-1


parm<-read.table(paste(tmpfile_path,path,"_parms",sep=""),header=TRUE) ##path from calling parameter
##parm<-read.table("./data0000001_parms",header=TRUE)
m<-c(parm["m"][1])
m<-m[["m"]]
n<-c(parm["n"][1])
n<-n[["n"]]
pp_value<-c(parm["p_value"][1])
pp_value<-pp_value[["p_value"]]

mafseq_path<-c(parm["mafseq_path"][1])
mafseq_path<-mafseq_path[["mafseq_path"]]

#cancerlist_path<-c(parm["cancerlist_path"][1])
#cancerlist_path<-cancerlist_path[["cancerlist_path"]]

#siteInfo_path<-c(parm["siteInfo_path"][1])
#siteInfo_path<-siteInfo_path[["siteInfo_path"]]

P<-matrix(scan(as.character(mafseq_path),0),ncol=m,nrow=5,byrow=TRUE)
Pcase<-P[1,]
Pcontrol<-P[2,]
Pall<-P[3,]


Wss<-rep(0,m)
Zs<-rep(0,m)
X<-rep(2,m)
for(i in 1:m)
{
  
  ##when mafcase <= mafcontrol, we let x=0 as a noncausal site
  if(Pcase[i]-Pcontrol[i]>0)
	Zs[i]<-2*(Pcase[i]-Pcontrol[i])/(sqrt(2/n)*sqrt((Pcase[i]+Pcontrol[i])*(2-Pcase[i]-Pcontrol[i])))
	else
	{##waiting modify..
	  Zs[i]<-0
    ##X[i]<-0
  }
}

Wss<-matrix(1,nc=m,nr=m)
Iss<-matrix(1,nc=m,nr=m)
tempZs<-numeric()
Omega<-rep(0,m)
Omega_mean<-rep(0,m)
Omega_31<-rep(0,m)
Omega_32<-rep(0,m)
temp_sort<-rep(0,m)
for(i in 1:m)
{
	for(j in 1:m)
	{
		if(Zs[i]!=0&&Zs[j]!=0)
		Wss[i,j]<-2*Zs[i]*Zs[j]/(Zs[i]^2+Zs[j]^2)
		else
		Wss[i,j]<-1
		
	}
	Wss[i,j]<-0
	if(omega==1)##mean
	{
	  ##Omega_mean[i]<-mean(Wss[i,])
    Omega[i]<-mean(Wss[i,])
	}
	if(omega==2)##1/3
	{
	  temp_sort<-sort(Wss[i,])
	  ##Omega_31[i]<-temp_sort[round(1/3*m,0)]
	  Omega[i]<-temp_sort[round(1/3*m,0)]
	}
	if(omega==3)##2/3
	{
	  temp_sort<-sort(Wss[i,])
	  ##Omega_31[i]<-temp_sort[round(1/3*m,0)]
	  Omega[i]<-temp_sort[round(2/3*m,0)]
	}
	
  
	
}




##Wss=0,neighbor
for(i in 1:m)
	for(j in 1:m)
	{
    if(Wss[i,j]<Omega[i])
        Wss[i,j]<-0
  }
		



##tmp path
write.table(Wss,paste(tmpfile_path,path,"_Wss.seq",sep=""))
write(Omega,paste(tmpfile_path,path,"_Omega.seq",sep=""),nc=m)
write(Zs,paste(tmpfile_path,path,"_Zs.seq",sep=""),nc=m)

##vs<-sqrt((1-Pall)/Pall)
vs<-rep(0,m)
for(i in 1:m)
{
  if(Pall[i]==0)
  {
    vs[i]<-0
    X[i]<-0  
  }
    
  else
    vs[i]<-sqrt((1-Pall[i])/Pall[i])
}
##vs<-sqrt((1-Pall)/Pall)

S<-numeric()
a<-0
b<-0
c<-0
for(i in 1:m)
{
  if(i==1)
  {
    if((vs[i]^2+vs[i+1]^2)!=0)
      S<-c(S,(vs[i]*Zs[i]+vs[i+1]*Zs[i+1])/sqrt(vs[i]^2+vs[i+1]^2))	
    else
      S<-c(S,999)
  }
  else if(i==m)
  {
    if((vs[i]^2+vs[i-1]^2)!=0)
      S<-c(S,(vs[i]*Zs[i]+vs[i-1]*Zs[i-1])/sqrt(vs[i]^2+vs[i-1]^2))
    else
      S<-c(S,999)
  }
  else
  {
    if((vs[i]^2+vs[i-1]^2+vs[i+1]^2)!=0)
      S<-c(S,(vs[i-1]*Zs[i-1]+vs[i]*Zs[i]+vs[i+1]*Zs[i+1])/sqrt(vs[i]^2+vs[i-1]^2+vs[i+1]^2))
    else
      S<-c(S,999)
  }
}


mean.test1<-function(x,mu=0,sigma=-1,n=1)
{
	source("r/P_value.r")
	
	if(sigma>0)
	{
		
		z<-(x-mu)/(sigma/sqrt(n))
		side=0
		P<-P_value(pnorm,z,side=side)
		return (P)
	}
	else
	{
		("sigma wrong")
	}
}

P<-rep(1,m)

Iss<-rep(1,m-1)
pre_x<-0;
count0<-0;
##waiting modify
for(i in 1:m)
{
	
  P[i]<-mean.test1(S[i],mu=0,sigma=1)
  if(X[i]==2)
  {
    ##print(i)
    ##print(P[i])
    if(P[i]>=0.001)
    {
      count0<-count0+1
      X[i]=0
    }
    else X[i]=1
  }
}

##
candidateCancerListExist = FALSE
siteInfoListExist=FALSE
if ( file.exists( cancerlist_path ) )
{
	cancerList<- read.table(cancerlist_path)
	cancerList <- t( cancerList )
	candidateCancerListExist=TRUE
}
if ( file.exists( siteInfo_path ) )
{ 
	siteInfo <- read.table(siteInfo_path)
	siteInfo <- t(siteInfo) # length(siteInfo) is equal to m
	siteInfoListExist=TRUE
}



iscontains<-function(x)
{
	
	for( i in 1: length(cancerList) )
	{
		if( cancerList[i] == x )
		{
			
			return (TRUE)
		}
	}

	return (FALSE)
}

#2. 
if( siteInfoListExist == TRUE && candidateCancerListExist == TRUE )
{
	for( i in 1:m)
	{
		if( iscontains( siteInfo[i] ) )
			X[i]<-1
	}
}


#3. 

for( i in 1:m)
{
	if( Pcontrol[i] != 0 && Pcase[i] != 0)
	{
		if( Pcase[i]/Pcontrol[i] >= 20)
			X[i] <- 1 
	}
}


##

xinfo<-data.frame(count0=count0,count1=m-count0)

write(X,paste(tmpfile_path,path,"_initX.seq",sep=""),nc=m)
write(X,paste(tmpfile_path,path,"_X.seq",sep=""),nc=m)
print("initialX1....")

