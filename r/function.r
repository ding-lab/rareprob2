I<-function(p1,p2)
{
	if(p1==p2)
	1
	else
	0
}
## X 
pie_thetarho<-function(x,p)
{
	obj<-(x^(p[1]+p[3]-1))*((1-x)^(p[2]+n-p[3]-1))
	obj
}
## R
pie_2<-function(x,p)
{
	obj<-(x^(p[1]+p[3]-1))*((1-x)^(p[2]+m-p[3]-1))
	obj
}
##log(gamma(n))
lngamma1<-function(n)
{	
		if(n<=171)
		return (log(gamma(n)))
	
	a<-round(n-171,0)
	b<-n-a
	result<-0
	temp<-n
    ##print(a)
	for(i in 1:a)
	{
		result<-result+log(temp-1)
		temp<-temp-1
	}
	result<-result+log(gamma(b))
	result
}
lngamma<-function(n)
{  
 
  a<-0
  b<-0
  for(i in 1:length(n))
  {
    if(n[i]<=171)
    { n[i]<-log(gamma(n[i]))}
    else
    {
      result<-0
      temp<-n[i]
      ##print(a)
      a<-round(n[i]-171,0)
      b<-n[i]-a
      for(j in 1:a)
      {
        result<-result+log(temp-1)
        temp<-temp-1
      }
      result<-result+log(gamma(b))
      n[i]<-result
    }
  }
    
  return (n)
}


#gamma(m)/gamma(n)
gammadiv<-function(m,n)
{
  if(m<=170&&n<=170)
    return (gamma(m)/gamma(n))
  return (exp(lngamma(m)-lngamma(n)))
}


##P(Y|X=1),L(X|R)
est2<-function(p,x1,x2,n1,n2)
{
P1=log(pgamma(p[1],p[2]))-log(gamma(p[1]))-log(gamma(p[2]))+lngamma(x1+p[1])+lngamma(n1-x1+p[2])-lngamma(n1+p[1]+p[2])
P2=log(pgamma(p[3],p[4]))-log(gamma(p[3]))-log(gamma(p[4]))+lngamma(x2+p[3])+lngamma(n2-x2+p[4])-lngamma(n2+p[3]+p[4])
-sum(P1+P2)
}
##P(Y|X=0) 
est3<-function(p,x)
{
P1=log(pgamma(p[1],p[2]))-log(gamma(p[1]))-log(gamma(p[2]))+lngamma(x+p[1])+lngamma(n-x+p[2])-lngamma(n+p[1]+p[2])
-P1
}
##L(X),L(R) ##unused
estTemp<-function(p)
{	
	if(!is.nan(1-p[1]))
	{
		if((1-p[1])>0)
		{	
			temp<-exp((1-p[1])*0.8+p[1]*P2)*0.8+exp((1-p[1])*0.2+p[1]*P2)*0.2
			PsXsXns<-(1-p[1])*sum(x*lxr)+p[1]*sum(P2)-log(temp)
			-sum(PsXsXns)
		}
		else
			return (Inf)
	}
	else return(Inf)
}
##unused
est4<-function(p)
{	
	if(!is.nan(1-p[1]))
	{
		if((p[2])>0)
		{	
			
			temp<-(1+exp(p[1]+p[2]*P2))
									
			PsXsXns<-p[1]*x+p[2]*P6-log(temp)
			sum(PsXsXns)
		}
		else
			return (Inf)
	}
	else return(Inf)
}
##unused
est6_1<-function(p)
{ 		 
	temp1<-(p[1]*sum(x)+p[2]*sum(P6))*lxr
	temp2<-0
	temp3<-0	
	for(i in 1:m)
	{
		temp3<-0
		for(j in 1:m)
		{
			if(Wss>Omega&&i!=j)
				temp3<-temp3+log(exp(lxr*p[1]+lxr*p[2]*Wss[i,j])*I(x[j],1)+exp(lxr*p[2]*Wss[i,j])*I(x[j],0))
		}
		temp2<-temp2+temp3;
	}
	obj<-temp2-temp1
	obj
}
##unused
est6<-function(p)
{ 		 
	temp1<-(p[1]*sum(x)+p[2]*sum(P6))*lxr
	
	temp2<-0
	temp3<-0	
	for(i in 1:m)
	{
		temp3<-sum(log(exp(lxr*p[1]+lxr*p[2]*Wss[i,])*I(x[j],1)+exp(lxr*p[2]*Wss[i,])*I(x[j],0)))
		temp2<-temp2+temp3;
	}
	obj<-temp2-temp1
	
	return (obj)
}
##unused
est6_2<-function(p)
{ 		 
	temp1<-(p[1]*sum(x)+p[2]*sum(P6))*lxr
	
	temp2<-0
	temp3<-numeric()
		
	for(i in 1:m)
	{
		
		temp2<-(exp(lxr*p[1]+lxr*p[2]*Wss[i,]-temp1)*I(x[j],1)+exp(lxr*p[2]*Wss[i,]-temp1)*I(x[j],0))
		temp3<-c(temp3,sum(temp2));
	}
	
	obj<-1/prod(temp3)
	
	-obj
}
##unused
est6_3<-function(p)
{ 		 
	temp1<-(p[1]*sum(x)+p[2]*sum(P6))*lxr	
	temp2<-print(lxr*p[1]*10000000000000+lxr*p[2]*Wss*10000000000000)
	temp4<-print(lxr*p[2]*Wss*10000000000000)
	temp3<-sum(log(exp(temp2)*Iss+exp(temp4)*Iss))
	obj<-temp3-m*temp1
	obj
}
##unused
est6_4<-function(p)
{
	temp1<-(p[1]*sum(x)+p[2]*sum(P6))
	temp2<-p[1]+p[2]*Wss
	temp4<-p[2]*Wss
	
	temp3<-sum(log(exp(temp2)*Iss+exp(temp4)*Iss))
	obj<-m*temp1-temp3
	obj
}
##unused
est6_5<-function(p)
{
	temp1<-sum(p[1]*x*0.8+p[2]*WI*0.2)
	print(temp1)

	temp2<-0
	temp3<-numeric()
	for(i in 1:m)
	{
		
		temp3<-0
		temp3<-temp3+log(exp(p[1]*0.8+p[2]*Wss[i,]*0.2)*I(x[i],1)+exp(p[2]*Wss[i,]*0.2)*I(x[i],0))
	}
	temp2<-sum(temp3)
	print(temp2)
	obj<-temp2-temp1
	obj
}
##use this function.
est6_6<-function(p,px0,px1)
{
	temp1<-sum(p[1]*x*px0+p[2]*WI*px1)
	

	temp2<-0
	temp3<-numeric()
  print(temp1)
	temp3<-log(exp(p[1]*px0+p[2]*Wss*px1)*Iss[i,j]+exp(p[2]*Wss[i,j]*px1)*Iss[i,j])
	temp2<-sum(temp3)
	print(temp2)
	obj<-temp2-temp1
  print("obj")
  print(obj)
	obj
}
est6_7<-function(p,px00,px01,px10,px11)
{
  
  
}
##"L(R)" used
est5<-function(p)
{
	##p<-c(1,1)
	temp1<-sum(p[1]*r+p[2]*WI)
	temp2<-0
	temp3<-0
	for(i in m)
	{
		temp3<-0
		temp3<-sum(log(exp(p[1]+p[2]*Wss[i,])*I(x[j],1)+exp(p[2]*Wss[i,])*I(x[j],0)))
		temp2<-temp2+temp3
	}
	obj<-temp2-temp1
  print("obj")
  print(obj)
	obj
}
##unused
est5_1<-function(p)
{
	PsRsRns<-sum(p[1]*r+p[2]*WI)
	temp1<-sum(p[1]*r+p[2]*WI)
	temp2<-0
	temp3<-0
	for(i in m)
	{
		temp3<-temp3+sum(log(exp(p[1]+p[2]*Wss[i,])*I(x[j],1)+exp(p[2]*Wss[i,])*I(x[j],0)))
		temp2<-temp2+temp3
		
	}
	obj<-temp2-temp1
	obj
}


