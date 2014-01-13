P_value<-function(cdf,x,paramet=numeric(0),side=0)
{
	print(".....\n")
	print("x")
	print(x)
	n<-length(paramet)
#	P <- cdf(x)	
	P<-switch(n+1,
			cdf(x),	
			cdf(x,paramet),
			cdf(x,paramet[1],	paramet[2]),
			cdf(x,paramet[1],paramet[2],paramet[3]),
			cat("n is out bound of length of paramet")
			)
	if( is.nan(P) )
	{
		P=0;
		print("P is nan")
	}
	print("P=")
	print(P)
	if(side<0)	P
	else if(side>0)	1-P
	else
		if(P<1/2) 2*P
		else 2*(1-P)
}
