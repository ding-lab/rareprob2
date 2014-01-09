source("./r/function.r")

parm<-read.table(paste(tmpfile_path,path,"_parms",sep=""),header=TRUE) ##path from calling parameter

m<-c(parm["m"][1])
m<-m[["m"]]
n<-c(parm["n"][1])
n<-n[["n"]]
na<-c(parm["na"][1])
na<-na[["na"]]
nc<-c(parm["nc"][1])
nc<-nc[["nc"]]


Wss<-read.table(paste(tmpfile_path,path,"_Wss.seq",sep=""))
par<-0.1
count<-0
P<-matrix(scan(paste(tmpfile_path,path,"_maf.seq",sep=""),0),ncol=m,nrow=5,byrow=TRUE)
maf_case<-P[1,]
maf_control<-P[2,]
maf_all<-P[3,]
caseCount<-P[4,]
controlCount<-P[5,]
allCount<-caseCount+controlCount

x<-c(scan(paste(tmpfile_path,path,"_X.seq",sep="")))
##x<-x-1

r<-c(scan(paste(tmpfile_path,path,"_initR.seq",sep="")))
ce<-0
cb<-0

cxe<-0
cxb<-0
("cx+,cx-,ce,cb")
for(i in 1:m)
{
	
	if(x[i]==0&&r[i]==0)
	{
		cb<-cb+1;
	}
	if(x[i]==0&&r[i]==1)
	{
		ce<-ce+1;
	}
	
	if(x[i]==1&&r[i]==0)
	{
		cb<-cb+1;
		cxb<-cxb+1;
	}
	if(x[i]==1&&r[i]==1)
	{
		ce<-ce+1;
		cxe<-cxe+1;
	}
}
##proportions of e region & b region
pe<-ce/m
pb<-cb/m
cxe+cxb

pre_pyx1<-numeric()
pre_pyx2<-numeric()
pre_pyx3<-numeric()
pre_pyx4<-numeric()
pre_pr<-numeric()
pre_px<-numeric()
pre_pxr<-numeric()

repeat{
	Iss<-numeric()
	for(i in 1:m)
	{
		Iss<-rbind(Iss,as.numeric(!xor(x[i]-1,x-1)))
	}
  ##WI for every site's neighbor
  WI<-rep(0,m)
  for(i in 1:m)
  {
    for(j in 1:m)
      if(Wss[i,j]!=0)
      {
        WI[i]<-WI[i]+Wss[i,j]*Iss[i,j]
      }
  }
	##WI<-Wss*Iss
	
	

	p0<-c(0.5,0.5,0.5,0.5)
	p1<-c(0.5,0.5)
	PYsXs<-numeric()
	
	pyx1<-numeric()
	pyx2<-numeric()
	pyx3<-numeric()
	pyx4<-numeric()
  #estimate
	for(i in 1:m)
	{
		if(x[i]==1)	
		{
			out<-nlminb(p0,est2,x1=caseCount[i],x2=controlCount[i], n1=na, n2=nc,  lower=0.00000001, upper=1)
			

			pyx1<-c(pyx1,out[["par"]][1])
			pyx2<-c(pyx2,out[["par"]][2])
			pyx3<-c(pyx3,out[["par"]][3])
			pyx4<-c(pyx4,out[["par"]][4])
		}
		else
		{
			out<-nlminb(p1,est3,x=allCount[i],lower=0.00000001,upper=1)

			pyx1<-c(pyx1,out[["par"]][1])##alpha_theta_s
			pyx2<-c(pyx2,out[["par"]][2])##beta_theta_s
			pyx3<-c(pyx3,out[["par"]][1])##alpha_rho_s
			pyx4<-c(pyx4,out[["par"]][2])##beta_rho_s
    }

	} 
	

  pysxs0<-rep(0,m)
  pysxs1<-rep(0,m)
  
  tmp1<-pgamma(pyx1,pyx2)*gamma(pyx1+caseCount)/gamma(pyx1)/gamma(pyx2)*gammadiv(na-caseCount+pyx2,pyx1+pyx2+na)
  tmp2<-pgamma(pyx3,pyx4)/gamma(pyx3)/gamma(pyx4)*gamma(pyx3+controlCount)*gammadiv(nc+pyx4-controlCount,pyx3+pyx4+nc)
  pysxs1<-tmp1*tmp2	
  pysxs0<-pgamma(pyx1,pyx2)*gamma(pyx1+allCount)/gamma(pyx1)/gamma(pyx2)*gammadiv(n-allCount+pyx2,pyx1+pyx2+n)
  ##print(PYsXs)

	
	theta_new<-numeric()
	rho_new<-numeric()
	for(i in 1:m)
	{
		p<-c(pyx1[i],pyx2[i],caseCount[i])
		out<-optimize(pie_thetarho,c(0,0.5),p=p,maximum=TRUE)
		theta_new<-c(theta_new,out[["maximum"]])
    ##print(i)##OK!
	}
	for(i in 1:m)
	{
		p<-c(pyx3[i],pyx4[i],controlCount[i])
		out<-optimize(pie_thetarho,c(0,0.5),p=p,maximum=TRUE)
		rho_new<-c(rho_new,out[["maximum"]])
		##print(i)##ok
	}
	
  
	pxr<-c(0.5,0.5,0.5,0.5)

	out<-optim(pxr,est2,x1=cxe,x2=cxb,n1=ce,n2=cb)
	
  lxr<-exp(-out[["value"]])
	
  pxr<-c(out[["par"]][1],out[["par"]][2],out[["par"]][3],out[["par"]][4])
	pxsrs<-c(0.95,0.5,0.05,0.5)
	
	
	##pxsrs[1]=P(Xs=0|Rs=0)
	##pxsrs[2]=P(Xs=0|Rs=1)
	##pxsrs[3]=P(Xs=1|Rs=0)
	##pxsrs[4]=P(Xs=1|Rs=1)
	out<-optimize(pie_2,c(0,0.5),p=c(pxr[1],pxr[2],cxe),maximum=TRUE)
  
	pxsrs[4]<-out[["maximum"]]
	
  
	out<-optimize(pie_2,c(0,0.5),p=c(pxr[3],pxr[4],cxb),maximum=TRUE)
	
	pxsrs[3]<-out[["maximum"]]
	pxsrs[1]<-(1-pxsrs[3])
	pxsrs[2]<-(1-pxsrs[4])
  
  

  ##total probability formula
  px0<-pxsrs[1]*pb+pxsrs[2]*pe
  px1<-pxsrs[3]*pb+pxsrs[4]*pe

	
  
  P2<-numeric()
	P1<-0
	P4<-numeric()
	P3<-0
  ##pxsrsx=P(Xs|Rs)*Xs    every site's wps= sum(Wss*P(Xs|Rs)*Xs)
  pxsrsx<-rep(0,m)
  wpx<-rep(0,m)
  temp_pxsrsx<-0
  temp_wpx<-0
  
  ##every site's wr= sum(Wss*Rs)
	wr<-rep(0,m)
  temp_wr<-0
  
	
	for(i in 1:m)
	{
	  temp_pxsrsx<-0
	  temp_wpx<-0
	  temp_wr<-0
		for(j in 1:m)
		{
			if(r[i]==0)
			{
        pxsrsx[i]<-x[i]*pxsrs[3]
			}
      else
      {
        pxsrsx[i]<-x[i]*pxsrs[4]        
      }
      if(i!=j)
			{
				P1<-(P1+Wss[i,j]*Iss[i,j]*x[j])
				P3<-(P3+Wss[i,j]*r[j]*Iss[i,j])
				
        temp_wr<-temp_wr+Wss[i,j]*r[i]
			
				##no need to plus
        if(x[i]==0&&r[i]==0)
				{
				  ##temp_pxsrsx<-temp_pxsrsx+pxsrs[0]*x
				}
        ##no need to plus
				if(x[i]==0&&r[i]==1)
				{
				  ##temp_pxsrsx<-temp_pxsrsx+pxsrs[0]*x
				}
				
				if(x[i]==1&&r[i]==0)
				{
				  temp_wpx<-temp_wpx+Wss[i,j]*pxsrs[3]
				}
				if(x[i]==1&&r[i]==1)
				{
				  temp_wpx<-temp_wpx+Wss[i,j]*pxsrs[4]
				}
			}
		}
    wpx[i]<-temp_wpx
    wr[i]<-temp_wr
		P2<-c(P2,P1)
		P4<-c(P4,P3)
	}
	
  ##
	px<-c(-1000,-1000)
	pr<-c(-1,-1)
	

  
  est6_7<-function(p)
  {
    temp1<-sum(p[1]*pxsrsx+p[2]*wpx)
    temp2<-sum(log(exp(p[2]*wpx)+exp(p[1]*pxsrsx+p[2]*wpx)))
    obj<-temp1-temp2
    -obj
  }
	
	out<-optim(c(0,0),est6_7,control=list(maxit=50))
  px<-c(out[["par"]][1],out[["par"]][2])
	
	PsXs1Xns<-exp(px[1]*pxsrsx+px[2]*wpx)/(exp(px[2]*wpx)+exp(px[1]*pxsrsx+px[2]*wpx))
	PsXs0Xns<-exp(px[2]*wpx)/(exp(px[2]*wpx)+exp(px[1]*pxsrsx+px[2]*wpx))
  
  ##print(PsXsXns)
 
  
  est5_2<-function(p)
  {
    temp1<-sum(p[1]*r+p[2]*wr)
    temp2<-sum(log(exp(p[2]*wr)+exp(p[1]*r+p[2]*wr)))
    obj<-temp1-temp2
    -obj
  }
  
	out<-optim(c(0,0),est5_2,control=list(maxit=50))
 
	pr<-c(out[["par"]][1],out[["par"]][2])
	PsRsRns<-1/(1+exp(pr[2]*wr-pr[1]*r-pr[2]*wr))   ##  1/(1+exp(b-a))

  
  PXs0Y<-pysxs0*PsXs0Xns
	PXs1Y<-pysxs1*PsXs1Xns

  
	
	for(i in 1:m)
	{
	  if(is.nan(PXs0Y[i])||is.nan(PXs1Y[i]))
	  {
      if(PsXs0Xns[i]<=PsXs1Xns[i])
          x[i]<-1
	  }
    else
    {
      if(PXs0Y[i]<=PXs1Y[i])
        x[i]<-1
      else x[i]<-0
    }
    
	}

  
  
  PRsX<-lxr*PsRsRns

	

	meanr<-mean(PRsX)

	for(i in 1:m)
	{

		if(PRsX[i]>=meanr)
			r[i]<-1
    else
			r[i]<-0
  }
	
	write(x,paste(tmpfile_path,path,"_newx.seq",sep=""),nc=m)
	write(r,paste(tmpfile_path,path,"_newr.seq",sep=""),nc=m)

	if(TRUE)
	{
	  break
    if(count==1)
		{
			con1<-( (abs(pre_pr[1]-pr[1]) / pre_pr[1]) < par)
			con2<-( (abs(pre_pr[2]-pr[2]) / pre_pr[2]) < par)
			con3<-( (abs(pre_px[1]-px[1]) / pre_px[1]) <par)
			con4<-( (abs(pre_px[2]-px[2]) / pre_px[2]) <par)

			con5<-( (abs(pre_pxr[1]-pxr[1])/pre_pxr[1]) <par)
			con6<-( (abs(pre_pxr[2]-pxr[2])/pre_pxr[2]) <par)
			con7<-( (abs(pre_pxr[3]-pxr[3])/pre_pxr[3]) <par)
			con8<-( (abs(pre_pxr[4]-pxr[4])/pre_pxr[4]) <par)
			
			if(con1&&con2&&con3&&con4&&con5&&con6&&con7&&con8)
			{	break}
			else
			{
				count<-count+1
				pre_pr<-pr
				pre_px<-px
				pre_pxr<-pxr
			}
		}
		else
		{
			count<-count+1
			pre_pr<-pr
			pre_px<-px
			pre_pxr<-pxr
		}
	}
  
	

}


source("./r/P_value.r")



temp1<-0
temp2<-0
for (i in 1:m)
{
	if(x[i]==1)
	{
		temp1<-temp1+(maf_case[i]-maf_control[i])/maf_all[i]
		temp2<-temp2+((1-maf_all[i])/maf_all[i])
	}
	
}

z<-2*temp1/(sqrt(2/n)*sqrt(temp2))

p_value<-P_value(pnorm,z,side=0)

write(x,paste(tmpfile_path,path,"_resultX.seq",sep=""),nc=m,sep="")
write(r,paste(tmpfile_path,path,"_resultR.seq",sep=""),nc=m,sep="")


write(p_value,paste(tmpfile_path,path,"_P_VALUE",sep=""))

write(z,paste(tmpfile_path,path,"_STATISTIC",sep=""))


