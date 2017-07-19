
q1=function(tao)
{
 p1=exp(tao)
 return(p1/sum(p1))
}

Q=function(tao,n,n1)
{
 p=q1(tao)
 return(sample(1:n,n1,prob=p,replace=TRUE))
}


GD=function(x,tao,beta0,f0,f1,f2,f3,f4,L0,L,epsilion,n,n1,Method='Constant')
{
  M=0
  i=1
  beta1=beta0
  L1=f0(beta1)
  T=0
  n2=length(tao)
  t1=Sys.time()
  if(Method=='AD')
  {
  while(T[i]<epsilion)
  {
   d=0
   d1=0
   d2=0
   beta0=beta1
   tao1=f4(tao)
   index=Q(tao1,n,n1)
   for(u in 1:n1)
   {
    u2=1/(n*f2(index[u],tao1))
    d2=d2+u2    
    u1=f1(x[[index[u]]],beta0)*u2
    d=d+u1
    d1=d1+sum(u1^2)*f3(index[u],tao)
   }
   d=d/n1
   d1=d1/n1
   d2=d2/n1
   beta1=beta0-d2/L*d
   tao=tao+1/L0*d1
   i=i+1
   L1[i]=f0(beta1)
   T[i]=Sys.time()-t1
  }
  }
 if(Method=='Constant')
  {
  while(T[i]<epsilion)
  {
   L0=0
   d=0
   d1=0
   d2=0
   beta0=beta1
   tao1=f4(tao)
   index=Q(tao1,n,n1)
   for(u in 1:n1)
   {    
    u1=f1(x[[index[u]]],beta0)
    d=d+u1
   }
   d=d/n1
   beta1=beta0-1/L*d
   tao=tao
   i=i+1
   L1[i]=f0(beta1)
   T[i]=Sys.time()-t1
  }
  }
 return(list(i-1,beta1,T,L1,tao))
}
    
#n1/n high dimension for tao
f0=function(beta)
 return(sum((Y-X%*%beta)^2)/n)

f_1=function(beta)
 return(-2*t(X)%*%(Y-X%*%beta)^2/n)


f1=function(x,beta)
 return(-2*x[-1]*(x[1]-sum(x[-1]*beta)))

f2=function(index,tao)
 return(q1(tao)[index])

f3=function(index,tao)
{
 m=q1(tao)
 m[index]=m[index]-1
 return(-m)
}

f4=function(tao)
 return(tao)

f_3=function(index,tao1)
{
 m1=exp(tao1)
 m=0
 m[1]=n2*m1[1]/(n2*m1[1]+n3*m1[2])
 m[2]=n3*m1[2]/(n2*m1[1]+n3*m1[2])
 if(index %in% index1)
  m[1]=m[1]-1
 else
  m[2]=m[2]-1
 return(-m)
}

 

f_4=function(tao1)
{
 tao=vector()
 length(tao)=n
 tao[index1]=tao1[1]
 tao[index2]=tao1[2]
 return(tao)
}

a=read.table('hdata.txt')
X=as.matrix(a[,1:13])
Y=as.matrix(a[,14])
for(i in 1:13)
 X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i])
Y=(Y-mean(Y))/sd(Y)

index1=Y>=median(Y)
index2=Y<median(Y)

n2=sum(index1)
n3=sum(index2)


n=dim(Y)[1]
x=list()
for(i in 1:n)
 x[[i]]=as.vector(c(Y[i],X[i,]))


beta0=rep(0,13)
L=2*max(eigen(t(X)%*%X/n)$values)
L0=200
L_0=10000
tao0=rep(0,2)
epsilion=1

tao=rep(0,n)
n1=90


 M1=vector()
 M2=vector()
 M3=vector()
 M4=vector()
 for(i in 1:1)
 {
 T=GD(x,tao0,beta0,f0,f1,f2,f_3,f_4,L_0,L,epsilion,n,n1,Method='AD')
 T1=GD(x,tao,beta0,f0,f1,f2,f3,f4,L0,L,epsilion,n,n1,Method='AD')
 M1=c(M1,T[[4]])
 M2=c(M2,T[[3]])
 M3=c(M3,T1[[4]])
 M4=c(M4,T1[[3]])

 }




DR=function(T1,T2,T3,T4)
{
 par(mfrow=c(1,3))
 plot(T1~T2,xlab='Time',ylab='Objective Function')
 plot(T3~T4,col='red',xlab='Time',ylab='Objective Function')
 plot(1, type="n",xlim=c(-0.01,1),ylim=c(0,1),xlab='Time',ylab='Objective Function')
 lines(stats::lowess(T1~T2,f=0.01,iter=1000),col='black')
 lines(stats::lowess(T3~T4,f=0.01,iter=1000),col='red')

}

DR(M1,M2,M3,M4)

