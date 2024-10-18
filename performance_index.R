index_band3=function(x, y, c, r, S, bse, bs)
{
  m=c%*%t(bse)
  lu=(m+r*S)
  lw=(m-r*S)
  
  luc=Data2fd(x, t(lu), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(luc-c)
    return(w.ff)
  }
  range=range(x)
  b=integrate(fint, lower=range[1], upper=range[2], stop.on.error=F)[[1]]
  
  lwc=Data2fd(x, t(lw), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(c - lwc)
    return(w.ff)
  }
  range=range(x)
  a=integrate(fint, lower=range[1], upper=range[2], stop.on.error=F)[[1]]
  Width=b+a
  
  B = eval.basis(x, bs)
  yf=y%*%t(B)
  a=lw-yf
  b=yf-lu
  #Cova=1-((length(a[a>0])!=0)||(length(b[b>0])!=0))
  Cova=covAG
  fint=function(x){
    #browser()
    B = eval.basis(x, bs)
    a=B%*%(luc-lwc)
    b=B%*%(lwc-y)
    b[b<0]=0
    c=B%*%(y-luc)
    c[c<0]=0
    return(a+(2/alpha)*b+(2/alpha)*c)
  }
  range=range(x)
  IntScore=integrate(fint, lower=range[1], upper=range[2], stop.on.error=F)[[1]]
  
  #browser()
  ### covalpha funzionale
  inout=matrix(0, ncol=1, nrow=length(x))
  a=lw-yf
  b=yf-lu
  inout[((1:length(x))[(a<0)&(b<0)]) ]=1
  Covat=mean(covAL)
  
  #####   indici che mi servono
  return(list(c(Cova, Width, IntScore), Covat))
  #
}

n=nrow(DOE0)
covt_buff=matrix(0, n, lt)
I=matrix(0, n, 3)
colnames(I)<-c("CovAlG %","Width","InterScore")
#bse=eval.basis(argavls, bse)
#Ps0C=Data2fd(argvals , Ps0s, bsb)$coef
for (i in 1:n) {
  buff1=index_band3(argvals, CfD0[,i], As[,i], raggi[i], Ss[,i], bse, b1.2)  
  I[i,]=buff1[[1]]
  covt_buff[i,]=buff1[[2]]
}


dfres_fr<-("CovAlG%"=median(I[,1]))
dfres_fr<-(cbind(dfres_fr,"Width" = median(I[,2])))
dfres_fr<-(cbind(dfres_fr,"Int"=median(I[,3])))

dfres_fr=cbind("CovAlL%"=median(rowMeans(covt_buff)), dfres_fr)
dfres_fr=cbind(dfres_fr,"TimeT"=sum(temp))
dfres_fr=cbind(dfres_fr,"TimeM"=mean(temp))

colnames(dfres_fr)=c("CovAlL%", "CovAlG%", "Width", "Int", "TimeT", "TimeM")

print(dfres_fr)