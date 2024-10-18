source('function.R')

data(maritimes.avg)
data(maritimes.coords)
data(maritimes.data)

alpha=0.05 
Functions0=maritimes.data
DOE0=as.data.frame(maritimes.coords)
n<-dim(Functions0)[1]
argvals<-seq(1,n, by=1)

n<-dim(Functions0)[2]
nb = 5:n
mdata = fdata(t(Functions0))
out<-optim.basis(mdata,lambda=0,numbasis=nb,type.basis="fourier")
# out <- optim.basis(mdata,
#                    lambda = 0,
#                    numbasis = nb,
#                    type.basis = "bspline")
k=65
b1.2 <- create.fourier.basis(rangeval =range(argvals) , nbasis=k)
#  fit the data without smoothing
fd1.2 <- Data2fd(argvals=argvals, Functions0, basisobj=b1.2)
CfD0=fd1.2$coefs
M=getbasispenalty(b1.2, 0)
L2ns0=l2.norm(n, fd1.2, M)
bse=eval.basis(argvals, b1.2)

fn=35
set.seed(2)
lt=length(argvals)
zz=dim(DOE0)[1]
verifica=matrix(0, 1, zz) #tutti8 zeri
verifica2l=matrix(0, 1, zz)
raggi=matrix(0, 1, zz) #semiampiezza
temp = matrix(0,nrow = fn ,ncol = 1)
Ps0s=matrix(0, lt, zz)
As=matrix(0, k, zz)
Ss=matrix(0, lt, zz)

pb <- progress_bar$new(
  format = paste0("conformal prediction [:bar] :percent eta: :eta"),
  total = fn,
  clear = FALSE,
  width = 60
)
for (z in 1:fn) {
  s0=z
  coord.cero <- matrix(DOE0[s0,], nrow=1, ncol=2)
  rd=Functions0[,s0]
  
  Functions=Functions0[,-s0]
  DOE=as.matrix(DOE0[-s0,])
  L2ns=L2ns0[-s0,-s0]
  CfD=CfD0[,-s0]
  
  n<-dim(Functions)[2]
  
  
  buffD=rbind(DOE0[s0,], DOE)
  buffD=as.matrix(dist(buffD))
  buffD=buffD[1,-1]
  train=which(buffD<median(buffD))
  #train = which(buffD<quantile(buffD,0.75))
  test=(1:n)[-train]
  l=length(test)
  lt= length(argvals)
  ## -- Conformal Prediction -  ##
  NCS=matrix(0, nrow = l, ncol = 1)
  
  pbuff=matrix(0, nrow = lt, ncol = l)
  start<-Sys.time()
  index_Crd=train
  index_X=index_Crd
  
  okfd.res<-okfd2(new.coords=coord.cero, coords=DOE[index_Crd,],
                  data=Functions[,index_X], nbasis=k, argvals=argvals,
                  fix.nugget=TRUE, L2norm=L2ns[index_X,index_X] , coefD=CfD[,index_X])
  Ps0 = okfd.res$krig.new.data
  plot(okfd.res)
  
  for (j in (1:l)) {
    index_Crd=c(train, test[j])
    index_X=index_Crd
    okfd.res<-okfd2(new.coords=coord.cero, coords=DOE[index_X,],
                    data=Functions[,index_X], nbasis=k, argvals=argvals,
                    fix.nugget=TRUE, L2norm=L2ns[index_X,index_X] , coefD=CfD[,index_X])
    pbuf=okfd.res$krig.new.data
    pbuff[,j]=bse%*%pbuf
#    cat(" - ",j)
  }
  S=matrix(0, nrow = lt, ncol = 1)
  A = Ps0
  Ps0=bse%*%Ps0
  for (j in 1:lt) {
    S[j]=max(abs((pbuff[j, ] - Ps0[j])), na.rm =T)
    
  }
#  S=sqrt(rowMeans((pbuff - matrix(Ps0, nrow = nrow(pbuff), ncol = ncol(pbuff)))^2))
  
  for (j in 1:l) {
    NCS[j]=max(abs((pbuff[,j ] - Ps0)/S), na.rm =T)
    
  }
  
  # for (j in 1:l) {
  #   NCS[j]=mean((pbuff[,j ] - Ps0)^2/S)
  #   NCS[j] = sqrt(NCS[j])
  # }
  end<-Sys.time()
#  cat("\n")
#  print(end-start)
  if(1){tempo_trasc=(end-start)}else{tempo_trasc=tempo_trasc+(end-start)}
#  print(paste0("estimated end at ", Sys.time()+(zz-z)*tempo_trasc))
  temp[z] = tempo_trasc
  alpha=0.05 
  r=(quantile(NCS, 1-alpha))
  low=Ps0 - r*S
  a=low - rd
  up=Ps0 + r*S
  b= rd - up
  verifica[z]=((length(a[a>0])!=0)||(length(b[b>0])!=0))
  raggi[z]=r
  As[,z] = A
  Ps0s[,z]=Ps0
  Ss[,z]=S
  verifica2l[z]=length(a[a>0]) + length(b[b>0])
  pb$tick()
}


###
### verifica2lp è il cov alfa funzionale!
###
covAG=mean(1- verifica)*100
covAL=(1-verifica2l/lt) * 100
covAG
covAL
hist(covAL)
r=(quantile(NCS, 1-alpha))
low=Ps0 - r*S
up=Ps0 + r*S
plot(argvals, Ps0, col=1, lwd=2, ylim = range(c(up,low, rd)),
     type="l", lty=1, main="Predictions", xlab="Day",
     ylab="Temperature (Degrees C)")
lines(argvals, up,col="red")
lines(argvals, low,col="blue")
lines(argvals, rd , type="p", pch=20,cex=0.5, col=2, lwd=1)

#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL.RData") sup sup
#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL_25.RData") #25
#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL_75.RData") #75

#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL2.RData") rad sup
#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL3.RData") rad rad
#save(Functions0,covAG,covAL,raggi,Ps0s,Ss,file = "LOL4.RData") #sup #rad
