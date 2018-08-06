#source three motifs and some functions
source("../three-motif.R")

#function to determine the outcomes of interaction: LC or DS
NM<-function(M,A,N=n,b,c){
  res<-matrix(NA,nrow=nrow(M),ncol=ncol(M))
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if(M[i,j]==0){res[i,j]<-A[i,j]}else{
        res[i,j]<-A[i,j]*(b*N[j]-c)
      }
    }
  }
  res
}



#main function to simulation local community dynamics
sim3<-function(M_l,b=1,c=runif(1,1,2),qu=c(0.4,0.8)){
  
  #intrsinc growth rate
  R<-rep(runif(1,0.5,1.5),3)
  #matrix indicating NM occuring
  
  #orginal competition matrix
  A<-CM(3,qu=qu)$Q
  
 
  EQ_l<-vector("list",length=16)
  n0<-runif(3,0,5)
  for (i in 1:16){
    parms<-list(R=R,A=A,M=M_l[[i]],b=b,c=c)
    model  <-  function(t,  n,  parms)  {
      with(parms,  {
        
        dn  <- n* (R - NM(M,A,n,b,c)%*%n)
        return(list(dn))
      })
    }
    library(deSolve)
    
    times<-seq(0,200,0.01)
    out  <- tryCatch(ode(y=n0, times=times, func=model, parms=parms,method="rk4"),error=function(e) NULL)
    EQ_l[[i]]<-tail(out[,-1],1)
    
  }
  
  
  return(list(R=R,n0=n0,q=CM(3,qu=qu)$q,b=b,c=c,EQ_l=do.call(rbind,EQ_l)))
}

library(parallel);library(pbapply)
cl <- makeCluster(47L)
clusterExport(cl, c("M_l","CM","comp","sim3"))
res_l<-pblapply(1:1000,function(x){sim3(M_l=M_l)}, cl =cl)
stopCluster(cl)

###########################meta competition community#########################
#meta community interaction matrix linking LC+DS community
MMMM_l<-vector("list",length=16)
for (i in 1:16){
  MM<-as.matrix(Matrix::bdiag(M_l[[1]],M_l[[i]]))
  MMMM_l[[i]]<-as.matrix(Matrix::bdiag(MM,MM))
}

#main function to simulation LC+DS community dynamics
sim33<-function(para){
  
  #intrsinc growth rate,replicate automotically
  R=para$R
  #matrix indicating NM occuring
  
  #orginal competition matrix
  qq<-c(para$q,para$q)
  QQ<-outer(qq,qq,comp)
  diag(QQ)<-1
  AAAA<-as.matrix(Matrix::bdiag(QQ,QQ))
  #threshold paramter 1
  
  #threshold paramter 2
  
  #dispersal
  D<-diag(runif(6,0.8,1.2))
  DD<-rbind(cbind(-D,D),cbind(D,-D))
  
  #the model
  model  <-  function(t,  n,  parms)  {
    with(parms,  {
      
      dn  <- n* (R - NM(M,A,n,b,c)%*%n)+DD%*%n
      return(list(dn))
    })
  }
  
  EQ_cl<-vector("list",length=16)
  for (i in 1:16){
    #i=1
    parms<-list(R=R,A=AAAA,M=MMMM_l[[i]],b=para$b,c=para$c)
    
    library(deSolve)
    times<-seq(0,200,0.01)
    n0<-c(para$EQ_l[1,],c(0,0,0),c(0,0,0),para$EQ_l[i,])
    out  <- tryCatch(ode(y=n0, times=times, func=model, parms=parms,method="rk4"),error=function(e) NULL)
    EQ_cl[[i]]<-tail(out[,-1],1)
  }
  
  return(list(EQ_cl=do.call(rbind,EQ_cl),DD=DD))
}


library(parallel);
cl <- makeCluster(48L)
clusterExport(cl, c("MMMM_l","NM","sim33","res_l","comp"))
res_cl<-pblapply(1:1000,function(x){sim33(res_l[[x]])}, cl =cl)
stopCluster(cl)

###linking any two DS communities

#main function to simulation DS+DS community dynamics
sim33_anytwo<-function(para,nets=c(1,2)){
  #combin two NM matrix
  MM<-as.matrix(Matrix::bdiag(M_l[[nets[1]]],M_l[[nets[2]]]))
  MMMM<-as.matrix(Matrix::bdiag(MM,MM))
  
  #intrsinc growth rate
  R=para$R
  #matrix indicating NM occuring
  
  #orginal competition matrix
  qq<-c(para$q,para$q)
  QQ<-outer(qq,qq,comp)
  diag(QQ)<-1
  AAAA<-as.matrix(Matrix::bdiag(QQ,QQ))
  
  #dispersal
  D<-diag(runif(6,0.8,1.2))
  DD<-rbind(cbind(-D,D),cbind(D,-D))
  
  #the model
  model  <-  function(t,  n,  parms)  {
    with(parms,  {
      
      dn  <- n* (R - NM(M,A,n,b,c)%*%n)+DD%*%n
      return(list(dn))
    })
  }
  
  parms<-list(R=R,A=AAAA,M=MMMM,b=para$b,c=para$c)
  
  library(deSolve)
  times<-seq(0,200,0.01)
  n0<-c(para$EQ_l[nets[1],],c(0,0,0),c(0,0,0),para$EQ_l[nets[2],])
  out  <- tryCatch(ode(y=n0, times=times, func=model, parms=parms,method="rk4"),error=function(e) NULL)
  
  #return(list(EQ_cl=tail(out[,-1],1),DD=DD))
  out
}

library(parallel);library(pbapply)
cl <- makeCluster(47L)
clusterExport(cl, c("NM","sim33_anytwo","res_l","comp","M_l"))
res_cl0203<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(2,3))}, cl =cl)
res_cl0403<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(4,3))}, cl =cl)
res_cl0503<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(5,3))}, cl =cl)
res_cl0603<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(6,3))}, cl =cl)
res_cl1014<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(10,14))}, cl =cl)
res_cl1114<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(11,14))}, cl =cl)
res_cl1214<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(12,14))}, cl =cl)
res_cl1314<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(13,14))}, cl =cl)
res_cl1516<-pblapply(1:1000,function(x){sim33_anytwo(res_l[[x]],c(15,16))}, cl =cl)
stopCluster(cl)


