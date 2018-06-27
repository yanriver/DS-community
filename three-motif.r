#3 species motif designed to indicating presence of DS interaction, 1 is DS, 0 is linear.
M_l<-list(
  M0=matrix(c(0,0,0,
              0,0,0,
              0,0,0),nrow=3, byrow = TRUE),
  M1.1=matrix(c(0,1,0,
                0,0,0,
                0,0,0),nrow=3, byrow = TRUE),
  M2.1=matrix(c(0,1,0,
                1,0,0,
                0,0,0),nrow=3, byrow = TRUE),
  M2.2=matrix(c(0,1,1,
                0,0,0,
                0,0,0),nrow=3, byrow = TRUE),
  M2.3=matrix(c(0,1,0,
                0,0,0,
                0,1,0),nrow=3, byrow = TRUE),
  M2.4=matrix(c(0,1,0,
                0,0,0,
                1,0,0),nrow=3, byrow = TRUE),
  M3.1=matrix(c(0,1,0,
                0,0,1,
                1,0,0),nrow=3, byrow = TRUE),
  M3.2=matrix(c(0,1,1,
                0,0,0,
                1,0,0),nrow=3, byrow = TRUE),
  M3.3=matrix(c(0,1,0,
                0,0,0,
                1,1,0),nrow=3, byrow = TRUE),
  M3.4=matrix(c(0,1,0,
                1,0,0,
                0,1,0),nrow=3, byrow = TRUE),
  
  M4.1=matrix(c(0,1,0,
                1,0,1,
                1,0,0),nrow=3, byrow = TRUE),
  M4.2=matrix(c(0,1,0,
                1,0,0,
                1,1,0),nrow=3, byrow = TRUE),
  M4.3=matrix(c(0,1,1,
                1,0,1,
                0,0,0),nrow=3, byrow = TRUE),
  M4.4=matrix(c(0,1,1,
                1,0,0,
                1,0,0),nrow=3, byrow = TRUE),
  M5.1=matrix(c(0,1,1,
                1,0,1,
                1,0,0),nrow=3, byrow = TRUE),
  M6.1=matrix(c(0,1,1,
                1,0,1,
                1,1,0),nrow=3, byrow = TRUE)
)

#function to get competition strength
comp<-function(x,y){
  y*y/(x+y)
}

#function to generate original linear competition,qu=range of species' ability,du=self feedback
CM<-function(S,qu=c(-0.5,-0.01),du=c(1,1)){
  #competition ability niche
  q=round(runif(S,qu[1],qu[2]),2)
  Q<-outer(q,q,FUN="comp")
  diag(Q)<-round(runif(S,du[1],du[2]),2)
  list(q=q,Q=Q)
}

