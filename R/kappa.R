
Kappa_cohen<-function(tab1, tab2,var,n,keep=TRUE,thr=0.75)
{
inter<-intersect(tab1[,var],tab2[,var])
for(g in inter){
  mat1<-tab1[tab1[,var] == g,]
  mat2<-tab2[tab2[,var] == g,]

  kc=matrix(rep(0,nrow(mat2)*nrow(mat1)),nrow(mat1),nrow(mat2))

  for(cpp in 1:nrow(mat1)){
    x=floor(c(mat1[cpp,]$start:mat1[cpp,]$end))
    
    for(rmt in 1:nrow(mat2)){
      y=floor(c(mat2[rmt,]$start:mat2[rmt,]$end))
      tot=unique(c(x,y))
      a1=length(intersect(x,y)) #++
      a2=length(setdiff(x,y)) #+-
      a3=length(setdiff(y,x)) #-+
      a4=n-length(tot) #--
      #
      xtab <- as.table(rbind(c(a1, a2),c(a3, a4)))
      kc[cpp,rmt]=((sum(diag(xtab))/sum(xtab))-sum((rowSums(xtab)/sum(xtab))*(colSums(xtab)/sum(xtab))))/(1-sum((rowSums(xtab)/sum(xtab))*(colSums(xtab)/sum(xtab))))
    }
  }
    for(rmt in 1:nrow(mat2)){
      if(max(abs(kc[,rmt]))>=thr & keep==TRUE) 
        tab2<-tab2[-which(mat2[rmt,]$start==tab2$start & mat2[rmt,]$end==tab2$end  & mat2[rmt,]$intensity==tab2$intensity ),]
    }
    
     if(keep==FALSE & nrow(mat2)>1) 
       {
       ikc=which.max(kc)
       for(j in (1:nrow(mat2))[-ikc])
       tab2<-tab2[-which(mat2[j,]$start==tab2$start & mat2[j,]$end==tab2$end),]
     }
  
  
}
return(tab2)
}
