validation<-function(est,tab_down, scale,k=k, pert,pert_underscale,h.int=h.int,threshold,thr_va=0.5,n=n)
{
  pert_underscale<-as.vector(pert_underscale)
  pert<-as.vector(pert)
  
  tab<-NULL
  for(g in pert)
  {
    submatrix<-subset(tab_down, get(scale[1],pos=tab_down)==g)
    intra=unique(get(scale[2],pos=submatrix))
    index= which(est$names_lev[[k+1]] %in% intra) 
    dd=est$med_lev[[k+1]][index,]
    outcome=apply(dd,1,moment.estimation,h.int=h.int,th=threshold[2],n=n)

    x=unlist(lapply(outcome,"[",1)) ##start
    #x1=unlist(lapply(outcome,"[",2)) ##start before (other strategy not considered yet)
    x2=unlist(lapply(outcome,"[",3)) ##end
    #x3=unlist(lapply(outcome,"[",4)) ##end of reaction (not considered yet)
    x4=unlist(lapply(outcome,"[",5)) ##intensity
    
    if (length(unique(round(x)))==1){
      pct2=1
      tab<-rbind(tab,data.frame(name=g,start=median(x),end=median(x2),intensity=median(x4)))
      
    } else {
      
      mc <- Mclust(round(x),verbose=FALSE)
      tmc=data.frame(table(mc$classification)) #mc table

      tmc$name<-rep(g,length(tmc$Var1))
      tmc$intensity<-tmc$end<-tmc$start<-rep(0,length(tmc$Var1))
      
      if (nrow(tmc)!=1){ ###Do not count 2 times a same disturbance having 2  starting points in a same cluster 
        for (cla in tmc$Var1){
          names_cla=names(mc$classification[mc$classification==cla]) ##starting point in class cla
          bb= matrix(unlist(strsplit(names_cla,".", fixed=TRUE)),2,length(names_cla))
          
          tmc$Freq[tmc$Var1==cla]=length(unique(bb[1,]))
          tmc$start[tmc$Var1==cla]<-median(x[mc$classification==cla])
          #tmc$M1[tmc$Var1==cla]<-median(x1[mc$classification==cla])
          tmc$end[tmc$Var1==cla]<-median(x2[mc$classification==cla])
          #tmc$M3[tmc$Var1==cla]<-median(x3[mc$classification==cla])
          tmc$intensity[tmc$Var1==cla]<-median(x4[mc$classification==cla])
          
        } #end for
        tmc$pct=tmc$Freq/length(intra)
      }else{
        tmc$pct=1
        tmc$start<-median(x)
        #tmc$M1<-median(x1)
        tmc$end<-median(x2)
        #tmc$M3<-median(x3)
        tmc$intensity<-median(x4)
        
      }
      tmc$pct2<-tmc$pct>=thr_va ##checked perturbations
      #if(sum(tmc$pct2)==0) cat(g,'\n')
    
    tmc<-tmc[tmc$pct2==TRUE,c("name","start","end","intensity")]
    if(nrow(tmc)==1) ####use the median median phenotype instead of median by class
    {
      index= which(est$names_lev[[k]] %in% g) 
      dd=est$med_lev[[k]][index,]
      #T=as.numeric(colnames(dd))

      outcome=apply(dd,1,moment.estimation,h.int=h.int,th=threshold[1],n=n)

      if(length(unlist(lapply(outcome,"[",1)))==1) tmc[2:4]<-unlist(outcome)[c(1,3,5)]
      if(length(unlist(lapply(outcome,"[",1)))>1) ###if multi disturbance in the median-median phenotype, use the nearest from the median by class
      {
        tempo<- rbind(unlist(lapply(outcome,"[",1)),unlist(lapply(outcome,"[",3)))
        tmc[2:3]<-tempo[,which.min(apply(tempo,2, function(x) {sum((x-tmc[2:3])^2)}))]
      }
    }
    
    tab<-rbind(tab,tmc)
    if(nrow(tmc)>=1){
      multip_pert<-matrix(unlist(strsplit(names(x),".", fixed=TRUE)),2,length(names(x)))
      nmp<-as.vector(names(which(table(multip_pert[1,])>1)))
      pert_underscale<-c(pert_underscale,nmp)
    }
    if(nrow(tmc)==0){
     rec_pert<-as.vector(intra[which(intra %in% est$mixing[[k+1]])])
     pert_underscale<-c(pert_underscale,rec_pert)
    }
  }#end else
  }#end pert
  
  return(list(tab=tab,pert_u=pert_underscale))
  
}
