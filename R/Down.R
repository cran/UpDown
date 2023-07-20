
#################
####DOWN#########
#################

Down<-function(est,h.int=10, kappa=NULL, thr_va=0.5) 
{
  n<-ncol(est$yp)
  levels<-est$levels
  nl<-length(levels)
  
  #initalizations
  tab_down<-est$tab_up
  
  ###Thresholds
  thr<-est$thr
  
  est_lvl<-list(); lvl<-list()

  for(k in 1:nl)
      est_lvl[[k]]<-unique(get(levels[k],pos=subset(tab_down,det==levels[k])))


  va<-list(); 
  
  ###validation level 1->level 2
  if(nl>1)
    va[[1]]=validation(est,tab_down,scale=levels[1:2],k=1,pert=est_lvl[[1]],pert_underscale=est_lvl[[2]],h.int=h.int,threshold = thr[1:2],thr_va=thr_va,n=n)
  
  ###validation level k->level k+1
  if(nl>2)
    for(k in 2:(nl-1))
      va[[k]]=validation(est,tab_down,scale=levels[k:(k+1)],k=k,pert=va[[k-1]]$pert_u,pert_underscale=est_lvl[[k+1]],h.int=h.int,threshold = thr[k:(k+1)],thr_va=thr_va,n=n)
  
  
  est_ind<-if(nl>1) va[[nl-1]]$pert_u else est_lvl[[nl]]
  
  tab_ind<-NULL
  
  ##moment estimations for individual disturbance
  for(g in est_ind)
  {
    ind<-get(levels[nl],pos=subset(tab_down, get(levels[nl],pos=tab_down)==g))
    index= which(est$names_lev[[nl]] %in% ind)   
    dd=est$yp[index,]

    outcome<-moment.estimation(dd,h.int=h.int,th=thr[nl],n=n)
    
    npoint=length(outcome$M0)
    x=matrix(unlist(outcome),nrow=npoint,ncol=5,byrow=FALSE)
    x<-data.frame(ind=rep(g,npoint),x[,1],x[,3],x[,5])
    colnames(x)<-c(levels[nl],"start","end","intensity")
    tab_ind<-rbind(tab_ind,x)
  }
  
  
  
  ##############################
  pert_down<-vector("list", nl)
  names(pert_down)<-levels
  
  if(nl>1)
    for(k in 1:(nl-1)){
      tt=subset(tab_down, get(levels[k],pos=tab_down)%in%va[[k]]$tab$name)
      pert_down[[k]]<-unique(data.frame(mget(levels[1:k],envir =as.environment(tt))))
      colnames(pert_down[[k]])<-levels[1:k]
      if(nrow(pert_down[[k]])>0)  
      {
      colnames(va[[k]]$tab)[1]<-levels[k]
      pert_down[[k]]<-merge(pert_down[[k]],va[[k]]$tab,by=levels[k])
      }
    }
  pert_down[[nl]]<-tab_ind
  tt=subset(tab_down, get(levels[nl],pos=tab_down)%in%tab_ind[,1])

    ########KAPPA concordance#####
  ##############################
  
 
  
  if(nl>1){
    pert_down[[nl]]<-merge(tt[,levels],tab_ind,by=levels[nl])
    
    if(!is.null(kappa))
      for(k in 1:(nl-1))
        for(l in (k+1):nl)
          if(!is.null(pert_down[[k]]) & !is.null(pert_down[[l]]))
            pert_down[[l]]<-Kappa_cohen(tab1=pert_down[[k]],tab2=pert_down[[l]],var=levels[k],n=n,thr=kappa)
  }
  
  return(list(pert_down=pert_down,levels=levels))
}

