########################
###########UP###########
########################


Up<-function(data,levels,obs,vtime, h.int=NULL,mixplot=FALSE, correction=NULL,
             options=list())
{
 
  ##replace colnames with space by '.'
  colnames(data)<-gsub(" ",".",colnames(data))
  levels<- gsub(" ",".",levels)
  nl<-length(levels)
  
  #default control parameters
  con <- list(epsilon = 1e-08, maxit = 500, maxrestarts=20,verb = FALSE, fast=FALSE, ECM = FALSE, arbmean = TRUE, arbvar = TRUE, seed=NULL, notest=FALSE, minobs=20)
  namc <- names(con)
  con[namc <- names(options)] <- options
  
  
  if(!con$notest){
  if(!obs %in% colnames(data)) stop("you need to specify a existing obs",call. = FALSE)
  if(!vtime %in% colnames(data)) stop("you need to specify a existing time variable",call. = FALSE)
  }
  for(i in 1:nl)
  {
      tt<-try(get(levels[i],pos=data),silent=TRUE)
      if(inherits(tt,"try-error")) stop("you need to specify a existing level",call. = FALSE)
      if(!inherits(tt,"try-error")) assign(paste(levels[i]),tt)
  }
  
  if(!con$notest){
  data %>% group_by(get(levels[nl])) %>%
    summarise(nb=n()) %>%
    unique() %>% data.frame() -> anim_out
  
  colnames(anim_out)[1]="ID"
  ao=anim_out[anim_out$nb<con$minobs,]$ID 

    if(length(ao)>=1)
    for(i in 1:length(ao)){
      warning(paste(levels[nl],ao[i],"has been removed because the trajectory contains less than ", con$minobs, " observations \n"))
    }
  data=data[!get(levels[nl],pos=data)%in%ao,]
  }
  

  
  if(!con$notest){
  if(!is.numeric(get(obs,pos=data))) stop(paste(obs,"is not numeric"))
  if(nl>1)
    for(i in 1:(nl-1)){
      if(length(unique(get(levels[i])))<4) stop(paste("Not enough elements in",levels[i],"level. Minimum allowed is 4",sep=" "),call. = FALSE)
      }
  

  if(nl>1)
    for(i in nl:2){
      if(nrow(unique(data.frame(mget(levels[1:i],as.environment(data))))) != length(unique(get(levels[i],pos=data))))
        stop(paste("wrong order for the hierarchical levels or duplicated names in",levels[i],"level",sep=" "),call. = FALSE)
    }
  }

  if(!is.null(correction)){
    for(j in 1:length(correction))
      if(!correction[j] %in% colnames(data))
        stop("you need to specify an existing 'correction' label",call. = FALSE)
    data<-trend_correction(data,obs,correction)
    
    
  }else  data$cobs<-get(obs,pos=data)
  ytest<-table(data[,c(levels[nl],vtime)])
  
  if(sum(apply(ytest[,-1],2,max)>1)) stop(sum(apply(ytest[,-1],2,max)>1)," ",levels[nl]," have duplicate times",call. = FALSE)
  y=spread( data[,c(levels[nl],"cobs",vtime)] ,get(vtime),"cobs")
  
  ind<-y[,1]
  y<-y[,-1]
  rownames(y)<-ind
  n<-max(apply(y,1,function(x) length(x[!is.na(x)])))
  
  if(is.null(h.int)) h.int<-sqrt(n)
  
  df.med_lev<-med_lev<-names_lev<-sc.x<-sc.y<-dY<-dX<-Vmin<-start<-classification<-mixmdl<-vector("list", nl)
  names(df.med_lev)<-names(med_lev)<-names(names_lev)<-names(sc.x)<-names(sc.y)<-names(dY)<-names(dX)<-names(Vmin)<-names(start)<-names(classification)<-names(mixmdl)<-levels
  med_levk<-data
  df.med_lev[[nl]]<-data
  
  if(nl>1)
    for(k in (nl-1):1){
      med_levk %>% group_by(get(levels[k],pos=med_levk),get(vtime)) %>%
        mutate(nb=n(),cobs=median(get("cobs"))) %>%
        summarise_at(c(levels[max(0,(k-1)):1],"nb","cobs"),unique)%>%
        data.frame() -> med_levk
        colnames(med_levk)[1:2]<-c(levels[k],vtime)

      if(sum(med_levk$nb<=(median(med_levk$nb)*0.5))>0) #if existing
      med_levk[med_levk$nb<=(median(med_levk$nb)*0.5),]$cobs<-NA
      
      med_lev[[k]]<-spread( med_levk[,c(levels[k],"cobs",vtime)] ,get(vtime),"cobs")
      
      names_lev[[k]]<-med_lev[[k]][,1]
      med_lev[[k]]<-med_lev[[k]][,-1]
      rownames(med_lev[[k]])<-names_lev[[k]]
      
      med_levk=med_levk[complete.cases(med_levk), ]##remove NA row
      
      df.med_lev[[k]]<-med_levk
    }

  names_lev[[nl]]<-ind
  med_lev[[nl]]<-y
  thr<-rep(0,nl)
  names(thr)<-levels

  for(k in 1:nl)
  {
    ri<-apply(apply(med_lev[[k]],1,complete.cases),2,sum)
    if(length(which(ri<2))>0)
    {    
    for(l in which(ri<2)){  
    if(k>1)  
    warning(paste(levels[k],names_lev[[k]][l],"has been removed due to complete NA in the median trajectory or too few observations"))
    if(k==nl)  
    warning(paste(levels[k],names_lev[[k]][l],"has been removed due to complete NA in the trajectory"))
    }
      names_lev[[k]]<-names_lev[[k]][-which(ri<2)]  
      med_lev[[k]]<-med_lev[[k]][-which(ri<2),]  
    }
  }
  
  for(k in 1:nl){
    N<-length(names_lev[[k]])
    sc.x[[k]] <- sc.y[[k]] <- matrix(rep(0,floor(h.int*n)*N),N,floor(h.int*n))
    dY[[k]]<-dX[[k]]<-matrix(rep(0,(floor(h.int*n)-1)*N),N,floor(h.int*n)-1)
    
    Vmin[[k]]<-start[[k]]<-rep(0,N)
    
    for(i in 1:N){
      t<-smoothing.curve(med_lev[[k]][i,],h.int=h.int,n=n)
      sc.y[[k]][i,]<-t$y
      sc.x[[k]][i,]<-t$x
      Vmin[[k]][i]<-t$Vmin
      start[[k]][i]<-t$start
      dY[[k]][i,]<-t$dY
      dX[[k]][i,]<-t$dX
    }
    
    if(!is.null(con$seed)) set.seed(con$seed)
    cat(paste("level ", levels[k],": ", "number of elements","=",length(Vmin[[k]]),"; ", sep=""))
    mixmdl[[k]] <- try(normalmixEM(Vmin[[k]], k=2, epsilon = con$epsilon, maxit = con$maxit, maxrestarts=con$maxrestarts,
                                   verb = con$verb, fast=con$fast, ECM = con$ECM, arbmean = con$arbmean, arbvar = con$arbvar), silent=TRUE)

    if(inherits(mixmdl[[k]],"try-error")) {classification[[k]]<-NA;  cat("mixmodel does not work \n"); stop(paste("The EM algorithm does not converge for level", levels[k], "try to remove this level"),call. = FALSE)}
    
    if(!inherits(mixmdl[[k]],"try-error")){
      gr.2=which.max(mixmdl[[k]]$mu) 
      gr.1=which.min(mixmdl[[k]]$mu)
      
      classification[[k]]<-names_lev[[k]][which(mixmdl[[k]]$posterior[,gr.2]<mixmdl[[k]]$posterior[,gr.1])]
      thr[k]<-quantile(Vmin[[k]],mean(mixmdl[[k]]$posterior[,gr.2]<mixmdl[[k]]$posterior[,gr.1]))
    }
    
  }

  if(mixplot==TRUE)
  {
    opar<-par(mfrow=c(1,nl))
    on.exit(par(opar))
    for(k in 1:nl){
      if(!inherits(mixmdl[[k]],"try-error")){
        plot(mixmdl[[k]], which = 2,breaks=30,main2=levels[k])
        abline(v=thr[k],lty=2,col="blue")}
    }
  }
  

  tab_up<-unique(data.frame(mget(levels)))
    
  tab_up$det<-0
  
    for(k in nl:1){
        tab_up$det[get(levels[k],pos=tab_up) %in% classification[[k]]]<-levels[k]
    }
  
  data<-data[order(get(levels[nl],pos=data)),]
  return(list(data=data,tab_up=tab_up,vmin=Vmin,yp=y,med_lev=med_lev, names_lev=names_lev,
              sc.x=sc.x,sc.y=sc.y,sc.dx=dX,sc.dy=dY,start=start,thr=thr,levels=levels,mixing=classification,df.med_lev=df.med_lev,mixmdl=mixmdl,obs=obs,vtime=vtime,h.int=h.int))
}

