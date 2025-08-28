smoothing.curve<-function(dat2i,h.int=10,n=n,method="Vmin")
{
  time=which(!is.na(dat2i))
  na.time=which(is.na(dat2i))
  
  dat2i=dat2i[!is.na(dat2i)]
  modelp <-ksmooth(time, dat2i, bandwidth = h.int,kernel="normal", n.points=floor(h.int*n)) #Nadaraya-Wattson
  
  
  dX=rowMeans(embed(modelp$x,2))
  df=diff(modelp$y)/diff(modelp$x)
  df[floor(dX) %in% na.time]<-NA
  
  if(method=="Vmin") {Vmin=min(df,na.rm=T);start=modelp$x[which.min(df)]}else{Vmin=max(df,na.rm=T);start=modelp$x[which.max(df)]}
  return(list(y=modelp$y,x=modelp$x,Vmin=Vmin,start=start,dY=df,dX=dX ))
}


moment.estimation<-function(dat2i,h.int=10,th,eps=5,n=n)
{
  time=which(!is.na(dat2i))
  na.time=which(is.na(dat2i))
  
  dd=dat2i[!is.na(dat2i)]
  #dat2i=med_batch[ibatch,]
  #plot(time,dd)
  ld=length(dd)
  #X <- data.frame(age=1:ld)
  #time=1:length(dat2i)
  modelp <-ksmooth(time, dd, bandwidth = h.int,kernel="normal", n.points=floor(h.int*n)) #Nadaraya-Wattson smoothing curve
  #lines(modelp$x,modelp$y,col="red")
  
  dX=rowMeans(embed(modelp$x,2))
  df=diff(modelp$y)/diff(modelp$x)
  df[floor(dX) %in% na.time]<-NA
  
  modelp2 <-ksmooth(dX, df, bandwidth = 1,kernel="normal", n.points=floor(h.int*n)) #Nadaraya-Wattson smoothing curve

  ddf=diff(modelp2$y)/diff(modelp2$x)
  ddX=rowMeans(embed(modelp2$x,2))
  #plot(ddX,ddf,type='l')

  dsp1<-data.frame(points=1:length(df),dX=dX,df=df,sign=c(diff(sign(df)),NA))
  dsp2<-data.frame(points=1:length(ddf),ddX=ddX,ddf=ddf,sign=c(diff(sign(ddf)),NA))


  #####start of the disturbance, solution 1###
  ############################################

  inflexion.decr<-as.vector(na.omit(dsp2$ddX[dsp2$sign>0]))

  tx=inflexion.decr[df[which(dsp2$sign>0)]<=th] ##inflexion points under the threshold thb

  #if no value
  if (length(tx)==0){
    tx=dX[which.min(df)]
  }
  M0=unique(tx)
  #plot(time,dd)
  #lines(modelp$x,modelp$y,col="red")
  #abline(v=M0)

  #####start of the disturbance, solution 2###
  ###########################################

  local.max<-as.vector(na.omit(dsp1$dX[dsp1$sign<0]))
  #plot(1:ld,dat2i)
  #lines(modelp$x,modelp$y,col="red")
  #abline(v=local.max)

  #if no-value, the start of the disturbance is the first time
  if (length(local.max)==0){
    local.max=time[1]
  }
  #mm1=unique(local.max)

  #####start of the disturbance, solution 3###
  ###########################################

  inflex<-as.vector(na.omit(dsp2$ddX[abs(dsp2$sign)>0]))

  if (length(local.max)==0){
    local.max=time[1]
  }

  mm1=unique( inflex[pmax(which(inflex %in% tx)-1,1)] )

  #####end of disturbance###
  ##########################

  local.min<-as.vector(na.omit(dsp1$dX[dsp1$sign>0]))
  #abline(v=local.min)

  #if no-value, the end of disturbance is the last time.
  if (length(local.min)==0){
    local.min=time[ld]
  }

  mm2=unique(local.min)

  #est$filtrexb[ibatch,mm2]

  #####end of reaction######
  ##########################

  inflexion.incr<-as.vector(na.omit(dsp2$ddX[dsp2$sign<0]))

  #if no-value, the end of reaction is the last time.
  if (length(inflexion.incr)==0){
    inflexion.incr=length(dd)
  }

  mm3<-unique(inflexion.incr)




  M1<-M2<-M3<-NULL
  for (txx in M0){
    if(length(mm1[mm1<txx])!=0){
      mmm1=max(mm1[mm1<txx])}else{
        mmm1<-time[1]
      }
    M1=c(M1,mmm1)

    if(length(mm2[mm2>txx])!=0){
      mmm2=min(mm2[mm2>txx])}else{
        mmm2<-time[ld]
      }
    M2=c(M2,mmm2)

    if(length(mm3[mm3>mmm2])!=0){
      mmm3=min(mm3[mm3>mmm2])} else{
        mmm3<-time[ld]}
    M3=c(M3,mmm3)
  }
  
   #intensity
  M4=rep(0,length(M0))
  for(k in 1:length(M0)){
  if(M2[k]==time[ld]) M4[k]=abs(tail(modelp$y,1)-modelp$y[which(dX %in% M0[k] | ddX %in% M0[k])])/(M2[k]-M0[k])
  if(M0[k]==time[1]) M4[k]=abs(modelp$y[which(dX %in% M2[k])]-modelp$y[1])/(M2[k]-M0[k]) 
  if(M2[k]!=time[ld] & M0[k]!=time[1])  M4[k]=abs(modelp$y[which(dX %in% M2[k] )]-modelp$y[which(dX %in% M0[k] | ddX %in% M0[k])])/(M2[k]-M0[k])
  } 

  return(list(M0=M0,M1=M1,M2=M2,M3=M3,M4=M4))
}
