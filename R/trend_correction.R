###correction of the trend

trend_correction<-function(dd,trait,var)
{
  
  
  dd$cobs<-get(trait,pos=dd)

    dd %>% group_by(data.frame(mget(var,envir=as.environment(dd)))) %>%
      summarise(nb=n(), med=median(get(trait),na.rm=TRUE)) %>%
      unique() %>% data.frame() ->medobs

    
    dd=merge(dd,medobs[,c(var,"med")],by=var)
    dd$cobs<-dd$cobs-dd$med
    
  
  return(dd)
  
}
