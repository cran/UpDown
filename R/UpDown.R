UpDown<-function(data,levels,obs,vtime, h.int=NULL,mixplot=FALSE, correction=NULL,kappa=NULL, thr_va=0.5,
                 options=list())
{

  if(!is.vector(levels)) stop("levels must be a vector of character strings",call. = FALSE)
  
  if(sum(is.na(get(vtime,pos=data)))>0) warning(paste(sum(is.na(get(vtime, pos=data))),"rows were removed of the original dataset due to missing times"))
  
    up=Up(data,levels=levels,vtime=vtime,h.int=h.int,mixplot = mixplot,correction=correction,obs=obs,options=options)
  
  down<-Down(up,h.int=up$h.int, kappa=kappa, thr_va=thr_va)

  return(list(data=up$data,levels=up$levels,trait=up$trait,mixmdl=rev(up$mixmdl),Up=up$tab_up,Down=rev(down$pert_down),obs=up$obs,vtime=up$vtime,sc.x_lev=rev(up$sc.x),sc.y_lev=rev(up$sc.y),sc.dx_lev=rev(up$sc.dx),sc.dy_lev=rev(up$sc.dy),med_lev=rev(up$med_lev), names_lev=rev(up$names_lev)))
}