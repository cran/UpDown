UpDownApp<-function(updown.out,obs=NULL,width=1000,height=1000){

levels=updown.out$levels
nl=length(levels)

pert_down<-rev(updown.out$Down)
vtime=updown.out$vtime
if(is.null(obs)) obs=updown.out$obs
names_lev<-rev(updown.out$names_lev)
med_lev<-rev(updown.out$med_lev)
sc.x_lev<-rev(updown.out$sc.x_lev)
sc.y_lev<-rev(updown.out$sc.y_lev)

data=updown.out$data




if(nl>5) stop("the maximum number of level is 5")


##############ui#################
#################################

ui<-navbarPage("Shiny disturbance",
            
             tabPanel("Plot",

                      # Application title
                      titlePanel(paste("Plot of detected group and individual disturbances")),
                      
                      # Sidebar with a slider input for number of bins
                      selectInput("phenotype", "observations",
                                  choices = obs),

                      
         
                      conditionalPanel( 
                        condition = "nl >1",
                        lapply(1, function(i) {if(nl>1)
                        selectInput("lev1", levels[1],
                                      choices = unique(get(levels[1],pos=data)))
                        })
                       
                      ),
                      conditionalPanel( 
                        condition = "nl >2",
                        lapply(1:(nl-1), function(i) { 
                          if(nl>(i+1)) uiOutput(paste("dynamic",i,sep=""))
                          })
                        ),
                      
                      
                      conditionalPanel(
                        condition = "nl ==1",
                        lapply(1, function(i) { if(nl==1)
                        selectInput("lev1", levels[nl], choices = unique(get(levels[nl],pos=data)))
                        }),
                        selectInput("correction", "Use corrected observations?",
                                    choices = c("Original observations"="No","Corrected observations"="Yes"))
                        
                        
                      ),


                      mainPanel(
                        plotOutput("distPlot",width=width,height=height)
                      )



             ),
             
             tabPanel("median plot",
  
                      titlePanel("median observations and smoothing curves"),
                      

                        
                        conditionalPanel(
                          condition = 'nl >1',
                          selectInput("ilevel", "level",
                                      choices = levels)),
                       
       
                         conditionalPanel( 
                           condition = "nl >1",
                              uiOutput("med1")

                          
                      ),
    
                      
                      

                      mainPanel(
                       plotOutput("distPlot2",width=width,height=height),
                      )
                      
                      
                      
             )
  )

##############server#############
#################################

server<-function(input, output) {


   output$med1 <- renderUI({


     for(k in levels)
       if(input$ilevel==k)  return(selectInput("element", k, choices = names_lev[[k]]))
     
     
   })

   
  ##for 3 hierarchical levels
  output$dynamic1 <- renderUI({
      tab<-subset(data, get(levels[1],pos=data) %in% input$lev1)
           selectInput("lev2", levels[2],
                       choices = unique(get(levels[2],pos=tab)))
  })

  ##for 4 hierarchical levels
  output$dynamic2 <- renderUI({
    tab<-subset(data, get(levels[2],pos=data) %in% input$lev2)
    selectInput("lev3", levels[3],
                choices = unique(get(levels[3],pos=tab)))
  })
  ##for 5 hierarchical levels
  output$dynamic3 <- renderUI({
    tab<-subset(data, get(levels[3],pos=data) %in% input$lev3)
    selectInput("lev4", levels[4],
                choices = unique(get(levels[4],pos=tab)))
  })
  
   
  
   tab<-reactive({
     tab<-data

     for(i in 1:max(nl-1,1))
     tab<-tab[get(levels[i],pos=tab) %in% input[[paste("lev",i,sep="")]], ] 
     
     tab<-merge(tab,pert_down[[nl]],by=levels[nl], all.x = TRUE)
     })
# 
# 
  output$distPlot <- renderPlot({

  
    tab<-tab()
    if(nrow(tab)==0){return()}else{
    p<-ggplot(tab,aes(x=get(vtime),y=if(input$correction=="No") get(input$phenotype,pos=tab) else tab$cobs))+ylab(paste(input$phenotype))+xlab(paste(vtime))+
      geom_point(size=1) +
      facet_wrap(~get(levels[nl]))
    
    

    if(nl>1)
    if(input$lev1 %in% unlist(pert_down[[1]][levels[1]])){
      dist<-pert_down[[1]][unlist(pert_down[[1]][levels[1]])== input$lev1 ,]
      for(k in 1:nrow(dist))
      p<-p+geom_vline(xintercept=as.numeric(dist[k,"start"]),linetype="dotted",color = "red")+geom_vline(xintercept=as.numeric(dist[k,"end"]),linetype="dotted",color = "blue")
      
    }
    if(nl>2)
      if(input$lev1 %in% unlist(pert_down[[2]][levels[1]]) & input$lev2 %in% unlist(pert_down[[2]][levels[2]])){
        dist<-pert_down[[2]][unlist(pert_down[[2]][levels[1]])== input$lev1 & unlist(pert_down[[2]][levels[2]])== input$lev2 ,]
        for(k in 1:nrow(dist))
          p<-p+geom_vline(xintercept=as.numeric(dist[k,"start"]),linetype="dotdash",color = "red")+geom_vline(xintercept=as.numeric(dist[k,"end"]),linetype="dotdash",color = "blue")
          
      }
    

    p<-p+geom_vline(aes(xintercept = start),linetype="longdash",color = "red")+geom_vline(aes(xintercept = end),linetype="longdash",color = "blue")

    suppressWarnings(print(p))
    }
  })
  

  temp<-reactive({
    temp<-list(temp1=input[["ilevel"]], temp2=input[["element"]])
 temp
  })
  # 
  # 
  
  output$distPlot2 <- renderPlot({
    
   
      temp<-temp()
      k<-which(input$ilevel ==levels)
      if(is.null(temp$temp2)){return()}
      else if(!temp$temp2 %in% names_lev[[k]]){return()}
      else{


      
    tt=input$element
    
    dat2i=med_lev[[k]][names_lev[[k]]==tt,]

    time=which(!is.na(dat2i))
    
        dist<-pert_down[[k]][pert_down[[k]]==tt,]
        
    main=if(nrow(dist)>0){paste(levels[k],"disturbance:", "start=",round(dist$start,2),"end=",round(dist$end,2),"intensity=",round(dist$intensity,2),sep=" ")}else{"no disturbance"}
    plot(time,dat2i[!is.na(dat2i)],xlab=paste(levels[k],tt,sep=" " ),main=main,ylab="",ylim=c(min(med_lev[[k]],na.rm=TRUE),max(med_lev[[k]],na.rm=TRUE)))
    lines(sc.x_lev[[k]][names_lev[[k]]==tt,],sc.y_lev[[k]][names_lev[[k]]==tt,],col="red")
    abline(v=dist$start,lty=2,col = "red")
    abline(v=dist$end,lty=2,col = "blue")

      }
  })

}


suppressWarnings(
  runApp(list(ui = ui, server = server))
)
}
