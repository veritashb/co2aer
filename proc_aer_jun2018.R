proc.co2v2<-function(data,pid="PID",co2="CO2",time="time",amb=400,output='summary_decay',strikes=0,sharp.peak="no",min.delta=100,back.to.background="yes"){
  #function(data,pid="id",pm="PM2.5",time="time",temp="temp",temp.out="temp.out",log.rate=5,type='summary',period=25,rolling=5,sensitivity=0.5,strikes=0)
  #load libraries
  require(gridExtra)
  require(grid)
  require(reshape2)
  require(ggplot2)
  require(Hmisc)
  require(plyr)
  require(useful)
  require(zoo)
  require(gdata)
  require(lubridate)
  require(scales)
  require(plotly)
  require(viridis)
  # is.POSIXct <- function(x) inherits(x, "POSIXct")
  
  #renaming vectors of interest
  # data$rolling=rolling
  # data$log.rate=log.rate
  # data$sensitivity=sensitivity
  # data$period=period
  #
  names(data)[names(data)==pid]="pid"
  names(data)[names(data)==co2]="value"
  # names(data)[names(data)==co2]="co2"
  names(data)[names(data)==time]="time"
  #names(data)[names(data)==amb]="amb"
  ##names(data)[names(data)==temp]="temp.in"
  ##names(data)[names(data)==temp.out]="temp.out"
  
# return(data)
# }
  #check basic initializing conditions are met
  stopifnot(is.data.frame(data),!is.na(data$value))
  
  
  #stopifnot(is.POSIXct(data$time))
  #creating new vector with vectors of interest exclusively
  #add into df vectors of ambient co2 concentration, logging rate and sensitivity
  
  # df<-data[,c("variable","time","value","temp.in","temp.out","log.rate","rolling","sensitivity","period")]
  #df<-data[,c("time","value")]
  data$value<-as.numeric(data$value)
  #data$amb<-as.numeric(data$amb)
  
  data<-data[order(data$pid,data$time),]
  # data$variable<-factor(data$variable)
  data$grp<-NA
  data$grp.buildup<-NA
  data$stks<-NA
  
  
  
#############################################  
 #custom funcs 
  pkfinder<-function(data,x,sharp.peak){
    require(Hmisc)
    if(sharp.peak=="yes"){
      flag.pk<-with(data,ifelse(x>Lag(x,-1)&x>Lag(x,1)&x>Lag(x,-2)&x>Lag(x,2),1,0))
      #flag.pk<-with(data,ifelse(x>Lag(x,-1)&x>Lag(x,1)&x>Lag(x,-2)&x>Lag(x,2)&x>3+Lag(x,-5)&x>10+Lag(x,-10)&x>Lag(x,3)&x>Lag(x,4),1,0))
    }
    else{
      flag.pk<-with(data,ifelse(x>Lag(x,-1)&x>Lag(x,1)&Lag(x,-1)>Lag(x,-2)&Lag(x,-2)>Lag(x,-3)&Lag(x,-4)>Lag(x,-5)&Lag(x,1)>Lag(x,2)&Lag(x,2)>Lag(x,3),1,0))
    }
    return(which(flag.pk==1))
  }
  
  lm.eqn = function(df){
    m = lm(ln.value~hourly, df);
    aer.result=round(abs(as.numeric(m$coefficients[2]))*3600,2);
    time.zero=as.POSIXct(df$time[1],origin="1970-01-01")
    time.length2=(as.numeric(df$time[nrow(df)])-as.numeric(df$time[1]))/60
    co2.start=as.numeric(format(df$value[1], digits=4))
    co2.end=as.numeric(format(df$value[nrow(df)], digits=4))
    delta.co2=as.numeric(format(df$value[1]-df$value[nrow(df)],digits=4))
    #temp.zero=as.numeric(format(df$temp.in[1], digits=4))
    #last.temp=as.numeric(format(df$temp.in[nrow(df)],digits=4))
    #temp.gradient=as.numeric(format(last.temp-temp.zero),digits=4)
    b.coef=as.numeric(round(m$coefficients[1],2))
    r2=as.numeric(format(summary(m)$r.squared, digits = 3))
    
    return(cbind(time.zero,co2.start,co2.end,delta.co2,time.length2,aer.result,r2)) 
    
  }
###############################################  
  data$svalue=rollapply(data$value,align="center",width=5,partial=TRUE,FUN=mean)
  peaks<-pkfinder(data,data$svalue,sharp.peak)
  
  #forward scan for decays  
  for(i in seq_along(peaks)){
    j=peaks[i]
    s=0
    # while(strikes<=s&!is.na(j<peaks[i+1]) ){
    while(j<nrow(data)&s<=strikes&j<ifelse(j<max(peaks),peaks[i+1],j+1)){
      
      if((0.5+data$svalue[j]>data$svalue[j+1])){
        data$grp[j]=i
        j=j+1
        data$stks[j]=s
      }
      else{
        if(s<=strikes&10+data$svalue[j]>data$svalue[j+1]){
          s=s+1
          data$grp[j]=i
          j=j+1
          data$stks[j]=s
          
        }
        else{
          break
        }
      }
    }
  }
  
  #backward scan for buildups  
  # for(i in seq_along(rev(peaks))){
  #   j=rev(peaks)[i]
  #   s=0
  #   while(j>1&s<=strikes&j>ifelse(j>min(peaks),rev(peaks)[i+1],j-1)){
  #     if(data$svalue[j]>5+data$svalue[j-1]){
  #       data$grp.buildup[j]=i
  #       j=j-1
  #     }
  #     else{
  #       if(s<=strikes&20+data$svalue[j]>data$svalue[j-1]){
  #         s=s+1
  #         data$grp.buildup[j]=i
  #         j=j-1
  #       }
  #       else{
  #         break
  #       }
  #     }
  #   }
  # }
  # 
  # if(output=="all_data"){
  #   return(data)
  # }
  # 
  if(output=="data"){
      data$ln.value<-log(data$value-amb)
      #data$ln.amb<-log(amb)
      data$hourly<-as.numeric(data$time)
      return(data)
  }
  
  if(output=="summary_decay"){
    data$ln.value<-log(data$value-amb)
    #data$ln.amb<-log(amb)
    data$hourly<-as.numeric(data$time)
    aer<-split(data,data$grp)
    reg.res<-lapply(aer,ddply,~grp,lm.eqn)
    reg.res<-ldply(reg.res)
    if(back.to.background==yes){
      return(subset(reg.res,!is.na(grp)&delta.co2>=min.delta&co2.end<amb+70))
    }
    else{
    return(subset(reg.res,!is.na(grp)&delta.co2>=min.delta))
    }
  }
  
  # if(output=="summary_buildup"){
  #   summary<-ddply(data,.(grp.buildup),summarise,
  #                  start.time=min(time,na.rm=TRUE),
  #                  end.time=max(time,na.rm=TRUE),
  #                  start.value=value[which(time==min(time,na.rm=TRUE))],
  #                  end.value=value[which(time==max(time,na.rm=TRUE))],
  #                  length.min=(as.numeric(max(time,na.rm=TRUE))-as.numeric(min(time,na.rm=TRUE)))/60,
  #                  delta=abs(round(value[which(time==min(time,na.rm=TRUE))]-value[which(time==max(time,na.rm=TRUE))],3)),
  #                  rate=abs(round((value[which(time==min(time,na.rm=TRUE))]-value[which(time==max(time,na.rm=TRUE))])/((as.numeric(max(time,na.rm=TRUE))-as.numeric(min(time,na.rm=TRUE)))/60),3))
  #   )
  #   return(subset(summary,!is.na(grp.buildup)&delta>=min.delta))
  # }
  
  if(output=="graph_decay_dark"){
    data$ln.value<-log(data$value-amb)
    #data$ln.amb<-log(amb)
      data$hourly<-as.numeric(data$time)
    aer.result<-ddply(data,.(grp),function(x) lm.eqn(x))
    if(back.to.background=="yes"){
      aer.result<-subset(aer.result, delta.co2>min.delta&co2.end<amb+70)
    }
    else{
    aer.result<-subset(aer.result, delta.co2>min.delta)
    }
    data<-merge(data,aer.result[,c("aer.result","grp")],by=c("grp"),all.x=TRUE)
    
    data<-data[order(data$pid,data$time),]
    data$aer.text<-data$aer.result
    aer<-split(data,data$grp)
    for(i in seq_along(aer)){if(!is.na(sum(aer[[i]]$aer.result))&nrow(aer[[i]])>1){aer[[i]]$aer.text[2:nrow(aer[[i]])]<-NA}}
    aer<-ldply(aer)
    aer<-aer[order(aer$pid,aer$time),]
    
      return(ggplot()+geom_point(data=data,aes(x=time,y=value),colour='lightgrey',alpha=0.5)+
                geom_line(data=subset(aer,!is.na(aer.result)),aes(x=time,y=value,group=grp,color=aer.result),size=1.1)+
                geom_text(data=aer,aes(x=time+15000,y=value,label=aer.text),size=3,color="white")+theme(legend.position="none")+
                theme_minimal(base_family="Gill Sans")+
               scale_color_viridis(name="ACH [1/h]")+
                theme(axis.text.x = element_text(angle = 45,hjust=1, color="white"),axis.text.y = element_text(color="white")) +
                scale_x_datetime(breaks=date_breaks('1 day'), minor_breaks=date_breaks('1 day'), labels=function(x) format(x,"%m/%d/%y"))+
             facet_wrap(~pid,scales="free_x")+labs(y="CO2 [ppm]",x="Date")+
               theme(
                 panel.background = element_rect(fill = "black",
                                                 colour = "black",
                                                 size = 0.5, linetype = "solid"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(size = 0.5, linetype = "solid",
                                          colour = "black"),
                 plot.background = element_rect(fill = "black"),
                 legend.title = element_text(colour="white", size=12),
                 legend.text = element_text(colour="white", size = 10))
             )
    
    }
    
  if(output=="graph_decay"){
    data$ln.value<-log(data$value-amb)
    #data$ln.amb<-log(amb)
    data$hourly<-as.numeric(data$time)
    aer.result<-ddply(data,.(grp),function(x) lm.eqn(x))
    if(back.to.background=="yes"){
      aer.result<-subset(aer.result, delta.co2>min.delta&co2.end<amb+70)
    }
    else{
      aer.result<-subset(aer.result, delta.co2>min.delta)
    }
    data<-merge(data,aer.result[,c("aer.result","grp")],by=c("grp"),all.x=TRUE)
    
    data<-data[order(data$pid,data$time),]
    data$aer.text<-data$aer.result
    aer<-split(data,data$grp)
    for(i in seq_along(aer)){if(!is.na(sum(aer[[i]]$aer.result))&nrow(aer[[i]])>1){aer[[i]]$aer.text[2:nrow(aer[[i]])]<-NA}}
    aer<-ldply(aer)
    aer<-aer[order(aer$pid,aer$time),]
    
    return(ggplot()+geom_point(data=data,aes(x=time,y=value),colour='grey',alpha=0.5)+
             geom_line(data=subset(aer,!is.na(aer.result)),aes(x=time,y=value,group=grp,color=aer.result),size=1.1)+
             geom_text(data=aer,aes(x=time+15000,y=value,label=aer.text),size=2)+theme(legend.position="none")+
             theme_minimal(base_family="Gill Sans")+
             scale_color_viridis(name="ACH [1/h]")+
             theme(axis.text.x = element_text(angle = 45,hjust=1, color="black"),axis.text.y = element_text(color="black")) +
             scale_x_datetime(breaks=date_breaks('1 day'), minor_breaks=date_breaks('1 day'), labels=function(x) format(x,"%m/%d/%y"))+
             facet_wrap(~pid,scales="free_x")+labs(y="CO2 [ppm]",x="Date")
             
    )
    
  }
  
  # if(output=="graph_buildup"){
  #   data<-ddply(data,.(grp.buildup),transform,delta=abs(round(value[which(time==min(time,na.rm=TRUE))]-value[which(time==max(time,na.rm=TRUE))],3)))
  #   data$grp.buildup<-with(data,ifelse(delta>min.delta,grp.buildup,NA))
  #   
  #   if(length(unique(data$pid))==1){
  #     return(  ggplotly(ggplot()+geom_line(data=subset(data,!is.na(grp.buildup)),aes(time,value,group=grp.buildup),color="blue")+
  #                         geom_point(data=data,aes(time,value),alpha=0.2)+facet_wrap(~pid))
  #     )
  #   }
  #   else
  #     return(
  #       ggplot()+geom_line(data=subset(data,!is.na(grp.buildup)),aes(time,value,group=grp.buildup),color="blue")+
  #         geom_point(data=data,aes(time,value),alpha=0.2)+facet_wrap(~pid, scales="free") 
  #     )
  # }
  
}
