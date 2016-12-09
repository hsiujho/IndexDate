
match_ex1_sqldf=function(){
  require(sqldf)
  on=Sys.time()
  y1=befmatch
  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]

  ctrl$case_id=""
  ncase=nrow(case)
  ps_threshold=0.1
  ctrl$ps_threshold=ps_threshold

  for(i in 1:ncase){
    ctrl$diff_ps=abs(ctrl$ps-case$ps[i])
    a1=sqldf("select ID from ctrl where diff_ps<=ps_threshold and case_id=''
             order by diff_ps,ID limit 1")
    if(nrow(a1)>0) ctrl[ctrl$ID==a1$ID[1],"case_id"]=case$ID[i]
    if(i%%100==0)cat(sprintf("i=%i\n",i))
  }
  Sys.time()-on
  df0=ctrl[ctrl$case_id!="",c("ID","case_id")]
  rownames(df0)=NULL
  print(Sys.time()-on)
  return(df0)
}

match_ex1_dplyr=function(){
  require(dplyr)
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]

  ctrl$case_id=""
  ncase=nrow(case)
  ps_threshold=0.1

  for(i in 1:ncase){
    a1=select(ctrl,ID,ps,case_id) %>%
      mutate(diff_ps=abs(case$ps[i]-ps)) %>%
      dplyr::filter(diff_ps<=ps_threshold&case_id=="") %>%
      top_n(1,desc(diff_ps)) %>%
      arrange(ID) #%>% select(ID)
    if(nrow(a1)>0) ctrl[ctrl$ID==a1$ID[1],"case_id"]=case$ID[i]
    if(i%%100==0)cat(sprintf("i=%i\n",i))
  }
  df1=ctrl[ctrl$case_id!="",c("ID","case_id")]
  rownames(df1)=NULL
  print(Sys.time()-on)
  print(sapply(df1,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df1)
}

match_ex1_data.table=function(){
  require(data.table)
  require(magrittr)
  require(dtplyr)
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,] %>% mutate(case_id="",diff_ps=NA) %>% data.table()
#   ctrl$case_id=""
# ctrl[,case_id:=""]
  setkey(ctrl,ID)
  ncase=nrow(case)
  ps_threshold=0.1

  for(i in 1:ncase){
    a1=ctrl[,c("ID","ps","case_id","diff_ps"),with=F]
#    a1[,`:=`(diff_ps=abs(case$ps[i]-ps))]
    a1$diff_ps=abs(case$ps[i]-a1$ps)
    a1=a1[diff_ps<=ps_threshold&case_id=="",]
    setkey(a1,diff_ps,ID)
    a1=a1[1,ID]
    if(length(a1)>0) ctrl[ctrl$ID==a1,"case_id"]=case$ID[i]
    if(i%%100==0)cat(sprintf("i=%i\n",i))
  }
  df2=ctrl[ctrl$case_id!="",c("ID","case_id"),with=F] %>% as.data.frame()
  rownames(df2)=NULL
  print(Sys.time()-on)
  print(sapply(df2,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df2)
}

match_ex1_rcpp=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match1(case,ctrl,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp2=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match2(case,ctrl,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp3=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match3(case,ctrl,"ps",.1,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp4=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match4(case,ctrl,"ps",.1,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp5=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match5(case,ctrl,"ps",.1,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp6=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match6(case,ctrl,"ps",.1,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}

match_ex1_rcpp7=function(){
  on=Sys.time()
  y1=befmatch

  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  aa<-Match7(case,ctrl,"ps",.1,NumMatch = 1,threads=8)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
  print(Sys.time()-on)
  print(sapply(df3,function(x){
    xx=befmatch[befmatch$ID%in%x,"ps"]
    return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
  }))
  return(df3)
}
