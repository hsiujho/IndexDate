
match_ex2_rcpp=function(){
  on=Sys.time()
  y$ID=as.vector(y$ID)
  df2=lapply(0:1,function(x){
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match1(case,ctrl,ps_var = "birn",ps_threshold = 365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}

match_ex2_rcpp2=function(){
  on=Sys.time()
  y$ID=as.vector(y$ID)
  df2=lapply(0:1,function(x){
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match2(case,ctrl,ps_var = "birn",ps_threshold = 365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}

match_ex2_rcpp_parallel_1=function(){
  require(parallel)

  on=Sys.time()
  mc_core=2
  cl <- makeCluster(mc_core)
  clusterEvalQ(cl,library(IndexDate))
  df2=parLapply(cl,0:1,function(x){
    y$ID=as.vector(y$ID)
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match1(case,ctrl,ps_var = "birn",ps_threshold = 365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  stopCluster(cl)
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}

match_ex2_rcpp2_parallel_1=function(){
  require(parallel)

  on=Sys.time()
  mc_core=2
  cl <- makeCluster(mc_core)
  clusterEvalQ(cl,library(IndexDate))
  df2=parLapply(cl,0:1,function(x){
    y$ID=as.vector(y$ID)
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match2(case,ctrl,ps_var = "birn",ps_threshold = 365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  stopCluster(cl)
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}

match_ex2_rcpp3_parallel_1=function(){
  require(parallel)

  on=Sys.time()
  mc_core=2
  cl <- makeCluster(mc_core)
  clusterEvalQ(cl,library(IndexDate))
  df2=parLapply(cl,0:1,function(x){
    y$ID=as.vector(y$ID)
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match3(case,ctrl,"birn",365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  stopCluster(cl)
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}

match_ex2_rcpp4_parallel_1=function(){
  require(parallel)

  on=Sys.time()
  mc_core=2
  cl <- makeCluster(mc_core)
  clusterEvalQ(cl,library(IndexDate))
  df2=parLapply(cl,0:1,function(x){
    y$ID=as.vector(y$ID)
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match4(case,ctrl,"birn",365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  stopCluster(cl)
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}


match_ex2_rcppparallel5_1=function(){

  on=Sys.time()
  df2=lapply(0:1,function(x){
    y$ID=as.vector(y$ID)
    case=y[y$gp==1&y$sex==x,]
    ctrl=y[y$gp==2&y$sex==x,]
    aa<-IndexDate::Match5(case,ctrl,"birn",365.0,NumMatch = 4)
    return(aa)
  })
  df3=do.call(rbind,df2)
  rownames(df3)=NULL
  df4=df3[df3$case_id!="",]
  print(Sys.time()-on)
  print(lapply(c("sex","birn"),function(colvar){
    sapply(df4,function(x){
      xx=y[y$ID%in%x,colvar]
      return(c(N=length(xx),mean=mean(xx),sd=sd(xx)))
    })
  }))
  return(df4)
}
