
index_date_ex_sqldf=function(){
  require(sqldf)
  on=Sys.time()
  y=MDD
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  z1=sqldf("select b.*,a.dg2_index_date as index_date, a.id as case_id
           ,abs(a.date_mdd-b.date_mdd)+abs(a.id_birn-b.id_birn) as diff
           from ctrl as b, case0 as a
           where a.id_sex=b.id_sex and a.dg2_index_date+56<b.ev_date
           group by b.id having diff=min(diff)
           order by b.id,diff")
  print(Sys.time()-on)
  return(z1)
}

index_date_ex_dplyr=function(){
  require(dplyr)
  on=Sys.time()
  y=MDD
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  z2=left_join(ctrl,select(case0,ID,ID_SEX,dg2_index_date,date_MDD,ID_BIRn),by="ID_SEX") %>%
    dplyr::filter(dg2_index_date.y+56<ev_date) %>%
    mutate(diff=abs(date_MDD.x-date_MDD.y)+abs(ID_BIRn.x-ID_BIRn.y),rannum=runif(length(ID_BIRn.x))) %>%
    group_by(ID.x) %>%
    dplyr::filter(diff==min(diff)) %>%
    dplyr::filter(rannum==min(rannum)) %>%
    arrange(ID.x)
  print(Sys.time()-on)
  return(z2)
}

index_date_ex_data.table=function(){
  require(data.table)
  on=Sys.time()
  y=MDD
  ctrl=y[y$cohort!=1,] %>% data.table(key="ID_SEX")
  case0=y[y$cohort==1,c("ID","ID_SEX","dg2_index_date","date_MDD","ID_BIRn")] %>% data.table(key="ID_SEX")

  z3=merge(ctrl,case0, allow.cartesian=TRUE)
  # allow.cartesian=TRUE for column of duplicate values
  z3=z3[dg2_index_date.y+56<ev_date,]
  z3[,c("diff","rannum"):=list(abs(date_MDD.x-date_MDD.y)+abs(ID_BIRn.x-ID_BIRn.y),runif(nrow(z3)))]
  z3.1=z3[,.(diff,rannum,ID.y,dg2_index_date.y,diff_min=min(diff)),by=ID.x]
  z3.2=z3.1[diff==diff_min,]
  z3.3=z3.2[,.(diff,rannum,ID.y,dg2_index_date.y,rannum_min=min(rannum)),by=ID.x]
  z3.4=z3.3[rannum==rannum_min,]
  setkey(z3.4,ID.x)
  print(Sys.time()-on)
  return(z3.4)
}

index_date_ex_rcpp=function(){
  on=Sys.time()
  y=MDD
  y$sexn=ifelse(y$ID_SEX=="M",1,0)
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  aa=IndexDate::index_date(case0,ctrl)
  a1=aa[aa$case_id!="",]
  a1$index_date=as.Date(a1$index_date,"1970-1-1")
  z5=a1[order(a1$ID),]
  print(Sys.time()-on)
  return(z5)
}

index_date_ex_rcpp_1=function(){
  on=Sys.time()
  y=MDD
  y$sexn=ifelse(y$ID_SEX=="M",1,0)
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  aa=IndexDate::index_date_1(case0,ctrl)
  a1=aa[aa$case_id!="",]
  a1$index_date=as.Date(a1$index_date,"1970-1-1")
  z5=a1[order(a1$ID),]
  print(Sys.time()-on)
  return(z5)
}

index_date_ex_rcpp_parallel=function(){
  require(parallel)
  on=Sys.time()
  y=MDD
  y$sexn=ifelse(y$ID_SEX=="M",1,0)
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  mc_core=2
  g0=nrow(ctrl)%/%mc_core
  g1=rep(1:mc_core,c(rep(g0,mc_core-1),nrow(ctrl)-g0*(mc_core-1)))
  ctrl_list=split(ctrl,g1)
  cl <- makeCluster(mc_core)
  clusterExport(cl,c("case0"),envir=environment())
  clusterEvalQ(cl,library(IndexDate))
  z6=parLapply(cl,ctrl_list,function(x){
    IndexDate::index_date(df_case=case0,df_ctrl=x)
  })
  stopCluster(cl)
  aa=do.call(rbind,z6)
  a1=aa[aa$case_id!="",]
  a1$index_date=as.Date(a1$index_date,"1970-1-1")
  z8=a1[order(a1$ID),]
  print(Sys.time()-on)
  return(z8)
}
