#
#
#
#
#

coxph.summary=function(data
                       ,ftime.var
                       ,fstatus.var
                       ,continuous.cov=NULL
                       ,discrete.cov=NULL){

  if(is.null(continuous.cov) & is.null(discrete.cov))
    stop("continuous.cov=NULL & discrete.cov=NULL")

  if(!is.null(continuous.cov) & is.null(discrete.cov)){
    continuous.cov=continuous.cov[continuous.cov%in%colnames(data)]
    for.cov=paste0(continuous.cov,collapse="+")
  }

  if(!is.null(discrete.cov)){
    a0=sapply(discrete.cov,function(x){
      any(length(unique(data[,x]))==1 #only one categorical
          ,any(table(data[data[,fstatus.var]%in%c(0,1),c(x,fstatus.var)])==0))  #any number of event is zero
    })
    discrete.cov=discrete.cov[!a0]
  }

  if(length(discrete.cov)>0){
    for(i in discrete.cov){
      data[[i]]=factor(data[[i]])
    }
  }

  if(is.null(continuous.cov) & !is.null(discrete.cov)){
    for.cov=paste0(discrete.cov,collapse="+")
  }

  if(!is.null(continuous.cov) & !is.null(discrete.cov)){
    continuous.cov=continuous.cov[continuous.cov%in%colnames(data)]
    for.cov= paste0(c(paste0(discrete.cov,collapse="+"),continuous.cov),collapse="+")
  }

  c1 <- coxph(as.formula(sprintf("Surv(%s,%s)~%s",ftime.var,fstatus.var,for.cov)),data = data)

  tmp=summary(c1)
  pval=2*(1-pnorm(abs(tmp$coef[,'z'])))
  tab1=data.frame(
    cov=rownames(tmp$conf.int)
    ,hr=tmp$conf.int[,1]
    ,hr.CI.lower=tmp$conf.int[,3]
    ,hr.CI.upper=tmp$conf.int[,4]
    ,stringsAsFactors = F
    ,check.names = F)
  tab1$HR=with(tab1,sprintf("%.2f",hr),'95% CI'=sprintf("%.2f-%.2f",hr.CI.lower,hr.CI.upper))
  tab1$'p-value'=with(tab1,ifelse(pval<.001,"<.001",ifelse(pval>.999,">.999",sprintf("%.3f",pval))))
  tab1$is.cont=tab1$cov%in%continuous.cov
  tab1$item=tab1$dis.var=""
  tab1$E.ref=tab1$N.ref=tab1$E=tab1$N=NA
  for(i in 1:NROW(tab1)){
    #i=1
    if(!tab1$is.cont[i]){
      k0=lapply(discrete.cov,function(x){
        regmatches(tab1$cov[i], regexpr(pattern=sprintf("^%s",x), tab1$cov[i]))
      })
      tab1$dis.var[i]=unlist(k0)
      tab1$item[i]=gsub(sprintf("^%s",tab1$dis.var[i]),"",tab1$cov[i])
      tab1$N[i]=sum(data[[tab1$dis.var[i]]]==tab1$item[i])
      tab1$E[i]=sum(data[[tab1$dis.var[i]]]==tab1$item[i]&data[[fstatus.var]]=="1")
      if(length(unique(data[[tab1$dis.var[i]]]))==2){
        tab1$N.ref[i]=sum(data[[tab1$dis.var[i]]]!=tab1$item[i])
        tab1$E.ref[i]=sum(data[[tab1$dis.var[i]]]!=tab1$item[i]&data[[fstatus.var]]=="1")
      }
    }
  }
  rownames(tab1)=NULL

  if(!is.null(discrete.cov)){
    b1=lapply(discrete.cov,function(x){
      #x=discrete.cov[1]
      c0=data.frame(table(data[[x]]))
      colnames(c0)=c("stat","N")
      c1=data.frame(table(data[[x]],data[[fstatus.var]]))
      c2=split(c1,c1$Var2)
      c2=lapply(c2,function(x){
        x[[sprintf("status_%s",as.character(x$Var2)[1])]]=x$Freq
        x$Var2=x$Freq=NULL
        return(x)
      })
      c3=Reduce(function(df1,df2)merge(df1,df2,by="Var1"),c2)
      c4=merge(c0,c3,by.x="stat", by.y="Var1")
      c4$stat=as.character(c4$stat)
      c4$variable=x
      return(c4)
    })
    b1=Reduce(rbind,b1)

    d2=tab1
    d2$variable=d2$dis.var
    d2$stat=d2$item
    d2$N=d2$E=d2$N.ref=d2$E.ref=NULL
    d2=merge(b1,d2,by=c("variable","stat"),all = T)
    d2$is.cont=d2$dis.var=d2$item=NULL
    rownames(d2)=NULL
    return(list(coxph=c1,summary=tab1,summary2=d2))
  } else {
    return(list(coxph=c1,summary=tab1))
  }

}
