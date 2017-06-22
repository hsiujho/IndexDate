#
#
#
#
#

coxph.summary=function(data,ftime.var,fstatus.var,continuous.cov=NULL,discrete.cov=NULL){

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

  if(is.null(continuous.cov) & !is.null(discrete.cov)){
    for.cov=paste0(sprintf("factor(%s)",discrete.cov),collapse="+")
  }

  if(!is.null(continuous.cov) & !is.null(discrete.cov)){
    continuous.cov=continuous.cov[continuous.cov%in%colnames(data)]
    for.cov=paste0(sprintf("factor(%s)",discrete.cov),collapse="+") %>%
      c(continuous.cov) %>%
      paste0(collapse="+")
  }

  covm = sprintf("~%s",for.cov) %>%
    as.formula() %>% model.matrix(data) %>% "["(,-1,drop=F)

  #  if(class(covm)!="matrix") covm=matrix(covm,ncol=1)
  # if(ncol(covm)==2){
  #   covm.colnames=colnames(covm)
  #   covm=matrix(covm[,-1],ncol=1)
  #   colnames(covm)=covm.colnames[2]
  # } else {
  #   covm=covm[,-1]
  # }

  c1 <- sprintf("Surv(%s,%s)~%s",ftime.var,fstatus.var,for.cov) %>%
    as.formula() %>% coxph(data = data)

  tmp=summary(c1)
  pval=2*(1-pnorm(abs(tmp$coef[,'z'])))
  tab1=data.frame(
    cov=rownames(tmp$conf.int)
    ,hr=tmp$conf.int[,1]
    ,hr.CI.lower=tmp$conf.int[,3]
    ,hr.CI.upper=tmp$conf.int[,4]
    ,stringsAsFactors = F
    ,check.names = F) %>%
    mutate(N=colSums(covm)
           ,E=colSums(subset(covm,data[,fstatus.var]==1))
           ,N.ref=c1$n-N
           ,E.ref=c1$nevent-E
           ,HR=sprintf("%.2f",hr),'95% CI'=sprintf("%.2f-%.2f",hr.CI.lower,hr.CI.upper)
           ,'p-value'=ifelse(pval<.001,"<.001",ifelse(pval>.999,">.999",sprintf("%.3f",pval)))
           )

  tab1[tab1$cov%in%continuous.cov,c("N","E","N.ref","E.ref")]=NA

  for(i in discrete.cov){
    if(length(unique(data[,i]))>2) {
      tab1[regexpr(i,tab1$cov)!=-1,c("N.ref","E.ref")]=NA
    }
  }

  if(!is.null(discrete.cov)){
    b1=lapply(discrete.cov,function(x){
      # x=discrete.cov[1]
      c0=group_by_(data,x) %>>% dplyr::summarise(N=n())
      attributes(c0[[x]])$label=NULL
      c1=group_by_(data,x,fstatus.var) %>>% dplyr::summarise(N=n()) %>>%
        dcast(formula(sprintf("%s~%s",x,fstatus.var)),value.var="N")
      left_join(c0,c1,by=x) %>>% rename_(.dots=c(stat=x)) %>>%
        mutate_at(funs(as.character),.vars="stat")
    }) %>>% setNames(discrete.cov) %>>% bind_rows(.id="variable")

    d0=tab1$cov %>>%
      (regmatches(., gregexpr("\\(.*?\\)", .))) %>>%
      (gsub("[\\(\\)]", "",.)) %>>%
      (ifelse(.=="character0","",.))
    d1=tab1$cov %>>% (gsub(".*\\)","",.))
    d2=mutate(tab1,variable=d0,stat=d1) %>>%
      dplyr::select(-N,-E,-N.ref,-E.ref) %>>%
      (full_join(b1,.,by=c("variable","stat")))

    return(list(coxph=c1,summary=tab1,summary2=d2))
  } else {
    return(list(coxph=c1,summary=tab1))
  }

}
