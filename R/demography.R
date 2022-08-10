pv_fn=function(x,digit=3){
  a=10^(-digit)
  tx1=sprintf("<.%s1",paste(rep(0,digit-1),sep="",collapse=""))
  tx2=sprintf(">.%s",paste(rep(9,digit),sep="",collapse=""))
#  tx3=sprintf("tx4=sprintf(\"%%.%if\",x)",digit)
#  eval(parse(text=tx3))
#  tx4=as.character(round(x,digit))
  tx4=sprintf(sprintf("%s.%if","%",digit),x)
  ifelse(x<a,tx1,ifelse(x>1-a,tx2,tx4))
}

NAR.fn = function(x,g,tp)
{
  #number at risk
  NAR=t(sapply(split(x,g),function(y){
    colSums(outer(y,tp,">="))
  }))
  colnames(NAR)=tp
  return(NAR)
}

IAT.fn=function(crr.obj,tp,alpha=0.05)
{
  #incidence at timepoints
  tt=timepoints(crr.obj, tp)
  aa = tt$est
  std = sqrt(tt$var)
  lower=100*(aa+qnorm(alpha/2)*std)
  lower=matrix(apply(cbind(c(lower),0),1,max),dim(lower)[1])
  upper=100*(aa+qnorm(1-alpha/2)*std)
  upper=matrix(apply(cbind(c(upper),100),1,min),dim(upper)[1])
  tab2=data.frame(rbind(matrix(sprintf("%.2f (%.2f-%.2f)",100*aa,lower,upper),nrow=nrow(aa))))
  colnames(tab2)=tp
  return(tab2)
}

demo_fn=function(y,group_var,disc_var=NULL,cont_var=NULL,...){
  var_names=colnames(y)
  if(all(var_names!=group_var)){stop(sprintf("%s is not in the data",group_var))}
  if(all(is.null(c(disc_var,cont_var)))){stop("Nothing to be done!")}
  if(!is.factor(y[,group_var])){y[,group_var]=factor(y[,group_var])}
  disc_var=disc_var[disc_var%in%var_names]
  cont_var=cont_var[cont_var%in%var_names]

  t0=c(table(y[,group_var]))
  ng=length(t0)
  k1=data.frame(y[,disc_var])
  colnames(k1)=disc_var
  disc_stat=do.call(rbind,lapply(k1,function(x){
    t1=table(y[,group_var],x)
    pv=try(chisq.test(x,y[,group_var])$p.value)
    t2=matrix(sprintf("%s (%.1f)",prettyNum(t1,...),t1/t0*100),ncol=ng,byrow=T)
    rownames(t2)=colnames(t1)
    colnames(t2)=rownames(t1)
    data.frame(t2,pv=ifelse("try-error" %in% class(pv),"",pv_fn(pv)),stringsAsFactors=F)
  }))

  k1=data.frame(y[,cont_var])
  colnames(k1)=cont_var
  cont_stat=do.call(rbind,lapply(k1,function(x){
    pv1=anova(lm(x~y[,group_var]))$'Pr(>F)'[1]
    if(ng==2){
      pv2=wilcox.test(x~y[,group_var])$p.value
    } else if (ng>2) {
      pv2=kruskal.test(x~y[,group_var])$p.value
    }
    t2=sapply(split(x,y[,group_var]),function(x){
      c(mean_SD=sprintf("%.1f\u00B1%.1f",mean(x, na.rm =T),sd(x, na.rm =T))
        ,median_q1q3=sprintf("%.1f (%.1f~%.1f)",median(x, na.rm =T),quantile(x,.25, na.rm =T),quantile(x,.75, na.rm =T)))
    })
    data.frame(t2,pv=c(pv_fn(pv1),pv_fn(pv2)),stringsAsFactors=F)
  }))
  return(rbind(cohort=c(prettyNum(table(y[,group_var]),...)," "),disc_stat,cont_stat))
}

#find string including specify substring
mymatch=function(x,txt) x[regexpr(txt,x)!=-1]

demo_fn_v2=function(y,group_var,disc_var,cont_var,big.mark=",",perc.sym=T,dig.disc=1,dig.cont=1){

  if(missing(group_var)) {
    group_var="Overall"
    y$Overall="Overall"
  }

  if(all(group_var!=colnames(y))) {
    group_var="Overall"
    y$Overall="Overall"
  }

  y[[group_var]]%<>%as.factor()

  var_names=colnames(y)
  l0=group_by_at(y,group_var) %>>% dplyr::summarise(gN=n())
  l3=mutate(l0,variable="cohort",stat="",gN=prettyNum(gN,big.mark=big.mark)) %>>%
    select_at(c("variable","stat","gN",group_var)) %>>%
    dcast(sprintf("variable+stat~%s",group_var),value.var = "gN")

  if(missing(disc_var)&missing(cont_var)){return(l3)}

  gn=nrow(l0)

  if(!missing(disc_var)) {
    disc_var=disc_var[disc_var%in%var_names]
    if(length(disc_var)==0) break()
    if(gn>1){
      k0=list()
      for(i in disc_var){
        tmp=try(chisq.test(y[[group_var]],y[[i]]))
        k0[[i]]=c(i,ifelse("try-error" %in% class(tmp),"",pv_fn(tmp$p.value)))
      }
      k1=do.call(rbind,k0) %>>% as.data.frame(stringsAsFactors =F) %>>% setNames(c("variable","pvalue"))
    }

    l1=mutate_at(y,.vars=disc_var,list(~as.character(.))) %>>%
      melt(id.vars = group_var,measure.vars = disc_var,value.name = "stat") %>>%
      group_by_at(c(group_var,"variable","stat")) %>>%
      dplyr::summarise(N=n(),.groups = "drop") %>>%
      left_join(l0,by=group_var) %>>%
      dplyr::mutate(lab=sprintf("%s (%s%s)",prettyNum(N,big.mark=big.mark),sprintf(sprintf("%%.%if",dig.disc),N/gN*100),ifelse(perc.sym,"%",""))) %>>%
      mutate_at(.vars="stat",list(~as.character(.))) %>>%
      dcast(sprintf("variable+stat~%s",group_var),value.var="lab",stringsAsFactors=F) %>>%
      mutate_at(.vars="variable",list(~as.character(.))) %>>%
      {
        if(gn>1){
          left_join(.,k1,by="variable")
        } else {
          .
        }
      }
    l3=bind_rows(l3,l1)
  }

  if(!missing(cont_var)) {
    cont_var=cont_var[cont_var%in%var_names]
    if(length(cont_var)==0) break()
    if(gn>1){
      k0=list()
      for(i in cont_var){
        pv1=try(anova(lm(y[[i]]~y[[group_var]])))
        if(gn==2){
          pv2=try(wilcox.test(y[[i]]~y[[group_var]]))
        } else if (gn>2) {
          pv2=try(kruskal.test(y[[i]]~y[[group_var]]))
        }
        pv11=ifelse("try-error" %in% class(pv1),"",pv_fn(pv1$'Pr(>F)'[1]))
        pv22=ifelse("try-error" %in% class(pv2),"",pv_fn(pv2$p.value))
        k0[[i]]=cbind(i,c("Mean\u00B1SD","Median (Q1~Q3)"),c(pv11,pv22))
      }
      k1=do.call(rbind,k0) %>>% as.data.frame(stringsAsFactors =F) %>>% setNames(c("variable","stat","pvalue"))
    }

    prefn=function(x,a=big.mark){
      a0=sprintf("%%.%if",dig.cont) %>>%
        sprintf(abs(x)) %>>%
        strsplit("\\.") %>>% "[["(1)
      as.numeric(a0[1]) %>>%
        prettyNum(big.mark=a) %>>%
        (paste0(c(.,a0[2]),collapse=".")) %>>%
        (ifelse(x<0,sprintf("-%s",.),.))
    }

    l2=mutate_at(y,.vars=cont_var,list(~as.numeric(.))) %>>%
      melt(id.vars = group_var,measure.vars = cont_var) %>>%
      group_by_at(c(group_var,"variable")) %>>%
      dplyr::summarise(Mean=prefn(mean(value,na.rm=T))
                       ,SD=prefn(sd(value,na.rm=T))
                       ,median=prefn(median(value,na.rm=T))
                       ,Q1=prefn(quantile(value,0.25,na.rm=T))
                       ,Q3=prefn(quantile(value,.75,na.rm=T))
                       ,.groups = "drop") %>>%
      dplyr::mutate("Mean\u00B1SD"=sprintf("%s\u00B1%s",Mean,SD)
                    ,`Median (Q1~Q3)`=sprintf("%s (%s~%s)",median,Q1,Q3)) %>>%
      melt(id.vars=c(group_var,"variable"),measure.vars=c("Mean\u00B1SD","Median (Q1~Q3)"),variable.name = "stat") %>>%
      mutate_at(.vars="stat",list(~as.character(.))) %>>%
      dcast(sprintf("variable+stat~%s",group_var),value.var="value") %>>%
      mutate_at(.vars="variable",list(~as.character(.))) %>>%
      {
        if(gn>1) {
          left_join(.,k1,by=c("variable","stat"))
        } else {
          .
        }
      }
    l3=bind_rows(l3,l2)
  }
  return(l3)
}
