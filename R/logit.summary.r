
# dist_var=c("ID_SEX","rf","hc","cva","ht","dm","acs","cirr"
#            ,outer(c("father","mother","siblings","offspring")
#                   ,c("pso","ra","ss2","sle"),function(x,y) sprintf("%s_%s_bin",x,y)) %>>%
#              t() %>>% as.vector() %>>% "["(-1)
#            ,"father_psatype"
# )
#
# data=mutate(y1,event=ifelse(cohort=="case",1,0))
# event.var="event"
# discrete.cov=dist_var
# continuous.cov="age"

#

logit.summary=function(data,event.var,discrete.cov,continuous.cov){

  f1=discrete.cov %>>%
    (sprintf("factor(%s)",.)) %>>%
    paste(collapse="+") %>>%
    (sprintf("%s~%s+%s",event.var,.,continuous.cov)) %>>%
    as.formula()

  m1=glm(f1,family=binomial(link='logit'),data=data)

  s1=lapply(discrete.cov,function(i){
    h0=group_by_(data,event.var,i) %>>% summarise(N=n()) %>>%
      dcast(as.formula(sprintf("%s~%s",event.var,i)),value.var="N")
    h1=lapply(3:ncol(h0),function(j){
      h0 %>>%
        (mutate(.,Num=sprintf("%i/%i","[["(.,j),"[["(.,2))
                ,Fact1=sprintf("factor(%s)%s",i,colnames(.)[j])
                ,Char=sprintf("%s (%s/%s)",i,colnames(.)[j],colnames(.)[2])))}) %>>%
      bind_rows() %>>%
      dcast(as.formula(sprintf("Fact1+Char~%s",event.var)),value.var="Num")
  }) %>>% bind_rows()

  s2=summary(m1)$coefficients %>>% data.frame() %>>%
    rownames_to_column("Fact1") %>>%
    mutate(Fact=Fact1 %>>%
             (regmatches(., gregexpr("\\(.*?\\)", .))) %>>%
             (gsub("[\\(\\)]", "",.)) %>>%
             (ifelse(.=="character0","",.)),OR=sprintf("%.3f",exp(Estimate))
           ,pvalue=`Pr...z..` %>>% (ifelse(.<.001,"<.001",ifelse(.>.999,">.999",sprintf("%.3f",.))))
    ) %>>% left_join(s1,by="Fact1") %>>%
    mutate(Char=ifelse(is.na(Char),Fact1,Char)) %>>%
    mutate_each_(funs(ifelse(is.na(.),"",.)),levels(data[[event.var]])) %>>%
#    select_(.dots=c("Char",levels(factor(data[[event.var]])),"OR","pvalue")) %>>%
    filter(Char!="(Intercept)")

  return(list(model=m1,summary=s2))
}

# s5=mutate(y1,cohort=factor(cohort,levels=c("ctrl","case"))) %>>%
#   logit.summary(event.var="cohort"
#                 ,discrete.cov=dist_var
#                 ,continuous.cov="age")
#
# s6=mutate(y1,event=ifelse(cohort=="case",1,0)) %>>%
#   logit.summary(event.var="event"
#                 ,discrete.cov=dist_var
#                 ,continuous.cov="age")

