# require(survival)
# #require(dplyr)
# require(pipeR)
# head(mgus)
#
# df=mgus
# df$age_lv=cut(df$age,breaks = c(-Inf,60,70,Inf))
# df$futime_yr=df$futime/365.25

subgroup_prop_v2=function(df,group_var
                       ,py_name
                       ,outcome_name,outcome_item="1")
{
  if(missing(group_var)){
    df$One_group="All"
    group_var="One_group"
  } else if (all(colnames(df)!=group_var)) {
    cat("\ngroup_var not in colnames(df)\n")
    df$One_group="All"
    group_var="One_group"
  }
  df$E=df[[outcome_name]]
  df$PY=df[[py_name]]
  df2=split(df[,c("E","PY",group_var)],df[[group_var]]) %>>%
    lapply(function(x){
      data.frame(subgroup=x[[group_var]][1]
                 ,N=NROW(x)
                 ,E=sum(x$E==outcome_item,na.rm = T)
                 ,PY=sum(x$PY,na.rm = T))
    }) %>>%
    (Reduce(rbind,.))

  df2$prop_est=df2$prop_l95ci=df2$prop_u95ci=NA
  df2$Cov=group_var

  for(i in 1:nrow(df2)){
    aa1=prop.test(x = df2$E[i],n = df2$PY[i])
    df2$prop_est[i]=aa1$estimate
    df2$prop_l95ci[i]=aa1$conf.int[1]
    df2$prop_u95ci[i]=aa1$conf.int[2]
  }

  if(nrow(df2)>1){
    aa2=prop.test(x = df2$E,n = df2$PY)
    df2$pvalue_test_for_equality=aa2$p.value
  } else {
    df2$pvalue_test_for_equality=NA
  }
  df3=df2[,c("Cov","subgroup","N","E","PY","prop_est","prop_l95ci","prop_u95ci","pvalue_test_for_equality")]
  return(df3)
}

# lapply(c("sex","age_lv"),function(x){
#   subgroup_prop_v2(df = df,group_var = x,py_name = "futime_yr"
#                 ,outcome_name = "death",outcome_item = 1)
# }) %>>%
#   (Reduce(rbind,.))
