#
#
#
#
#


subgroup_prop=function(df,group_var
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
  df1=group_by_at(df,.vars = group_var)

  df2=merge(summarise(df1,N=n())
            ,summarise_at(df1,.vars = setNames(outcome_name,"E"),funs(sum(.==outcome_item,na.rm = T)))
            ,by=group_var) %>>%
    merge(summarise_at(df1,.vars = setNames(py_name,"PY"),funs(sum(.,na.rm = T))),by=group_var)

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
  df2=rename_at(df2,setNames(group_var,"subgroup"),funs(names)) %>>%
    mutate(subgroup=as.character(subgroup)) %>>%
    select_at(c("Cov","subgroup","N","E","PY","prop_est","prop_l95ci","prop_u95ci","pvalue_test_for_equality"))
  return(df2)
}


