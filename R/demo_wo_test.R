











demo_wo_test=function(y,group_var,disc_var,cont_var){
  a0=lapply(disc_var,function(x){
    mutate_each_(y,funs(as.character),x) %>%
      melt(id.vars=group_var,measure.vars = x) %>%
      group_by_(group_var,"variable","value") %>% summarise(n=n()) %>%
      ungroup() %>%
      mutate_each(funs(as.character),variable)
  }) %>% bind_rows() %>%
    left_join(group_by_(y,group_var)%>%summarise(n=n())%>%mutate_each_(funs(as.character),group_var),by=group_var) %>%
    mutate(n=sprintf("%i (%.2f%s)",n.x,n.x/n.y*100,"%")) %>%
    dcast(formula(sprintf("variable+value~%s",group_var)),value.var="n")

  a1=lapply(cont_var,function(x){
    melt(y,id.vars=group_var,measure.vars = x) %>%
      group_by_(group_var,"variable") %>%
      summarise(mean=mean(value),sd=sd(value),median=median(value),q1=quantile(value,prob=0.25),q3=quantile(value,prob=0.75)) %>%
      mutate(meansd=sprintf("%.2f\u00B1%.2f",mean,sd),medianIQR=sprintf("%.2f (%.2f-%.2f)",median,q1,q3)) %>%
      melt(id.vars=c(group_var,"variable"),measure.vars=c("meansd","medianIQR"),variable.name = "value",value.name="n") %>%
      mutate_each(funs(as.character))
  }) %>% bind_rows() %>%
    dcast(formula(sprintf("variable+value~%s",group_var)),value.var="n")
  return(bind_rows(a0,a1))
}
