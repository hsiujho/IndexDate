#
#
#
#
#
#
#


tofactor=function(x,lev,rev=T) {
  if(rev==T) {
    x=factor(x,levels = rev(lev))
  } else {
    x=factor(x,levels = lev)
  }
  return(x)
}

plot_CI=function(
  df
  ,point_est_col
  ,lower_CI_col
  ,upper_CI_col
  ,x_breaks
  ,x_limits
  ,y_var
  ,y_levels
  ,strata_vars
  ,strata_levels
  ,arrow_length
  ,y_expand_add
){
  if(missing(x_breaks)) {
    x_breaks=pretty(unlist(c(df[,c(point_est_col,lower_CI_col,upper_CI_col)])))
  }
  if(missing(x_limits)){
    x_limits=range(x_breaks)
  }
  if(length(strata_vars)>0&missing(strata_levels)){
    strata_levels=list()
    strata_levels[[strata_vars[1]]]=unique(df[[strata_vars[1]]])
  }
  if(missing(y_levels)){
    if(class(df[[y_var]])=="factor"){
      y_levels=levels(df[[y_var]])
    } else {
      y_levels=df[[y_var]] %>>% (.[!duplicated(.)])
    }
  }
  x0=df %>>%
    mutate_at(y_var,list(~tofactor(.,y_levels,rev=T))) %>>%
    mutate_at(strata_vars[1],list(~tofactor(.,strata_levels[[strata_vars[1]]],rev=F)))
  x0$LCI=x0[[lower_CI_col]] %>>%
    (ifelse(.<x_limits[1],x_limits[1],.))
  x0$UCI=x0[[upper_CI_col]] %>>%
    (ifelse(.>x_limits[2],x_limits[2],.))
  x0$left_angle=ifelse(x0[[lower_CI_col]]<x_limits[1],45,90)
  x0$right_angle=ifelse(x0[[upper_CI_col]]>x_limits[2],45,90)
  x0$Psudo_strip=""
  x0$seq=1:NROW(x0)

  f2=paste0(strata_vars,collapse = "+") %>>%
    (sprintf("%s~Psudo_strip",.)) %>>%
    as.formula()
  p2=ggplot()+
    geom_segment(aes_string(x="LCI",xend=point_est_col,y="seq",yend="seq"),data = filter(x0,left_angle==90)
                 , lineend = "butt", linejoin ="round"
                 , inherit.aes = F
                 ,arrow = arrow(angle = 90,ends = "first",length = arrow_length))+
    geom_segment(aes_string(x=point_est_col,xend="UCI",y="seq",yend="seq"),data = filter(x0,right_angle==90)
                 , lineend = "butt", linejoin ="round"
                 , inherit.aes = F
                 ,arrow = arrow(angle = 90,ends = "last",length = arrow_length))+
    geom_segment(aes_string(x="LCI",xend=point_est_col,y="seq",yend="seq"),data = filter(x0,left_angle==45)
                 , lineend = "butt", linejoin ="round"
                 , inherit.aes = F
                 ,arrow = arrow(angle = 45,ends = "first",length = arrow_length))+
    geom_segment(aes_string(x=point_est_col,xend="UCI",y="seq",yend="seq"),data = filter(x0,right_angle==45)
                 , lineend = "butt", linejoin ="round"
                 , inherit.aes = F
                 ,arrow = arrow(angle = 45,ends = "last",length = arrow_length))+
    geom_point(aes_string(x=point_est_col,y="seq"),data=x0)+
    facet_grid(f2,scales = "free_y")+
    scale_x_continuous(breaks = x_breaks,trans = "log"
                       ,limits = x_limits,expand = c(0,0.05)#,position = "top"
                       ,sec.axis = dup_axis())+
    scale_y_reverse(breaks=x0$seq,labels=x0[[y_var]]
                    ,minor_breaks=NULL
                    ,expand = expand_scale(add = y_expand_add))

  return(p2)
}


plot_table=function(
  df
  ,y_var
  ,y_levels
  ,strata_vars
  ,strata_levels
  ,cols
  ,cols_labels
  ,cols_group
  ,cols_breaks
  ,y_expand_add
){
  if(length(strata_vars)>0&missing(strata_levels)){
    strata_levels=list()
    strata_levels[[strata_vars[1]]]=unique(df[[strata_vars[1]]])
  }
  if(missing(y_levels)){
    if(class(df[[y_var]])=="factor"){
      y_levels=levels(df[[y_var]])
    } else {
      y_levels=df[[y_var]] %>>% (.[!duplicated(.)])
    }
  }
  if(missing(cols_group)){
    cols_group=cols
  }
  if(missing(cols_labels)){
    cols_labels=cols
  }
  if(missing(cols_breaks)){
    cols_breaks=1:length(cols)
  }
  x0=df %>>%
    mutate_at(y_var,list(~tofactor(.,y_levels,rev=T))) %>>%
    mutate_at(strata_vars[1],list(~tofactor(.,strata_levels[[strata_vars[1]]],rev=F)))
  x0$seq=1:NROW(x0)
  x1=melt(x0,id.vars = c(strata_vars,"seq"),measure.vars = cols
          ,value.name = "text",variable.name = "Size") %>>%
    dplyr::mutate(top.strip=plyr::mapvalues(Size,from = cols,to=cols_group)
                  ,Size=factor(Size,levels = cols)
                  ,x=plyr::mapvalues(Size,from=cols,to=cols_breaks) %>>%
                    as.character() %>>% as.numeric()
    )
  f1=paste0(strata_vars,collapse = "+") %>>%
    (sprintf("%s~top.strip",.)) %>>%
    as.formula()
  p1=ggplot(x1,aes_string(y="seq",x="x",label="text"))+
    geom_text(size=4)+
    facet_grid(f1,scales = "free", switch="y")+
    scale_x_continuous(breaks =cols_breaks
                       ,labels = cols_labels
                       ,sec.axis = dup_axis()
                       ,expand = c(0,0.5)
    )+
    scale_y_reverse(breaks=x0$seq,labels=x0[[y_var]]
                    ,minor_breaks=NULL
                    ,expand = expand_scale(add = y_expand_add))
  return(p1)
}


plot_forest=function(
  df
  ,point_est_col
  ,lower_CI_col
  ,upper_CI_col
  ,x_breaks
  ,x_limits
  ,y_var
  ,y_levels
  ,y_size
  ,strata_vars
  ,strata_levels
  ,strata_size
  ,cols
  ,cols_labels
  ,cols_group
  ,cols_group_widths
  ,cols_breaks
  ,cols_limits
  ,panel_widths
  ,arrow_length
  ,y_expand_add
){
  if(missing(x_breaks)) {
    x_breaks=pretty(unlist(c(df[,c(point_est_col,lower_CI_col,upper_CI_col)])))
  }
  if(missing(x_limits)){
    x_limits=range(x_breaks)
  }
  if(length(strata_vars)>0&missing(strata_levels)){
    strata_levels=list()
    strata_levels[[strata_vars[1]]]=unique(df[[strata_vars[1]]])
  }
  if(missing(y_levels)){
    if(class(df[[y_var]])=="factor"){
      y_levels=levels(df[[y_var]])
    } else {
      y_levels=df[[y_var]] %>>% (.[!duplicated(.)])
    }
  }
  if(missing(cols_group)){
    cols_group=cols
  }
  if(missing(cols_group_widths)){
    cols_group_widths=rep(1,length(cols_group))
  }
  if(missing(cols_labels)){
    cols_labels=cols
  }
  if(missing(cols_breaks)){
    cols_breaks=1:length(cols)
  }
  if(missing(cols_limits)){
    cols_limits=c(head(cols_breaks,1),tail(cols_breaks,1))+c(-0.5,1)
  }
  if(missing(panel_widths)){
    panel_widths=c(2,1)
  }
  if(missing(y_size)){
    y_size=15
  }
  if(missing(strata_size)){
    strata_size=20
  }
  if(missing(arrow_length)){
    arrow_length=unit(0.2,units = "cm")
  }
  if(missing(y_expand_add)){
    y_expand_add=0.95
  }
  cip1=plot_CI(df=df
               ,point_est_col=point_est_col
               ,lower_CI_col=lower_CI_col
               ,upper_CI_col=upper_CI_col
               ,y_var = y_var
               ,y_levels = y_levels
               ,x_breaks=x_breaks
               ,strata_vars=strata_vars
               ,arrow_length=arrow_length
               ,y_expand_add=y_expand_add
  )
  cip2=cip1+
    geom_vline(xintercept = 1,linetype = "dashed")+
    theme(axis.title = element_blank()
          ,axis.ticks.x.top = element_blank()
          ,axis.text.x.top = element_blank()
          ,axis.text.y = element_blank()
          ,axis.ticks.y = element_blank()
          ,strip.placement = "outside"
          ,strip.background = element_blank()
          ,strip.text = element_blank()
          ,panel.spacing.y = unit(0.2,"lines")
          # ,panel.background = element_blank()
    )
  # if(length(strata_vars)==2){
  #   cip2=cip2+scale_y_discrete(limits=rev(y_levels))
  # }
  tabp1=plot_table(df=df
                   ,y_var = y_var
                   ,y_levels = y_levels
                   ,strata_vars=strata_vars
                   ,cols=cols
                   ,cols_labels=cols_labels
                   ,cols_group = cols_group
                   ,y_expand_add=y_expand_add
  )
  tabp2=tabp1+
    theme(axis.title = element_blank()
          ,axis.ticks.x.bottom = element_blank()
          ,axis.text.x.bottom = element_blank()
          ,axis.text=element_text(size=y_size)
          ,strip.background = element_blank()
          ,strip.placement = "outside"
          ,strip.text.x = element_text(size=strata_size)
          ,strip.text.y = element_text(angle = 180,size=strata_size)
          ,panel.background = element_blank()
          ,axis.ticks = element_blank()
          ,panel.spacing = unit(0,"lines")
    )
  if(all(cols==cols_group)){
    tabp2=tabp2+theme(strip.text.x = element_blank())
  }
  # Get maximum widths and heights
  gt1 <- ggplot_gtable(ggplot_build(tabp2))
  panel_posi=which(as.character(gt1$widths)=="1null")
  gt1$widths[panel_posi]=unit(cols_group_widths,units = "null")
  gt2 <- ggplot_gtable(ggplot_build(cip2))
  maxWidth <- unit.pmax(gt1$heights, gt2$heights)
  if(length(strata_vars)==1){
    panel_posi=which(as.character(maxWidth)=="max(1null, 1null)")
    maxWidth[panel_posi]=split(df[[y_var]],df[[strata_vars[1]]]) %>>%
      sapply(length) %>>%
      (data.frame(n=.)) %>>%
      rownames_to_column(strata_vars[1]) %>>%
      mutate_at(strata_vars[1],list(~tofactor(.,strata_levels[[strata_vars[1]]],rev=F))) %>>%
      arrange_at(strata_vars[1]) %$%
      unit(n,units = "null")
  }
  #  Set the maximums in the gtables for gt1 and gt2
  gt1$heights <- as.list(maxWidth)
  gt2$heights <- as.list(maxWidth)
  gt0=grid.arrange(gt1,gt2,nrow=1,ncol=2,widths=panel_widths)
  return(list(ci=cip2,tab=tabp2,bind=gt0))
}
