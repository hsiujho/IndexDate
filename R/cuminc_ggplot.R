#
#
#
#
#
#
#

cuminc_ggplot=function(
  cuminc_obj
  ,nat_tbl
  ,x_breaks
  ,xlim
  ,line_labs
  ,line_color
  ,panel_heights
  ,upper_anno
  ,upper_labs
  ,upper_theme
  ,lower_theme
) {
  # 把 nat_tbl 轉成 data.frame 格式 ####

  if(length(x_breaks)!=NCOL(nat_tbl)) stop("length(x_breaks)!=NCOL(nat_tbl)")
  if(missing(xlim)){
    xlim=range(x_breaks)
    x_expand=c(0.1,0)
  } else {
    x_expand=c(0,0)
  }
  if(missing(line_labs)){
    line_labs=names(cuminc_obj) %>>%
      (grep(" 1$",.,value = T))
  }
  if(missing(panel_heights)){
    panel_heights=c(4,1)
  }
  #rownames(nat_tbl)=gsub("[[:punct:]]|[[:space:]]","_",rownames(nat_tbl))
  rownames(nat_tbl)=sprintf("lab%i",1:NROW(nat_tbl))
  nat_df=rbind(nat_tbl,x_breaks) %>>%
    t() %>>%
    data.frame() %>>%
    melt(id.vars="x_breaks",measure.vars=rownames(nat_tbl),value.name = "NAT") %>>%
    mutate(variable=plyr::mapvalues(variable,from = rownames(nat_tbl),to=line_labs))

  # 把 cuminc_obj 轉成 data.frame 格式 ####

  pos=names(cuminc_obj) %>>%
    (grepl(" 1$",.)) %>>%
    which()

  ylim=lapply(cuminc_obj[pos],function(x) if(class(x)=="list"){
    return(max(x$est))
  } else {
    return()
  }) %>>% unlist() %>>% max() %>>% (c(0,.))

  data=NULL
  for(i in pos){
    data=rbind(data,cbind(line=i,time=cuminc_obj[[i]]$time,est=cuminc_obj[[i]]$est))
  }
  data=data.frame(data)
  data$line=factor(data$line, levels=pos,labels = line_labs)

  np1<-ggplot(data = data, aes(x = time, y = est, group = line, color = factor(line))) +
    geom_line(size=1.5)+
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_y_continuous(labels=scales::percent) +
    scale_x_continuous(expand=x_expand,breaks=x_breaks,labels=x_breaks) +
    theme(legend.position=c(0.95,0.05)
          ,legend.justification=c(0.95,0.05)
          #        ,axis.text = element_text(size=15)
          #        ,axis.title = element_text(size=20)
    )

  if(missing(upper_anno)){
    if(any(names(cuminc_obj)=="Tests")) {
      pv=cuminc_obj$Tests[1,2]
      if(pv<.001 | pv>.999){
        pvtext=ifelse(pv<.001,"<0.001",">0.999")
      } else {
        pvtext=sprintf("=%.3f",pv)
      }
      np1=np1+annotate("text",x=xlim[2],y=0.3,vjust=0,hjust=1,label=sprintf("Log-rank P%s",pvtext),size=8)
    }
  } else {
    np1=np1+upper_anno
  }

  if(missing(upper_labs)){
    np1=np1+labs(x="Follow-up",y="Cumulative incidence",colour ="Cohort")
  } else {
    np1=np1+upper_labs
  }

  if(missing(upper_theme)){
    np1=np1+theme(axis.text = element_text(size=15)
                  ,axis.title = element_text(size=20)
    )
  } else {
    np1=np1+upper_theme
  }

  if(!missing(line_color)){
    np1=np1+scale_color_manual(values = line_color)
  }

  np2 <- ggplot(nat_df, aes(x=x_breaks, y=variable , label=NAT)) +
    geom_text(hjust=1,vjust=.5,size=4) +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand=x_expand)+
    coord_cartesian(xlim = xlim,ylim=c(0.5,NROW(nat_tbl)+1.5)) +
    theme(axis.title = element_blank()
          ,axis.ticks=element_blank()
          ,axis.text.y = element_text(size=13)
          ,axis.text.x = element_blank()
          ,panel.background = element_blank()
          ,legend.position = "none") +
    annotate("text", x = xlim[1], y = NROW(nat_tbl)+1
             , label = "Number at risk",hjust=0,vjust=0.5,size=6)
  if(!missing(lower_theme)){
    np2=np2+lower_theme
  }
  # Get maximum widths and heights
  gt1 <- ggplot_gtable(ggplot_build(np1))
  gt2 <- ggplot_gtable(ggplot_build(np2))
  rep_pos=2:4
  maxWidth <- unit.pmax(gt1$widths[rep_pos], gt2$widths[rep_pos])
  # Set the maximums in the gtables for gt1 and gt2
  gt1$widths[rep_pos] <- as.list(maxWidth)
  gt2$widths[rep_pos] <- as.list(maxWidth)
  outp=grid.arrange(gt1,gt2,nrow=2,ncol=1,heights=panel_heights)
  grid.newpage()
  grid.draw(outp)
  return(outp)
}
