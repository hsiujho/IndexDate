#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//配對條件一般化, 可輸入欄位名稱, 及每個欄位差異的範圍, 數量不限制

Function Order("order");

inline IntegerVector ordernmea(int nmea, NumericMatrix diff1){
  switch(nmea){
  case 1:
    return Order(diff1(_,0));
  case 2:
    return Order(diff1(_,0),diff1(_,1));
  case 3:
    return Order(diff1(_,0),diff1(_,1),diff1(_,2));
  case 4:
    return Order(diff1(_,0),diff1(_,1),diff1(_,2),diff1(_,3));
  case 5:
    return Order(diff1(_,0),diff1(_,1),diff1(_,2),diff1(_,3),diff1(_,4));
  }
  return 0;
}


// [[Rcpp::export]]

DataFrame Match3(DataFrame df_case, DataFrame df_ctrl, CharacterVector meaname
                   , NumericVector threshold, int NumMatch=1){
  int ncase=df_case.nrows();
  int nctrl=df_ctrl.nrows();
  int nmea=meaname.size();
  CharacterVector id_ctrl=df_ctrl["ID"];
  CharacterVector id_case=df_case["ID"];
  IntegerVector ctrl_seq=seq_len(nctrl)-1;
  CharacterVector case_id(nctrl);
  NumericMatrix diff(nctrl,nmea);
  NumericMatrix x_ctrl(nctrl,nmea);
  NumericMatrix x_case(ncase,nmea);
//從data.frame抽取資料放進矩陣
  for(int j=0;j<nmea;j++){
    String v1=meaname[j];
    x_ctrl(_,j)=as<NumericVector>(df_ctrl[v1]);
    x_case(_,j)=as<NumericVector>(df_case[v1]);
  }
//  int i=0;
  for(int i=0;i<ncase;i++){
//計算差異
    for(int j1=0;j1<nmea;j1++){
      diff(_,j1)=abs(x_ctrl(_,j1)-x_case(i,j1));
    }
  //差異有無同時在範圍內且還未配對
    LogicalVector InRange(nctrl);
    for(int j2=0;j2<nctrl;j2++){
      InRange[j2]=is_true(all(diff(j2,_)<=threshold)) && case_id[j2]=="";
    }
  //符合配對條件的control id的位置
    IntegerVector v1=ctrl_seq[InRange];

    if(v1.size()>=NumMatch){
  //取得submatrix
      NumericMatrix diff1(v1.size(),nmea);
      for(int j3=0;j3<v1.size();j3++){
        diff1(j3,_)=diff(v1[j3],_);
      }
  //order//呼叫R函數order, 輸出記得減一
//      IntegerVector v2=Order(diff1(_,0),diff1(_,1));
      IntegerVector v2=ordernmea(nmea,diff1);
      //配對到control的ID位置指派case的ID
      for(int j4=0;j4<NumMatch;j4++){
        case_id[v1[v2[j4]-1]]=id_case[i];
      }
    }
  }

  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["stringsAsFactors"] = false);
//  return List::create(_["x_case"]=case_id);
}

// [[Rcpp::export]]
arma::uvec myorder_arma_2(arma::mat x, arma::uvec posi=0, unsigned int head=0) {
  //uvec 為positive integer向量, vec為double向量
  vec y=x.col(0);  //抽出第i行
  vec sortedy=sort(y);
  unsigned int ny=y.n_rows;
  unsigned int nx=x.n_cols;
  if(head==0) head=ny;
  if(posi(0)==0) posi=linspace<uvec>(1,ny,ny);
  unsigned int nh=std::min(ny,head);
  uvec orderx=zeros<uvec>(nh);

  unsigned int m=0;
  while (m<nh) {
    uvec x0=find(y==sortedy(m));
    uvec x1=posi(x0);
    unsigned int nx1=x1.n_rows;
    if(nx1==1){
      orderx(m)=x1(0);
      m++;
    } else {
      uvec x2=find(find(x1>=0)<head-m);
      if(nx==1){
        orderx(m+x2)=x1(x2);
      } else {
        //submatrix
        uvec xcs=linspace<uvec>(1,nx-1,nx-1);
        mat y1=x.submat(x0,xcs);
//        Rcout << "before recusive function" << std::endl;
        orderx(m+x2)=myorder_arma_2(y1,x1,head-m);
      }
      m=m+x2.n_elem;
    }
  }
  return orderx;
}

/*** R
function(){
  require(IndexDate)
  require(magrittr)
  y$ID=as.vector(y$ID)
  case=y[y$gp==1,]
  ctrl=y[y$gp==2,]
  outp=Match3(case,ctrl,c("sex","birn"),c(0,365),4)
  sum(outp[[1]]!="")
  table(table(outp[[1]]!=""))
#  head(outp[[1]])
#  d0=sweep(ctrl[,c("birn","sex")],2,unlist(case[1,c("birn","sex")]),"-") %>%
#    abs() %>% sweep(2,c(365,0),"<=") %>% apply(1,all) %>% sum()
  require(IndexDate)
  require(magrittr)
  y1=befmatch
  set.seed(123456)
  y1$unif=runif(nrow(y1))
  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]
  outp=Match3(case,ctrl,c("ps","unif"),c(0.1,1),NumMatch=2)
  sum(outp$case_id!="")

  outp4=Match4(case,ctrl,c("ps","unif"),c(0.1,1),NumMatch=2)
#  head(outp4)
#  head(abs(sweep(ctrl[,c("ps","unif")],2,unlist(case[1,c("ps","unif")]),"-")))
  sum(outp$case_id!="")

}


function(){
  require(microbenchmark)
  microbenchmark(outp=Match3(case,ctrl,c("ps","unif"),c(0.1,1),NumMatch=2)
                 ,outp4=Match4(case,ctrl,c("ps","unif"),c(0.1,1),NumMatch=2)
                   ,times=10)

}


  */
