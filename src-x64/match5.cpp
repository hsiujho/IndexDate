
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;


inline arma::uvec myorder2(arma::mat x, arma::uvec posi=0, unsigned int head=0) {
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
        orderx(m+x2)=myorder2(y1,x1,head-m);
      }
      m=m+x2.n_elem;
    }
  }
  return orderx;
}

struct check: public Worker {
//  CharacterVector case_id;
//  uword b;
  vec threshold;
  mat diff;
  uvec& InRange;
  check(vec threshold, mat diff, uvec& InRange) : threshold(threshold), diff(diff), InRange(InRange) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) {
//      if(a(i)!="") InRange(i)=0;
//      unsigned int q=0;
//      while(InRange(i)==1 && q<b){
//        if(d(i,q)>c(q)) InRange(i)=0;
//        q++;
//      }

      if(InRange(i)==1 && any(diff.row(i)>threshold)) InRange(i)=0;

    }
  }
};



// [[Rcpp::export]]

DataFrame Match5(DataFrame df_case, DataFrame df_ctrl, CharacterVector meaname
                   , arma::vec threshold, unsigned int NumMatch=1){
  unsigned int ncase=df_case.nrows();
  unsigned int nctrl=df_ctrl.nrows();
  unsigned int nmea=meaname.size();
  CharacterVector id_ctrl=df_ctrl["ID"];
  CharacterVector id_case=df_case["ID"];
//  IntegerVector ctrl_seq=seq_len(nctrl)-1;
  CharacterVector case_id(nctrl);
  NumericMatrix x_ctrlr(nctrl,nmea);
  NumericMatrix x_caser(ncase,nmea);
  //從data.frame抽取資料放進矩陣

  mat x_ctrl(nctrl,nmea);
  mat x_case(ncase,nmea);
  for(unsigned int j=0;j<nmea;j++){
    String v1=meaname[j];
    x_ctrl.col(j)=as<vec>(df_ctrl[v1]);
    x_case.col(j)=as<vec>(df_case[v1]);
  }

  LogicalVector aa(nctrl,1); //記錄case_id是否還未被配對, 1(true)為還未配對, 0(false)為已配對
  uvec InRange(nctrl,fill::ones); //記錄差異是否在範圍內, 1(true)為在範圍內, 0(false)為在範圍外

//  Rcout << "before definition of myorder" << std::endl;
  uvec v3(1,fill::zeros);

  for(unsigned int i=0;i<ncase;i++){
    //計算差異
    mat diff=abs(x_ctrl-repmat(x_case.row(i),nctrl,1));
    //差異有無同時在範圍內且還未配對
//    uvec InRange=ones<uvec>(nctrl);
//    InRange(case_id!="")=0;
//    for(uword p=0;p<nctrl;p++){
//      if(case_id(p)!="") InRange(p)=0;
//    }
    InRange=as<arma::uvec>(aa);

//    for(unsigned int p=0;p<nctrl;p++){
//      if(case_id(p)!="") InRange(p)=0;
//      unsigned int q=0;
//      while(InRange(p)==1 && q<nmea){
//        if(diff(p,q)>threshold(q)) InRange(p)=0;
//        q++;
//      }
//    }

    // create the worker
    check Check(threshold, diff, InRange);
    // call it with parallelFor
    parallelFor(0, nctrl, Check);


    //符合配對條件的control id的位置
    uvec v1=find(InRange==1);
    if(v1.n_elem>=NumMatch){
      //取得submatrix
      mat diff1=diff.rows(v1);
      uvec v2=myorder2(diff1,v3,NumMatch)-1;
      //配對到control的ID位置指派case的ID
      //      case_id(v1(v2-1))=id_case[i];
      for(uword j4=0;j4<NumMatch;j4++){
        case_id[v1[v2[j4]]]=id_case[i];
        aa[v1[v2[j4]]]=0;
      }
    }
  }

  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["stringsAsFactors"] = false);
}
