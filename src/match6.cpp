
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]

#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <progress.hpp>
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
  const NumericMatrix diff2;
  const NumericVector thr;
  uvec& InRange;
  check(const NumericMatrix diff2, const NumericVector thr, uvec& InRange) : thr(thr), diff2(diff2), InRange(InRange) {}
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) {
      if(InRange(i)==1 && Rcpp::any(diff2(i,_)>thr).is_true()) InRange(i)=0;
    }
  }
};

// [[Rcpp::export]]

DataFrame Match6(DataFrame df_case, DataFrame df_ctrl, CharacterVector meaname
                   , NumericVector threshold, unsigned int NumMatch=1){
  unsigned int ncase=df_case.nrows();
  unsigned int nctrl=df_ctrl.nrows();
  unsigned int nmea=meaname.size();
  CharacterVector id_ctrl=df_ctrl["ID"];
  CharacterVector id_case=df_case["ID"];
  IntegerVector ctrl_seq=seq_len(nctrl)-1;
  CharacterVector case_id(nctrl);
  NumericMatrix x_ctrlr(nctrl,nmea);
  NumericMatrix x_caser(ncase,nmea);
  //從data.frame抽取資料放進矩陣
  for(unsigned int j=0;j<nmea;j++){
    String v1=meaname[j];
    x_ctrlr(_,j)=as<NumericVector>(df_ctrl[v1]);
    x_caser(_,j)=as<NumericVector>(df_case[v1]);
  }
  mat x_ctrl(x_ctrlr.begin(), x_ctrlr.nrow(), x_ctrlr.ncol(), false);
  mat x_case(x_caser.begin(), x_caser.nrow(), x_caser.ncol(), false);
//  Rcout << "before definition of myorder" << std::endl;

  uvec aa(nctrl,fill::ones); //記錄case_id是否還未被配對, 1(true)為還未配對, 0(false)為已配對
//  uvec InRange(nctrl,fill::ones); //記錄差異是否在範圍內, 1(true)為在範圍內, 0(false)為在範圍外
  uvec v3(1,fill::zeros);
//  NumericVector thr(threshold.begin(),threshold.end());

  //  mat diff(nctrl,nmea); 因為不會使用到之前的訊息, 迴圈內宣告效率會比較好, 可能不用找被定義過的記憶體位址.

  //  unsigned int i=0;
  Progress p(ncase, 1);
  for(unsigned int i=0;i<ncase;i++){
    if (Progress::check_abort() )
      return -1.0;
    p.increment(); // update progress
    //計算差異
//    mat diff=abs(x_ctrl-repmat(x_case.row(i),nctrl,1));

    //差異有無同時在範圍內且還未配對
    uvec InRange=aa;

    NumericMatrix diff2(nctrl,nmea);
    for(uword j=0;j<nmea;j++){
      diff2(_,j)=abs(x_ctrlr(_,j)-x_caser(i,j));
    }

for(uword j=0;j<nctrl;j++){
  if(InRange(j)==1 && Rcpp::any(diff2(j,_)>threshold).is_true()) InRange(j)=0;
    }

//    for(uword j=0;j<nctrl;j++){
//      if(InRange(j)==1){
//        diff2(j,_)=abs(x_ctrlr(j,_)-x_caser(i,_));
//        if(Rcpp::any(diff2(j,_)>threshold).is_true()) InRange(j)=0;
//      }
//    }

    mat diff(diff2.begin(),diff2.nrow(),diff2.ncol(),false);



//    for(unsigned int p=0;p<nctrl;p++){
//      uword q=0;
//      while(InRange(p)==1 && q<nmea){
//        if(diff2(p,q)>threshold(q)) InRange(p)=0;
//        q++;
//      }
//    }

//          Rcout << "here" << std::endl;
    // create the worker
//    check Check(thr, diff2, InRange);
    // call it with parallelFor
//    parallelFor(0, nctrl, Check);

    //符合配對條件的control id的位置
    uvec v1=find(InRange==1);
    IntegerVector v11(v1.begin(),v1.end());
    if(v1.n_elem>=NumMatch){
      //取得submatrix
      mat diff1=diff.rows(v1);
      //order
      uvec v2=myorder2(diff1,v3,NumMatch)-1;
      //配對到control的ID位置指派case的ID
      for(unsigned int j4=0;j4<NumMatch;j4++){
        case_id[v1[v2[j4]]]=id_case[i];
        aa[v1[v2[j4]]]=0;
      }
    }
  }

  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["stringsAsFactors"] = false);
}
