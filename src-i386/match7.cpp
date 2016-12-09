// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

//#ifdef _OPENMP
#include <omp.h>
//#endif
#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace Rcpp;
using namespace arma;

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

// [[Rcpp::export]]

DataFrame Match7(DataFrame df_case, DataFrame df_ctrl, CharacterVector meaname
                   , arma::vec threshold, unsigned int NumMatch=1, int threads=1){
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

  uvec aa(nctrl,fill::ones); //記錄case_id是否還未被配對, 1(true)為還未配對, 0(false)為已配對
  uvec v3(1,fill::zeros);

#ifdef _OPENMP
  if ( threads > 0 )
    omp_set_num_threads( threads );
  REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif

  Progress p(ncase, 1);
  for(unsigned int i=0;i<ncase;i++){
    if (Progress::check_abort() )
      return -1.0;
    p.increment(); // update progress
    //計算差異
    mat diff=abs(x_ctrl-repmat(x_case.row(i),nctrl,1));
    //差異有無同時在範圍內且還未配對
    uvec InRange=aa;//記錄差異是否在範圍內, 1(true)為在範圍內, 0(false)為在範圍外
    uword q=0;
#pragma omp parallel for private(q)
    for(unsigned int j=0;j<nctrl;j++){
      q=0;
      while(InRange(j)==1 && q<nmea){
        if(diff(j,q)>threshold(q)) InRange(j)=0;
        q++;
      }
//      if(InRange(j)==1 && any(diff.row(j)>threshold)) InRange(j)=0;
    }

    //符合配對條件的control id的位置
    uvec v1=find(InRange==1);
    if(v1.n_elem>=NumMatch){
      //取得submatrix
      mat diff1=diff.rows(v1);
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
