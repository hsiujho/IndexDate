#include <Rcpp.h>
using namespace Rcpp;

//Function Order1("order");

/*
typedef std::pair<int, double> paired;

bool cmp_second2(const paired & left, const paired & right) {
  return left.second < right.second;
}

inline Rcpp::IntegerVector Order(const Rcpp::NumericVector & x) {
  const size_t n = x.size();
  std::vector<paired> pairs; pairs.reserve(n);

  for(size_t i = 0; i < n; i++)
  {
//    if(R_IsNA(x(i))) x(i)=R_PosInf;
    pairs.push_back(std::make_pair(i, x(i)));
  }
  std::sort(pairs.begin(), pairs.end(), cmp_second2);

  Rcpp::IntegerVector result(n);
  for(size_t i = 0; i < n; i++)
    result(i) = pairs[i].first;
  return result;
}*/

// [[Rcpp::export]]

DataFrame Match2(DataFrame df_case, DataFrame df_ctrl, String ps_var="ps", double ps_threshold=0.1, int NumMatch=1){
  int ncase=df_case.nrows();
  int nctrl=df_ctrl.nrows();
  NumericVector ps_case=df_case[ps_var];
  NumericVector ps_ctrl=df_ctrl[ps_var];
  CharacterVector id_ctrl=df_ctrl["ID"];
  CharacterVector id_case=df_case["ID"];
  IntegerVector ctrl_seq=seq_len(nctrl)-1;
  CharacterVector case_id(nctrl);
  NumericVector diff_ps(nctrl);
//  ncase=10;
  for(int i=0;i<ncase;i++){
    diff_ps=abs(ps_ctrl-ps_case[i]);
    for(int j=0;j<nctrl;j++){
//      diff_ps[j]=std::abs(ps_ctrl[j]-ps_case[i]);
      if (diff_ps[j]>ps_threshold || case_id[j]!="") diff_ps[j]=-1;
    }
    IntegerVector d2=ctrl_seq[diff_ps>=0];// position of non-negative value
    //IntegerVector d2=ctrl_seq[!is_na(diff_ps)];
    //Rcpp::Rcout << "Number of non-negative value: " << d2.size()  << std::endl;
    if(d2.size()>=NumMatch){
//      NumericVector d3=diff_ps[d2];
/*      IntegerVector d4=Order(diff_ps[d2]);
      for(int k=0;k<NumMatch;k++){
        case_id[d2[d4[k]]]=id_case[i];
      }*/
      NumericVector d3=diff_ps[d2];
      NumericVector d4=clone(d3).sort();
      int l=0;
      while(l<NumMatch){
        IntegerVector d5=d2[d3==d4[l]];
        int m=0;
        while(m<d5.size() && l<NumMatch ){
          case_id[d5[m]]=id_case[i];
          m++;
          l++;
        }
      }
    }
  }
  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["stringsAsFactors"] = false);
}

/*** R
function(){
  y1=befmatch
  case=y1[y1$cohort==1,]
  ctrl=y1[y1$cohort==0,]

  on=Sys.time()
  aa<-Match2(case,ctrl,NumMatch = 1)
  Sys.time()-on
  df4=aa[aa$case_id!="",]
  rownames(df4)=NULL
  sapply(df4,function(x){
    xx=y1[y1$ID%in%x,"ps"]
    return(c(mean=mean(xx),sd=sd(xx)))
  })
}
*/
