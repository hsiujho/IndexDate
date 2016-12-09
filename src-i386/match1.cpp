#include <Rcpp.h>
using namespace Rcpp;

typedef std::pair<int, double> paired;

bool cmp_second(const paired & left, const paired & right) {
  return left.second < right.second;
}

inline Rcpp::IntegerVector Order(const Rcpp::NumericVector & x) {
  const size_t n = x.size();
  std::vector<paired> pairs; pairs.reserve(n);

  for(size_t i = 0; i < n; i++)
    pairs.push_back(std::make_pair(i, x(i)));

  std::sort(pairs.begin(), pairs.end(), cmp_second);

  Rcpp::IntegerVector result(n);
  for(size_t i = 0; i < n; i++)
    result(i) = pairs[i].first;
  return result;
}

// [[Rcpp::export]]

DataFrame Match1(DataFrame df_case, DataFrame df_ctrl, String ps_var="ps", double ps_threshold=0.1, int NumMatch=1){
  int ncase=df_case.nrows();
  int nctrl=df_ctrl.nrows();
//  double ps_threshold=0.1;
//  NumericVector ps_case=df_case["ps"];
//  NumericVector ps_ctrl=df_ctrl["ps"];
  NumericVector ps_case=df_case[ps_var];
  NumericVector ps_ctrl=df_ctrl[ps_var];
  CharacterVector id_ctrl=df_ctrl["ID"];
  CharacterVector id_case=df_case["ID"];
  CharacterVector case_id(nctrl);
  NumericVector diff_ps(nctrl);

  for(int i=0;i<ncase;i++){
    diff_ps=abs(ps_ctrl-ps_case[i]);
    for(int j=0;j<nctrl;j++){
//      diff_ps[j]=std::abs(ps_ctrl[j]-ps_case[i]);
      if (diff_ps[j]>ps_threshold || case_id[j]!="") diff_ps[j]=9.0;
    }
    IntegerVector ps_order=Order(diff_ps);

    if(diff_ps[ps_order[NumMatch-1]]!=9.0){
      for(int k=0;k<NumMatch;k++){
        case_id[ps_order[k]]=id_case[i];
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
  aa<-Match1(case,ctrl,NumMatch = 1)
  Sys.time()-on
  df3=aa[aa$case_id!="",]
  rownames(df3)=NULL
}
*/
