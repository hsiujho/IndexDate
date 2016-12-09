#include <Rcpp.h>
using namespace Rcpp;

Function order("order");

// [[Rcpp::export]]

DataFrame index_date(DataFrame df_case, DataFrame df_ctrl){
  RNGScope scope; // ensure RNG gets set/reset
  int ncase=df_case.nrows()
    , nctrl=df_ctrl.nrows(), k;
  CharacterVector id_case=df_case["ID"]
  , id_ctrl=df_ctrl["ID"]
  , case_id(nctrl);
  IntegerVector sex_case=df_case["sexn"]
  , sex_ctrl=df_ctrl["sexn"]
  , order_case(ncase)
    ,equalmin(nctrl);
  NumericVector mdd_case=df_case["date_MDD"]
  , mdd_ctrl=df_ctrl["date_MDD"]
  , bir_case=df_case["ID_BIRn"]
  , bir_ctrl=df_ctrl["ID_BIRn"]
  , ev_ctrl=df_ctrl["ev_date"]
  , index_date=df_case["dg2_index_date"]
  , index_date56=index_date+56.0
  , min_diff(nctrl)
    , index_date_ctrl(nctrl)
    , diff(ncase)
    , rannum(ncase);
  for(int i=0;i<nctrl;i++){
    diff=ifelse((index_date56>=ev_ctrl(i))+(sex_case!=sex_ctrl(i))>0
                  ,999999.0
                  ,abs(mdd_case-mdd_ctrl(i))+abs(bir_case-bir_ctrl(i)));
    equalmin[i]=sum(diff==min(diff));
    if(equalmin[i]==1) {
      order_case=order(diff);
    } else {
      rannum=runif(ncase);
      order_case=order(diff,rannum);
    }
    k=order_case(0)-1;
    if(diff(k)!=999999.0){
      case_id(i)=id_case(k);
      min_diff(i)=diff(k);
      index_date_ctrl(i)=index_date(k);
    }
  }
  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["min_diff"]=min_diff
                             ,_["index_date"]=index_date_ctrl
                             ,_["stringsAsFactors"] = false);
}
