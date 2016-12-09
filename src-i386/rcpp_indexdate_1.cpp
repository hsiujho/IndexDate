#include <Rcpp.h>
using namespace Rcpp;

inline IntegerVector find_pos(double y, NumericVector x, int n=1){
  IntegerVector z(n);
  int i=0, j=0;
  while(i<x.size() && j<n){
    if (x[i]==y){
      z[j]=i;
      j++;
    }
    i++;
  }
  return z;
}

// [[Rcpp::export]]

DataFrame index_date_1(DataFrame df_case, DataFrame df_ctrl){
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
    double o1=min(diff);
    equalmin[i]=sum(diff==o1);
    if(equalmin[i]==1) {
      NumericVector o2(1);
      o2[0]=o1;
      IntegerVector o3=match(o2,diff);
      k=o3(0)-1;
    } else {
      IntegerVector o3=find_pos(o1,diff,equalmin[i]);
      NumericVector o4=runif(equalmin[i]);
      NumericVector o6(1);
      o6[0]=min(o4);
      IntegerVector o5=match(o6,o4);
      k=o3(o5[0]-1);
    }
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

/*** R
function(){
  load('E:/soft/R/IndexDate/data/MDD.rda')
  y=MDD
  y$sexn=ifelse(y$ID_SEX=="M",1,0)
  ctrl=y[y$cohort!=1,]
  case0=y[y$cohort==1,]
  on=Sys.time()
  aa=index_date_1(case0,ctrl)
  Sys.time()-on
  sum(aa$case_id!="")
}
*/
