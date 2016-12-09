#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//可以配對兩組以上的群組

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

DataFrame Match11(DataFrame df_case
                  , DataFrame df_ctrl
                  , CharacterVector meaname
                  , arma::vec threshold
                  , String ID_var
                  , String group_var
                  , CharacterVector lab_ctrl
                  , arma::uvec NumMatch){
  unsigned int ncase=df_case.nrows();
  unsigned int nctrl=df_ctrl.nrows();
  unsigned int nmea=meaname.size();
  unsigned int gctrl=lab_ctrl.size();
  uword TN=sum(NumMatch);
  CharacterVector ctrl_group=df_ctrl[group_var];
  CharacterVector id_ctrl=df_ctrl[ID_var];
  CharacterVector id_case=df_case[ID_var];
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

//  mat x_ctrl(nctrl,nmea);
//  mat x_case(ncase,nmea);
//  for(unsigned int j=0;j<nmea;j++){
//    String v1=meaname[j];
//    x_ctrl.col(j)=as<vec>(df_ctrl[v1]);
//    x_case.col(j)=as<vec>(df_case[v1]);
//  }



//  Rcout << "before definition of myorder" << std::endl;

  //uvec aa(nctrl,fill::ones); //記錄case_id是否還未被配對, 1(true)為還未配對, 0(false)為已配對

  umat aa(nctrl,gctrl,fill::zeros);
  for(unsigned int j=0;j<gctrl;j++){
//    aa.col(j)= as<vec>(ctrl_group==lab_ctrl(j));
    for(unsigned int k=0;k<nctrl;k++){
      aa(k,j)= ctrl_group(k)==lab_ctrl(j);
    }
  }

  uvec InRange(nctrl,fill::ones); //記錄差異是否在範圍內, 1(true)為在範圍內, 0(false)為在範圍外
  //  uword q=0;
  uvec v3(1,fill::zeros);


  //  unsigned int i=0;
  for(unsigned int i=0;i<ncase;i++){
    //計算差異
    mat diff=abs(x_ctrl-repmat(x_case.row(i),nctrl,1));
    //暫時存放配對到的case_id的位置, 待其他cohort都有配對完成, 再移置完成的物件
    uvec v4(nctrl,fill::zeros);
    for(uword j=0;j<gctrl;j++){
    //差異有無同時在範圍內且還未配對
    //    uvec InRange=ones<uvec>(nctrl);
      InRange=aa.col(j);
      for(unsigned int p=0;p<nctrl;p++){
        //      if(case_id(p)!="") InRange(p)=0;
        uword q=0;
        while(InRange(p)==1 && q<nmea){
          if(diff(p,q)>threshold(q)) InRange(p)=0;
          q++;
        }
      }
      //符合配對條件的control id的位置
      uvec v1=find(InRange==1);
      if(v1.n_elem>=NumMatch(j)){
        //取得submatrix
        mat diff1=diff.rows(v1);
        uvec v2=myorder2(diff1,v3,NumMatch(j))-1;
        //配對到control的ID位置指派case的ID
        //      case_id(v1(v2-1))=id_case[i];
        for(unsigned int j4=0;j4<NumMatch(j);j4++){
          v4[v1[v2[j4]]]=1;
        }
      }

    }
    if(sum(v4)==TN){
      uvec v5=find(v4==1);
      for(uword j4=0;j4<TN;j4++){
//        v4[v1[v2[j4]]]=1;
//        case_id[v1[v2[j4]]]=id_case[i];
//        aa.row(v5[j4])=zeros<uvec>(gctrl,);
        case_id[v5[j4]]=id_case[i];
        for(uword k=0;k<gctrl;k++){
          aa(v5[j4],k)=0;
        }
      }

    }
  }


//  return Rcpp::List::create(_["aa"]= aa);
  return DataFrame::create(_["ID"]=id_ctrl,_["case_id"]=case_id
                             ,_["stringsAsFactors"] = false);
}

/*** R
#
# 這是 R 的程式碼

#dbpath="G:/20160826 GCA+HP/R data"
#y=readRDS(file.path(dbpath,"gca_res_20161027.rds"))

#  y0=dplyr::filter(y,!is.na(res_date)&!is.na(ID_SEX)&!is.na(ID_BIRn)) %>%
#    mutate(age=as.numeric(index_date-ID_BIRn)/365.25
#           ,pre_sto_op=ifelse(is.na(res2_date),0,ifelse(res2_date<res_date,1,0))
#           ,pre_HP_era=ifelse(is.na(hp_fst_date),0,ifelse(is.na(hp_date_after_res)&hp_fst_date<res_date,1,0)))

#  y1=dplyr::filter(y0,fu_yr>1&non151can==0&pre_sto_op==0&pre_HP_era==0) %>>% select(ID,cohort,age,ID_SEX,pu)
#  GCAHP=y1
#save(GCAHP,file="data/GCAHP.rdata")

load("data/GCAHP.rdata")

table(GCAHP$cohort)

a0=Match11(df_case=dplyr::filter(GCAHP,cohort==">1y")
        , df_ctrl=dplyr::filter(GCAHP,cohort!=">1y")
        , meaname=c("age","ID_SEX","pu")
        , threshold=c(3,0,0)
        , ID_var="ID"
        , group_var="cohort"
        , lab_ctrl=c("<1y","No")
        , NumMatch=c(1,4))

a1=dplyr::filter(a0,case_id!="")
table(table(a1$case_id))
a2=unique(do.call(c,a1))

a3=dplyr::filter(GCAHP,ID%in%a2)
table(a3$cohort)
*/
