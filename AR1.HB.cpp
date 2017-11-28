//#include <Rcpp.h>
// [[Rcpp::depends(Matrix, RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List HBUpdate(cube HB_current_arr, vec alpha_current, vec sigma2_current,
               vec sigma2_update_HB,  mat X_phi, cube y_arr, cube O_arr,
               cube XtHB_arr, int d, int S, int tp, vec HBacc)
{  
  //create new objects
  double HB_prior_mean, HB_prior_var, prior_top, prior_bottom, top_likelihood, bottom_likelihood;
  mat   top, bottom, r, XtHB_prop_slice, XtHB_prop;
  cube HB_proposal_arr;
  vec y, O_prop, O_curr;
  
  HB_proposal_arr = HB_current_arr;
  
  for(int j=0; j<d; j++){
    for(int s=0; s<S; s++){
      for(int t=0; t<tp; t++){
        //Calculate the mean and variance for the prior distribution based
        //on which time point we're at.
        if(t==0){
          HB_prior_mean=HB_current_arr((t+1),s,j)/alpha_current[j];
          HB_prior_var=sigma2_current[j]/pow(alpha_current[j],2);
        }
        else if(t==(tp-1)){
          HB_prior_mean = alpha_current[j]*HB_current_arr((t-1),s,j);
          HB_prior_var = sigma2_current[j];
        }
        else{
          HB_prior_mean=(alpha_current[j]*(HB_current_arr((t-1),s,j) + HB_current_arr((t+1),s,j)))/(1 + pow(alpha_current[j],2));
          HB_prior_var=sigma2_current[j]/(1 + pow(alpha_current[j],2)); 
        }
        
        //draw proposal value
        HB_proposal_arr(t,s,j)=rnorm(1,HB_current_arr(t,s,j), sqrt(sigma2_update_HB[j]))[0];
        
        //create new XtHB for proposal values when calculating densities
        //use reshape to turn cube form of HB_proposal_arr into a matrix to allow for 
        //correct matrix multiplication
;
        XtHB_prop_slice = HB_proposal_arr.slice(j);
        XtHB_prop=X_phi * trans(XtHB_prop_slice.row(t));
        
        //calculate prior densities 
        prior_top = 0.5 * pow((HB_proposal_arr(t,s,j) - HB_prior_mean),2) / HB_prior_var;
        prior_bottom = 0.5 * pow((HB_current_arr(t,s,j) - HB_prior_mean),2) / HB_prior_var;

        //calculate acceptance ratio

        O_prop = exp(O_arr.slice(j).col(t) + XtHB_prop);
        O_curr = exp(O_arr.slice(j).col(t) + XtHB_arr.slice(j).col(t));
        
        top_likelihood = sum((y_arr.slice(j).col(t)%log(O_prop) - O_prop));
        bottom_likelihood = sum((y_arr.slice(j).col(t)%log(O_curr) - O_curr));
        
        top = top_likelihood - prior_top;
        bottom = bottom_likelihood - prior_bottom;
        
        r = top - bottom;
        
        //accept/reject proposal value
        if (log(runif(1)[0]) < r[0]){
          HB_current_arr(t,s,j) = HB_proposal_arr(t,s,j);
          XtHB_arr.slice(j).col(t) = XtHB_prop;
          //count the number of accepted draws for this block
          HBacc[j] += 1;
        }
        else{
          HB_proposal_arr(t,s,j) = HB_current_arr(t,s,j);
          XtHB_prop = XtHB_arr.slice(j).col(t);
        }

        
        
        
        
      }
    }
  }
  
  
  
  List out(2);
  out[0] = HB_current_arr;
  out[1] = HBacc;
  
  return out; 
  
}