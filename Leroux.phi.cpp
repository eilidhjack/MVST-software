//#include <Rcpp.h>
  // [[Rcpp::depends(Matrix, RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List PhiUpdate(mat phi_current, double rho_current, double proposal_sd_phi, int n_phi, 
               vec n_neighbours, mat W_triplet, vec W_start, 
               cube y,
               cube O, mat sigma_inv, int d, int t, vec phi_acc)
{  
  //create new objects
  int track=0, row=0, ele=0; 
  double top_likelihood_sum, bottom_likelihood_sum;
  vec  new_phi(d),prior_mean(d),top_likelihood(d), bottom_likelihood(d),top_prior, bottom_prior,r;
  mat phi_proposal(n_phi,d), prior_var(d,d), inv_mat(d,d),new_lambda(d,t), current_lambda(d,t);
  
  phi_proposal=phi_current;
  //update each phi
  for(int j=0; j<n_phi; j++){
    
    //Initiate starting points to keep track of which rows we are working with
    track = n_neighbours[j] ;
    row = W_start[j];
    
    //create empty vectors for value of W matrix and phi for neighbours
    vec Wvalue(track);
    mat phi_value(track,d);
    
    
    //Now store the value from W and phi for neighbouring areas of phi[j]
    for(int i=0; i<d; i++){
      
      for(int k=0; k<track; k++){
        Wvalue[k] = W_triplet((row-1),2);
        ele = W_triplet((row-1),1)-1;
        phi_value(k,i) = phi_current(ele,i);
        row += 1;
      }
      row = W_start[j];
    }
    
    //Propose new value for phi.j
    for(int i=0; i<d; i++){
      new_phi[i] = rnorm(1,phi_current(j,i), proposal_sd_phi)[0];
    }
    
    //Calculate Acceptance ratio, r
    //Calculate quantities needed for acceptance ratio

    for(int i=0; i<d; i++){
      mat O1=O.slice(i);
      new_lambda.row(i) = exp(O1.row(j)  + new_phi[i]);
      current_lambda.row(i) = exp(O1.row(j)  + phi_current(j,i));
    }
    //Calculate prior mean
    for(int i=0; i<d; i++){
      prior_mean[i] = (sum(rho_current*trans(Wvalue)*phi_value.col(i)))/ (rho_current * n_neighbours[j] + (1-rho_current));
    }
    
    //Calculate prior variance
    inv_mat = sigma_inv * (rho_current*n_neighbours[j] + (1-rho_current));

    
    //vec top_likelihood(d);
    for(int i=0; i<d; i++){
      mat y1=y.slice(i);
      vec y11=trans(y1.row(j));
      vec ln_lambda_new=log(trans(new_lambda.row(i)));
      top_likelihood(i)=sum((y11%ln_lambda_new)-trans(new_lambda.row(i)));
    }
    
    //vec bottom_likelihood(d);
    for(int i=0; i<d; i++){
      mat y1=y.slice(i);
      vec y11=trans(y1.row(j));
      vec ln_lambda_curr=log(trans(current_lambda.row(i)));
      bottom_likelihood(i)=sum((y11%ln_lambda_curr)-trans(current_lambda.row(i)));
    }
    
    top_likelihood_sum = sum(top_likelihood);
    bottom_likelihood_sum = sum(bottom_likelihood);
    
    
    top_prior = -((trans(new_phi-prior_mean))*inv_mat*(new_phi-prior_mean))/2;
    
    bottom_prior = -((trans(trans(phi_current.row(j))-prior_mean)*inv_mat*(trans(phi_current.row(j))-prior_mean)))/2;
  
    r = (top_likelihood_sum + top_prior) - (bottom_likelihood_sum  + bottom_prior) ;
    
    //set phi_current = phi_proposal with probability r
    
    if (log(runif(1)[0]) < r[0]){
      phi_proposal.row(j) = trans(new_phi);
      phi_acc[0] += 1;
    }else{
    }
    phi_acc[1] += 1;
    
    
  }
  
  List out(2);
  out[0] = phi_proposal;
  out[1] = phi_acc;
  
  
  
  
  
  return out;
}






