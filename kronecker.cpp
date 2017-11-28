//#include <Rcpp.h>
// [[Rcpp::depends(Matrix, RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List KroneckerUpdate(mat phi_current, double rho_current, double rho_proposal,
               int n_phi, int triplet_length, vec n_neighbours, mat W_triplet, 
                mat sigma_inv)
{  
  //create new objects
  double nondiagonalQ_current, diagonalsum_current=0, nondiagonalsum_current=0, kronecker_current, nondiagonalQ_proposal, diagonalsum_proposal=0, nondiagonalsum_proposal=0, kronecker_proposal;
  vec  diagonalQ_current(n_neighbours), diagonalQ_proposal(n_neighbours);
 
 //********************************************************// 
 //*************Kronecker term for current rho*************//
 //********************************************************// 
 
 //Calculate the Q coefficient for diagonal elements
  for(int j=0; j<n_phi; j++){
    diagonalQ_current[j]=rho_current*n_neighbours[j]+(1-rho_current);
  }
  
  //Calculate the Q coefficient for the non diagonal elements
  nondiagonalQ_current= -rho_current*W_triplet(0,2);
  
  //Compute the quadratic form for the diagonal elements and multiply by Q(i,i)
  for(int j=0; j<n_phi; j++){
    vec quad = phi_current.row(j) * sigma_inv * (trans(phi_current.row(j)));
    diagonalsum_current +=  diagonalQ_current[j] * (quad[0]);
  }
  
  //Compute the quadratic form for the non diagonal elements and multiply by Q(i,j)
  for(int j=0; j<triplet_length; j++){
    int row1 = W_triplet(j,0)-1;
    int col1 = W_triplet(j,1)-1;
    vec quad2 = phi_current.row(row1) * sigma_inv * trans(phi_current.row(col1));
    nondiagonalsum_current += nondiagonalQ_current * quad2[0]; 
  }
  
  //Compute the sum
  kronecker_current = 0.5*sum(diagonalsum_current+nondiagonalsum_current);
  
  
  //********************************************************// 
  //*************Kronecker term for proposal rho*************//
  //********************************************************// 
  
  //Calculate the Q coefficient for diagonal elements
  for(int j=0; j<n_phi; j++){
    diagonalQ_proposal[j]=rho_proposal*n_neighbours[j]+(1-rho_proposal);
  }
  
  //Calculate the Q coefficient for the non diagonal elements
  nondiagonalQ_proposal= -rho_proposal*W_triplet(0,2);
  
  //Compute the quadratic form for the diagonal elements and multiply by Q(i,i)
  for(int j=0; j<n_phi; j++){
    vec quad = phi_current.row(j) * sigma_inv * (trans(phi_current.row(j)));
    diagonalsum_proposal +=  diagonalQ_proposal[j] * (quad[0]);
  }
  
  //Compute the quadratic form for the non diagonal elements and multiply by Q(i,j)
  for(int j=0; j<triplet_length; j++){
    int row1 = W_triplet(j,0)-1;
    int col1 = W_triplet(j,1)-1;
    vec quad2 = phi_current.row(row1) * sigma_inv * trans(phi_current.row(col1));
    nondiagonalsum_proposal += nondiagonalQ_proposal * quad2[0]; 
  }
  
  //Compute the sum
  kronecker_proposal = 0.5*sum(diagonalsum_proposal+nondiagonalsum_proposal);
  
  List out(2);
  out[0] = kronecker_current;
  out[1] = kronecker_proposal;
  
  
  
  
  
  return out;
}

