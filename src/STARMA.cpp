
#include <TMB.hpp>                                // Links in the TMB libraries
#include <math.h>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA ARRAY/MATRIX
  //DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_ARRAY(W);
  //DATA_VECTOR(init);
  
  
  // DATA_MATRIX(w); //neighbour matrix size m x m
  //DATA_SCALAR(h0);
  
  // PARAMETERS
  PARAMETER(mu);
  PARAMETER_MATRIX(phi);
  PARAMETER_MATRIX(theta);
  
  PARAMETER(sigma);
  //PARAMETER(omega);
  //PARAMETER_MATRIX(alpha);
  //PARAMETER_MATRIX(beta);
  
  

  
// dimension of x
  array<Type>yhat(y.rows(),y.cols()); 
  array<Type>x(y.rows(),y.cols()); 
  //array<Type>x2(y_dim(0),y_dim(1)); 
  //array<Type>sigma(y_dim(0),y_dim(1)); 
  
   // for(int k=0;k<y_dim(0);k++){
  //  x(k,0)=0.0;
  //  x2(k,0)=0.0;
  //  sigma(k,0)=init(k);
  //  }
  // Declare the "objective function" (neg. log. likelihood)
  Type f;
  f=0;
  for (int t = 0; t < phi.cols();t++) {
  for(int k=0; k < y.rows(); k++) {
    yhat(k,t)=mu;
  }
  }
  for (int t = phi.cols(); t < y.cols()-1; t++) {
    for(int k=0; k<y.rows(); k++){
      yhat(k,t+1) = mu;
      for(int i = 0; i < phi.cols(); i++){
        for(int j = 0; j < phi.rows(); j++){
          for(int u =0; u < y.rows(); u++) {
            yhat(k,t+1) +=  phi(j,i)*W(k,u,j) * (y(u,t-i)-mu);
          }
        }
      }
      x(k,t+1)= y(k,t+1)-yhat(k,t+1);
      for(int i = 0; i < theta.cols(); i++){
        for(int j = 0; j < theta.rows(); j++){
          for(int u =0; u < y.rows(); u++) {
            x(k,t+1) -=  theta(j,i)*W(k,u,j) * x(u,t-i);
          }
        }
      }
      f -= dnorm(x(k,t+1),Type(0.0),sigma,true);
    //f -= dnorm(y(k,t+1),muy(k,t+1),sigma,true);
    }
  }
  
  //REPORT(x);
  //REPORT(sigma);
  REPORT(yhat);
  return f;
}


