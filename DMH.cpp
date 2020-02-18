//
// This file contains our implementation of the Double Metropolis Hastings algorithm (Liang, 2010)
// applied for estimation of spatial point model parameters.
// Follow the comments for explanaitions and refer to the report for the underlying theory. 
// The functions are sourced into R and used in the runALL.R file. 
// 
// The code is currently not fully debugged, it stops at some point during the cycle of iterations. 
//
//

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// ==== INDICATOR FN

// [[Rcpp::export]]
int Indicator_CPP(int cell, int mark){
  int result;
  if (mark == cell) {result = 1;}
  else      {result = 0;}
  return (result);
}


// ==== DOUBLE INDICATOR FN

// [[Rcpp::export]]
int Indicator2_CPP(int cell, int mark, int cell2, int mark2){
  if (mark == cell & mark2 == cell2 ) {return 1;}
  else      {return 0;}; 
}


// ==== CONDITIONAL PROB (5) ==== vector probability (mark1,mark2,mark3)

// [[Rcpp::export]]
arma::vec conditionalZi_CPP(NumericVector z, NumericVector x, NumericVector y, 
                            NumericVector omega, NumericMatrix theta,
                            double lambda, NumericMatrix dmatrix,
                            int i, double c, int Q, int n) { 
  
  NumericVector brackets(Q);
  int q, qq, ii;
  double qqsum_theta, sum_ii, e, sumexp;
  
  for (q = 0; q < Q; q++) { 
    brackets(q) = - omega(q);
    qqsum_theta = 0;
    for (qq = 0; qq < Q; qq++)  {          // q' 
      sum_ii = 0;
      for (ii = 0; ii < n; ii++) {         // i' 
        if ((dmatrix(i-1,ii) <= c) && ((i-1) != ii)) 
        {
          e = exp(-lambda*dmatrix(i-1,ii));
          sum_ii = sum_ii + e*Indicator_CPP(z(ii),qq+1);
        }
      }
      qqsum_theta = qqsum_theta + theta(q,qq)*sum_ii;
    }
    brackets(q) = brackets(q) - qqsum_theta;
  }
  sumexp = sum(exp(brackets));
  return {exp(brackets)/sumexp}; // 
}


// FULL LOG LIKELIHOOD (4) with partition function cancelled = (ENERGY FUNCTION) (3) ====

// [[Rcpp::export]]
double Energy_CPP(NumericVector z, NumericVector x, NumericVector y,
                  NumericVector omega, NumericMatrix theta,
                  double lambda, NumericMatrix dmatrix, 
                  double c, int Q, int n) {
  
  int q, qq, i, ii, j, w; 
  double e, sum_omega, sum_indicator, qqsum_theta, sum_ii, sum_theta;
  
  sum_omega = 0;
  for (w = 0; w < Q; w++){
    sum_indicator = 0;
    for (j = 0; j < n; j++) {
      sum_indicator = sum_indicator + Indicator_CPP(z(j),(w+1));
    }
    sum_omega = sum_omega + (omega(w) * sum_indicator);
  }
  
  sum_theta = 0;
  for (q = 0; q < Q; q++) {
    qqsum_theta = 0;
    for(qq = 0; qq < Q; qq++) {                        //q'
      sum_ii=0;
      for (i = 0; i < n; i++) {                        //i
        for (ii = 0; ii < n; ii++) {                   //i'
          
          if ((dmatrix(i,ii) <= c) && (ii>=i)) 
          { 
            e = exp(-lambda*dmatrix(i,ii));
            sum_ii = sum_ii + e*Indicator2_CPP(z(i),(q+1), z(ii),(qq+1)) ;
          }
        }
      }
      qqsum_theta = qqsum_theta + theta(q,qq)*sum_ii;
    }
    sum_theta = sum_theta + qqsum_theta;
  }
  return(-sum_omega - sum_theta);
}



// ==== DMH ===

// [[Rcpp::export]]
Rcpp::List DMH(NumericVector omega, NumericVector omega_start, 
               NumericMatrix theta,
               NumericVector x, NumericVector y,
               NumericVector z, arma::vec  marks,
               NumericMatrix dmatrix,
               int iter, int m, int n, int Q, double c,  
               double lambda, double mu_omega, double sigma_omega, double tau2,
               double mu_theta, double sigma_theta, double lambda_a, double lambda_b)  {
  
  double z_star = 0;
  double LogRatio = 0;
  NumericVector omega_star = omega_start; 
  NumericMatrix theta_star(Q,Q);
  double la_star = 0;
  arma::vec probabilityZi(Q); 
  arma::mat omega_matrix(iter, Q);
  arma::mat zstar_mat(m,n);
  NumericVector zstar(n); 
  NumericVector lambda_vec(iter);
  int accept_omega = 0;
  int accept_theta = 0;
  int accept_la = 0;
  //int itercount = 0;
  double u;
  int t, q, i, k, qq, s, f;
  
  
  // ====external DMH cycle
  
     for (t = 0; t < iter; t++)  { 
    
    // === omega update == 
    // propose new omega -> omega_star vector by elements
    
    for (q = 0; q < Q-1; q++)  {                       // could be Q, but we add dfs
      for (s = 0; s < Q; s++){
        omega_star(s) = omega(s);}
      for (s = 0; s < n; s++){
        zstar(s) = z(s);}
      
        omega_star(q) = ::rnorm(1, omega(q),tau2) (0); 
      
    
   // generate vector z* with conditional (5), M times gibbs sampler
      for (k = 0; k < m; k++) {
        for (i = 0; i < n; i++) {
          probabilityZi = conditionalZi_CPP(zstar, x, y, omega_star, theta, lambda, dmatrix, i, c, Q, n);
          z_star = RcppArmadillo::sample(marks, 1, true, probabilityZi) (0);           //length (marks) should be = Q number
          zstar(i) = z_star;
          
          if ( t == iter-1) {zstar_mat(k,i) = zstar(i);}
        }
      }
      
      // calculate log of ratio, accept/reject omega_star    
      LogRatio = Energy_CPP(zstar, x, y, omega, theta, lambda, dmatrix, c, Q, n) +
        Energy_CPP(z, x, y, omega_star, theta, lambda, dmatrix, c, Q, n)+
        ((omega_star(q)-mu_omega)*(omega_star(q)-mu_omega) /2/sigma_omega/sigma_omega) -
        Energy_CPP(z, x, y, omega, theta, lambda, dmatrix, c, Q, n) -
        Energy_CPP(zstar, x, y, omega_star, theta, lambda, dmatrix, c, Q, n) -
        ((omega(q)-mu_omega)*(omega(q)-mu_omega)/2/sigma_omega/sigma_omega);
      
      // no accept_prob = min(1, ratio) as not nessesary 
      u = arma::randu();
      if (log(u) < LogRatio) {
        omega(q) = omega_star(q);
        omega_matrix(t,q) = omega_star(q);
        accept_omega = accept_omega + 1;
      }
      else { omega_matrix(t,q) = omega(q); }
    }
    
    // === theta update == 
    // propose new theta -> theta_star matrix by elements  
    
    for (q = 0; q < Q-1; q++)  {                     // could be Q, but we add dfs
      for (qq = 0; qq < Q; qq++) {
        
        if (q<=qq) {
          for (s = 0; s < Q; s++){
            for (f = 0; f < Q; f++){
              theta_star(s,f) = theta(s,f); 
              theta_star(f,s) = theta(f,s);}
          }
          
          for (s = 0; s < n; s++){
            zstar(s) = z(s);}
          
          theta_star(q,qq) = ::rnorm(1, theta(q,qq), tau2) (0); 
          
          // generate set z* with conditional eq (5), M times gibbs sampler
          for (k = 0; k < m; k++) {
            for (i = 0; i < n; i++) {
              probabilityZi = conditionalZi_CPP(zstar, x, y, omega, theta_star, lambda, dmatrix, i, c, Q, n);
              z_star = RcppArmadillo::sample(marks, 1, true, probabilityZi) (0);      //length (marks) should be = Q number
              zstar(i) = z_star;
            }
          }
          
          // calculate log of ratio, accept/reject theta_star    
          LogRatio = Energy_CPP(zstar, x, y, omega, theta, lambda, dmatrix, c, Q, n) +
            Energy_CPP(z, x, y, omega, theta_star, lambda, dmatrix, c, Q, n)+
            ((theta_star(q,qq)-mu_theta)*(theta_star(q,qq)-mu_theta)/2/sigma_theta/sigma_theta) -
            Energy_CPP(z, x, y, omega, theta, lambda, dmatrix, c, Q, n) -
            Energy_CPP(zstar, x, y, omega, theta_star, lambda, dmatrix, c, Q, n) -
            ((theta(q,qq)-mu_theta)*(theta(q,qq)-mu_theta) /2/sigma_theta/sigma_theta);
          
          // no accept_prob 
          u = arma::randu();
          if (log(u) < LogRatio) {
            theta(q,qq) = theta_star(q,qq);
            theta(qq,q) = theta_star(q,qq); // symmetrical
            accept_theta = accept_theta + 1;
          }
        }
      }
    }
    
    
    // === lambda update == 
    // propose new lambda -> la_star  
    
      la_star = lambda; 
      for (s = 0; s < n; s++){
        zstar(s) = z(s);}
      do { 
        la_star = rnorm(1, lambda,0.05) (0);
      } while ( la_star<=0);
      
      // generate vector z* with conditional (5), M times gibbs sampler
      for (k = 0; k < m; k++) {
        for (i = 0; i < n; i++) {
          probabilityZi = conditionalZi_CPP(zstar, x, y, omega, theta, la_star, dmatrix, i, c, Q, n);
          z_star = RcppArmadillo::sample(marks, 1, true, probabilityZi) (0); 
          //length (marks) should be = Q number
          zstar(i) = z_star;
        }
      }
      
      // calculate log of ratio, accept/reject la_star    
      LogRatio = Energy_CPP(zstar, x, y, omega, theta, lambda, dmatrix, c, Q, n) +
        Energy_CPP(z, x, y, omega, theta, la_star, dmatrix, c, Q, n)+
        (lambda_a - 1)*(log(la_star) - log(lambda)) + lambda_b*(la_star - lambda) -
        Energy_CPP(z, x, y, omega, theta, lambda, dmatrix, c, Q, n)-
        Energy_CPP(zstar, x, y, omega, theta, la_star, dmatrix, c, Q, n);
      
      // no accept_prob 
      u = arma::randu();
      if (log(u) < LogRatio) {
        lambda = la_star;
        lambda_vec(t) = la_star;
        accept_la = accept_la + 1;
      }
      else { lambda_vec(t) = lambda; }
    
    
    // ========= iter report
    
  if ((t)%500==0) {Rcout << "iterations done = " << t+1 << std::endl;}
    // ========= 
  }
  
  
  return Rcpp::List::create(Rcpp::Named("omega_matrix") = omega_matrix, 
                            Rcpp::Named("t") = t,
                            Rcpp::Named("theta_matrix") = theta,
                            Rcpp::Named("lambda_vec") = lambda_vec,
                            Rcpp::Named("accept_omega") = accept_omega,
                            Rcpp::Named("accept_theta") = accept_theta,
                            Rcpp::Named("accept_lambda") = accept_la,
                            Rcpp::Named("zstar_mat") = zstar_mat);
}
