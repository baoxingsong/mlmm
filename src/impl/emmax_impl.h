//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_EMMAX_IMPL_H
#define MLMM_CPP_EMMAX_IMPL_H

#include "../util/util.h"
#include "../model/model.h"
#include "My_linear_algebra_impl.h"
#include <iomanip>
Eigen_result _get_eigen_L_(const My_matrix<double> & k);

Eigen_result _get_eigen_R_(const My_matrix<double> & x, const My_matrix<double> & k);

// log-likelihoods (eq. 6 from paper), this is the emma.delta.ML.LL.wo.Z function in R emma
double _ll_(const double & delta, const Eigen_result & eigen_R, const Eigen_result & eigen_L, const My_Vector<double> & sq_etas);

//diffrentiated log-likelihoods (eq. 8 from paper) the emma.delta.ML.dLL.wo.Z function in R emma
double _dll_( const double & delta, const Eigen_result & eigen_R, const Eigen_result & eigen_L, const My_Vector<double> & sq_etas);

//newton_ll( new_opt_delta, eps, maxiter, n, eig_R_values, eig_L_values, sq_etas, p); 
int newton_ll( double* new_opt_delta, const double & eps, const int & js, const int & n, const Eigen_result & eig_R, const Eigen_result & eig_L, const My_Vector<double> & sq_etas);

//log-likelihoods (eq. 7 from paper)
double _rell_(const double & delta, const Eigen_result & eigen_R, const My_Vector<double> & sq_etas);

// diffrentiated log-likelihoods (*2) (eq. 9 from paper)
double _redll_( const double & delta, const Eigen_result & eig_R, const My_Vector<double> & sq_etas );

//int newton_reml( double  new_opt_delta, const double & eps, const int & js, Eigen_result & eig_R, My_Vector<double> & sq_etas );
//int newton_reml( double* new_opt_delta, const double & eps, const int & js, double delta, Eigen_result & eig_R, My_Vector<double> & sq_etas );
Emma_result emma_estimates ( const My_Vector<double> & y, const My_matrix<double> & k, const My_matrix<double> & x,
                             const My_matrix<double> & xs, const Eigen_result & eigen_L, const int & ngrids,
                             const double & llim, const double &ulim, const double & eps,
                            const std::string & method, const int & maxiter );
Emma_result emma_estimates ( const My_Vector<double> & y, const My_matrix<double> & k, const My_matrix<double> & x,
                             const My_matrix<double> & xs, const Eigen_result & eigen_L, const int & ngrids,
                             const double & llim, const double &ulim, const double & eps,
                             const std::string & method, const int & maxiter, double & delta);
#endif //MLMM_CPP_EMMAX_IMPL_H
