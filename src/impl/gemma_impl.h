//
// Created by song on 5/20/18.
//

#ifndef MLMM_CPP_GEMMA_IMPL_H
#define MLMM_CPP_GEMMA_IMPL_H
#include "../util/util.h"
#include "../model/model.h"
#include "My_linear_algebra_impl.h"
//void _get_UtW_( const Eigen_result & eigen_L, const My_matrix<double> & W, const int & n, const int & q, My_matrix<double> & UtW);
//void _get_Uty_( const Eigen_result & eigen_L, const My_Vector<double> & y, const size_t & n, My_Vector<double> & Uty);
double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty, const My_matrix<double> & x,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                         const std::string & method, const int & maxiter );
double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const  double &ulim, const double & eps,
                         const std::string & method, const int & maxiter );
double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                         const std::string & method, const int & maxiter, double & lambda );
double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty, const My_matrix<double> & x,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                         const std::string & method, const int & maxiter, double & lambda );

void CalcLmmVgVe(const My_Vector<double> & y, const Eigen_result & eigen_L, const My_matrix<double> &UtW,
                 const My_Vector<double> & Uty, const double & lambda, double &vg,
                 double &ve);

#endif //MLMM_CPP_GEMMA_IMPL_H
