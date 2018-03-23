//
// Created by Baoxing song on 12.03.18.
//

#ifndef MLMM_CPP_AVERAGE_INFORMATION_IMPL_H
#define MLMM_CPP_AVERAGE_INFORMATION_IMPL_H

#include "../util/util.h"
#include "../model/model.h"
#include "My_linear_algebra_impl.h"
#include "emmax_impl.h"


class Average_information_impl {
private:
    My_Vector<double> y;
    std::vector<My_matrix<double>> ks;
    My_matrix<double> x;
    My_matrix<double> xs;
    double eps;
    int maxiter;
    double _ll_();
    double _dll_( );
    double get_var_y();
    My_matrix<double> V;
    void update_V( const My_Vector<double> & etas );
public:
    Average_information_impl(const My_Vector<double> & _y, const std::vector<My_matrix<double>> & _ks,
                             const My_matrix<double> & _x, const My_matrix<double> & _xs, double & _eps,
                             const int & _maxiter);
    void get_ai_estimates();

};

#endif //MLMM_CPP_AVERAGE_INFORMATION_IMPL_H
