//
// Created by baoxing on 2/24/18.
//

#include "Qr_decomposition_result.h"
Qr_decomposition_result::Qr_decomposition_result( const My_matrix<double> & _r, const My_matrix<double> & _q){
    this->r=My_matrix<double>(_r);
    this->q=My_matrix<double>(_q);
}
const My_matrix<double> &Qr_decomposition_result::get_r() const{
    return r;
}
const My_matrix<double> &Qr_decomposition_result::get_q() const{
    return q;
}
