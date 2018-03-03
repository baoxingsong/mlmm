//
// Created by baoxing on 2/27/18.
//

#include "Emma_result.h"

Emma_result::Emma_result(const My_Vector<double> & _Y_t, const My_matrix<double> & _H_sqrt_inv, const double & mahalanobis_rss){
    this->Y_t = My_Vector<double>(_Y_t);
    this->H_sqrt_inv = My_matrix<double> (_H_sqrt_inv);
    this->mahalanobis_rss=double(mahalanobis_rss);
}

Emma_result::Emma_result(const My_Vector<double> & _Y_t, const My_matrix<double> & _H_sqrt_inv, const double & mahalanobis_rss, const double & _pvalue){
    this->Y_t = My_Vector<double>(_Y_t);
    this->H_sqrt_inv = My_matrix<double> (_H_sqrt_inv);
    this->mahalanobis_rss=double(mahalanobis_rss);
    this->pvalue=double(_pvalue);
}

const My_Vector<double> &Emma_result::getY_t() const {
    return this->Y_t;
}

const My_matrix<double> &Emma_result::getH_sqrt_inv() const {
    return this->H_sqrt_inv;
}

double Emma_result::getMahalanobis_rss() const {
    return this->mahalanobis_rss;
}

double Emma_result::getPvalue() const{
    return this->pvalue;
}


