//
// Created by baoxing on 2/27/18.
//

#ifndef MLMM_CPP_EMMA_RESULT_H
#define MLMM_CPP_EMMA_RESULT_H
#include "./My_Vector.h"
#include "./My_matrix.h"

class Emma_result {
    private:
        My_Vector<double> Y_t;
        My_matrix<double> H_sqrt_inv;
        double mahalanobis_rss;
        double pvalue;
    public:
        Emma_result(const My_Vector<double> & _Y_t, const My_matrix<double> & _H_sqrt_inv, const double & mahalanobis_rss);
        Emma_result(const My_Vector<double> & _Y_t, const My_matrix<double> & _H_sqrt_inv, const double & mahalanobis_rss, const double & _pvalue);
        const My_Vector<double> &getY_t() const;
        const My_matrix<double> &getH_sqrt_inv() const;
        double getMahalanobis_rss() const;
        double getPvalue() const;
};


#endif //MLMM_CPP_EMMA_RESULT_H
