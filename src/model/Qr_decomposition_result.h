//
// Created by baoxing on 2/24/18.
//

#ifndef MLMM_CPP_QR_DECOMPOSITION_RESULT_H
#define MLMM_CPP_QR_DECOMPOSITION_RESULT_H

#include "./My_matrix.h"

class Qr_decomposition_result {
    private:
        My_matrix<double> r;
        My_matrix<double> q; //each column is a eigen vector
    public:
        Qr_decomposition_result( const My_matrix<double> & _r, const My_matrix<double> & _q);
        const My_matrix<double> &get_r() const;
        const My_matrix<double> &get_q() const;
};


#endif //MLMM_CPP_QR_DECOMPOSITION_RESULT_H
