//
// Created by baoxing on 2/21/18.
//

#ifndef MLMM_CPP_MY_LINEAR_ALGEBRA_IMPL_H
#define MLMM_CPP_MY_LINEAR_ALGEBRA_IMPL_H

#include "../model/model.h"
#include <vector>
#include <algorithm>
#include <set>
#include <vector>
template <typename T>
T sum_of_a_vector( const My_Vector<T> & a ){
    double dsum=0;
    for (int i=0; i<a.get_length(); i++){
        dsum+=a.get_array()[i];
    }
    return dsum;
}


template <typename T>
T sum_of_powed_vector(const My_Vector<T> & a, const int & ex){
    double ReSum = 0;
    for (int i=0; i<a.get_length(); i++){
        ReSum += pow(a.get_array()[i], ex);
    }
    return ReSum;
}


template <typename T1, typename T2, typename T3>
void sum_of_two_vectors( const My_Vector<T1> & a, const My_Vector<T2> & b, My_Vector<T3> & c ){
    assert(a.get_length() == b.get_length() && c.get_length()==a.get_length());
    for (int i=0; i<a.get_length(); i++){
        c.get_array()[i] = a.get_array()[i]+b.get_array()[i];
    }
}

template <typename T1, typename T2, typename T3>
void sum_of_a_vector_a_number(const My_Vector<T1> & a, const T2 & z, My_Vector<T3> & c ){
    assert( c.get_length()==a.get_length());
    for (int i=0; i<a.get_length(); i++){
        c.get_array()[i]=a.get_array()[i]+z;
    }
}

template <typename T>
void production_of_two_vectors( const My_Vector<T> & a, const My_Vector<T> & b, My_Vector<T> & c ){
    assert(a.get_length() == b.get_length()&& c.get_length()==a.get_length());
    for (int i=0; i<a.get_length(); i++){
        c.get_array()[i] = a.get_array()[i] * b.get_array()[i];
    }
}

template <typename T1, typename T2, typename T3>
void quotient_of_two_vectors(const My_Vector<T1> & a, const My_Vector<T2> & b,  My_Vector<T2> & c ){
    assert(a.get_length() == b.get_length()&& c.get_length()==a.get_length());
    for (int i=0; i<a.get_length(); i++){
        c.get_array()[i] = a.get_array()[i] / b.get_array()[i];
    }
}

template <typename T>
void append_two_matrix(const My_matrix<T> & a, const My_matrix<T> & b, My_matrix<T> & c){
    assert( a.get_num_row() == b.get_num_row() && a.get_num_row()==c.get_num_row() && a.get_num_column()+b.get_num_column()==c.get_num_column());
    int i, j;
    for( i=0; i<a.get_num_row(); ++i ){
        for( j=0; j<a.get_num_column(); ++j ){
            c.get_matrix()[i][j] =a.get_matrix()[i][j];
        }
        for ( j=0; j<b.get_num_column(); ++j ){
            c.get_matrix()[i][a.get_num_column()+j] =b.get_matrix()[i][j];
        }
    }
}

template <typename T>
void trmul(const My_matrix<T> & a, const My_matrix<T> & b, My_matrix<T> & c){
    assert( a.get_num_column() == b.get_num_row() && a.get_num_row()==c.get_num_row() && b.get_num_column()==c.get_num_column());
    int i, j, l;
    for (i = 0; i < a.get_num_row(); ++i) {
        for (j = 0; j < b.get_num_column(); ++j) {
            c.get_matrix()[i][j] = 0.0;
            for (l = 0; l < a.get_num_column(); ++l) {
                c.get_matrix()[i][j] = c.get_matrix()[i][j] + a.get_matrix()[i][l] * b.get_matrix()[l][j];
            }
        }
    }
}

template <typename T>
void T_matrix( const My_matrix<T> & a, My_matrix<T> & b){
    assert(a.get_num_row()==b.get_num_column() && a.get_num_column()==b.get_num_row());
    int i, j;
    for( i=0; i<a.get_num_row(); ++i ){
        for ( j=0; j<a.get_num_column(); ++j ){
            b.get_matrix()[j][i]=a.get_matrix()[i][j];
        }
    }
}

void inverse_matrix(My_matrix<double> & c);
Qr_decomposition_result qr_decomposition( const My_matrix<double> & _a);
void lsq(  const My_matrix<double> & x,  const My_Vector<double> & y, My_Vector<double> & b );
void strq(const My_matrix<double> & a, double ** q, double *b, double * c);
int sstq(int & n, double * b, double *c, double **q, double eps, int l);
Eigen_result eigen( const My_matrix<double> & a, const double & eps, const int & l, const int & keep);
Eigen_result eigen( const My_matrix<double> & _a, const int & keep);
Eigen_result eigen( const My_matrix<double> & _a, const double & eps);
Eigen_result eigen( const My_matrix<double> & _a);
#endif //MLMM_CPP_MY_LINEAR_ALGEBRA_IMPL_H
// here is a good document for PCA analysis
// http://www.sthda.com/english/wiki/print.php?id=206
