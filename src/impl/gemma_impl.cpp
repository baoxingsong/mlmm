// =====================================================================================
//
//       Filename:  gemma_impl.cpp
//
//    Description:
//
//        Version:  1.0
//        Created:  05/20/2018 05:53:23 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//        Company:  unemployed
//
// =====================================================================================
/*************************************************************************
 *
 *
 *  This is a re impelmentation of the GEMMA method
 *  Xiang Zhou and Matthew Stephens (2012).
 *  Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821–824.
 *
 *  I ever read the source code of GEMMA, but I did copy any piece of code
 *
************************************************************************/


//
// Created by song on 5/20/18.
// This file implemented the algorithm in the pater ogf gemma published on nature genetics
//

#include "gemma_impl.h"

// the euqations have nothing with lambda value begin
//void _get_UtW_( const Eigen_result & eigen_L, const My_matrix<double> & W, const int & n, const int & q, My_matrix<double> & UtW){
//    int i, j, k;
//    for ( i=0; i<n; ++i ) {
//        for (j = 0; j < q; ++j) {
//            UtW.get_matrix()[i][j] = 0.0;
//        }
//    }
//    for ( i=0; i<n; ++i ) {
//        for (j = 0; j < q; ++j) {
//            for ( k=0; k<n; ++k ){
//                UtW.get_matrix()[i][j] += (eigen_L.get_eigen_vectors().get_matrix()[k][i] * W.get_matrix()[k][j]);
//            }
//        }
//    }
//}
//
//void _get_Uty_( const Eigen_result & eigen_L, const My_Vector<double> & y, const size_t & n, My_Vector<double> & Uty){
//    // this function has been compared with gemma, it should be correct
//    int i, j;
//    for ( i=0; i<n; ++i ){
//        Uty.get_array()[i]=0.0;
//        for ( j=0; j<n; ++j ){
//            Uty.get_array()[i]+=eigen_L.get_eigen_vectors().get_matrix()[j][i]*y.get_array()[j];
//        }
//    }
//}

// the equations have nothing with lambda value end
void get_lambda_theta_p_1_( const double & lambda, const Eigen_result & eigen_L, My_Vector<double> & lambda_theta_p_1, const int & n){
    for ( int i=0; i<n; ++i ){
        lambda_theta_p_1.get_array()[i]=(lambda*eigen_L.get_eigen_values().get_array()[i]+1.0);
    }
}

// the equation on the second page of online method
void _get_yTH_1_y_s_(const My_Vector<double> & Uty, const int & n,
                     My_Vector<double> & lambda_theta_p_1, double & yTH_1_y, double & yTH_1H_1H_1_y ){
    yTH_1_y=0.0;
    yTH_1H_1H_1_y=0.0;
    double temp1;
    for( int i=0; i<n; ++i ){
        temp1 = Uty.get_array()[i]*Uty.get_array()[i];
        yTH_1_y += (temp1/lambda_theta_p_1.get_array()[i]);
        yTH_1H_1H_1_y += (temp1/pow(lambda_theta_p_1.get_array()[i], 3.0));
    }
}

// _get_yTPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxy, ypw, wpw )
void _get_yTPxy_(const int & n, const My_Vector<double> & y, const My_matrix<double> & UtW, const My_Vector<double> & Uty,
                 const Eigen_result & eigen_L, const double & lambda, const int & q, const My_Vector<double> & lambda_theta_p_1,
                 double & yTPxy, My_matrix<double> & ypw, double *** wpw ) {
    yTPxy=0.0;
    int i, j, k, l;
    //for aTPiPib begin
    // for the wpw term in the function begin

    k=-1;// put k=0 part out of for loop is good for auto vectorization
            // here i is the index of a, j is the index of b and k is the index of P
    for( i=0; i<q; ++i ){ // since the index could not be -1, so here we use k+1 as the third index, attention, this only affect the third index
        for( j=0; j<q; ++j ){
            wpw[i][j][k+1] = 0.0;
            for( l=0; l<n; ++l ) {
                wpw[i][j][k+1] += UtW.get_matrix()[l][i]*UtW.get_matrix()[l][j]/lambda_theta_p_1.get_array()[l]; // todo check, there maybe something wrong with the indeX here
            }
        }
    }
    for( k=0; k<q; ++k ){
        for( i=0; i<q; ++i ){
            for( j=0; j<q; ++j ){
                wpw[i][j][k+1] = wpw[i][j][k+1-1] - wpw[i][k][k+1-1]*wpw[j][k][k+1-1]/wpw[k][k][k+1-1];
            }
        }
    }
    // for the wpw term in the function end

    // this here we have a==b==y,
    // for the ypw term in the function begin
    k=-1; // j is the index of b and k is the index of P
    for( j=0; j<q; ++j ){
        ypw.get_matrix()[j][k+1] = 0.0;
        for( l=0; l<n; ++l ) {
            ypw.get_matrix()[j][k+1] += Uty.get_array()[l]*UtW.get_matrix()[l][j]/lambda_theta_p_1.get_array()[l];
        }
    }
    for( k=0; k < q; ++k ){
        for( j=0; j<q; ++j ){
            ypw.get_matrix()[j][k+1] = ypw.get_matrix()[j][k+1-1] - ypw.get_matrix()[k][k+1-1]*wpw[j][k][k+1-1]/wpw[k][k][k+1-1];
        }
    }
    // for the ypw term in the function end
    k=-1;
    for( l=0; l<n; ++l ) {
        yTPxy += Uty.get_array()[l] *  Uty.get_array()[l] / lambda_theta_p_1.get_array()[l];
    }
    for(  k=0; k < q; ++k ){
        yTPxy = yTPxy - ypw.get_matrix()[k][k+1-1]*ypw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1];
    }
}

void _get_yTPxPxy_(const int & n, const My_Vector<double> & y, const My_matrix<double> & UtW, const My_Vector<double> & Uty,
                 const Eigen_result & eigen_L, const double & lambda, const int & q, const My_Vector<double> & lambda_theta_p_1,
                 double & yTPxPxy, const My_matrix<double> & ypw, My_matrix<double> & yppw, double *** wpw, double *** wppw ) {
    int i, j, k, l;
    yTPxPxy=0;
    //for aTPiPib begin
    // for the wppw term in the function begin

    k=-1;// put k=0 part out of for loop is good for auto vectorization
    // here i is the index of a, j is the index of b and k is the index of P
    for( i=0; i<q; ++i ){
        for( j=0; j<q; ++j ){
            wppw[i][j][k+1] = 0.0;
            for( l=0; l<n; ++l ) {
                wppw[i][j][k+1] += UtW.get_matrix()[l][i]*UtW.get_matrix()[l][j]/pow(lambda_theta_p_1.get_array()[l], 2);
            }
        }
    }

    for( k=0; k < q; ++k ){
        for( i=0; i<q; ++i ){
            for( j=0; j<q; ++j ){
                wppw[i][j][k+1] = wppw[i][j][k+1-1] + wpw[i][k][k+1-1]*wpw[j][k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1],2)-
                                  wpw[i][k][k+1-1]*wppw[j][k][k+1-1]/wpw[k][k][k+1-1]-
                                  wpw[j][k][k+1-1]*wppw[i][k][k+1-1]/wpw[k][k][k+1-1];
            }
        }
    }
    // for the wppw term in the function end

    // this here we have a==b==y
    // for the yppw term in the function begin
    k=-1; // j is the index of b and k is the index of P
    for( j=0; j<q; ++j ){
        yppw.get_matrix()[j][k+1] = 0.0;
        for( l=0; l<n; ++l ) {
            yppw.get_matrix()[j][k+1] += Uty.get_array()[l]*UtW.get_matrix()[l][j]/pow(lambda_theta_p_1.get_array()[l], 2);
        }
    }
    for( k=0; k < q; ++k ){
        for( j=0; j<q; ++j ){
            yppw.get_matrix()[j][k+1] = yppw.get_matrix()[j][k+1-1] + ypw.get_matrix()[k][k+1-1]*wpw[j][k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1],2)-
                                        ypw.get_matrix()[k][k+1-1]*wppw[j][k][k+1-1]/wpw[k][k][k+1-1]-
                                        wpw[j][k][k+1-1]*yppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1];
        }
    }
    // for the yppw term in the function end
    k=-1;
    for( j=0; j<n; ++j ) {
        yTPxPxy += Uty.get_array()[j] *  Uty.get_array()[j] / pow(lambda_theta_p_1.get_array()[j], 2);
    }
    for(  k=0; k < q; ++k ){
        yTPxPxy = yTPxPxy + ypw.get_matrix()[k][k+1-1]*ypw.get_matrix()[k][k+1-1]*wppw[k][k][k+1-1]/pow((wpw[k][k][k+1-1]), 2) -
                  2*(ypw.get_matrix()[k][k+1-1]*yppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1]);
    }
    //for aTPiPib end
}

void _get_yTPxPxPxy_(const int & n, const My_Vector<double> & y, const My_matrix<double> & UtW,
                 const My_Vector<double> & Uty, const Eigen_result & eigen_L, const double & lambda, const int & q, const My_Vector<double> & lambda_theta_p_1,
                 double & yTPxPxPxy, const My_matrix<double> & ypw, const My_matrix<double> & yppw,
                 My_matrix<double> & ypppw, double *** wpw, double *** wppw, double *** wpppw) {
    yTPxPxPxy=0.0;
    int i, j, k, l;

    //for aTPiPiPib begin
    // for the wpppw term in the function begin

    k=-1;// put k=0 part out of for loop is good for auto vectorization
    // here i is the index of a, j is the index of b and k is the index of P
    for( i=0; i<q; ++i ){
        for( j=0; j<q; ++j ){
            wpppw[i][j][k+1] = 0.0;
            for( l=0; l<n; ++l ) {
                wpppw[i][j][k+1] += UtW.get_matrix()[l][i]*UtW.get_matrix()[l][j]/pow(lambda_theta_p_1.get_array()[l], 3);
            }
        }
    }
    for( k=0; k < q; ++k ){
        for( i=0; i<q; ++i ){
            for( j=0; j<q; ++j ){
                wpppw[i][j][k+1] = wpppw[i][j][k+1-1] - wpw[i][k][k+1-1]*wpw[j][k][k+1-1]*pow(wppw[k][k][k+1-1],2)/pow(wpw[k][k][k+1-1],3)-
                                   wpw[i][k][k+1-1]*wpppw[j][k][k+1-1]/wpw[k][k][k+1-1]-
                                   wpw[j][k][k+1-1]*wpppw[i][k][k+1-1]/wpw[k][k][k+1-1]-
                                   wppw[i][k][k+1-1]*wppw[j][k][k+1-1]/wpw[k][k][k+1-1]+
                                   wpw[i][k][k+1-1]*wppw[j][k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                                   wpw[j][k][k+1-1]*wppw[i][k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                                   wpw[i][k][k+1-1]*wpw[j][k][k+1-1]*wpppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2);
            }
        }
    }
    // for the wpppw term in the function end

    // this here we have a==b==y
    // for the ypppw term in the function begin
    k=-1; // j is the index of b and k is the index of P
    for( j=0; j<q; ++j ){
        ypppw.get_matrix()[j][k+1] = 0.0;
        for( l=0; l<n; ++l ) {
            ypppw.get_matrix()[j][k+1] += Uty.get_array()[l]*UtW.get_matrix()[l][j]/pow(lambda_theta_p_1.get_array()[l], 3);
        }
    }
    for( k=0; k < q; ++k ){
        for( j=0; j<q; ++j ){
            ypppw.get_matrix()[j][k+1] = ypppw.get_matrix()[j][k+1-1] - ypw.get_matrix()[k][k+1-1]*wpw[j][k][k+1-1]*pow(wppw[k][k][k+1-1], 2)/pow(wpw[k][k][k+1-1], 3)-
                                         ypw.get_matrix()[k][k+1-1]*wpppw[j][k][k+1-1]/wpw[k][k][k+1-1]-
                                         wpw[j][k][k+1-1]*ypppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1]-
                                         yppw.get_matrix()[k][k+1-1]*wppw[j][k][k+1-1]/wpw[k][k][k+1-1]+
                                         ypw.get_matrix()[k][k+1-1]*wppw[j][k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                                         wpw[j][k][k+1-1]*yppw.get_matrix()[k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                                         ypw.get_matrix()[k][k+1-1]*wpw[j][k][k+1-1]*wpppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2);
        }
    }
    // for the ypppw term in the function end
    k=-1;
    for( j=0; j<n; ++j ) {
        yTPxPxPxy += Uty.get_array()[j] *  Uty.get_array()[j] / pow(lambda_theta_p_1.get_array()[j], 3);
    }
    for(  k=0; k < q; ++k ){
        yTPxPxPxy = yTPxPxPxy - pow(ypw.get_matrix()[k][k+1-1], 2) * pow(wppw[k][k][k+1-1], 2)/pow(wpw[k][k][k+1-1], 3)-
                    ypw.get_matrix()[k][k+1-1]*ypppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1]-
                    ypw.get_matrix()[k][k+1-1]*ypppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1]-
                    yppw.get_matrix()[k][k+1-1]*yppw.get_matrix()[k][k+1-1]/wpw[k][k][k+1-1]+
                    ypw.get_matrix()[k][k+1-1]*yppw.get_matrix()[k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                    ypw.get_matrix()[k][k+1-1]*yppw.get_matrix()[k][k+1-1]*wppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2)+
                    ypw.get_matrix()[k][k+1-1]*ypw.get_matrix()[k][k+1-1]*wpppw[k][k][k+1-1]/pow(wpw[k][k][k+1-1], 2);
    }
    //for aTPiPib end
}

// the euqation on the second page of online method
double _get_log_det_H_(const My_Vector<double> & lambda_theta_p_1, const int & n){
    double result = 0.0;
    for( int i=0; i<n; ++i ){
        if( lambda_theta_p_1.get_array()[i]<0 ){
            result += log(-lambda_theta_p_1.get_array()[i]); //try to avoid to use if for avx
        }else if( lambda_theta_p_1.get_array()[i]>0 ) {
            result += log(lambda_theta_p_1.get_array()[i]); //try to avoid to use if for avx
        }
    }
    return result;
}

// the eqaution on the second page of online method
void _get_trace_H_1_s_(const My_Vector<double> & lambda_theta_p_1, const int & n, double & trace_H_1, double & trace_H_1_H_1){
    trace_H_1 = 0.0;
    trace_H_1_H_1 = 0.0;
    for( int i=0; i<n; ++i ){
        trace_H_1 += 1.0/lambda_theta_p_1.get_array()[i];
        trace_H_1_H_1 += 1.0/pow(lambda_theta_p_1.get_array()[i], 2);
    }
}

// the first equation of 3.1.4 on page 5 of supplementary document
// this one rely on the result from _get_trace_H_1_
double _get_trace_H_1_G_( const double & lambda, const int & n, const double & trace_H_1){
    return (n-trace_H_1)/lambda;
}

// the second equation of 3.1.4 on page 5 of supplementary document
// this one rely on the result from _get_trace_H_1_ and _get_trace_H_1_H_1_
double _get_trace_H_1_G_H_1_G_(const double & lambda, const int & n, const double & trace_H_1, const double & trace_H_1_H_1){
    return (n+trace_H_1_H_1-2.0*trace_H_1)/pow(lambda, 2);
}

void _get_trace_Px_s_( const double & trace_H_1, const double & trace_H_1_H_1, const My_Vector<double> & lambda_theta_p_1, const My_matrix<double> & UtW,
                       const int & n, const int & q, double *** wpw, double *** wppw, double *** wpppw, double & trace_P, double & trace_PP){
    int k=0;
    trace_P  = trace_H_1 -     wppw[k][k][k+1-1]/wpw[k][k][k+1-1];
    trace_PP = trace_H_1_H_1 + pow(wppw[k][k][k+1-1], 2)/pow(wpw[k][k][k+1-1], 2)-2*wpppw[k][k][k+1-1]/wpw[k][k][k+1-1];
    for(  k=1; k < q; ++k ){
        trace_P  = trace_P  -  wppw[k][k][k+1-1]/wpw[k][k][k+1-1];
        trace_PP = trace_PP +  pow(wppw[k][k][k+1-1], 2)/pow(wpw[k][k][k+1-1], 2)-2*wpppw[k][k][k+1-1]/wpw[k][k][k+1-1];
    }
}

// the third equation of 3.1.4 on page 5 of supplementary document
double _get_trace_Px_G_(const double & lambda, const int & p, const double & trace_Px){
    return (p-trace_Px)/lambda;
}

// the fourth equation of 3.1.4 on page 5 of supplementary document
double _get_trace_Px_G_Px_G_(const double & lambda, const int & p, const double & trace_Px, const double & trace_Px_Px){
    return (p+trace_Px_Px-2*trace_Px)/pow(lambda,2);
}

double get_det_WTH_1W_( const int & n, const int & q, const My_matrix<double> & UtW, double *** wpw,
                        const My_Vector<double> & lambda_theta_p_1){ // todo check, there maybe something wrong
    double WTH_1W = 0.0;
    int i, j = 0;
    for( i=0; i<n; ++i ){
        WTH_1W += UtW.get_matrix()[i][0]*UtW.get_matrix()[i][0]/lambda_theta_p_1.get_array()[i];
    }
    WTH_1W = WTH_1W * wpw[j][j][j+1-1];
    for( j=1; j<q; ++j ){
        WTH_1W = WTH_1W * wpw[j][j][j+1-1];
    }
    return WTH_1W;
}

// the fifth equation of 3.1.4 on page 5 of supplementary document
// rely one the function get_ytPxy and get_ytPxPxy
double get_ytPxGPxy( const double & ytPxy, const double & ytPxPxy, const double & lambda){
    return (ytPxy-ytPxPxy)/lambda;
}
// the sixth equation of 3.1.4 on page 5 of supplementary document
// rely one the function get_ytPxy, get_ytPxPxy and get_ytPxPxPxy
double get_yTPxGPxGPxy( const double & ytPxy, const double & ytPxPxPxy, const double & ytPxPxy, const double & lambda ){
    return (ytPxy+ytPxPxPxy-2*ytPxPxy)/pow(lambda, 2);
}

// simplified log-likelihoods (eq. 3 from gemma paper)
double _sll_( const int & n, const double & log_det_H, const double & ytPxy){
    return 0.5*n*log((double)n/(2*M_PI)) - 0.5*n - 0.5*log_det_H-0.5*n*log(ytPxy);
}

// simplified log-restricted likelihoods (eq. 4 from gemma paper) // could not be used, since the function to calculate det_WTW is not implemented
double _srll_( const int & p, const double & log_det_H, const double & det_WTW, const double & det_WTH_1_W, const double & ytPxy){
    return 0.5*(p)*log((double)p/(2*M_PI))-0.5*p+0.5*log(det_WTW)-0.5*log_det_H-0.5*log(det_WTH_1_W) - 0.5*p*log(ytPxy);
}

// simplified log-restricted likelihoods (eq. 4 from gemma paper) // this is a simplified function, since det_WTH_1_W is not affected by lambda value. It could be only used for newton optimization
double _srll_simple_( const int & p, const double & det_h, const double & det_WTH_1_W, const double & ytPxy){
    return -0.5*log(det_h)-0.5*log(det_WTH_1_W) - 0.5*p*log(ytPxy);
}

// (eq. 5 from gemma paper)
double _dsll_(const double & trace_H_1_G, const int & n, const double & ytPxGPxy, const double & ytPxy){
    return -0.5*trace_H_1_G+0.5*n*ytPxGPxy/ytPxy;
}

// (eq. 7 from gemma paper)
double _dsrll_(const double & trace_Px_G, const int & p, const double & ytPxGPxy, const double & ytPxy){
    return -0.5/trace_Px_G + 0.5*p*ytPxGPxy/ytPxy;
}

// (eq. 6 from gemma paper)
double _ddsll_( const double & trace_H_1_G_H_1_G, const int & n, const double & yTPxGPxGPxy, const double & ytPxy, const double & ytPxGPxy ){
    return 0.5*trace_H_1_G_H_1_G-0.5*n*(2*yTPxGPxGPxy*ytPxy-pow(ytPxGPxy,2))/pow(ytPxy, 2);
}

// (eq. 8 from gemma paper)
double _ddsrll_(const double & trace_Px_G_Px_G, const int & p, const double & yTPxGPxGPxy, const double & ytPxy, const double & ytPxGPxy){
    return 0.5*trace_Px_G_Px_G - 0.5*p*(2*yTPxGPxGPxy*ytPxy-pow(ytPxGPxy, 2)) / pow(ytPxy,2);
}

// this is an implementation of the euqation 2.6 from publication:
//Generalized Newton Raphson’s method free from second derivative, Journal of Nonlinear Science and Applications 9(5):2823-2831 · April 2016 DOI: 10.22436/jnsa.009.05.77
// however this is not the original publication, so do not cite this paper for publication purpose
int newton_raphson_ll( double & new_lambda, const double & eps, const int & js, const int & n, const int & q, const Eigen_result & eigen_L,
                       const My_Vector<double> & y, const My_matrix<double> & UtW, const My_Vector<double> & Uty,
                       double *** wpw, double *** wppw, double *** wpppw, My_matrix<double> & ypw, My_matrix<double> & yppw, My_matrix<double> & ypppw ){ // js maximum number of iteration
    int k=1, l=js;
    double y0, y1, y2, d, u, x0=new_lambda, x1=new_lambda;

    My_Vector<double> lambda_theta_p_1(n);
    double trace_H_1;
    double trace_H_1_H_1;
    get_lambda_theta_p_1_( x1, eigen_L, lambda_theta_p_1, n);
    double log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
    double yTH_1_y;
    double yTH_1H_1H_1_y;

    _get_yTH_1_y_s_(Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );
    double yTPxy;
    double yTPxPxy;
    double yTPxPxPxy;
    _get_yTPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxy, ypw, wpw ); //there maybe something wrong with this function. Because yTPxy should not be negative
    //std::cout << "line 404 _get_yTPxy_ " << yTPxy << std::endl;

    _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
    _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );

    _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);
    double trace_H_1_G = _get_trace_H_1_G_( x1, n, trace_H_1);
    double trace_H_1_G_H_1_G = _get_trace_H_1_G_H_1_G_(x1, n, trace_H_1, trace_H_1_H_1);
    double ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxPxy, x1);
    double yTPxGPxGPxy = get_yTPxGPxGPxy(yTPxy, yTPxPxPxy, yTPxPxy, x1);
    y0=_sll_(n, log_det_H, yTPxy);
    y1=_dsll_(trace_H_1_G, n, ytPxGPxy, yTPxy);
    y2=_ddsll_(trace_H_1_G_H_1_G, n, yTPxGPxGPxy, yTPxy, ytPxGPxy);
    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
        //std::cout << "eps " << eps << " " << l << " " << x1 << " " << x0 << " " << y0 << " " << y1 << " " << y2 << std::endl;
        if( fabs(y1)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        //x1 = x0 - (y1-sqrt(pow(y1, 2) - 2*y0*y2))/y2;
        //x1 = x0 - (exp(y1)-sqrt(pow(exp(y1), 2) - 2*exp(y0)*exp(y2)))/exp(y2);
        x1 = x0 - y1/y2;
        //x1 = x0 - exp(y0-y1);
//        x1 = x0 - (2*y0*y1)/(2*y1*y1 -y1*y2);
//        x1 = x0 - (y1-sqrt(pow(y1, 2) - 2*y0*y2))/y2;
        get_lambda_theta_p_1_( x1, eigen_L, lambda_theta_p_1, n);
        log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
        _get_yTH_1_y_s_(Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );
        _get_yTPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxy, ypw, wpw );
        _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
        _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );
        _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);

        trace_H_1_G = _get_trace_H_1_G_( x1, n, trace_H_1);
        trace_H_1_G_H_1_G = _get_trace_H_1_G_H_1_G_(x1, n, trace_H_1, trace_H_1_H_1);
        ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxPxy, x1);
        yTPxGPxGPxy = get_yTPxGPxGPxy(yTPxy, yTPxPxPxy, yTPxPxy, x1);

        y0=_sll_(n, log_det_H, yTPxy);
        //std::cout << "n " << n << " log_det_H " << log_det_H << " yTPxy " << yTPxy  << std::endl;
        y1=_dsll_(trace_H_1_G, n, ytPxGPxy, yTPxy);
        y2=_ddsll_(trace_H_1_G_H_1_G, n, yTPxGPxGPxy, yTPxy, ytPxGPxy);
        //std::cout << "eps " << eps << " " << l << " " << x1 << " " << x0 << " " << y0 << " " << y1 << " " << y2 << std::endl;
        //std::cout << "d " << d << " l " << l << std::endl;

        d=fabs(x1-x0);
//        u=fabs(y0);
//        if( u>d ){
//            d=u;
//        }
        if( x1 < 10e-5 ){
            x1 = 10e-5;
        }
        if( x1 > 10e5 ){
            x1 = 10e5;
        }
        x0=x1;
        --l;
        //std::cout << "d " << d << " l " << l << std::endl;
    }
    new_lambda=x1;
    k=js-l;
    return k;
}

int newton_raphson_reml( double & new_lambda, const double & eps, const int & js, const int & n, const Eigen_result & eigen_L,
                         const My_Vector<double> & y,
                         const int & q, const int & p, double *** wpw, double *** wppw, double *** wpppw, const My_matrix<double> & UtW,
                         const My_Vector<double> & Uty, My_matrix<double> & ypw, My_matrix<double> & yppw, My_matrix<double> & ypppw){ // js maximum number of iteration
    int k=1,l=js;
    double y0, y1, y2, d,  x0, x1=new_lambda, u;
    x0=new_lambda;

    My_Vector<double> lambda_theta_p_1(n);
    double trace_H_1;
    double trace_H_1_H_1;
    get_lambda_theta_p_1_( x1, eigen_L, lambda_theta_p_1, n);
    double log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);

    double yTH_1_y;
    double yTH_1H_1H_1_y;
    _get_yTH_1_y_s_(Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );

    double yTPxy;
    double yTPxPxy;
    double yTPxPxPxy;
    _get_yTPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxy, ypw, wpw );
    _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
    _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );

    _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);

    double ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxPxy, x1);
    double yTPxGPxGPxy = get_yTPxGPxGPxy(yTPxy, yTPxPxPxy, yTPxPxy, x1);
    double trace_P;
    double trace_PP;
    _get_trace_Px_s_( trace_H_1, trace_H_1_H_1, lambda_theta_p_1, UtW, n, q, wpw, wppw, wpppw, trace_P, trace_PP);


    double trace_Px_G = _get_trace_Px_G_(x1, p, trace_P);
    double trace_Px_G_Px_G = _get_trace_Px_G_Px_G_(x1, p, trace_P, trace_PP);

    double det_WTH_1_W = get_det_WTH_1W_( n, q, UtW, wpw, lambda_theta_p_1 );

    y0=_srll_simple_(p, log_det_H, det_WTH_1_W, yTPxy);
    y1=_dsrll_(trace_Px_G, p, ytPxGPxy, yTPxy);
    y2=_ddsrll_(trace_Px_G_Px_G, p, yTPxGPxGPxy, yTPxy, ytPxGPxy);
    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
        if( fabs(y2)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        x1 = x0 - y1/y2;

        get_lambda_theta_p_1_( x1, eigen_L, lambda_theta_p_1, n);
        log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
        _get_yTH_1_y_s_(Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );

        _get_yTPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxy, ypw, wpw );
        _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
        _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x1, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );
        _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);

        ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxy, x1);
        yTPxGPxGPxy = get_yTPxGPxGPxy(yTPxy, yTPxPxPxy, yTPxPxy, x1);
        _get_trace_Px_s_( trace_H_1, trace_H_1_H_1, lambda_theta_p_1, UtW, n, q, wpw, wppw, wpppw, trace_P, trace_PP);
        trace_Px_G = _get_trace_Px_G_(x1, p, trace_P);
        trace_Px_G_Px_G = _get_trace_Px_G_Px_G_(x1, p, trace_P, trace_PP);
        det_WTH_1_W = get_det_WTH_1W_( n, q, UtW, wpw, lambda_theta_p_1 );
        y0=_srll_simple_(p, log_det_H, det_WTH_1_W, yTPxy);
        y1=_dsrll_(trace_Px_G, p, ytPxGPxy, yTPxy);
        y2=_ddsrll_(trace_Px_G_Px_G, p, yTPxGPxGPxy, yTPxy, ytPxGPxy);
        d=fabs(x1-x0);
        u=fabs(y0);
        if( u>d ){
            d=u;
        }
        x0=x1;
        --l;
    }
    new_lambda=x1;
    k=js-l;
    return k;
}

double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty, const My_matrix<double> & x,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                         const std::string & method, const int & maxiter ){
    double lambda;
    return gemma_estimates ( y, Uty, x, w, eigen_L, ngrids, llim, ulim, eps, method, maxiter, lambda );
}

double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty,
                         const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                         const std::string & method, const int & maxiter ){
    My_matrix<double> x(y.get_length(), 0);
    double lambda;
    return gemma_estimates ( y, Uty, x, w, eigen_L, ngrids, llim, ulim, eps, method, maxiter, lambda );

}

double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty,
                  const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                         const double & llim, const double &ulim, const double & eps,
                  const std::string & method, const int & maxiter, double & lambda ){
    My_matrix<double> x(y.get_length(), 0);
    return gemma_estimates ( y, Uty, x, w, eigen_L, ngrids, llim, ulim, eps, method, maxiter, lambda );
}

double gemma_estimates ( const My_Vector<double> & y, const My_Vector<double> & Uty, const My_matrix<double> & x,
                       const My_matrix<double> & w, const Eigen_result & eigen_L, const int & ngrids,
                       const double & llim, const double &ulim, const double & eps,
                            const std::string & method, const int & maxiter, double & lambda ){

    assert( y.get_length()==w.get_num_row() );

    int i, j;
    int n = y.get_length();
    int c = w.get_num_column();
    int cx = x.get_num_column();
    int q = c + cx;
    int p = n - q; // in the gemma paper the cx == 1, so in their function they have n-c-1 in the equation (4)
    int m = ngrids+1;

    My_matrix<double> W(x.get_num_row(), q);
    if( x.get_num_column()>0 ){
        append_two_matrix( w, x, W);
    } else {
        W.value_copy(w);
    }
    std::cout << "q " << q << std::endl;
    for( i=0; i<n; ++i ){
        std::cout << "y " << y.get_array()[i] << std::endl;
        std::cout << "W " << W.get_matrix()[i][0] << std::endl;
    }

    My_matrix<double> UtW(n, q);
    trmul( eigen_L.get_eigen_vectors(), W, UtW);

    My_Vector<double> log_lambdas(m); // the space for deltas to search
    My_Vector<double> lambdas(m);
    for ( i=0; i < m; ++i ){
        log_lambdas.get_array()[i] = (double(i) / ngrids)*(ulim - llim) + llim;
        lambdas.get_array()[i] = pow(10, log_lambdas.get_array()[i]);
    }
    double x0;
    double log_det_H;
    double trace_H_1;
    double trace_H_1_H_1;
    double trace_H_1_G;
    double yTH_1_y;
    double yTPxPxPxy;
    double yTH_1H_1H_1_y;
    double yTPxy;
    double yTPxPxy;
    double ytPxGPxy;
    double trace_P;
    double trace_PP;
    double trace_Px_G;
    double det_WTH_1_W;

    My_matrix<double> P(n, n);
    My_Vector<double> lls(m);
    My_Vector<double> dlls(m);

    //matrix for lambda estimation begin
    double *** wpw = new double** [q];
    for( i=0; i<q; ++i ){
        wpw[i]= new double* [q];
        for( j=0; j<q; ++j ){
            wpw[i][j] = new double[q+1];
        }
    }

    double *** wppw = new double** [q];
    for( i=0; i<q; ++i ){
        wppw[i]= new double* [q];
        for( j=0; j<q; ++j ){
            wppw[i][j] = new double[q+1];
        }
    }

    double *** wpppw = new double** [q];
    for( i=0; i<q; ++i ){
        wpppw[i]= new double* [q];
        for( j=0; j<q; ++j ){
            wpppw[i][j] = new double[q+1];
        }
    }
    My_matrix<double> ypw(q, q+1);
    My_matrix<double> yppw(q, q+1);
    My_matrix<double> ypppw(q, q+1);
    //matrix for lambda estimation end
    My_Vector<double> lambda_theta_p_1(n);

    if (method.compare("REML") == 0){
        for( i=0; i<m; ++i ){
            x0 = lambdas.get_array()[i];
            get_lambda_theta_p_1_( x0, eigen_L, lambda_theta_p_1, n);
            log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
            _get_yTH_1_y_s_( Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );
            _get_yTPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxy, ypw, wpw );
            _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
            _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );

            _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);
            ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxy, x0);
            _get_trace_Px_s_( trace_H_1, trace_H_1_H_1, lambda_theta_p_1, UtW, n, q, wpw, wppw, wpppw, trace_P, trace_PP);
            trace_Px_G = _get_trace_Px_G_(x0, p, trace_P);
            det_WTH_1_W = get_det_WTH_1W_( n, q, UtW, wpw, lambda_theta_p_1 );


            lls.get_array()[i] =_srll_simple_(p, log_det_H, det_WTH_1_W, yTPxy);
            dlls.get_array()[i]=_dsrll_(trace_Px_G, p, ytPxGPxy, yTPxy);
//            std::cout << "lls.get_array()[i] "<< lls.get_array()[i] << " dlls.get_array()[i] " << dlls.get_array()[i] <<std::endl;
//            printf("inter loop log_det_H:%10.20f\tyTPxy:%10.20f\tlls.get_array()[i]:%10.20f\tdlls.get_array()[i]:%10.20f\n", log_det_H, yTPxy, lls.get_array()[i], dlls.get_array()[i]);
        }
    }else if (method.compare("ML") == 0){
        for( i=0; i<m; ++i ){
            x0 = lambdas.get_array()[i];
            get_lambda_theta_p_1_( x0, eigen_L, lambda_theta_p_1, n);
            log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
            _get_yTH_1_y_s_(Uty, n, lambda_theta_p_1, yTH_1_y, yTH_1H_1H_1_y );
            _get_yTPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxy, ypw, wpw );
            _get_yTPxPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxPxy, ypw, yppw, wpw, wppw );
            _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);
            trace_H_1_G = _get_trace_H_1_G_( x0, n, trace_H_1);
            ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxPxy, x0);

            _get_yTPxPxPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxPxPxy, ypw, yppw, ypppw, wpw, wppw, wpppw );

            _get_trace_H_1_s_(lambda_theta_p_1, n, trace_H_1, trace_H_1_H_1);
            double trace_H_1_G = _get_trace_H_1_G_( x0, n, trace_H_1);
            double trace_H_1_G_H_1_G = _get_trace_H_1_G_H_1_G_(x0, n, trace_H_1, trace_H_1_H_1);
            double ytPxGPxy = get_ytPxGPxy( yTPxy, yTPxPxy, x0);
            double yTPxGPxGPxy = get_yTPxGPxGPxy(yTPxy, yTPxPxPxy, yTPxPxy, x0);

            lls.get_array()[i] = _sll_( n, log_det_H, yTPxy);
            dlls.get_array()[i] = _dsll_( trace_H_1_G, n, ytPxGPxy, yTPxy);
        }
    }
    int max_ll_i = argmax(lls);
    double max_ll = lls.get_array()[max_ll_i];
    double last_dll = dlls.get_array()[0];

    std::map<int, double> zero_intervals;
    for ( i =0; i<m; ++i ){
        if ( dlls.get_array()[i]<0 && last_dll>0 ){
            zero_intervals[i] = lls.get_array()[i];
        }
//        std::cout << i << " " << zero_intervals[i] << " " << last_dll << std::endl;
        last_dll = dlls.get_array()[i];
    }

    double opt_lambda;
    double opt_ll;

    if (zero_intervals.size() > 0){
        opt_ll = -DBL_MAX;
        int opt_i;
        for(std::map<int, double>::iterator it=zero_intervals.begin();it!=zero_intervals.end();++it){
            if( it->second > opt_ll ){
                opt_ll = it->second;
                opt_i = it->first;
            }
        }
//        if( (opt_i) <0 ){
            opt_lambda = lambdas.get_array()[opt_i];
//        }else{
//            opt_lambda = lambdas.get_array()[opt_i - 1];
//        }

        double new_opt_lambda = opt_lambda;
        if (method.compare("REML") == 0) {
            newton_raphson_reml( new_opt_lambda , eps, maxiter, n, eigen_L, y, q, p, wpw, wppw, wpppw, UtW, Uty, ypw, yppw, ypppw);
        }else if (method.compare("ML") == 0){
            newton_raphson_ll( new_opt_lambda, eps, maxiter, n, q, eigen_L, y, UtW, Uty, wpw, wppw, wpppw, ypw, yppw, ypppw);
        }

        if (opt_i > 1 && lambdas.get_array()[opt_i - 1] - eps < new_opt_lambda && new_opt_lambda < lambdas.get_array()[opt_i] + eps){
            opt_lambda = new_opt_lambda;
        }else if (opt_i == 1 && 0.0 < new_opt_lambda && new_opt_lambda < lambdas.get_array()[opt_i] + eps){
            opt_lambda = new_opt_lambda;
        }else if(opt_i == m - 1 && new_opt_lambda > lambdas.get_array()[opt_i - 1] - eps ){
            opt_lambda = new_opt_lambda;
        }

        get_lambda_theta_p_1_( opt_lambda, eigen_L, lambda_theta_p_1, n);
        log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
        _get_yTPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxy, ypw, wpw );

        opt_ll = _sll_( n, log_det_H, yTPxy);
        std::cout << "759 opt_lambda " << opt_lambda << std::endl;
        if (opt_ll < max_ll || std::isnan(opt_ll) ){
            opt_lambda = lambdas.get_array()[max_ll_i];
        }
        std::cout << "763 opt_lambda " << opt_lambda << std::endl;
    } else {
        std::cout << " 762 opt_lambda " << opt_lambda << std::endl;
        opt_lambda = lambdas.get_array()[max_ll_i];
        std::cout << " 764 opt_lambda " << opt_lambda << std::endl;
    }
    // no matter use ml or reml to optimize the parameters, always use the ml function to get the likelyhood

    get_lambda_theta_p_1_( opt_lambda, eigen_L, lambda_theta_p_1, n);
    log_det_H = _get_log_det_H_(lambda_theta_p_1,  n);
    _get_yTPxy_( n, y, UtW, Uty, eigen_L, x0, q, lambda_theta_p_1, yTPxy, ypw, wpw );

    opt_ll = _sll_( n, log_det_H, yTPxy);

    for ( i=0; i<q; ++i ){
        for ( j=0; j<q; ++j ){
            delete [] wpw[i][j];
            delete [] wppw[i][j];
            delete [] wpppw[i][j];
        }
        delete [] wpw[i];
        delete [] wppw[i];
        delete [] wpppw[i];
    }
    delete [] wpw;
    delete [] wppw;
    delete [] wpppw;

    lambda = opt_lambda;
    std::cout << "792 opt_lambda " << opt_lambda << std::endl;
    return opt_ll;
    // this function should return beta and likelihood for p-value computing
}

void CalcLmmVgVe(const My_Vector<double> & y, const Eigen_result & eigen_L, const My_matrix<double> &UtW,
                 const My_Vector<double> & Uty, const double & lambda, double &vg,
                     double &ve) {

    int q = UtW.get_num_column();
    int n = UtW.get_num_row();
    double yTPxy;
    My_matrix<double> ypw(q, q+1);

    int i, j;
    double *** wpw = new double** [q];
    for( i=0; i<q; ++i ){
        wpw[i]= new double* [q];
        for( j=0; j<q; ++j ){
            wpw[i][j] = new double[q+1];
        }
    }

    My_Vector<double> lambda_theta_p_1(n);

    get_lambda_theta_p_1_( lambda, eigen_L, lambda_theta_p_1, n);
    _get_yTPxy_( n, y, UtW, Uty, eigen_L, lambda, q, lambda_theta_p_1, yTPxy, ypw, wpw );

    ve = yTPxy / (double)(n - q);
    vg = ve * lambda;
}

/*
 * get_lambda_theta_p_1_
 * _get_yTH_1_y_s_
 *
 * get_lambda_theta_p_1_
 * _get_yTPxy_s_
 *
 * get_lambda_theta_p_1_
 * _get_det_H_
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * double _get_trace_H_1_G_
 *
 * get_lambda_theta_p_1_
 * _get_yTPxy_s_
 * double get_det_WTH_1W_
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * double _get_trace_H_1_G_H_1_G_
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * _get_yTPxy_s_
 * void _get_trace_Px_s_
 * double _get_trace_Px_G_
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * _get_yTPxy_s_
 * void _get_trace_Px_s_
 * double _get_trace_Px_G_Px_G_
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * _get_yTPxy_s_
 * double get_ytPxGPxy
 *
 * get_lambda_theta_p_1_
 * void _get_trace_H_1_s_
 * _get_yTPxy_s_
 * double get_yTPxGPxGPxy
 * */



/*My_Vector<double> & get_a_t_h_1_b( My_matrix<double> & a, My_matrix<double> & b, const Eigen_result & eigen_L, const double & lamda, My_Vector<double> & result ){
    My_matrix UTa(a.get_num_row(), a.get_num_column());
    trmul(a, eigen_L.get_eigen_vectors(), UTa);

    My_matrix UTb(b.get_num_row(), b.get_num_column());
    trmul(b, eigen_L.get_eigen_vectors(), UTb);
    int j;
    for( j=0; j<a.get_num_column(); ++j ){
        result.get_array()[j]=0;
    }
    double lamda_p_theta_1;
    for( int i=0; i<a.get_num_row(); ++i ){
        lamda_p_theta_1 = 1/(lamda*eigen_L.get_eigen_values().get_array()[i]+1);
        for( j=0; j<a.get_num_column(); ++j ){
            result.get_array()[j] = UTa.get_matrix()[i][j]*UTb.get_matrix()[i][j]*lamda_p_theta_1;
        }
    }
    return result;
}

My_Vector<double> & get_a_t_h_1_h_1_b( My_matrix<double> & a, My_matrix<double> & b, const Eigen_result & eigen_L, const double & lamda, My_Vector<double> & result ){
    My_matrix UTa(a.get_num_row(), a.get_num_column());
    trmul(a, eigen_L.get_eigen_vectors(), UTa);

    My_matrix UTb(b.get_num_row(), b.get_num_column());
    trmul(b, eigen_L.get_eigen_vectors(), UTb);
    int j;
    for( j=0; j<a.get_num_column(); ++j ){
        result.get_array()[j]=0;
    }
    double lamda_p_theta_2;
    for( int i=0; i<a.get_num_row(); ++i ){
        lamda_p_theta_2 = 2/(lamda*eigen_L.get_eigen_values().get_array()[i]+1);
        for( j=0; j<a.get_num_column(); ++j ){
            result.get_array()[j] = UTa.get_matrix()[i][j]*UTb.get_matrix()[i][j]*lamda_p_theta_2;
        }
    }
    return result;
}

My_Vector<double> & get_a_t_h_1_h_1_h_1_b( My_matrix<double> & a, My_matrix<double> & b, const Eigen_result & eigen_L, const double & lamda, My_Vector<double> & result ){
    My_matrix UTa(a.get_num_row(), a.get_num_column());
    trmul(a, eigen_L.get_eigen_vectors(), UTa);

    My_matrix UTb(b.get_num_row(), b.get_num_column());
    trmul(b, eigen_L.get_eigen_vectors(), UTb);
    int j;
    for( j=0; j<a.get_num_column(); ++j ){
        result.get_array()[j]=0;
    }
    double lamda_p_theta_3;
    for( int i=0; i<a.get_num_row(); ++i ){
        lamda_p_theta_3 = 3/(lamda*eigen_L.get_eigen_values().get_array()[i]+1);
        for( j=0; j<a.get_num_column(); ++j ){
            result.get_array()[j] = UTa.get_matrix()[i][j]*UTb.get_matrix()[i][j]*lamda_p_theta_3;
        }
    }
    return result;
}


 My_matrix<double> & get_Pi(int i, My_matrix<double> & H1, const My_matrix<double> & W, const My_matrix<double> & Wt, My_matrix<double> & Pi, const int & n){ // i = #covariants + independent variants
    if( -1 == i ){
        return H1;
    }else{
        My_matrix<double> Pi_1(n, n);
        get_Pi(i-1, H1, W, Wt, Pi_1, n);

        int j, z, is;
        My_Vector<double> wiTPi_1(n);
        for( j=0; j<n; ++j ){
            wiTPi_1.get_array()[j]=0;
            for( z=0; z<n; ++z ){
                wiTPi_1.get_array()[j] +=  Wt.get_matrix()[i][z]*Pi_1.get_matrix()[j][z];
            }
        }

        double wiTPi_1_wi_1_1 = 0.0;
        for( j=0; j<n; ++j ){
            wiTPi_1_wi_1_1 += wiTPi_1.get_array()[j]*W.get_matrix()[j][i];
        }
        wiTPi_1_wi_1_1=1/wiTPi_1_wi_1_1;

        My_Vector<double> Pi_1_wi(W.get_num_row());
        for( j=0; j<n; ++j ){
            Pi_1_wi.get_array()[j]=0;
            for( z=0; z<n; ++z ){
                Pi_1_wi.get_array()[j] += Pi_1.get_matrix()[j][z]*W.get_matrix()[j][i];
            }
        }
        for( j=0; j<n; ++j ){
            for( z=0; z<n; ++z ){
                Pi.get_matrix()[j][z] = Pi_1.get_matrix()[j][z] - Pi_1_wi.get_array()[j] * wiTPi_1_wi_1_1 * wiTPi_1.get_array()[z];
            }
        }
        return Pi;
    }
}

// this function is somehow complex, should be careful for debuging
double get_det_W_T_h_1_W( int i, My_matrix<double> & W, const My_matrix<double> & H1, const int & n, const My_matrix<double> * P ){
    if( i==1 ){
        double witpi_1wi=0.0;
        My_Vector witpi_1(n);
        int j, z;
        for(j=0; j<n; ++j){
            witpi_1.get_array()[j]=0;
            for(z=0; z<n; ++z){
                witpi_1.get_array()[j] += W.get_matrix()[z][i] * P[i-1].get_matrix()[j][z];
            }
            witpi_1wi += witpi_1.get_array()[j] * W.get_matrix()[j][i];
        }
        My_Vector<double> wi1t_h1(n);
        for( j=0; j<n; ++j ){
            wi1t_h1.get_array()[j] = 0;
            for(z=0; z<n; ++z){
                wi1t_h1.get_array()[j] += W.get_matrix()[z][0] * H1.get_matrix()[j][z];
            }
        }
        double wi1t_h1_wi1(0.0);
        for( j=0; j<W.get_num_row(); ++j ){
            wi1t_h1_wi1 += wi1t_h1.get_array()[j]*W.get_matrix()[j][0];
        }
        return witpi_1wi * wi1t_h1_wi1;
    }else{
        double witpi_1wi=0.0;
        My_Vector witpi_1(n);
        int j, z;
        for(j=0; j<n; ++j){
            witpi_1.get_array()[j]=0;
            for(z=0; z<n; ++z){
                witpi_1.get_array()[j] += W.get_matrix()[z][i] * P[i-1].get_matrix()[j][z];
            }
            witpi_1wi += witpi_1.get_array()[j] * W.get_matrix()[j][i];
        }
        int next_i = i - 1;
        return witpi_1wi * get_det_W_T_h_1_W( next_i, W, H1, n, P );
    }
}
*/


/*
// the euqation on the online method part of calculation of the basic quantities
double _get_aTh_1b( const My_Vector<double> & va, const My_Vector<double> & vb, const int & n, const Eigen_result & eigen_L, const double & lambda ){
    double result = 0.0;
    for( int i=0; i<n; ++i ){
        result += (va.get_array()[i]*vb.get_array()[i])/(lambda*eigen_L.get_eigen_values().get_array()[i]+1);
    }
    return result;
}

// the euqation on the online method part of calculation of the basic quantities
double _get_aTh_1h_1b( const My_Vector<double> & va, const My_Vector<double> & vb, const int & n, const Eigen_result & eigen_L, const double & lambda ){
    double result = 0.0;
    for( int i=0; i<n; ++i ){
        result += (va.get_array()[i]*vb.get_array()[i])/pow((lambda*eigen_L.get_eigen_values().get_array()[i]+1),2);
    }
    return result;
}

// the euqation on the online method part of calculation of the basic quantities
double _get_aTh_1h_1h_1b( const My_Vector<double> & va, const My_Vector<double> & vb, const int & n, const Eigen_result & eigen_L, const double & lambda ){
    double result = 0.0;
    for( int i=0; i<n; ++i ){
        result += (va.get_array()[i]*vb.get_array()[i])/pow((lambda*eigen_L.get_eigen_values().get_array()[i]+1),3);
    }
    return result;
}
*/
/*// the third equation of 3.1.5 on page 6 of supplementary document
void _set_aTPibs( const int & q, const double & aTh_1b, const double & aTh_1Wi, const double & bTh_1Wi, const double & WiTh_1Wi, const My_Vector<double> & aTPi_1Wi, const My_Vector<double> & bTPi_1Wi, const My_Vector<double> & WiTPi_1Wi, My_Vector<double> & aTPib){//todo
    aTPib.get_array()[0] = aTh_1b - aTh_1Wi*bTh_1Wi*WiTh_1Wi;
    for( int i=1; i<q; ++i ){
        aTPib.get_array()[i] = aTPi_1Wi.get_array()[i-1] - aTPi_1Wi.get_array()[i]*bTPi_1Wi.get_array()[i]/WiTPi_1Wi.get_array()[i];
    }
}


// the forth equation of 3.1.5 on page 6 of supplementary document
void _set_atPiPib( const int & q, const double & aTh_1h_1b, const double & aTH_1Wi, const double & bTH_1Wi,
                  const double & WiTh_1Wi, const double & WiTH_1H_1Wi, const double & bTH_1H_1Wi, const double & WiTH_1Wi,
                  const double & aTH_1H_1Wi, const My_Vector<double> & aTPi_1Wi,
                const My_Vector<double> & bTPi_1Wi, const My_Vector<double> & WiTPi_1Wi,
            const My_Vector<double> & aTPi_1Pi_1Wi, const My_Vector<double> & bTPi_1Pi_1Wi,
        const My_Vector<double> & WiTPi_1Pi_1Wi, My_Vector<double> & aTPiPib){//todo

    aTPiPib.get_array()[0] = aTh_1h_1b + aTH_1Wi*bTH_1Wi*WiTH_1H_1Wi/pow(WiTh_1Wi, 2)-
            aTH_1Wi*bTH_1H_1Wi/WiTH_1Wi-
            bTH_1Wi*aTH_1H_1Wi/WiTH_1Wi;
    for( int i=1; i<q; ++i ){
        aTPiPib.get_array()[i] = aTPiPib.get_array()[i-1] + aTPi_1Wi.get_array()[i-1]*bTPi_1Wi.get_array()[i-1]*
                                 WiTPi_1Pi_1Wi.get_array()[i-1]/pow(WiTPi_1Wi.get_array()[i-1], 2)-
                                 aTPi_1Wi.get_array()[i-1]*bTPi_1Pi_1Wi.get_array()[i-1]/WiTPi_1Wi.get_array()[i-1]-
                                 bTPi_1Wi.get_array()[i-1]*aTPi_1Pi_1Wi.get_array()[i-1]/WiTPi_1Wi.get_array()[i-1];
    }
}

void _set_atPiPiPib( const int & q, const double & aTh_1h_1h_1b, const double & bTH_1H_1H_1Wi, const double & aTH_1H_1H_1Wi,
                  const double & aTH_1Wi, const double & bTH_1Wi, const double & WiTH_1H_1H_1Wi,
                  const double & WiTh_1Wi, const double & WiTH_1H_1Wi, const double & bTH_1H_1Wi, const double & WiTH_1Wi,
                  const double & aTH_1H_1Wi, const My_Vector<double> & aTPi_1Wi,
                  const My_Vector<double> & bTPi_1Wi, const My_Vector<double> & WiTPi_1Wi,
                  const My_Vector<double> & aTPi_1Pi_1Wi, const My_Vector<double> & bTPi_1Pi_1Wi,
                  const My_Vector<double> & WiTPi_1Pi_1Wi, const My_Vector<double> & bTPi_1Pi_1Pi_1Wi,
                  const My_Vector<double> & aTPi_1Pi_1Pi_1Wi, const My_Vector<double> & WiTPi_1Pi_1Pi_1Wi, My_Vector<double> & atPiPiPib){//todo
    atPiPiPib.get_array()[0] = aTh_1h_1h_1b
                               - aTH_1Wi*bTH_1Wi*pow(WiTH_1H_1Wi,2)/pow(WiTh_1Wi, 3)
                               - aTH_1Wi*bTH_1H_1H_1Wi/WiTH_1Wi
                               - bTH_1Wi*aTH_1H_1H_1Wi/WiTH_1Wi
                               - aTH_1H_1Wi*bTH_1H_1Wi/WiTH_1Wi
                               + aTH_1Wi*bTH_1H_1Wi*WiTH_1H_1Wi/pow(WiTH_1Wi, 2)
                               + bTH_1Wi*aTH_1H_1Wi*WiTH_1H_1Wi/pow(WiTH_1Wi, 2)
                               + aTH_1Wi*bTH_1Wi*WiTH_1H_1H_1Wi/pow(WiTH_1Wi, 2);
    for( int i=1; i<q; ++i ){
        atPiPiPib.get_array()[i] = atPiPiPib.get_array()[i-1]
        - aTPi_1Wi.get_array()[i]*bTPi_1Wi.get_array()[i]*pow(WiTPi_1Pi_1Wi.get_array()[i],2)/pow(WiTPi_1Wi.get_array()[i], 3)
        - aTPi_1Wi.get_array()[i]*bTPi_1Pi_1Pi_1Wi.get_array()[i]/WiTPi_1Wi.get_array()[i]
        - bTPi_1Wi.get_array()[i]*aTPi_1Pi_1Pi_1Wi.get_array()[i]/WiTPi_1Wi.get_array()[i]
        - aTPi_1Pi_1Wi.get_array()[i]*bTPi_1Pi_1Wi.get_array()[i]/WiTPi_1Wi.get_array()[i]
        + aTPi_1Wi.get_array()[i]*bTPi_1Pi_1Wi.get_array()[i]*WiTPi_1Pi_1Wi.get_array()[i]/pow(WiTPi_1Wi.get_array()[i], 2)
        + bTPi_1Wi.get_array()[i]*aTPi_1Pi_1Wi.get_array()[i]*WiTPi_1Pi_1Wi.get_array()[i]/pow(WiTPi_1Wi.get_array()[i], 2)
        + aTPi_1Wi.get_array()[i]*bTPi_1Wi.get_array()[i]*WiTPi_1Pi_1Pi_1Wi.get_array()[i]/pow(WiTPi_1Wi.get_array()[i], 2);
    }
}



// the first equation of 3.1.5 on page 6 of supplementary document
void _setTracePis( const int & q, const double & trace_H_1, My_Vector<double> & TracePis, My_Vector<double> & WiTPi_1wi, My_Vector<double> & WiTPi_1Pi_1wi){ //todo set WiTPi_1Pi_1wi and WiTPi_1wi
    TracePis.get_array()[0] = trace_H_1 - WiTPi_1Pi_1wi.get_array()[0]/WiTPi_1wi.get_array()[0];
    for( int i=1; i<q; ++i ){
        TracePis.get_array()[i]=0;
        TracePis.get_array()[i] = TracePis.get_array()[i-1]-WiTPi_1Pi_1wi.get_array()[i]/WiTPi_1wi.get_array()[i];
    }
}

// the second equation of 3.1.5 on page 6 of supplementary document
void _setTracePiPis( const int & q, const double & trace_H_1_H_1, My_Vector<double> & TracePiPis, My_Vector<double> & WiTPi_1wi, My_Vector<double> & WiTPi_1Pi_1wi, My_Vector<double> & WiTPi_1Pi_1Pi_1wi){ //todo set WiTPi_1Pi_1wi and WiTPi_1wi
    double PP0 = trace_H_1_H_1;
    TracePiPis.get_array()[0] = PP0 + pow(WiTPi_1Pi_1wi.get_array()[0], 2)/pow(WiTPi_1wi.get_array()[0],2)
            -2.0*WiTPi_1Pi_1Pi_1wi.get_array()[0]/WiTPi_1wi.get_array()[0];
    for( int i=1; i<q; ++i ){
        TracePiPis.get_array()[i]=0;
        TracePiPis.get_array()[i] = TracePiPis.get_array()[i-1]+pow(WiTPi_1Pi_1wi.get_array()[i], 2)/pow(WiTPi_1wi.get_array()[i],2)
                  -2.0*WiTPi_1Pi_1Pi_1wi.get_array()[i]/WiTPi_1wi.get_array()[i];
    }
}
*/

//double _get_trace_PxH_(const int & n, const int & c){
//    return double(n-c-1);
//}



//
//// todo is there any special way to to this????
//double get_trace_Pi( const My_matrix<double> & Pi, const int & n ){
//    double result=0.0;
//    for( int i=0; i<n; ++i ){
//        result += Pi.get_matrix()[i][i];
//    }
//    return result;
//}
//
//// todo is there any special way to to this????
//double get_trace_Pi_Pi( const My_matrix<double> & Pi, const int & n ){
//    double result=0.0;
//    for( int i=0; i<n; ++i ){
//        result += Pi.get_matrix()[i][i]*Pi.get_matrix()[i][i];
//    }
//    return result;
//}

//// this function will generate the result of ytPx
//double _get_ytPxy_( const My_Vector<double> & y, My_Vector<double> & ytPx, const My_matrix<double> & Pi, const int & n ){
//    int i, j;
//    for( i=0; i<n; ++i ){
//        ytPx.get_array()[i]=0;
//        for( j=0; j<n; ++j ) {
//            ytPx.get_array()[i] += y.get_array()[j] * Pi.get_matrix()[i][j];
//        }
//    }
//    double result=0.0;
//    for( i=0; i<n; ++i ){
//        result += ytPx.get_array()[i] * y.get_array()[i];
//    }
//    return result;
//}
//
//// this function rely on get_ytPxy, you should run get_ytPxy to get the value of ytPx
//double get_ytPxPxy(const My_Vector<double> & y, const My_Vector<double> & ytPx, const My_matrix<double> & P, const int & n, My_Vector<double> & ytPxPx ){
//    int i, j;
//    for( i=0; i<n; ++i ){
//        ytPxPx.get_array()[i]=0;
//        for( j=0; j<n; ++j ) {
//            ytPxPx.get_array()[i] += ytPx.get_array()[j] * P.get_matrix()[i][j];
//        }
//    }
//
//    double result=0.0;
//    for( i=0; i<n; ++i ){
//        result += ytPxPx.get_array()[i] * y.get_array()[i];
//    }
//    return result;
//}
//
//// this function rely on get_ytPxPxy, you shoudl run get_ytPxPxy to get the value of ytPxPx
//double get_ytPxPxPxy(const My_Vector<double> & y, const My_matrix<double> & P, const int & n, const My_Vector<double> & ytPxPx ){
//    My_Vector<double> ytPxPxPx(n);
//    int i, j;
//    for( i=0; i<n; ++i ){
//        ytPxPxPx.get_array()[i]=0;
//        for( j=0; j<n; ++j ) {
//            ytPxPxPx.get_array()[i] += ytPxPx.get_array()[j] * P.get_matrix()[i][j];
//        }
//    }
//
//    double result=0.0;
//    for( i=0; i<n; ++i ){
//        result += ytPxPxPx.get_array()[i] * y.get_array()[i];
//    }
//    return result;
//}




