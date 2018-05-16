// =====================================================================================
//
//       Filename:  emmax_impl.cpp
//
//    Description:
//
//        Version:  1.0
//        Created:  02/27/2018 05:53:23 PM
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
 *  The R source code of emma and the intel version c code of emmax is referred
 *  The GEMMA paper is referred for the NR algorithm
 *  And the SVS gives a lot of details could not be found any where else
 *  http://doc.goldenhelix.com/SVS/latest/svsmanual/mixedModelMethods/overview.html#overviewofmixedlinearmodels
 *
 *
 *
************************************************************************/


#include "emmax_impl.h"


Eigen_result _get_eigen_L_(const My_matrix<double> & k){
    return eigen(k);
}

Eigen_result _get_eigen_R_(const My_matrix<double> & x, const My_matrix<double> & k){
    int i, j;

    My_matrix<double> xt(x.get_num_column(), x.get_num_row());
    T_matrix(x, xt);

    My_matrix<double> x_t_x(xt.get_num_row(), x.get_num_column());
    trmul(xt, x, x_t_x);

    My_matrix<double> x_squared_inverse(x_t_x);
    inverse_matrix(x_squared_inverse);

    My_matrix<double> x_x_squared_inverse(x.get_num_row(), x_squared_inverse.get_num_column());
    trmul(x, x_squared_inverse, x_x_squared_inverse);

    My_matrix<double> hat_matrix (x_x_squared_inverse.get_num_row(), xt.get_num_column());
    trmul(x_x_squared_inverse, xt, hat_matrix);

    for( i=0; i<x.get_num_row(); i++ ){
        hat_matrix.get_matrix()[i][i] = 1 - hat_matrix.get_matrix()[i][i];
        for( j=0; j<i; j++ ){
            hat_matrix.get_matrix()[i][j]=-hat_matrix.get_matrix()[i][j];
            hat_matrix.get_matrix()[j][i]=-hat_matrix.get_matrix()[j][i];
        }
    }
    My_matrix<double> s_k(hat_matrix.get_num_row(), k.get_num_column());
    trmul(hat_matrix, k, s_k);
    My_matrix<double> s_k_s(s_k.get_num_row(), hat_matrix.get_num_column());
    trmul(s_k, hat_matrix, s_k_s);
    int keep = k.get_num_column()-x.get_num_column();
    Eigen_result eigen_result1 = eigen(s_k_s, keep);
//
//    int size = k.get_num_column() - x.get_num_column();
//    My_Vector eigen_value( size );
//    My_matrix eigen_vector(size, k.get_num_column());
//    for( i=0; i< size; ++i ){
//        eigen_value.get_array()[i]=eigen_result1.get_eigen_values().get_array()[x.get_num_column()+i];
//        memcpy(eigen_vector.get_matrix()[i], eigen_result1.get_eigen_vectors().get_matrix()[x.get_num_column()+i], k.get_num_column()*sizeof(double) );
//    }
//    Eigen_result eigen_result(eigen_value, eigen_vector);
    return eigen_result1;
}

// log-likelihoods (eq. 6 from paper), this is the emma.delta.ML.LL.wo.Z function in R emma
double _ll_(const double & delta, const Eigen_result & eigen_R, const Eigen_result & eigen_L, const My_Vector<double> & sq_etas){
    int n = eigen_L.get_eigen_values().get_length();
    double c_1 = 0.5 * n * (log(n / (2.0 * M_PI)) - 1);
    My_Vector<double> v1(eigen_R.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_R.get_eigen_values(), delta, v1);

    My_Vector<double> v2(eigen_L.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_L.get_eigen_values(), delta, v2);

    My_Vector<double> sq_etas_d_v1(sq_etas.get_length());
    quotient_of_two_vectors<double, double, double>(sq_etas, v1, sq_etas_d_v1);

    double sum_log_v2=0.0;
    for( int i=0; i<n; ++i ){
        sum_log_v2 += log(v2.get_array()[i]);
    }
    return c_1 - 0.5*(n*log(sum_of_a_vector(sq_etas_d_v1))) + sum_log_v2;
}

//diffrentiated log-likelihoods (eq. 8 from paper) the emma.delta.ML.dLL.wo.Z function in R emma
double _dll_( const double & delta, const Eigen_result & eigen_R, const Eigen_result & eigen_L, const My_Vector<double> & sq_etas){
    int n = eigen_L.get_eigen_values().get_length();
    My_Vector<double> v1(eigen_R.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_R.get_eigen_values(), delta, v1);

    My_Vector<double> v2(sq_etas.get_length());
    quotient_of_two_vectors<double, double, double>(sq_etas, v1, v2);

    My_Vector<double> v3(eigen_L.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_L.get_eigen_values(), delta, v3);

    My_Vector<double> v2_2_v1(v2.get_length());
    quotient_of_two_vectors<double, double, double>(v2, v1, v2_2_v1);

    for( int i=0; i<n ; ++i ){
        v3.get_array()[i] = 1.0/v3.get_array()[i];
    }
    double res = (n * sum_of_a_vector(v2_2_v1) / sum_of_a_vector(v2) - sum_of_a_vector(v3));
    return res;
}


//newton_ll( new_opt_delta, eps, maxiter, n, eig_R_values, eig_L_values, sq_etas, p);
int newton_ll( double* new_opt_delta, const double & eps, const int & js, const int & n, const Eigen_result & eig_R, const Eigen_result & eig_L, const My_Vector<double> & sq_etas){ // js maximum number of iteration
    int k=1,l=js;
    double y1, y2, d, p, x0, x1=*new_opt_delta;
    x0=*new_opt_delta;
    y1=_ll_(x0, eig_R, eig_L, sq_etas);
    y2=_dll_(x0, eig_R, eig_L, sq_etas);
    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
        if( fabs(y2)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        x1 = x0 - y1/y2;
        y1=_ll_( x1, eig_R, eig_L, sq_etas);
        y2=_dll_( x1, eig_R, eig_L, sq_etas);
        d=fabs(x1-x0);
        p=fabs(y1);
        if( p>d ){
            d=p;
        }
        x0=x1;
        l=l-1;
    }
    *new_opt_delta=x1;
    k=js-l;
    return k;
}

//log-likelihoods (eq. 7 from paper)
double _rell_(const double & delta, const Eigen_result & eigen_R, const My_Vector<double> & sq_etas){
    double c_1 = 0.5 * eigen_R.get_eigen_values().get_length() * (log(eigen_R.get_eigen_values().get_length() / (2.0 * M_PI)) - 1);

    My_Vector<double> v1(eigen_R.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_R.get_eigen_values(), delta, v1);

    My_Vector<double> sq_etas_d_v1 (sq_etas.get_length());
    quotient_of_two_vectors<double, double, double>(sq_etas, v1, sq_etas_d_v1);

    double res = 0.0;
    for ( int i=0; i< eigen_R.get_eigen_values().get_length(); ++i ){
        res += log(v1.get_array()[i]);
    }
    res = c_1 - 0.5 * (eigen_R.get_eigen_values().get_length() * log(sum_of_a_vector(sq_etas_d_v1)) + res);
    return res;
}

// diffrentiated log-likelihoods (*2) (eq. 9 from paper)
double _redll_( const double & delta, const Eigen_result & eig_R, const My_Vector<double> & sq_etas ){
    My_Vector<double> v1(eig_R.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eig_R.get_eigen_values(), delta, v1);

    My_Vector<double> v2(eig_R.get_eigen_values().get_length());
    quotient_of_two_vectors<double, double, double>(sq_etas, v1, v2);

    My_Vector<double> v2_2_v1( eig_R.get_eigen_values().get_length() );
    quotient_of_two_vectors<double, double, double>(v2, v1, v2_2_v1);

    double res=0.0;
    for ( int i=0; i<eig_R.get_eigen_values().get_length(); ++i){
        res += 1.0 / v1.get_array()[i];
    }
    res = (eig_R.get_eigen_values().get_length() * sum_of_a_vector(v2_2_v1) / sum_of_a_vector(v2) -res );
    return res;
}

int newton_reml( double* new_opt_delta, const double & eps, const int & js, double delta, Eigen_result & eig_R, My_Vector<double> & sq_etas ){ // js maximum number of iteration
    int k=1,l=js;
    double y1, y2, d, p, x0, x1=*new_opt_delta;
    x0=*new_opt_delta;
    y1=_rell_( x0, eig_R, sq_etas );
    y2=_redll_( x0, eig_R, sq_etas);
    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
//        std::cout << "y2 " << y2 << std::endl;
        if( fabs(y2)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        x1 = x0 - y1/y2;
        y1=_rell_( x1, eig_R, sq_etas );
        y2=_redll_( x1, eig_R, sq_etas );
        d=fabs(x1-x0);
        p=fabs(y1);
        if( p>d ){
            d=p;
        }
        x0=x1;
        l -= 1;
    }
    *new_opt_delta=x1;
    k=js-l;
    return k;
}

Emma_result get_estimates ( const My_Vector<double> & y, const My_matrix<double> & k, const My_matrix<double> & x, const My_matrix<double> & xs, const Eigen_result & eigen_L, int & ngrids, double & llim, double &ulim, double & eps,
                     const std::string & method, const int & maxiter ){
    //std::cout << y.get_length() << "\t" << k.get_num_column() << "\t" << k.get_num_row() << "\t" << x.get_num_row() << std::endl;
    assert(y.get_length() == k.get_num_column() && y.get_length()==k.get_num_row() && y.get_length()==x.get_num_row() );

    int i, j;
    int n = y.get_length();
    int q = x.get_num_column()+xs.get_num_column();
    int p = n - q;
    int m = ngrids+1;

    My_matrix<double> full_x(x.get_num_row(), x.get_num_column()+xs.get_num_column());
    if( xs.get_num_column()>0 ){
        append_two_matrix( x, xs, full_x);
    }else{
        full_x=x;
    }

//    std::cout << "line 210 " << currentDateTime() << std::endl;
    Eigen_result eigen_R = _get_eigen_R_(full_x, k);
//    std::cout << "line 212 " << currentDateTime() << std::endl;
    //Eigen_result eigen_L = _get_eigen_L_(k);
//    std::cout << "line 213 " << currentDateTime() << std::endl;
//
//    std::cout << "eigen_L_V " << std::endl;
//    for (i=0; i<n; ++i) {
//        std::cout << eigen_L.get_eigen_values().get_array()[i]<< "    ";
//    }
//    std::cout <<  std::endl;
//    for (i=0; i<n; ++i){
//        std::cout << "eigen_L " << std::endl;
//        for( j=0; j < n; ++j ){
//            std::cout << eigen_L.get_eigen_vectors().get_matrix()[i][j] << "    ";
//        }
//        std::cout <<  std::endl;
//    }
//    std::cout <<  std::endl;
//    std::cout <<  std::endl;
//    for (i=0; i<m; ++i){
//        std::cout << "eigen_R " << std::endl;
//        for( j=0; j < n; ++j ){
//            std::cout << eigen_R.get_eigen_vectors().get_matrix()[i][j] << "    ";
//        }
//       std::cout <<  std::endl;
//    }
//    for( j=0; j < n; ++j ){
//        std::cout << "phenotype " << y.get_array()[j] << std::endl;
//    }


    My_Vector<double> etas(p);
    double t;
    for( i=0; i < p; ++i ){
        t = 0.0;
        for( j=0; j < n; ++j ){
            t = t + eigen_R.get_eigen_vectors().get_matrix()[i][j]*y.get_array()[j];
        }
        etas.get_array()[i]=t;
//        std::cout << "etas " << etas.get_array()[i] << std::endl;
    }
//    for( i=0; i < p; ++i ){
//        std::cout << "etas " << i << " " << etas.get_array()[i] << std::endl;
//    }
    My_Vector<double> sq_etas(p);
    for( i=0; i < p; ++i ){
        sq_etas.get_array()[i] = etas.get_array()[i] * etas.get_array()[i];
//        std::cout << "sq_etas " << sq_etas.get_array()[i] << std::endl;
    }
//    i=0;
//    for( j=0; j < n; ++j ){
//        std::cout << "eigen " << eigen_R.get_eigen_vectors().get_matrix()[i][j] << std::endl;
//    }
//    for( j=0; j < n; ++j ){
//        std::cout << "phenotype " << y.get_array()[j] << std::endl;
//    }

    My_Vector<double> log_deltas(m); // the space for deltas to search
    My_Vector<double> deltas(m);
    for ( i=0; i < m; ++i ){
        log_deltas.get_array()[i] = (double(i) / ngrids)*(ulim - llim) + llim;
        deltas.get_array()[i] = exp(log_deltas.get_array()[i]);
    }

    My_matrix<double> lambdas(p, m);
    for ( i=0; i < p; ++i ){
        for ( j=0; j < m; ++j ){
            lambdas.get_matrix()[i][j] = eigen_R.get_eigen_values().get_array()[i] + deltas.get_array()[j];
        }
    }

    My_Vector<double> s1(m);
    for ( i=0; i < m; ++i ){
        s1.get_array()[i] = 0.0;
        for ( j=0; j < p; ++j ){
            s1.get_array()[i] += sq_etas.get_array()[j]/lambdas.get_matrix()[j][i];
        }
//        std::cout << "line 204 " << i << " " << s1.get_array()[i] << std::endl;
    }

    My_Vector<double> s3(m);
    for( i=0; i<m; ++i ){
        s3.get_array()[i] = 0.0;
        for ( j=0; j<p; ++j ){
            s3.get_array()[i]+=sq_etas.get_array()[j]/(lambdas.get_matrix()[j][i]*lambdas.get_matrix()[j][i]);
        }
//        std::cout << "s3  " << i << " " << s3.get_array()[i] << std::endl;
    }

    My_Vector<double> lls(m);
    My_Vector<double> dlls(m);
    if (method.compare("REML") == 0){
        My_Vector<double> s2(m);
        for( i=0; i<m; ++i ){
            s2.get_array()[i]=0.0;
            for( j=0; j<p; ++j  ){
                s2.get_array()[i] += log(lambdas.get_matrix()[j][i]);
            }
//            std::cout << "line 284 " << s2.get_array()[i] << std::endl;
        }
        for( i=0; i<m; ++i ){
//            std::cout << "PI " << M_PI << " " << s1.get_array()[i] << " " << log(s1.get_array()[i]) << " " << s2.get_array()[i] << std::endl;
            lls.get_array()[i] = 0.5 * (p * (log((p) / (2.0 * M_PI)) - 1 - log(s1.get_array()[i])) - s2.get_array()[i]);
//            std::cout << "line 288 " << lls.get_array()[i] << std::endl;
        }
        My_Vector<double> s4(m);
        for( i=0; i<m; ++i ){
            s4.get_array()[i] = 0.0;
            for( j=0; j<p; ++j ){
                s4.get_array()[i] += 1.0/lambdas.get_matrix()[j][i];
            }
//            std::cout << "s4 " << s4.get_array()[i] << std::endl;
        }

        for( i=0; i<m; ++i ){
            dlls.get_array()[i] = 0.5 * (p * s3.get_array()[i] / s1.get_array()[i] - s4.get_array()[i]);
//            std::cout << "dlls " << dlls.get_array()[i] << " p " << p << " s3.get_array()[i] " << s3.get_array()[i] << "  s1.get_array()[i] " <<  s1.get_array()[i] << " s4.get_array()[i] " << s4.get_array()[i] << std::endl;
        }
    }else if (method.compare("ML") == 0){
        My_matrix<double> xis(n, m);
        for ( i=0; i < n; ++i ){
            for ( j=0; j < m; ++j ){
                xis.get_matrix()[i][j] = eigen_L.get_eigen_values().get_array()[i] + deltas.get_array()[j];
            }
        }
        My_Vector<double> s2(m);
        for( i=0; i<n; ++i ){
            s2.get_array()[i]=0.0;
            for( j=0; j<m; ++j  ){
                s2.get_array()[j] += log(xis.get_matrix()[i][j]);
            }
        }
        for( i=0; i<m; ++i ){
            lls.get_array()[i] = 0.5 * (p * (log((p) / (2.0 * M_PI)) - 1 - log(s1.get_array()[i])) - s2.get_array()[i]);
//            std::cout << lls.get_array()[i] << std::endl;
        }
        My_Vector<double> s4(m);
        for( i=0; i<n; ++i ){
            s4.get_array()[i] = 0.0;
            for( j=0; j<m; ++j ){
                s4.get_array()[j]=1.0/xis.get_matrix()[i][j];
            }
        }
        for( i=0; i<m; ++i ){
            dlls.get_array()[i] = 0.5 * (n*s3.get_array()[i]/s1.get_array()[i]-s4.get_array()[i]);
        }
    }

    int max_ll_i = argmax(lls);
    double max_ll = lls.get_array()[max_ll_i];
    double last_dll = dlls.get_array()[0];
    double  last_ll = lls.get_array()[0];

//    std::cout << "max_ll_i " << max_ll_i << std::endl;
//    std::cout << "max_ll " << max_ll << std::endl;
//    std::cout << "last_dll " << last_dll << std::endl;
//    std::cout << "last_ll " << last_ll << std::endl;

    std::map<int, double> zero_intervals;
    for (i =0; i<m; ++i){
        if (dlls.get_array()[i] < 0 && last_dll > 0){
            zero_intervals[i]=(lls.get_array()[i] + last_ll) * 0.5;
        }
        last_ll = lls.get_array()[i];
        last_dll = dlls.get_array()[i];
    }

    double opt_delta;
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
        opt_delta = 0.5 * (deltas.get_array()[opt_i - 1] + deltas.get_array()[opt_i]);
        std::cout << "line 344 opt_delta " << opt_delta << std::endl;
        double * new_opt_delta;
        *new_opt_delta = opt_delta;
        if (method.compare("REML") == 0) {
            newton_reml( new_opt_delta, eps, maxiter, opt_delta, eigen_R, sq_etas);
        }else if (method.compare("ML") == 0){
            newton_ll( new_opt_delta, eps, maxiter, n, eigen_R, eigen_L, sq_etas);
        }
        if (opt_i > 1 && deltas.get_array()[opt_i - 1] - eps < *new_opt_delta && *new_opt_delta < deltas.get_array()[opt_i] + eps){
            opt_delta = *new_opt_delta;
            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
        }else if (opt_i == 1 && 0.0 < *new_opt_delta && *new_opt_delta < deltas.get_array()[opt_i] + eps){
            opt_delta = *new_opt_delta;
            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
        }else if(opt_i == m - 1 && *new_opt_delta > deltas.get_array()[opt_i - 1] - eps ){
            opt_delta = *new_opt_delta;
            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
        }
        if (method.compare("REML") == 0) {
            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
        }else if (method.compare("ML") == 0) {
            opt_ll = _ll_(opt_delta, eigen_R, eigen_L, sq_etas);
        }
        if (opt_ll < max_ll){
            opt_delta = deltas.get_array()[max_ll_i];
        }
    }else{
        opt_delta = deltas.get_array()[max_ll_i];
        opt_ll = max_ll;
    }
//    std::cout << "line 418 opt_delta " << opt_delta << std::endl;
//    std::cout << "line 419 opt_ll " << opt_ll << std::endl;


    My_matrix<double> R(p, p);
    for ( j=0; j<p; ++j ){
        for ( i=0; i<p; ++i ) {
            R.get_matrix()[i][j] = sq_etas.get_array()[i] / (eigen_R.get_eigen_values().get_array()[j] + opt_delta);

        }

//        std::cout << i << " " << "sq_etas.get_array()[i] " << sq_etas.get_array()[i]
//                  << " eigen_R.get_eigen_values().get_array()[i] "
//                  << eigen_R.get_eigen_values().get_array()[i] << " opt_delta " << opt_delta << " R.get_array()[i] "
//                  << R.get_array()[i] << std::endl;
    }
    double opt_vg = 0.0;
    for ( i=0; i<p; ++i ){
        for ( j=0; j<p; ++j ) {
//            std::cout << i << " " << j << " " << R.get_matrix()[i][j] << " " << std::endl;
            opt_vg += R.get_matrix()[i][j];
        }
//        std::cout << std::endl;
    }
    opt_vg=opt_vg/p; //This is the REML estimation. the ML estimation is np.sum(R) / n
//    std::cout << "opt_vg " << opt_vg << std::endl;
    double opt_ve = opt_vg * opt_delta;
//    std::cout << "opt_ve " << opt_ve << std::endl;

    My_matrix<double> H_sqrt_inv (n, n);
    for( i=0; i<n; i++ ){
        double this_value = 1.0/sqrt(eigen_L.get_eigen_values().get_array()[i]+opt_delta);
        for( j=0; j<n; j++ ){
            H_sqrt_inv.get_matrix()[i][j] = this_value*eigen_L.get_eigen_vectors().get_matrix()[i][j];
        }
    }
//    for( i=0; i<n; i++ ){
//        for( j=0; j<n; j++ ){
//            std::cout << H_sqrt_inv.get_matrix()[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
    My_matrix<double> X_t(n, q);
    trmul(H_sqrt_inv, full_x, X_t);
    My_Vector<double> Y_t (n);
    for(i=0; i<n; ++i){
        Y_t.get_array()[i] = 0;
        for(j=0; j<n; ++j) {
            Y_t.get_array()[i] += H_sqrt_inv.get_matrix()[i][j] * y.get_array()[j];
        }
//        std::cout << " Y_t " << Y_t.get_array()[i] << std::endl;
    }
//    std::cout << "472" << std::endl;
    My_Vector<double> beta_est(q);
//    std::cout << "474" << std::endl;
    lsq(X_t, Y_t, beta_est);
//    std::cout << "475" << std::endl;
    My_Vector<double> x_beta(n);
    for( i=0; i<n; ++i ){
        x_beta.get_array()[i]=0;
        for( j=0; j<q; ++j ){
            x_beta.get_array()[i] += X_t.get_matrix()[i][j]*beta_est.get_array()[j];
        }
//        std::cout << " x_beta " << x_beta.get_array()[i] << " Y_t " << Y_t.get_array()[i] << " Y_t - x_beta " << Y_t.get_array()[i]-x_beta.get_array()[i] << std::endl;
    }
    double mahalanobis_rss = 0.0;
    for ( i=0; i<n; ++i ){
        mahalanobis_rss += pow( Y_t.get_array()[i]-x_beta.get_array()[i], 2);
//        std::cout << Y_t.get_array()[i] << "  " << x_beta.get_array()[i] << std::endl;
    }
//    std::cout << " mahalanobis_rss " << mahalanobis_rss << std::endl;
    //get the likelyhood value for LRT test
//    opt_ll = _ll_(opt_delta, eigen_R, eigen_L, sq_etas);
//    std::cout << "opt_ll " << opt_ll << std::endl;

//    double opt_dll = _dll_(opt_delta, eigen_R, eigen_L, sq_etas);
//    std::cout << "opt_dll " << opt_dll << std::endl;

//    double opt_rell = _rell_(opt_delta, eigen_R, sq_etas);
//    std::cout << "opt_rell " << opt_rell << std::endl;

//    double opt_redll = _redll_(opt_delta, eigen_R, sq_etas);
//    std::cout << "opt_redll " << opt_redll << std::endl;


    if( xs.get_num_column() > 0 ){
        My_matrix<double> h0_X (n, x.get_num_column());
        trmul(H_sqrt_inv, x, h0_X);
//        std::cout << "line 434" << std::endl;
        My_Vector<double> h0_beta_est(h0_X.get_num_column());
        lsq( h0_X, Y_t, h0_beta_est);

        My_Vector<double> h0_x_beta(n);
        for( i=0; i<n; ++i ){
            h0_x_beta.get_array()[i]=0;
            for( j=0; j<x.get_num_column(); ++j ){
                h0_x_beta.get_array()[i] += h0_X.get_matrix()[i][j]*h0_beta_est.get_array()[j];
            }
        }

        double h0_rss = 0.0;
        for ( i=0; i<n; ++i ){
            h0_rss += pow(Y_t.get_array()[i]-h0_x_beta.get_array()[i], 2);
        }
//        std::cout << " h0_rss " << h0_rss << std::endl;
        double f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.get_num_column();
//        std::cout << "line 444 " << h0_rss << " " << mahalanobis_rss << std::endl;
        int freedome1 = xs.get_num_column();
        double p_val = sf(f_stat, freedome1, p);
        Emma_result r_re(Y_t, H_sqrt_inv, mahalanobis_rss, p_val);
        return r_re;
    }
    Emma_result r_re(Y_t, H_sqrt_inv, mahalanobis_rss);
    return r_re;
}
