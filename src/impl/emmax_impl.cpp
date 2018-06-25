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

    for( i=0; i<x.get_num_row(); ++i ){
        hat_matrix.get_matrix()[i][i] = 1 - hat_matrix.get_matrix()[i][i];
        for( j=0; j<i; ++j ){
            hat_matrix.get_matrix()[i][j]=-hat_matrix.get_matrix()[i][j];
            hat_matrix.get_matrix()[j][i]=-hat_matrix.get_matrix()[j][i];
        }
    }
    My_matrix<double> k_p_1(k.get_num_row(), k.get_num_column());
    k_p_1.value_copy(k);
    for( i=0; i<k.get_num_row(); ++i ){
        k_p_1.get_matrix()[i][i] += 1;
    }
    My_matrix<double> s_k(hat_matrix.get_num_row(), k.get_num_column());
    trmul(hat_matrix, k_p_1, s_k);
    My_matrix<double> s_k_s(s_k.get_num_row(), hat_matrix.get_num_column());
    trmul(s_k, hat_matrix, s_k_s);
    int keep = k.get_num_column()-x.get_num_column();
    Eigen_result e_r = eigen(s_k_s, keep);
    for( i=0; i<keep; ++i ){
        e_r.get_eigen_values().get_array()[i] -= 1;
    }
    return e_r;
}
// log-likelihoods (eq. 6 from paper), this is the emma.delta.ML.LL.wo.Z function in R emma
double _ll_(const double & delta, const Eigen_result & eigen_R, const Eigen_result & eigen_L, const My_Vector<double> & sq_etas){
//    std::cout << "line 77 eigen_R.get_eigen_values().get_length() " << eigen_R.get_eigen_values().get_length() << std::endl;
    int n = eigen_L.get_eigen_values().get_length();
    My_Vector<double> v1(eigen_R.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_R.get_eigen_values(), delta, v1);
//    std::cout << "line 81" << std::endl;
    My_Vector<double> v2(eigen_L.get_eigen_values().get_length());
    sum_of_a_vector_a_number<double, double, double>(eigen_L.get_eigen_values(), delta, v2);
//    std::cout << "line 84" << std::endl;
    My_Vector<double> sq_etas_d_v1(sq_etas.get_length());
    quotient_of_two_vectors<double, double, double>(sq_etas, v1, sq_etas_d_v1);
//    std::cout << "line 87" << std::endl;
    double sum_log_v2=0.0;
    for( int i=0; i<n; ++i ){
        sum_log_v2 += log(v2.get_array()[i]);
    }
//    std::cout << "line 92" << std::endl;
    return  0.5*(n*(  log(n / (2.0 * M_PI)) - 1- log(sum_of_a_vector(sq_etas_d_v1))) - sum_log_v2);
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
    //return 0.5 * delta * (n * sum_of_a_vector(v2_2_v1) / sum_of_a_vector(v2) - sum_of_a_vector(v3)); // this funtion used in the R code of emma at the emma.MLE function, but maybe it is wrong
    return 0.5 * (n * sum_of_a_vector(v2_2_v1) / sum_of_a_vector(v2) - sum_of_a_vector(v3));
}

//newton_ll( new_opt_delta, eps, maxiter, n, eig_R_values, eig_L_values, sq_etas, p);
int newton_ll( double & new_opt_delta, const double & eps, const int & js, const int & n, const Eigen_result & eig_R,
               const Eigen_result & eig_L, const My_Vector<double> & sq_etas){ // js maximum number of iteration
    int k, l=js;
    double y2, y3, d, p, x0, x1=new_opt_delta;
    double small_number = 0.00001;
    x0=new_opt_delta;
    y2=_dll_(x0, eig_R, eig_L, sq_etas);
    y3=_dll_(x0+small_number, eig_R, eig_L, sq_etas);

    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
        x0=x1;
        if( fabs(y2)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        x1 = x0 - y2/((y3-y2)/small_number);
        y2=_dll_(x1, eig_R, eig_L, sq_etas);
        y3=_dll_(x1+small_number, eig_R, eig_L, sq_etas);
        d=fabs(x1-x0);
        --l;
    }
    new_opt_delta=0;
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

int newton_reml( double & new_opt_delta, const double & eps, const int & js, Eigen_result & eig_R, My_Vector<double> & sq_etas ){ // js maximum number of iteration
    int k, l=js;
    double y2, y3, d, p, x0, x1=new_opt_delta;

    double small_number = 0.00001;

    x0=new_opt_delta;
    y2=_redll_( x0, eig_R, sq_etas);
    y3=_redll_( x0+small_number, eig_R, sq_etas);
    d=eps+1.0;
    while( (d>=eps) && (l!=0) ){
        x0=x1;
        if( fabs(y2)+1.0 == 1.0 ){
            std::cout << "error" << std::endl;
            return -1;
        }
        x1 = x0 - y2/((y3-y2)/small_number);
        y2=_redll_( x1, eig_R, sq_etas );
        y3=_redll_( x1+small_number, eig_R, sq_etas);
        d=fabs(x1-x0);
        --l;
    }
    new_opt_delta=x0;
    k=js-l;
    return k;
}

Emma_result emma_estimates ( const My_Vector<double> & y, const My_matrix<double> & k, const My_matrix<double> & x,
                             const My_matrix<double> & xs, const Eigen_result & eigen_L, const int & ngrids,
                             const double & llim, const double &ulim, const double & eps,
                             const std::string & method, const int & maxiter ){
    double delta;
    return emma_estimates (  y, k, x, xs, eigen_L, ngrids, llim, ulim, eps, method, maxiter, delta );
}

Emma_result emma_estimates ( const My_Vector<double> & y, const My_matrix<double> & k, const My_matrix<double> & x,
                             const My_matrix<double> & xs, const Eigen_result & eigen_L, const int & ngrids,
                             const double & llim, const double &ulim, const double & eps,
                            const std::string & method, const int & maxiter, double & delta ){
//    std::cout << "  " << y.get_length() << "    " << k.get_num_column() << "    " << k.get_num_row() << "   " << x.get_num_row() << std::endl;
//    assert(y.get_length() == k.get_num_column() && y.get_length()==k.get_num_row() && y.get_length()==x.get_num_row() );
    std::cout << std::setprecision(16);
    int i, j;
    int n = y.get_length();
    int q = x.get_num_column()+xs.get_num_column();
    int p = n - q;
    int m = ngrids+1;

    My_matrix<double> full_x(x.get_num_row(), x.get_num_column()+xs.get_num_column());
    if( xs.get_num_column()>0 ){
        append_two_matrix( x, xs, full_x);
    }else{
        full_x.value_copy(x);
    }
    Eigen_result eigen_R = _get_eigen_R_(full_x, k);
    My_Vector<double> etas(p);
    double t;
    for( i=0; i < p; ++i ){
        t = 0.0;
        for( j=0; j < n; ++j ){
            t = t + eigen_R.get_eigen_vectors().get_matrix()[i][j]*y.get_array()[j];
        }
        etas.get_array()[i]=t;
    }

    My_Vector<double> sq_etas(p);
    for( i=0; i < p; ++i ){
        sq_etas.get_array()[i] = etas.get_array()[i] * etas.get_array()[i];
    }

    My_Vector<double> log_deltas(m); // the space for deltas to search
    My_Vector<double> deltas(m);
    for ( i=0; i < m; ++i ){
        log_deltas.get_array()[i] = (double(i) / ngrids)*(ulim - llim) + llim;
        deltas.get_array()[i] = exp(log_deltas.get_array()[i]);
    }

    My_Vector<double> lls(m);
    My_Vector<double> dlls(m);
    if (method.compare("REML") == 0){
        for( i=0; i<m; ++i ){
            lls.get_array()[i] = _rell_(deltas.get_array()[i], eigen_R, sq_etas);
            dlls.get_array()[i] = _redll_( deltas.get_array()[i], eigen_R, sq_etas );
        }
    }else if (method.compare("ML") == 0){
        for( i=0; i<m; ++i ){
            lls.get_array()[i] = _ll_(deltas.get_array()[i], eigen_R, eigen_L, sq_etas);
            dlls.get_array()[i] = _dll_( deltas.get_array()[i], eigen_R, eigen_L, sq_etas );
        }
    }
    int max_ll_i = argmax(lls);
    double max_ll = lls.get_array()[max_ll_i];
    double last_dll = dlls.get_array()[0];
    double  last_ll = lls.get_array()[0];

    std::map<int, double> zero_intervals;
    for (i =0; i<m; ++i){
        if (dlls.get_array()[i] < 0 && last_dll > 0){
            zero_intervals[i]=last_ll;
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
        opt_delta = deltas.get_array()[opt_i - 1];
        double new_opt_delta = opt_delta;
        if (method.compare("REML") == 0) {
//            std::cout << "line 314 " << new_opt_delta << std::endl;
            newton_reml( new_opt_delta, eps, maxiter, eigen_R, sq_etas);
//            std::cout << "line 316 " << new_opt_delta << std::endl;
        }else if (method.compare("ML") == 0){
//            std::cout << "line 336 " << new_opt_delta << std::endl;
            newton_ll( new_opt_delta, eps, maxiter, n, eigen_R, eigen_L, sq_etas);
//            std::cout << "line 338 " << new_opt_delta << std::endl;
        }
        if (opt_i > 0 && deltas.get_array()[opt_i - 1] - eps < new_opt_delta && new_opt_delta < deltas.get_array()[opt_i] + eps){
            opt_delta = new_opt_delta;
        }else if (opt_i == 0 && 0.0 < new_opt_delta && new_opt_delta < deltas.get_array()[opt_i] + eps){
            opt_delta = new_opt_delta;
        }else if(opt_i == m - 1 && new_opt_delta > deltas.get_array()[opt_i - 1] - eps ){
            opt_delta = new_opt_delta;
        }
        if (method.compare("REML") == 0) {
            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
        }else if (method.compare("ML") == 0) {
            opt_ll = _ll_(opt_delta, eigen_R, eigen_L, sq_etas);
        }
//        std::cout << "line 352 " << opt_delta << std::endl;
        if (opt_ll < max_ll){
            opt_delta = deltas.get_array()[max_ll_i];
            if (method.compare("REML") == 0) {
                opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
            }else if (method.compare("ML") == 0) {
                opt_ll = _ll_(opt_delta, eigen_R, eigen_L, sq_etas);
            }
        }
    }else{
        opt_delta = deltas.get_array()[max_ll_i];
        opt_ll = max_ll;
    }
    delta = opt_delta;
    std::cout << "line 404 " << delta << " " << opt_ll << std::endl;

    My_Vector<double> R(p);
    for ( i=0; i<p; ++i ) {
        R.get_array()[i] = sq_etas.get_array()[i] / (eigen_R.get_eigen_values().get_array()[i] + opt_delta);
    }
    double opt_vg = 0.0;
    for ( i=0; i<p; ++i ){
        opt_vg += R.get_array()[i];
    }
    opt_vg=opt_vg/p; //This is the REML estimation. the ML estimation is np.sum(R) / n


    My_matrix<double> H_sqrt_inv_pre (n, n);
    H_sqrt_inv_pre.set_values_zero();
    for( i=0; i<n; ++i ){
        double temp = 1.0/sqrt(eigen_L.get_eigen_values().get_array()[i]+opt_delta);
        H_sqrt_inv_pre.get_matrix()[i][i] = temp; //this is from the R code of emma
    }
    My_matrix<double> H_sqrt_inv (n, n);
    My_matrix<double> temp(n, n);
    T_matrix(eigen_L.get_eigen_vectors(), temp);
    trmul( H_sqrt_inv_pre, eigen_L.get_eigen_vectors(), H_sqrt_inv );

    My_matrix<double> X_t(n, q);
    trmul(H_sqrt_inv, full_x, X_t);
    My_Vector<double> Y_t (n);
    trmul(H_sqrt_inv, y, Y_t);

    My_Vector<double> beta_est(q);
    lsq(X_t, Y_t, beta_est);
    My_Vector<double> x_beta(n);
    for( i=0; i<n; ++i ){
        x_beta.get_array()[i]=0;
        for( j=0; j<q; ++j ){
            x_beta.get_array()[i] += X_t.get_matrix()[i][j]*beta_est.get_array()[j];
        }
    }
    double mahalanobis_rss = 0.0;
    for ( i=0; i<n; ++i ){
        mahalanobis_rss += pow( Y_t.get_array()[i]-x_beta.get_array()[i], 2);
    }
//    std::cout << "mahalanobis_rss " << mahalanobis_rss << std::endl;
    if( xs.get_num_column() > 0 ){
        My_matrix<double> h0_X (n, x.get_num_column());
        trmul(H_sqrt_inv, x, h0_X);
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
        double f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.get_num_column();
        int freedome1 = xs.get_num_column();
        double p_val = sf(f_stat, freedome1, p);
        Emma_result r_re(Y_t, H_sqrt_inv, mahalanobis_rss, p_val);
        return r_re;
    }
    Emma_result r_re(Y_t, H_sqrt_inv, mahalanobis_rss);
    return r_re;
}

