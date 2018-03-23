//
// Created by Baoxing song on 12.03.18.
//

#include "Average_information_impl.h"

Average_information_impl::Average_information_impl(const My_Vector<double> & _y, const std::vector<My_matrix<double>> & _ks,
                         const My_matrix<double> & _x, const My_matrix<double> & _xs, double & _eps,
                         const int & _maxiter){
    My_Vector<double> y = _y;
    std::vector<My_matrix<double>> ks= _ks;
    My_matrix<double> x = _x;
    My_matrix<double> xs = _xs;
    double eps = _eps;
    int maxiter = _maxiter;
    V = My_matrix<double>(_y.get_length(), _y.get_length());
}

double Average_information_impl::_ll_(){

}

double Average_information_impl::_dll_( ){

}

double Average_information_impl::get_var_y(){
    double sum_y=0.0;
    int i;
    for( i=0; i<y.get_length(); ++i ){
        sum_y += y.get_array()[i];
    }
    double mean_y = sum_y / y.get_length();
    double var_y=0.0;
    for( i=0; i<y.get_length(); ++i ){
        var_y += pow(y.get_array()[i]-mean_y, 2 );
    }
    return var_y;
}

void Average_information_impl::update_V( const My_Vector<double> & etas ){
    int i, j, k;
    for( i = 0; i<y.get_length(); ++i ){
        for( j = 0; j<y.get_length(); ++j  ){
            V.get_matrix()[i][j] = 0;
            for( k=0; k<ks.size(); ++k ){
                V.get_matrix()[i][j] += etas.get_array()[k] * ks[k].get_matrix()[i][j];
            }
        }
    }
    for( i = 0; i<y.get_length(); ++i ){
        V.get_matrix()[i][i] += etas.get_array()[ks.size()];
    }
}


void Average_information_impl::get_ai_estimates ( ){


    int i, j;
    int n = y.get_length();
    int q = x.get_num_column()+xs.get_num_column();
    int p = n - q;

    My_matrix<double> full_x(x.get_num_row(), x.get_num_column()+xs.get_num_column());
    if( xs.get_num_column()>0 ){
        append_two_matrix( x, xs, full_x);
    }else{
        full_x=x;
    }
//    Eigen_result* eigen_Rs = new Eigen_result[ks.size()];
//    for ( i=0; i<ks.size(); ++i ){
//        eigen_Rs[i] = _get_eigen_R_(full_x, ks[i]);
//    }
    double eta_p=this->get_var_y();

    My_Vector<double> etas(ks.size()+1);
    double initial_etas = eta_p/(ks.size()+2);// here we think fixed terms explain a peace of phenotypic variance,
    // which is different with the  SVS document http://doc.goldenhelix.com/SVS/latest/svsmanual/mixedModelMethods/overview.html#overviewofmixedlinearmodels
    // the GCTA paper did not explain how to initialize the error term https://www.sciencedirect.com/science/article/pii/S0002929710005987?via%3Dihub
    // but it does not matter, since the following iterats will update that
    for( i=0; i<ks.size(); ++i ){
        etas.get_array()[i]=initial_etas;
    }


//|V| means det(V)

//
//    //todo break
//
//    double t;
//    for( i=0; i < p; ++i ){
//        t = 0.0;
//        for( j=0; j < n; ++j ){
//            t = t + eigen_R.get_eigen_vectors().get_matrix()[i][j]*y.get_array()[j];
//        }
//        etas.get_array()[i]=t;
////        std::cout << "etas " << etas.get_array()[i] << std::endl;
//    }
////    for( i=0; i < p; ++i ){
////        std::cout << "etas " << i << " " << etas.get_array()[i] << std::endl;
////    }
//    My_Vector<double> sq_etas(p);
//    for( i=0; i < p; ++i ){
//        sq_etas.get_array()[i] = etas.get_array()[i] * etas.get_array()[i];
////        std::cout << "sq_etas " << sq_etas.get_array()[i] << std::endl;
//    }
////    i=0;
////    for( j=0; j < n; ++j ){
////        std::cout << "eigen " << eigen_R.get_eigen_vectors().get_matrix()[i][j] << std::endl;
////    }
////    for( j=0; j < n; ++j ){
////        std::cout << "phenotype " << y.get_array()[j] << std::endl;
////    }
//
//    My_Vector<double> log_deltas(m); // the space for deltas to search
//    My_Vector<double> deltas(m);
//    for ( i=0; i < m; ++i ){
//        log_deltas.get_array()[i] = (double(i) / ngrids)*(ulim - llim) + llim;
//        deltas.get_array()[i] = exp(log_deltas.get_array()[i]);
//    }
//
//    My_matrix<double> lambdas(p, m);
//    for ( i=0; i < p; ++i ){
//        for ( j=0; j < m; ++j ){
//            lambdas.get_matrix()[i][j] = eigen_R.get_eigen_values().get_array()[i] + deltas.get_array()[j];
//        }
//    }
//
//    My_Vector<double> s1(m);
//    for ( i=0; i < m; ++i ){
//        s1.get_array()[i] = 0.0;
//        for ( j=0; j < p; ++j ){
//            s1.get_array()[i] += sq_etas.get_array()[j]/lambdas.get_matrix()[j][i];
//        }
////        std::cout << "line 204 " << i << " " << s1.get_array()[i] << std::endl;
//    }
//
//    My_Vector<double> s3(m);
//    for( i=0; i<m; ++i ){
//        s3.get_array()[i] = 0.0;
//        for ( j=0; j<p; ++j ){
//            s3.get_array()[i]+=sq_etas.get_array()[j]/(lambdas.get_matrix()[j][i]*lambdas.get_matrix()[j][i]);
//        }
////        std::cout << "s3  " << i << " " << s3.get_array()[i] << std::endl;
//    }
//
//    My_Vector<double> lls(m);
//    My_Vector<double> dlls(m);
//    if (method.compare("REML") == 0){
//        My_Vector<double> s2(m);
//        for( i=0; i<m; ++i ){
//            s2.get_array()[i]=0.0;
//            for( j=0; j<p; ++j  ){
//                s2.get_array()[i] += log(lambdas.get_matrix()[j][i]);
//            }
////            std::cout << "line 284 " << s2.get_array()[i] << std::endl;
//        }
//        for( i=0; i<m; ++i ){
////            std::cout << "PI " << M_PI << " " << s1.get_array()[i] << " " << log(s1.get_array()[i]) << " " << s2.get_array()[i] << std::endl;
//            lls.get_array()[i] = 0.5 * (p * (log((p) / (2.0 * M_PI)) - 1 - log(s1.get_array()[i])) - s2.get_array()[i]);
////            std::cout << "line 288 " << lls.get_array()[i] << std::endl;
//        }
//        My_Vector<double> s4(m);
//        for( i=0; i<m; ++i ){
//            s4.get_array()[i] = 0.0;
//            for( j=0; j<p; ++j ){
//                s4.get_array()[i] += 1.0/lambdas.get_matrix()[j][i];
//            }
////            std::cout << "s4 " << s4.get_array()[i] << std::endl;
//        }
//
//        for( i=0; i<m; ++i ){
//            dlls.get_array()[i] = 0.5 * (p * s3.get_array()[i] / s1.get_array()[i] - s4.get_array()[i]);
////            std::cout << "dlls " << dlls.get_array()[i] << " p " << p << " s3.get_array()[i] " << s3.get_array()[i] << "  s1.get_array()[i] " <<  s1.get_array()[i] << " s4.get_array()[i] " << s4.get_array()[i] << std::endl;
//        }
//    }else if (method.compare("ML") == 0){
//        My_matrix<double> xis(n, m);
//        for ( i=0; i < n; ++i ){
//            for ( j=0; j < m; ++j ){
//                xis.get_matrix()[i][j] = eigen_L.get_eigen_values().get_array()[i] + deltas.get_array()[j];
//            }
//        }
//        My_Vector<double> s2(m);
//        for( i=0; i<n; ++i ){
//            s2.get_array()[i]=0.0;
//            for( j=0; j<m; ++j  ){
//                s2.get_array()[j] += log(xis.get_matrix()[i][j]);
//            }
//        }
//        for( i=0; i<m; ++i ){
//            lls.get_array()[i] = 0.5 * (p * (log((p) / (2.0 * M_PI)) - 1 - log(s1.get_array()[i])) - s2.get_array()[i]);
////            std::cout << lls.get_array()[i] << std::endl;
//        }
//        My_Vector<double> s4(m);
//        for( i=0; i<n; ++i ){
//            s4.get_array()[i] = 0.0;
//            for( j=0; j<m; ++j ){
//                s4.get_array()[j]=1.0/xis.get_matrix()[i][j];
//            }
//        }
//        for( i=0; i<m; ++i ){
//            dlls.get_array()[i] = 0.5 * (n*s3.get_array()[i]/s1.get_array()[i]-s4.get_array()[i]);
//        }
//    }
//
//    int max_ll_i = argmax(lls);
//    double max_ll = lls.get_array()[max_ll_i];
//    double last_dll = dlls.get_array()[0];
//    double  last_ll = lls.get_array()[0];
//
////    std::cout << "max_ll_i " << max_ll_i << std::endl;
////    std::cout << "max_ll " << max_ll << std::endl;
////    std::cout << "last_dll " << last_dll << std::endl;
////    std::cout << "last_ll " << last_ll << std::endl;
//
//    std::map<int, double> zero_intervals;
//    for (i =0; i<m; ++i){
//        if (dlls.get_array()[i] < 0 && last_dll > 0){
//            zero_intervals[i]=(lls.get_array()[i] + last_ll) * 0.5;
//        }
//        last_ll = lls.get_array()[i];
//        last_dll = dlls.get_array()[i];
//    }
//
//    double opt_delta;
//    double opt_ll;
//    if (zero_intervals.size() > 0){
//        opt_ll = -DBL_MAX;
//        int opt_i;
//        for(std::map<int, double>::iterator it=zero_intervals.begin();it!=zero_intervals.end();++it){
//            if( it->second > opt_ll ){
//                opt_ll = it->second;
//                opt_i = it->first;
//            }
//        }
//        opt_delta = 0.5 * (deltas.get_array()[opt_i - 1] + deltas.get_array()[opt_i]);
//        std::cout << "line 344 opt_delta " << opt_delta << std::endl;
//        double * new_opt_delta;
//        *new_opt_delta = opt_delta;
//        if (method.compare("REML") == 0) {
//            newton_reml( new_opt_delta, eps, maxiter, opt_delta, eigen_R, sq_etas);
//        }else if (method.compare("ML") == 0){
//            newton_ll( new_opt_delta, eps, maxiter, n, eigen_R, eigen_L, sq_etas);
//        }
//        if (opt_i > 1 && deltas.get_array()[opt_i - 1] - eps < *new_opt_delta && *new_opt_delta < deltas.get_array()[opt_i] + eps){
//            opt_delta = *new_opt_delta;
//            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
//        }else if (opt_i == 1 && 0.0 < *new_opt_delta && *new_opt_delta < deltas.get_array()[opt_i] + eps){
//            opt_delta = *new_opt_delta;
//            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
//        }else if(opt_i == m - 1 && *new_opt_delta > deltas.get_array()[opt_i - 1] - eps ){
//            opt_delta = *new_opt_delta;
//            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
//        }
//        if (method.compare("REML") == 0) {
//            opt_ll = _rell_( opt_delta, eigen_R, sq_etas);
//        }else if (method.compare("ML") == 0) {
//            opt_ll = _ll_(opt_delta, eigen_R, eigen_L, sq_etas);
//        }
//        if (opt_ll < max_ll){
//            opt_delta = deltas.get_array()[max_ll_i];
//        }
//    }else{
//        opt_delta = deltas.get_array()[max_ll_i];
//        opt_ll = max_ll;
//    }
////    std::cout << "line 418 opt_delta " << opt_delta << std::endl;
////    std::cout << "line 419 opt_ll " << opt_ll << std::endl;
//
//
//    My_matrix<double> R(p, p);
//    for ( j=0; j<p; ++j ){
//        for ( i=0; i<p; ++i ) {
//            R.get_matrix()[i][j] = sq_etas.get_array()[i] / (eigen_R.get_eigen_values().get_array()[j] + opt_delta);
//
//        }
//
////        std::cout << i << " " << "sq_etas.get_array()[i] " << sq_etas.get_array()[i]
////                  << " eigen_R.get_eigen_values().get_array()[i] "
////                  << eigen_R.get_eigen_values().get_array()[i] << " opt_delta " << opt_delta << " R.get_array()[i] "
////                  << R.get_array()[i] << std::endl;
//    }
//    double opt_vg = 0.0;
//    for ( i=0; i<p; ++i ){
//        for ( j=0; j<p; ++j ) {
////            std::cout << i << " " << j << " " << R.get_matrix()[i][j] << " " << std::endl;
//            opt_vg += R.get_matrix()[i][j];
//        }
////        std::cout << std::endl;
//    }
//    opt_vg=opt_vg/p; //This is the REML estimation. the ML estimation is np.sum(R) / n
////    std::cout << "opt_vg " << opt_vg << std::endl;
//    double opt_ve = opt_vg * opt_delta;
////    std::cout << "opt_ve " << opt_ve << std::endl;
//
//    My_matrix<double> H_sqrt_inv (n, n);
//    for( i=0; i<n; i++ ){
//        double this_value = 1.0/sqrt(eigen_L.get_eigen_values().get_array()[i]+opt_delta);
//        for( j=0; j<n; j++ ){
//            H_sqrt_inv.get_matrix()[i][j] = this_value*eigen_L.get_eigen_vectors().get_matrix()[i][j];
//        }
//    }
////    for( i=0; i<n; i++ ){
////        for( j=0; j<n; j++ ){
////            std::cout << H_sqrt_inv.get_matrix()[i][j] << " ";
////        }
////        std::cout << std::endl;
////    }
//    My_matrix<double> X_t(n, q);
//    trmul(H_sqrt_inv, full_x, X_t);
//    My_Vector<double> Y_t (n);
//    for(i=0; i<n; ++i){
//        Y_t.get_array()[i] = 0;
//        for(j=0; j<n; ++j) {
//            Y_t.get_array()[i] += H_sqrt_inv.get_matrix()[i][j] * y.get_array()[j];
//        }
////        std::cout << " Y_t " << Y_t.get_array()[i] << std::endl;
//    }
////    std::cout << "472" << std::endl;
//    My_Vector<double> beta_est(q);
////    std::cout << "474" << std::endl;
//    lsq(X_t, Y_t, beta_est);
////    std::cout << "475" << std::endl;
//    My_Vector<double> x_beta(n);
//    for( i=0; i<n; ++i ){
//        x_beta.get_array()[i]=0;
//        for( j=0; j<q; ++j ){
//            x_beta.get_array()[i] += X_t.get_matrix()[i][j]*beta_est.get_array()[j];
//        }
////        std::cout << " x_beta " << x_beta.get_array()[i] << " Y_t " << Y_t.get_array()[i] << " Y_t - x_beta " << Y_t.get_array()[i]-x_beta.get_array()[i] << std::endl;
//    }
//    double mahalanobis_rss = 0.0;
//    for ( i=0; i<n; ++i ){
//        mahalanobis_rss += pow( Y_t.get_array()[i]-x_beta.get_array()[i], 2);
////        std::cout << Y_t.get_array()[i] << "  " << x_beta.get_array()[i] << std::endl;
//    }
}