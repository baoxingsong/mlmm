//
// Created by Baoxing song on 06.06.18.
//

#include "mt_gemma_impl.h"

//this function related with the
double EigenProc(const My_matrix<double> & V_g, const My_matrix<double> & V_e, const size_t & d_size,
                 My_Vector<double> & D_l,
                 My_matrix<double> & UltVeh, My_matrix<double> & UltVehi) { // here h mean half and i means -1
    double d, logdet_Ve = 0.0;

    // Eigen decomposition of V_e.

    My_matrix<double> V_e_h(d_size, d_size);
    My_matrix<double> V_e_hi(d_size, d_size);

    Eigen_result eigen_l = eigen(V_e);

    // the original paper did not give the mathematical theory underlying this part.
    // i.e. A = U*D*U^T
    // let B = U*D^(1/2)*U^T
    // then B*B = U*D^(1/2)*U^T * U*D^(1/2)*U^T =  U*D^(1/2)*D^(1/2)*U^T = U*D^(1/2)*U^T using the feature that U^T*U=I

    // Calculate V_e_h and V_e_hi.
    V_e_h.set_values_zero();
    V_e_hi.set_values_zero();
    size_t i, j, k;
    for (i = 0; i < d_size; ++i) {
        d = eigen_l.get_eigen_values().get_array()[i];
        if (d <= 0) { //why not test negative value, // we have some problem caused from here
            continue;
        }
        logdet_Ve += log(d);
        d = sqrt(d);
        //gsl_vector_view U_col = gsl_matrix_column(U_l, i);
        //gsl_blas_dsyr(CblasUpper, d, &U_col.vector, V_e_h);
        //gsl_blas_dsyr(CblasUpper, d, &U_col.vector, V_e_hi);
        for ( j = 0; j < d_size; ++j) {
            for( k = 0; k < d_size; ++k ){
                V_e_h.get_matrix()[j][k] += eigen_l.get_eigen_vectors().get_matrix()[i][j]*
                                            eigen_l.get_eigen_vectors().get_matrix()[i][k]*d;
            }
            // since my function return U^t, So the index should be adapted
        }
//        d = 1.0/d;
        for ( j = 0; j < d_size; ++j) {
            for( k = 0; k < d_size; ++k ){
                V_e_hi.get_matrix()[j][k] += eigen_l.get_eigen_vectors().get_matrix()[i][j]*
                                             eigen_l.get_eigen_vectors().get_matrix()[i][k]/d;
            }
            // since my function return U^t, So the index should be adapted
        }
    }
//
//    for ( i = 0; i < d_size; ++i) {
//        for( j = 0; j < d_size; ++j ){
//            std::cout << "i: " << i << " j: " << j << " V_e_h.get_matrix()[i][j] " << V_e_h.get_matrix()[i][j] << std::endl;
//            std::cout << "i: " << i << " j: " << j << " V_e_hi.get_matrix()[i][j] " << V_e_hi.get_matrix()[i][j] << std::endl;
//        }
//    }

    My_matrix<double> Lambda(d_size, d_size);
    My_matrix<double> VgVehi(d_size, d_size);
//    My_matrix<double> U_l(d_size, d_size);
    // Calculate Lambda=V_ehi V_g V_ehi.
    trmul(V_g, V_e_hi, VgVehi);
    trmul(V_e_hi, VgVehi, Lambda);
    // Eigen decomposition of Lambda.
    eigen_l = eigen(Lambda);
    //T_matrix(eigen_l.get_eigen_vectors(), U_l);
    // Calculate UltVeh and UltVehi.
//    T_matrix(eigen_l.get_eigen_vectors(), U_l);
    trmul(eigen_l.get_eigen_vectors(), V_e_h, UltVeh);
    trmul(eigen_l.get_eigen_vectors(), V_e_hi, UltVehi);
    for( i=0; i<d_size; ++i ){
        D_l.get_array()[i] = eigen_l.get_eigen_values().get_array()[i];
    }
    return logdet_Ve;
}

// Qi=(\sum_{k=1}^n x_kx_k^T\otimes(delta_k*Dl+I)^{-1} )^{-1}.
double CalcQi(const Eigen_result & eigen_r, const My_Vector<double> & D_l, const size_t & n_size, const size_t & d_size,
              const size_t & c_size,
              const My_matrix<double> & X, My_matrix<double> & Qi) {

    double delta, dl, d1, d2, d, logdet_Q;

    My_matrix<double> Q(Qi.get_num_row(), Qi.get_num_row());
    Q.set_values_zero();
    size_t i, j, k, l;
    for (i = 0; i < c_size; ++i) {
        for (j = 0; j < c_size; ++j) {
            for (l = 0; l < d_size; ++l) {
                dl = D_l.get_array()[l];
                if (j < i) {
                    d = Q.get_matrix()[j * d_size + l][i * d_size + l];
                } else {
                    d = 0.0;
                    for ( k = 0; k < n_size; ++k) {
                        d1 = X.get_matrix()[i][k];
                        d2 = X.get_matrix()[j][k];
                        delta = eigen_r.get_eigen_values().get_array()[k];
                        d += d1 * d2 / (dl * delta + 1.0);
                        // this function was explained on Xiang Zhou and Matthew Stephens (2012).
                        // Genome-wide efficient mixed-model analysis for association studies. Nature Genetics 44, 821–824.
                    }
                }
                Q.get_matrix()[i * d_size + l][j * d_size + l]=d;
            }
        }
    }
    // Calculate LU decomposition of Q, and invert Q and calculate |Q|.
    logdet_Q = log(determinant(Q));
    Qi.value_copy(Q);
    inverse_matrix(Qi);
    return logdet_Q;
}

// xHiy=\sum_{k=1}^n x_k\otimes ((delta_k*Dl+I)^{-1}Ul^TVe^{-1/2}y.
// this function was explained in Xiang Zhou and Matthew Stephens (2012). Genome-wide efficient mixed-model analysis for
// association studies. Nature Genetics 44, 821–824.
void CalcXHiY(const Eigen_result & eigen_r, const My_Vector<double> & D_l,
              const My_matrix<double> & X, const My_matrix<double> & UltVehiY, const size_t & n_size,
              const size_t & d_size, const size_t & c_size,
              My_Vector<double> & xHiy) {

    xHiy.set_values_zero();
    double x, delta, dl, y, d;
    size_t i, j, k;
    for (i = 0; i < d_size; ++i) {
        dl = D_l.get_array()[i];
        for (j = 0; j < c_size; ++j) {
            d = 0.0;
            for ( k = 0; k < n_size; ++k) {
                x = X.get_matrix()[j][k];
                y = UltVehiY.get_matrix()[i][k];
                delta = eigen_r.get_eigen_values().get_array()[k];
                d += x * y / (delta * dl + 1.0);
            }
            xHiy.get_array()[j * d_size + i]=d;
        }
    }
}

// OmegaU=D_l/(delta Dl+I)^{-1}
// OmegaE=delta D_l/(delta Dl+I)^{-1}
void CalcOmega(const Eigen_result & eigen_r, const My_Vector<double> & D_l, const size_t & n_size, const size_t & d_size,
               My_matrix<double> &OmegaU, My_matrix<double> &OmegaE) {
    double delta, dl, d_u, d_e;
    size_t k, i;
    for (k = 0; k < n_size; ++k) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        for (i = 0; i < d_size; ++i) {
            dl = D_l.get_array()[i];

            d_u = dl / (delta * dl + 1.0);  // @@
            d_e = delta * d_u;
            OmegaU.get_matrix()[i][k] = d_u;
            OmegaE.get_matrix()[i][k] = d_e;
        }
    }
}

void UpdateU(const My_matrix<double> &OmegaE, const My_matrix<double> & UltVehiY, const size_t & n_size,
             const size_t & d_size,
             const My_matrix<double> &UltVehiBX, My_matrix<double> &UltVehiU) {
    size_t i, j;
    for( i=0; i<d_size; ++i ){
        for( j=0; j<n_size; ++j ){
            UltVehiU.get_matrix()[i][j] =
                    (UltVehiY.get_matrix()[i][j] - UltVehiBX.get_matrix()[i][j])*OmegaE.get_matrix()[i][j];
        }
    }
}

void UpdateE(const My_matrix<double> &UltVehiY, const My_matrix<double> &UltVehiBX, const size_t & n_size,
             const size_t & d_size, const My_matrix<double> &UltVehiU, My_matrix<double> &UltVehiE) {

    size_t i, j;
    for( i=0; i<d_size; ++i ){
        for( j=0; j<n_size; ++j ){
            UltVehiE.get_matrix()[i][j] = (UltVehiY.get_matrix()[i][j] -
                                          UltVehiBX.get_matrix()[i][j]) - UltVehiU.get_matrix()[i][j];
        }
    }
}

void UpdateL_B(const My_matrix<double> &X, const My_matrix<double> &XXti, const size_t & n_size, const size_t & d_size,
               const size_t & c_size,
               const My_matrix<double> &UltVehiY, const My_matrix<double> &UltVehiU,
               My_matrix<double> &UltVehiBX, My_matrix<double> &UltVehiB) {  // (Y-G)/X = beta
    My_matrix<double> YUX(d_size, c_size);

    size_t i, j;
    for( i=0; i<d_size; ++i ){
        for( j=0; j<n_size; ++j ){
            UltVehiBX.get_matrix()[i][j] = UltVehiY.get_matrix()[i][j] - UltVehiU.get_matrix()[i][j];
        }
    }
    My_matrix<double> X_t(X.get_num_column(), X.get_num_row());
    T_matrix(X, X_t);
    trmul(UltVehiBX, X_t, YUX);
    trmul(YUX, XXti, UltVehiB);
}

void UpdateRL_B(const My_Vector<double> & xHiy, const My_matrix<double> & Qi, const size_t & n_size,
                const size_t & d_size, const size_t & c_size, const size_t & dc_size,
                My_matrix<double> & UltVehiB) { //QI with size dc * dc  xHiy with size dc

    My_Vector<double> b(dc_size);
    size_t i, j;
    for( i=0; i<dc_size; ++i ){
        b.get_array()[i] = 0;
        for( j=0; j<dc_size; ++j ){
            b.get_array()[i] += Qi.get_matrix()[i][j]*xHiy.get_array()[j];
        }
    }

    // Copy b to UltVehiB.
    for (i = 0; i < c_size; ++i) {
        for( j=0; j<d_size; ++j ){
            UltVehiB.get_matrix()[j][i] = b.get_array()[i * d_size + j];
        }
    }
}

void UpdateV(const Eigen_result & eigen_r, const My_matrix<double> & U, const My_matrix<double> & E,
             const size_t & n_size, const size_t & d_size, const size_t & c_size,
             const My_matrix<double> & Sigma_uu, const My_matrix<double> & Sigma_ee,
             My_matrix<double> & V_g, My_matrix<double> & V_e) {

    V_g.set_values_zero();
    V_e.set_values_zero();
    double delta;
    size_t k, i, j;
    // Calculate the first part: UD^{-1}U^T and EE^T.
    for (k = 0; k < n_size; ++k) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        if (delta == 0) {
            continue;
        }
       // gsl_vector_const_view U_col = gsl_matrix_const_column(U, k);
        //gsl_blas_dsyr(CblasUpper, 1.0 / delta, &U_col.vector, V_g);
        for( i=0; i<d_size; ++i ){
            for( j=0; j<d_size; ++j ){
                V_g.get_matrix()[i][j] += U.get_matrix()[i][k]*U.get_matrix()[j][k]/delta;
            }
        }
    }
    My_matrix<double> E_t(E.get_num_column(), E.get_num_row());
    T_matrix(E, E_t);
    trmul(E, E_t, V_e);

    // Add Sigma.
    for( i=0; i<d_size; ++i ){
        for( j=0; j<d_size; ++j ){
            V_g.get_matrix()[i][j] += Sigma_uu.get_matrix()[i][j];
            V_e.get_matrix()[i][j] += Sigma_ee.get_matrix()[i][j];
        }
    }
    for( i=0; i<d_size; ++i ) {
        for (j = 0; j < d_size; ++j) {
            V_g.get_matrix()[i][j] = V_g.get_matrix()[i][j]/(double)n_size;
            V_e.get_matrix()[i][j] = V_e.get_matrix()[i][j]/(double)n_size;
        }
    }
}

void CalcSigma(const char func_name, const Eigen_result & eigen_r,
               const size_t & n_size, const size_t & d_size, const size_t & c_size, const size_t & dc_size,
               const My_Vector<double> & D_l, const My_matrix<double> & X,
               const My_matrix<double> & OmegaU, const My_matrix<double> & OmegaE,
               const My_matrix<double> & UltVeh, const My_matrix<double> & Qi,
               My_matrix<double> & Sigma_uu, My_matrix<double> & Sigma_ee) {

    Sigma_uu.set_values_zero();
    Sigma_ee.set_values_zero();

    double delta, dl, x, d;

    // Calculate the first diagonal term.
    size_t k, i, j;
    for (k = 0; k < n_size; ++k) {
        for (i = 0; i < d_size; ++i) {
            Sigma_uu.get_matrix()[i][i] += OmegaU.get_matrix()[i][k];
            Sigma_ee.get_matrix()[i][i] += OmegaE.get_matrix()[i][k];
        }
    }
    // Calculate the second term for REML.
    if (func_name == 'R' || func_name == 'r') {
        My_matrix<double> M_u(dc_size, d_size);
        My_matrix<double> M_u_t(d_size, dc_size);

        My_matrix<double> M_e(dc_size, d_size);
        My_matrix<double> M_e_T(d_size, dc_size);
        My_matrix<double> QiM(dc_size, d_size);
        My_matrix<double> M_u_t_QiM(d_size, d_size);
        My_matrix<double> M_e_QiM(d_size, d_size);

        M_u.set_values_zero();
        M_e.set_values_zero();
        for (k = 0; k < n_size; ++k) {
            delta = eigen_r.get_eigen_values().get_array()[k];//gsl_vector_get(eval, k);

            for ( i = 0; i < d_size; ++i) {
                dl = D_l.get_array()[i];//gsl_vector_get(D_l, i);
                for ( j = 0; j < c_size; ++j) {
                    x = X.get_matrix()[j][k];//gsl_matrix_get(X, j, k);
                    d = x / (delta * dl + 1.0);
                    M_e.get_matrix()[j * d_size + i][i]=d;
                    M_u.get_matrix()[j * d_size + i][i]=d * dl;
                }
            }
            trmul(Qi, M_u, QiM);

            T_matrix(M_u, M_u_t);
            trmul(M_u_t, QiM, M_u_t_QiM);
            for( i = 0; i < d_size; ++i ){
                for ( j = 0; j < d_size; ++j) {
                    Sigma_uu.get_matrix()[i][j] += delta*M_u_t_QiM.get_matrix()[i][j];
                }
            }
            trmul(Qi, M_e, QiM);
            T_matrix(M_e, M_e_T);
            trmul(M_e_T, QiM, M_e_QiM);
            for( i = 0; i < d_size; ++i ){
                for ( j = 0; j < d_size; ++j) {
                    Sigma_ee.get_matrix()[i][j] += M_e_QiM.get_matrix()[i][j];
                }
            }
        }
    }

    // Multiply both sides by VehUl.
    My_matrix<double> M(d_size, d_size);

    trmul(Sigma_uu, UltVeh, M);
    My_matrix<double> UltVeh_t(UltVeh.get_num_column(), UltVeh.get_num_row());
    T_matrix(UltVeh, UltVeh_t);
    trmul(UltVeh_t, M, Sigma_uu);
    trmul(Sigma_ee, UltVeh, M);
    trmul(UltVeh_t, M, Sigma_ee);
}

// 'R' for restricted likelihood and 'L' for likelihood.
// 'R' update B and 'L' don't.
// only calculate -0.5*\sum_{k=1}^n|H_k|-0.5yPxy.
double MphCalcLogL(const Eigen_result & eigen_r, const My_Vector<double> & xHiy,
                   const size_t & n_size, const size_t & d_size, const size_t & c_size, const size_t & dc_size,
                   const My_Vector<double> & D_l, const My_matrix<double> & UltVehiY,
                   const My_matrix<double> & Qi) {

    double logl = 0.0, delta, dl, y, d;

    // Calculate yHiy+log|H_k|.
    size_t k;
    size_t i;
    for ( k = 0; k < n_size; ++k) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        for ( i = 0; i < d_size; ++i) {
            y = UltVehiY.get_matrix()[i][k];
            dl = D_l.get_array()[i];
            d = delta * dl + 1.0;
            logl += y * y / d + log(d);
            //std::cout << "k: " << k << " i: " <<  i << " y: " << y << " dl: " << dl << " d: " << d << " logl: " << logl << std::endl;
        }
    }
    // Calculate the rest of yPxy.
    My_Vector<double> Qiv(dc_size);
    Qiv.set_values_zero();
    for( k=0; k<dc_size; ++k ){
        for( i=0; i<dc_size; ++i ){
            Qiv.get_array()[k] += Qi.get_matrix()[k][i]*xHiy.get_array()[i];
        }
    }

    d=0.0;
    for( k=0; k<dc_size; ++k ){
        d += xHiy.get_array()[k]*Qiv.get_array()[k];
    }
    logl -= d;
    return -0.5 * logl;
}

// Y is a dxn matrix, X is a cxn matrix, B is a dxc matrix, V_g is a
// dxd matrix, V_e is a dxd matrix, eval is a size n vector
//'R' for restricted likelihood and 'L' for likelihood.
double MphEM(const char func_name, const size_t & max_iter, const double & max_prec,
             const Eigen_result & eigen_r, const My_matrix<double> & X, const My_matrix<double> & Y,
             const size_t & n_size, const size_t & d_size, const size_t & c_size,
             My_matrix<double> & U_hat, My_matrix<double> & E_hat, My_matrix<double> & OmegaU,
             My_matrix<double> & OmegaE, My_matrix<double> & UltVehiY, My_matrix<double> & UltVehiBX,
             My_matrix<double> & UltVehiU, My_matrix<double> & UltVehiE, My_matrix<double> & V_g,
             My_matrix<double> & V_e, My_matrix<double> & B) {
    const size_t dc_size = d_size * c_size;

    My_matrix<double> Xt(n_size, c_size);
    My_matrix<double> XXt(c_size, c_size);
    My_Vector<double> D_l(d_size);
    My_matrix<double> UltVeh(d_size, d_size);
    My_matrix<double> UltVehi(d_size, d_size);
    My_matrix<double> UltVehiB(d_size, c_size);
    My_matrix<double> Qi(dc_size, dc_size);
    My_matrix<double> Sigma_uu(d_size, d_size);
    My_matrix<double> Sigma_ee(d_size, d_size);
    My_Vector<double> xHiy(dc_size);

    double logl_const = 0.0, logl_old = 0.0, logl_new = 0.0;
    double logdet_Q, logdet_Ve;

    // Calculate |XXt| and (XXt)^{-1}.
    T_matrix(X, Xt);
    trmul(X, Xt, XXt);
    My_matrix<double> XXti(XXt);
    inverse_matrix(XXti);
//    std::cout << "line 415" << std::endl;
//
//    for( int i=0; i<c_size; ++i ){
//        for( int j=0; j<c_size; ++j ){
//            std::cout << "i " << i << " j " << j << " XXt " << XXt.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " XXti " << XXti.get_matrix()[i][j] << std::endl;
//        }
//    }

    // Calculate the constant for logl.
    if (func_name == 'R' || func_name == 'r') {
        logl_const =
                -0.5 * (double)(n_size - c_size) * (double)d_size * log(2.0 * M_PI) +
                0.5 * (double)d_size * log(fabs(determinant(XXt)));//LULndet(XXt);
    } else {
        logl_const = -0.5 * (double)n_size * (double)d_size * log(2.0 * M_PI);
    }
//    std::cout << "logl_const " << logl_const << std::endl;
    // Start EM.
    for (size_t t = 0; t < max_iter; ++t) {
//        std::cout << " max_iter " << max_iter << " max_prec " << max_prec << " t " << t << std::endl;
        logdet_Ve = EigenProc( V_g, V_e, d_size, D_l, UltVeh, UltVehi);
//        std::cout << "line 428 logdet_Ve: " << logdet_Ve << std::endl;
//        for( int song_t=0; song_t<D_l.get_length(); ++song_t ){
//            std::cout << "line 446 D_l: " << D_l.get_array()[song_t] << std::endl;
//        }
//        std::cout << "line 441" << std::endl;
        logdet_Q = CalcQi(eigen_r, D_l, n_size, d_size, c_size, X, Qi);
//        std::cout << "line 430 logdet_Q: " << logdet_Q << std::endl;
        trmul(UltVehi, Y, UltVehiY);
        CalcXHiY( eigen_r, D_l, X, UltVehiY, n_size, d_size, c_size, xHiy);
//        std::cout << "line 446" << std::endl;
        // Calculate log likelihood/restricted likelihood value, and
        // terminate if change is small.
        double temp_value = MphCalcLogL( eigen_r, xHiy, n_size, d_size, c_size, dc_size, D_l, UltVehiY, Qi);
//        std::cout << "line 446 temp_value: " << temp_value << std::endl;
        logl_new = logl_const + temp_value -
                   0.5 * (double)n_size * logdet_Ve;
//        std::cout << "line 444 logl_new: " << logl_new << std::endl;
        if (func_name == 'R' || func_name == 'r') {
            logl_new += -0.5 * (logdet_Q - (double)c_size * logdet_Ve);
        }
//        std::cout << "line 457" << std::endl;
        if (t != 0 && fabs(logl_new - logl_old) < max_prec) {
            break;
        }
        logl_old = logl_new;
//        std::cout << "line 462" << std::endl;
        CalcOmega( eigen_r,  D_l,  n_size,  d_size, OmegaU, OmegaE);
//        std::cout << "line 464" << std::endl;
//        std::cout << " logl_new " << logl_new << std::endl;
        /*
        for( int ti=0; ti<d_size; ++ti ){
            for( int tj=0; tj<n_size; ++tj ){
                std::cout << "ti " << ti << " tj " << tj << " OmegaU " << OmegaU.get_matrix()[ti][tj] << std::endl;
                std::cout << "ti " << ti << " tj " << tj << " OmegaE " << OmegaE.get_matrix()[ti][tj] << std::endl;
            }
        }*/
        // Update UltVehiB, UltVehiU.
        if (func_name == 'R' || func_name == 'r') {
            UpdateRL_B( xHiy, Qi, n_size, d_size, c_size, dc_size, UltVehiB);
            trmul(UltVehiB, X, UltVehiBX);
        } else if (t == 0) {
            trmul(UltVehi, B, UltVehiB);
            trmul(UltVehiB, X, UltVehiBX);
        }
//        std::cout << "line 480" << std::endl;
///*
//        for( int ti=0; ti<d_size; ++ti ){
//            for( int tj=0; tj<n_size; ++tj ){
//                std::cout << "ti " << ti << " tj " << tj << " UltVehiBX " << UltVehiBX.get_matrix()[ti][tj] << std::endl;
//                std::cout << "ti " << ti << " tj " << tj << " OmegaE " << OmegaE.get_matrix()[ti][tj] << std::endl;
//                std::cout << "ti " << ti << " tj " << tj << " UltVehiY " << UltVehiY.get_matrix()[ti][tj] << std::endl;
//            }
//        }
//*/
        UpdateU(OmegaE, UltVehiY, n_size, d_size, UltVehiBX, UltVehiU);
//        std::cout << "line 492" << std::endl;
        //        for( int ti=0; ti<d_size; ++ti ){
//            for( int tj=0; tj<n_size; ++tj ){
////                std::cout << "ti " << ti << " tj " << tj << " UltVehiBX " << UltVehiBX.get_matrix()[ti][tj] << std::endl;
////                std::cout << "ti " << ti << " tj " << tj << " OmegaE " << OmegaE.get_matrix()[ti][tj] << std::endl;
////                std::cout << "ti " << ti << " tj " << tj << " UltVehiY " << UltVehiY.get_matrix()[ti][tj] << std::endl;
////                std::cout << "ti " << ti << " tj " << tj << " UltVehiU " << UltVehiU.get_matrix()[ti][tj] << std::endl;
//            }
//        }

//
//        for( int ti=0; ti<d_size; ++ti ){
//            for( int tj=0; tj<n_size; ++tj ){
//              std::cout << "ti " << ti << " tj " << tj << " UltVehiU " << UltVehiU.get_matrix()[ti][tj] << std::endl;
//            }
//        }

        if (func_name == 'L' || func_name == 'l') {
            // UltVehiBX is destroyed here.
            UpdateL_B(X, XXti, n_size,  d_size,  c_size, UltVehiY, UltVehiU, UltVehiBX, UltVehiB);
//            std::cout << "line 512" << std::endl;
            trmul(UltVehiB, X,UltVehiBX);
        }
/*
        for( int ti=0; ti<d_size; ++ti ){
            for( int tj=0; tj<n_size; ++tj ){
                std::cout << "ti " << ti << " tj " << tj << " UltVehiBX " << UltVehiBX.get_matrix()[ti][tj] << std::endl;
            }
        }
*/
//        std::cout << "line 521" << std::endl;
        UpdateE(UltVehiY, UltVehiBX, n_size, d_size, UltVehiU, UltVehiE);

//        for( int ti=0; ti<d_size; ++ti ){
//            for( int tj=0; tj<n_size; ++tj ){
//                std::cout << "ti " << ti << " tj " << tj << " UltVehiE " << UltVehiE.get_matrix()[ti][tj] << std::endl;
//            }
//        }

        // Calculate U_hat, E_hat and B.
        My_matrix<double> UltVeh_t(UltVeh.get_num_column(), UltVeh.get_num_row());
        T_matrix(UltVeh, UltVeh_t);
        trmul(UltVeh_t, UltVehiU, U_hat);
        trmul(UltVeh_t, UltVehiE, E_hat);
        trmul(UltVeh_t, UltVehiB, B);

        // Calculate Sigma_uu and Sigma_ee.
        CalcSigma(func_name, eigen_r, n_size, d_size, c_size, dc_size, D_l, X,  OmegaU, OmegaE,
                  UltVeh, Qi, Sigma_uu, Sigma_ee);
/*
        for( int i_t=0; i_t<d_size; ++i_t ){
            for( int j_t=0; j_t<d_size; ++j_t ){
                std::cout << " i_t " << i_t << " j_t " << j_t << " Sigma_uu " << Sigma_uu.get_matrix()[i_t][j_t] << std::endl;
                std::cout << " i_t " << i_t << " j_t " << j_t << " Sigma_ee " << Sigma_ee.get_matrix()[i_t][j_t] << std::endl;
            }
        }
*/
        // Update V_g and V_e.
//        std::cout << "line 548" << std::endl;
        UpdateV( eigen_r, U_hat, E_hat,  n_size,  d_size, c_size, Sigma_uu, Sigma_ee, V_g, V_e);
/*
        for( int i_t=0; i_t<d_size; ++i_t ){
            for( int j_t=0; j_t<d_size; ++j_t ){
                std::cout << " i_t " << i_t << " j_t " << j_t << " V_g " << V_g.get_matrix()[i_t][j_t] << std::endl;
                std::cout << " i_t " << i_t << " j_t " << j_t << " V_e " << V_e.get_matrix()[i_t][j_t] << std::endl;
            }
        }
        */
//
//        for( int i_t=0; i_t<d_size; ++i_t ){
//            for( int j_t=0; j_t<n_size; ++j_t ){
//                std::cout << " i_t " << i_t << " j_t " << j_t << " E_hat " << E_hat.get_matrix()[i_t][j_t] << std::endl;
//            }
//        }

    }
    return logl_new;
}

// Calculate p-value, beta (d by 1 vector) and V(beta).
// todo it sames this is to get the wald-test pvalue
// I should change it to log likelyhood ratio test

double MphCalcP(const Eigen_result & eigen_r, const My_Vector<double> & x_vec,
                const size_t & n_size, const size_t & d_size, const size_t & c_size,
                const My_matrix<double> & W, const My_matrix<double> & Y, const My_matrix<double> & V_g,
                const My_matrix<double> &V_e,  My_matrix<double> & UltVehiY, My_Vector<double> & beta,
                My_matrix<double> &Vbeta) {

    double delta, dl, d, d1, d2, dy, dx, dw; //  logdet_Ve, logdet_Q, p_value;
    size_t dc_size = d_size*c_size;
    My_Vector<double> D_l(d_size);
    My_matrix<double> UltVeh(d_size, d_size);
    My_matrix<double> UltVehi(d_size, d_size);
    My_matrix<double> Qi(dc_size, dc_size);
    My_matrix<double> WHix(dc_size, d_size);
    My_matrix<double> QiWHix(dc_size, d_size);
    My_matrix<double> xPx(d_size, d_size);
    My_Vector<double> xPy(d_size);
    My_Vector<double> WHiy(dc_size);
    xPx.set_values_zero();
    WHix.set_values_zero();
    xPy.set_values_zero();
    WHiy.set_values_zero();

    // Eigen decomposition and calculate log|Ve|.
    EigenProc( V_g, V_e, d_size, D_l, UltVeh, UltVehi);
    // Calculate Qi and log|Q|.
    CalcQi( eigen_r, D_l, n_size, d_size, c_size, W,  Qi);
    // Calculate UltVehiY.
    trmul(UltVehi, Y, UltVehiY);
    // Calculate WHix, WHiy, xHiy, xHix.
    size_t i, j, k;
    for ( i = 0; i < d_size; ++i) {
        dl = D_l.get_array()[i];

        d1 = 0.0;
        d2 = 0.0;
        for ( k = 0; k < n_size; ++k) {
            delta = eigen_r.get_eigen_values().get_array()[k];
            dx = x_vec.get_array()[k];
            dy = UltVehiY.get_matrix()[i][k];
            d1 += dx * dy / (delta * dl + 1.0);
            d2 += dx * dx / (delta * dl + 1.0);
        }
        xPy.get_array()[i]=d1;
        xPx.get_matrix()[i][i]=d2;

        for ( j = 0; j < c_size; ++j) {
            d1 = 0.0;
            d2 = 0.0;
            for ( k = 0; k < n_size; ++k) {
                delta = eigen_r.get_eigen_values().get_array()[k];
                dx = x_vec.get_array()[k];
                dw = W.get_matrix()[j][k];
                dy = UltVehiY.get_matrix()[i][j];

                d1 += dx * dw / (delta * dl + 1.0);
                d2 += dy * dw / (delta * dl + 1.0);
            }
            WHix.get_matrix()[j * d_size + i][i]=d1;
            WHiy.get_array()[j * d_size + i]=d2;
        }
    }

    trmul(Qi, WHix, QiWHix);
    My_matrix<double> WHix_t(WHix.get_num_column(), WHix.get_num_row());
    T_matrix(WHix, WHix_t);
    trmul(WHix_t, QiWHix, xPx);
    for( i=0; i<xPx.get_num_row(); ++i ){
        for( j=0; j<xPx.get_num_column(); ++j ){
            xPx.get_matrix()[i][j] = -1 * xPx.get_matrix()[i][j];
        }
    }
    for( i=0; i<d_size; ++i ){
        for( j=0; j<d_size; ++j ){
            xPx.get_matrix()[i][j] = -1*xPx.get_matrix()[i][j];
        }
    }
    My_matrix<double> QiWHix_t(d_size, dc_size);
    T_matrix(QiWHix, QiWHix_t);
    My_Vector<double> QiWHix_tWHiy(d_size);
    trmul(QiWHix_t, WHiy, QiWHix_tWHiy);
    for( i=0; i <d_size; ++i ){
        xPy.get_array()[i] -= QiWHix_tWHiy.get_array()[i];
    }

    Vbeta.value_copy(xPx);
    inverse_matrix(Vbeta);
    // Need to multiply UltVehi on both sides or one side.
    My_matrix<double> UltVeh_t(UltVeh.get_num_column(), UltVeh.get_num_row());
    T_matrix(UltVeh, UltVeh_t);
    trmul(UltVeh_t, D_l, beta);
    trmul( Vbeta, UltVeh, xPx);
    trmul(UltVeh_t, xPx, Vbeta);

    d=0.0;
    for( i=0; i<d_size; ++i ){
        d += D_l.get_array()[i]*xPy.get_array()[i];
    }
    return 1-chii(d, d_size);
}


// Calculate B and its standard error (which is a matrix of the same
// dimension as B).
void MphCalcBeta(const Eigen_result & eigen_r, const My_matrix<double> & W,
                 const size_t & n_size, const size_t & d_size, const size_t & c_size,
                 const  My_matrix<double> &Y, const  My_matrix<double> &V_g,
                 const  My_matrix<double> &V_e,  My_matrix<double> &UltVehiY,  My_matrix<double> &B,
                 My_matrix<double> &se_B) {
    double delta, dl, d, dy, dw; // , logdet_Ve, logdet_Q;

    size_t dc_size = d_size * c_size;

    My_Vector<double> D_l(d_size);
    My_matrix<double> UltVeh(d_size, d_size);
    My_matrix<double> UltVehi(d_size, d_size);
    My_matrix<double> Qi(dc_size, dc_size);
    My_matrix<double> Qi_temp(dc_size, dc_size);
    My_Vector<double> WHiy(dc_size);
    My_Vector<double> QiWHiy(dc_size);
    My_Vector<double> beta(dc_size);
    My_matrix<double> Vbeta(dc_size, dc_size);

    WHiy.set_values_zero();

    // Eigen decomposition and calculate log|Ve|.
    EigenProc( V_g, V_e, d_size, D_l, UltVeh, UltVehi);

    // Calculate Qi and log|Q|.
    CalcQi(eigen_r, D_l, n_size, d_size, c_size, W, Qi);

    // Calculate UltVehiY.
    trmul(UltVehi, Y, UltVehiY);

    // Calculate WHiy.
    size_t i, j, k, l, m;
    for (i = 0; i < d_size; ++i) {
        dl = D_l.get_array()[i];
        for (j = 0; j < c_size; ++j) {
            d = 0.0;
            for ( k = 0; k < n_size; ++k) {
                delta = eigen_r.get_eigen_values().get_array()[k];
                dw = W.get_matrix()[j][k];
                dy = UltVehiY.get_matrix()[i][k];
                d += dy * dw / (delta * dl + 1.0);
            }
            WHiy.get_array()[j * d_size + i]=d;
        }
    }
    trmul(Qi, WHiy, QiWHiy);
    // Need to multiply I_c\otimes UltVehi on both sides or one side.
    for (i = 0; i < c_size; ++i) {
        for ( j = 0; j < d_size; ++j) {
            beta.get_array()[i * d_size+j]=0;
            for ( k = 0; k < d_size; ++k) {
                beta.get_array()[i * d_size+j] += UltVeh.get_matrix()[k][j]*QiWHiy.get_array()[i*d_size+k];
            }
        }

        for ( j = 0; j < c_size; ++j) {
            if (j < i) {
                My_matrix<double> Vbeta_sym(dc_size, dc_size);
                for( k=0;k<d_size;++k ) {
                    for (l = 0; l < d_size; ++l) {
                        Vbeta_sym.get_matrix()[k][l] = Vbeta.get_matrix()[j * d_size+k][i * d_size+l];
                    }
                }

                for( k=0;k<d_size;++k ){
                    for( l=0;l<d_size;++l ){
                        Vbeta.get_matrix()[i*d_size+k][j * d_size+l] = Vbeta_sym.get_matrix()[l][k];
                    }
                }

            } else {
                for( k=0;k<d_size;++k ){
                    for( l=0;l<d_size;++l ){
                        Qi_temp.get_matrix()[i*d_size+k][j*d_size+l] = 0;
                        for( m=0;m<d_size;++m ){
                            Qi_temp.get_matrix()[i*d_size+k][j*d_size+l] +=
                                    Qi.get_matrix()[i*d_size+k][j*d_size+m]*UltVeh.get_matrix()[m][l];
                        }
                    }
                }
                for( k=0;k<d_size;++k ) {
                    for (l = 0; l < d_size; ++l) {
                        Vbeta.get_matrix()[i*d_size+k][j*d_size+l] = 0;
                        for( m=0;m<d_size;++m ){
                            Vbeta.get_matrix()[i*d_size+k][j*d_size+l] +=
                                    UltVeh.get_matrix()[m][k]*Qi_temp.get_matrix()[i * d_size+m][j * d_size+l];
                        }
                    }
                }
            }
        }
    }
    // Copy beta to B, and Vbeta to se_B.
    for (j = 0; j < c_size; ++j) {
        for (i = 0; i < d_size; ++i) {
            B.get_matrix()[i][j] = beta.get_array()[j * d_size + i];
            se_B.get_matrix()[i][j] = sqrt(Vbeta.get_matrix()[j * d_size + i][j * d_size + i]);
        }
    }
}

// Below are functions for Newton-Raphson's algorithm.

// Calculate all Hi and return logdet_H=\sum_{k=1}^{n}log|H_k|
// and calculate Qi and return logdet_Q
// and calculate yPy.
void CalcHiQi(const Eigen_result & eigen_r, const My_matrix<double> & X,
              const size_t & n_size, const size_t & d_size, const size_t & c_size, const size_t & dc_size,
              const My_matrix<double> & V_g, const My_matrix<double> & V_e, My_matrix<double> & Hi_all,
              My_matrix<double> & Qi, double &logdet_H, double &logdet_Q) {
    Hi_all.set_values_zero();
    Qi.set_values_zero();
    logdet_H = 0.0;
    logdet_Q = 0.0;

    double logdet_Ve = 0.0, delta, dl, d;

    My_matrix<double> mat_dd(d_size, d_size);
    My_matrix<double> UltVeh(d_size, d_size);
    My_matrix<double> UltVeh_t(d_size, d_size);
    My_matrix<double> UltVehi(d_size, d_size);
    My_Vector<double> D_l(d_size);

    // Calculate D_l, UltVeh and UltVehi.
    logdet_Ve = EigenProc(V_g, V_e, d_size, D_l, UltVeh, UltVehi);

    size_t k, i, j, l;
    // Calculate each Hi and log|H_k|.
    logdet_H = (double)n_size * logdet_Ve;
    for (k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        mat_dd.value_copy(UltVehi);
        for (i = 0; i < d_size; ++i) {
            dl = D_l.get_array()[i];
            d = delta * dl + 1.0;

            for( j=0; j <d_size; ++j ){
                mat_dd.get_matrix()[i][j] = mat_dd.get_matrix()[i][j]/d;
            }
            logdet_H += log(d);
        }
        for( i=0; i<d_size; ++i ){
            for( j=0; j<d_size; ++j ){
                Hi_all.get_matrix()[i][k*d_size+j] = 0;
                for( l=0; l<d_size; ++l ){
                    Hi_all.get_matrix()[i][k*d_size+j] += UltVehi.get_matrix()[l][i]*mat_dd.get_matrix()[l][j];
                }
            }
        }
    }

    // Calculate Qi, and multiply I\o times UtVeh on both side and
    // calculate logdet_Q, don't forget to substract
    // c_size*logdet_Ve.
    logdet_Q = CalcQi(eigen_r, D_l, n_size, d_size, c_size, X, Qi) - (double)c_size * logdet_Ve;
    //std::cout << " line 814 logdet_Q " << logdet_Q << std::endl;
    for ( i = 0; i < c_size; ++i) {
        for ( j = 0; j < c_size; ++j) {
            if (j < i) {
                My_matrix<double> Qi_sub(d_size, d_size);
                My_matrix<double> Qi_sym(d_size, d_size);
                for( k=0; k<d_size; ++k ) {
                    for (l = 0; l < d_size; ++l) {
                        Qi_sym.get_matrix()[l][k] = Qi.get_matrix()[j * d_size+l][i * d_size+k];
                    }
                }
                T_matrix(Qi_sym, Qi_sub);
                for( k=0; k<d_size; ++k ){
                    for( l=0; l<d_size; ++l ){
                        Qi.get_matrix()[i * d_size+l][j * d_size+k] = Qi_sub.get_matrix()[l][k];
                    }
                }
            } else {
                My_matrix<double> Qi_sub(d_size, d_size);
                for( k=0; k<d_size; ++k ) {
                    for (l = 0; l < d_size; ++l) {
                        Qi_sub.get_matrix()[l][k] = Qi.get_matrix()[i * d_size+l][j * d_size+k];
                    }
                }
                trmul(Qi_sub, UltVeh, mat_dd);
                T_matrix(UltVeh, UltVeh_t);
                trmul(UltVeh_t, mat_dd, Qi_sub);
                for( k=0; k<d_size; ++k ){
                    for( l=0; l<d_size; ++l ){
                        Qi.get_matrix()[i * d_size+l][j * d_size+k] = Qi_sub.get_matrix()[l][k];
                    }
                }
            }
        }
    }
}

// Calculate all Hiy.
void Calc_Hiy_all(const My_matrix<double> & Y, const My_matrix<double> & Hi_all,
                  const size_t & n_size, const size_t & d_size,
                  My_matrix<double> & Hiy_all) {
    Hiy_all.set_values_zero();
    size_t i, j, k;
    for ( k= 0; k < n_size; ++k) {
        for( i=0; i<d_size; ++i ){
            Hiy_all.get_matrix()[i][k]=0;
            for( j=0; j<d_size; ++j ){
                Hiy_all.get_matrix()[i][k] += Hi_all.get_matrix()[i][k*d_size+j] * Y.get_matrix()[j][k];
            }
        }
    }
}

// Calculate all xHi.
void Calc_xHi_all(const My_matrix<double> & X, const My_matrix<double> &Hi_all,
                  const size_t & n_size, const size_t & d_size, const size_t & c_size,
                  My_matrix<double> &xHi_all) {
    xHi_all.set_values_zero();
    double d;
    size_t i, j, k, l;
    for (k = 0; k < n_size; ++k) {
        for ( i = 0; i < c_size; ++i) {
            d = X.get_matrix()[i][k];
            for( j=0; j<d_size; ++j ){
                for( l=0; l<d_size; ++l ){
                    xHi_all.get_matrix()[i * d_size+j][k * d_size+l] = Hi_all.get_matrix()[j][k * d_size+l];
                    xHi_all.get_matrix()[i*d_size+j][k*d_size+l] = xHi_all.get_matrix()[i*d_size+j][k*d_size+l]*d;
                }
            }
        }
    }
}

// Calculate scalar yHiy.
double Calc_yHiy(const My_matrix<double> & Y, const My_matrix<double> & Hiy_all, const size_t & n_size,
                 const size_t & d_size) {
    double yHiy = 0.0;
    size_t i, j;//k
    for ( i=0; i<n_size; ++i) {
        for( j=0; j<d_size; ++j ){
            yHiy += Hiy_all.get_matrix()[j][i]*Y.get_matrix()[j][i];
        }
    }
    return yHiy;
}

// Calculate the vector xHiy.
void Calc_xHiy(const My_matrix<double> & Y, const My_matrix<double> & xHi,
               const size_t & n_size, const size_t & d_size, const size_t & dc_size, My_Vector<double> & xHiy) {
    xHiy.set_values_zero();
    size_t i, j, k;
    My_Vector<double> xHi_k_y_K(dc_size);
    for ( k = 0; k < n_size; ++k) {
        xHi_k_y_K.set_values_zero();
        for( i=0; i<dc_size; ++i ){
            for( j=0; j<d_size; ++j ){
                xHi_k_y_K.get_array()[i] += xHi.get_matrix()[i][k*d_size+j]*Y.get_matrix()[j][k];
            }
            xHiy.get_array()[i] += xHi_k_y_K.get_array()[i];
        }
    }
}

// 0<=i,j<d_size
size_t GetIndex(const size_t & i, const size_t & j, const size_t & d_size) {
    size_t s, l;
    if (j < i) {
        s = j;
        l = i;
    } else {
        s = i;
        l = j;
    }
    return (2 * d_size - s + 1) * s / 2 + l - s;
}

void Calc_yHiDHiy(const Eigen_result & eigen_r, const My_matrix<double> & Hiy, const size_t & i,
                  const size_t & j, const size_t & n_size, double & yHiDHiy_g, double & yHiDHiy_e) {
    yHiDHiy_g = 0.0;
    yHiDHiy_e = 0.0;
    double delta, d1, d2;
    for (size_t k = 0; k < n_size; ++k) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        d1 = Hiy.get_matrix()[i][k];
        d2 = Hiy.get_matrix()[j][k];
        if (i == j) {
            yHiDHiy_g += delta * d1 * d2;
            yHiDHiy_e += d1 * d2;
        } else {
            yHiDHiy_g += delta * d1 * d2 * 2.0;
            yHiDHiy_e += d1 * d2 * 2.0;
        }
    }
}

void Calc_xHiDHiy(const Eigen_result & eigen_r, const My_matrix<double> & xHi,
                  const size_t & n_size, const size_t & d_size,
                  const My_matrix<double> & Hiy, const size_t & i, const size_t & j,
                  My_Vector<double> & xHiDHiy_g, My_Vector<double> & xHiDHiy_e) {
    xHiDHiy_g.set_values_zero();
    xHiDHiy_e.set_values_zero();

    double delta, d;
    size_t k, l;
    if (i != j) {
        for (k = 0; k < n_size; ++k) {
            delta = eigen_r.get_eigen_values().get_array()[k];
            d = Hiy.get_matrix()[j][k];
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_g.get_array()[l] += d * delta * xHi.get_matrix()[l][k * d_size + i];
            }
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_e.get_array()[l] += d * xHi.get_matrix()[l][k * d_size + i];
            }

            d = Hiy.get_matrix()[i][k];
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_g.get_array()[l] += d * delta * xHi.get_matrix()[l][k * d_size + j];
            }
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_e.get_array()[l] += d * xHi.get_matrix()[l][k * d_size + j];
            }
        }
    }else{
        for (k = 0; k < n_size; ++k) {
            delta = eigen_r.get_eigen_values().get_array()[k];
            d = Hiy.get_matrix()[j][k];
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_g.get_array()[l] += d * delta * xHi.get_matrix()[l][k * d_size + i];
            }
            for (l = 0; l < xHi.get_num_row(); ++l) {
                xHiDHiy_e.get_array()[l] += d * xHi.get_matrix()[l][k * d_size + i];
            }
        }
    }
}

void Calc_xHiDHix(const Eigen_result & eigen_r, const My_matrix<double> & xHi, const size_t & ii,
                  const size_t & jj,
                  const size_t & n_size, const size_t & d_size, const size_t & dc_size, My_matrix<double> & xHiDHix_g,
                  My_matrix<double> & xHiDHix_e) {
    xHiDHix_g.set_values_zero();
    xHiDHix_e.set_values_zero();

    double delta;

    My_matrix<double> mat_dcdc(dc_size, dc_size);
    My_matrix<double> mat_dcdc_t(dc_size, dc_size);
    size_t i, j;
    for (size_t k = 0; k < n_size; ++k) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        mat_dcdc.set_values_zero();
        for( i=0; i<xHi.get_num_row(); ++i ){
            for( j=0; j<xHi.get_num_row(); ++j ){
                mat_dcdc.get_matrix()[i][j] = xHi.get_matrix()[i][k * d_size + ii]*xHi.get_matrix()[j][k * d_size + jj];
            }
        }
        T_matrix(mat_dcdc, mat_dcdc_t);

        for( i=0; i<dc_size; ++i ){
            for( j=0; j<dc_size; ++j ){
                xHiDHix_e.get_matrix()[i][j] += mat_dcdc.get_matrix()[i][j];
            }
        }
        for( i=0; i<dc_size; ++i ){
            for( j=0; j<dc_size; ++j ){
                mat_dcdc.get_matrix()[i][j] = mat_dcdc.get_matrix()[i][j] * delta;
            }
        }
        for( i=0; i<dc_size; ++i ){
            for( j=0; j<dc_size; ++j ){
                xHiDHix_g.get_matrix()[i][j] += mat_dcdc.get_matrix()[i][j];
            }
        }

        if (ii != jj) {
            for( i=0; i<dc_size; ++i ) {
                for (j = 0; j < dc_size; ++j) {
                    xHiDHix_e.get_matrix()[i][j] += mat_dcdc_t.get_matrix()[i][j];
                }
            }
            for( i=0; i<dc_size; ++i ){
                for( j=0; j<dc_size; ++j ){
                    mat_dcdc_t.get_matrix()[i][j] = mat_dcdc_t.get_matrix()[i][j] * delta;
                }
            }
            for( i=0; i<dc_size; ++i ){
                for( j=0; j<dc_size; ++j ){
                    xHiDHix_g.get_matrix()[i][j] += mat_dcdc_t.get_matrix()[i][j];
                }
            }
        }
    }
}

void Calc_yHiDHiDHiy(const Eigen_result & eigen_r, const My_matrix<double> & Hi,
                     const My_matrix<double> & Hiy, const size_t & i1, const size_t & j1,
                     const size_t & i2, const size_t & j2,
                     const size_t & n_size, const size_t & d_size,
                     double &yHiDHiDHiy_gg, double &yHiDHiDHiy_ee, double &yHiDHiDHiy_ge) {
    yHiDHiDHiy_gg = 0.0;
    yHiDHiDHiy_ee = 0.0;
    yHiDHiDHiy_ge = 0.0;

    double delta, d_Hiy_i1, d_Hiy_j1, d_Hiy_i2, d_Hiy_j2;
    double d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;

    for (size_t k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];

        d_Hiy_i1 = Hiy.get_matrix()[i1][k];
        d_Hiy_j1 = Hiy.get_matrix()[j1][k];
        d_Hiy_i2 = Hiy.get_matrix()[i2][k];
        d_Hiy_j2 = Hiy.get_matrix()[j2][k];

        d_Hi_i1i2 = Hi.get_matrix()[i1][k * d_size + i2];
        d_Hi_i1j2 = Hi.get_matrix()[i1][k * d_size + j2];
        d_Hi_j1i2 = Hi.get_matrix()[j1][k * d_size + i2];
        d_Hi_j1j2 = Hi.get_matrix()[j1][k * d_size + j2];

        if (i1 == j1) {
            yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);
            yHiDHiDHiy_ee += (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);
            yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2);

            if (i2 != j2) {
                yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
                yHiDHiDHiy_ee += (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
                yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2);
            }
        } else {
            yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 +
                                              d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);
            yHiDHiDHiy_ee +=
                    (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 + d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);
            yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1i2 * d_Hiy_j2 +
                                      d_Hiy_j1 * d_Hi_i1i2 * d_Hiy_j2);

            if (i2 != j2) {
                yHiDHiDHiy_gg += delta * delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 +
                                                  d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
                yHiDHiDHiy_ee +=
                        (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 + d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
                yHiDHiDHiy_ge += delta * (d_Hiy_i1 * d_Hi_j1j2 * d_Hiy_i2 +
                                          d_Hiy_j1 * d_Hi_i1j2 * d_Hiy_i2);
            }
        }
    }
}

void Calc_xHiDHiDHiy(const Eigen_result & eigen_r, const My_matrix<double> & Hi,
                     const My_matrix<double> &xHi, const My_matrix<double> &Hiy,
                     const size_t & i1, const size_t & j1, const size_t & i2,
                     const size_t & j2, const size_t & n_size, const size_t & d_size,
                     My_Vector<double> & xHiDHiDHiy_gg, My_Vector<double> & xHiDHiDHiy_ee,
                     My_Vector<double> & xHiDHiDHiy_ge) {
    xHiDHiDHiy_gg.set_values_zero();
    xHiDHiDHiy_ee.set_values_zero();
    xHiDHiDHiy_ge.set_values_zero();

    double delta, d_Hiy_i, d_Hiy_j, d_Hi_i1i2, d_Hi_i1j2;
    double d_Hi_j1i2, d_Hi_j1j2;
    size_t k, i;
    for (k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];

        d_Hiy_i = Hiy.get_matrix()[i2][k];
        d_Hiy_j = Hiy.get_matrix()[j2][k];

        d_Hi_i1i2 = Hi.get_matrix()[i1][k * d_size + i2];
        d_Hi_i1j2 = Hi.get_matrix()[i1][k * d_size + j2];
        d_Hi_j1i2 = Hi.get_matrix()[j1][k * d_size + i2];
        d_Hi_j1j2 = Hi.get_matrix()[j1][k * d_size + j2];

        if (i1 == j1) {
            for( i=0; i< xHiDHiDHiy_gg.get_length(); ++i){
                xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_j1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + i1];
                xHiDHiDHiy_ee.get_array()[i] += d_Hi_j1i2 * d_Hiy_j * xHi.get_matrix()[i][k * d_size + i1];
                xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_j1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + i1];
            }

            if (i2 != j2) {
                for( i=0; i< xHiDHiDHiy_gg.get_length(); ++i){
                    xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_j1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + i1];
                    xHiDHiDHiy_ee.get_array()[i] += d_Hi_j1j2 * d_Hiy_i * xHi.get_matrix()[i][k * d_size + i1];
                    xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_j1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + i1];
                }
            }
        } else {

            for( i=0; i< xHiDHiDHiy_gg.get_length(); ++i){
                xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_j1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + i1];
                xHiDHiDHiy_ee.get_array()[i] += d_Hi_j1i2 * d_Hiy_j * xHi.get_matrix()[i][k * d_size + i1];
                xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_j1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + i1];

                xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_i1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + j1];
                xHiDHiDHiy_ee.get_array()[i] += d_Hi_i1i2 * d_Hiy_j * xHi.get_matrix()[i][k * d_size + j1];
                xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_i1i2 * d_Hiy_j*xHi.get_matrix()[i][k * d_size + j1];
            }

            if (i2 != j2) {
                for( i=0; i< xHiDHiDHiy_gg.get_length(); ++i){
                    xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_j1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + i1];
                    xHiDHiDHiy_ee.get_array()[i] += d_Hi_j1j2 * d_Hiy_i * xHi.get_matrix()[i][k * d_size + i1];
                    xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_j1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + i1];

                    xHiDHiDHiy_gg.get_array()[i] += delta * delta * d_Hi_i1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + j1];
                    xHiDHiDHiy_ee.get_array()[i] += d_Hi_i1j2 * d_Hiy_i * xHi.get_matrix()[i][k * d_size + j1];
                    xHiDHiDHiy_ge.get_array()[i] += delta * d_Hi_i1j2 * d_Hiy_i*xHi.get_matrix()[i][k * d_size + j1];
                }
            }
        }
    }
}

void Calc_xHiDHiDHix(const Eigen_result & eigen_r, const My_matrix<double> & Hi,
                     const My_matrix<double> & xHi, const size_t & i1, const size_t & j1,
                     const size_t & i2, const size_t & j2,
                     const size_t & n_size, const size_t & d_size, const size_t &  dc_size,
                     My_matrix<double> & xHiDHiDHix_gg,  My_matrix<double> & xHiDHiDHix_ee,
                     My_matrix<double> & xHiDHiDHix_ge) {
    xHiDHiDHix_gg.set_values_zero();
    xHiDHiDHix_ee.set_values_zero();
    xHiDHiDHix_ge.set_values_zero();
    double delta, d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;
    My_matrix<double> mat_dcdc(dc_size, dc_size);
    size_t l, m;
    for (size_t k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];
        d_Hi_i1i2 = Hi.get_matrix()[i1][k * d_size + i2];//gsl_matrix_get(Hi, i1, k * d_size + i2);
        d_Hi_i1j2 = Hi.get_matrix()[i1][k * d_size + j2];//gsl_matrix_get(Hi, i1, k * d_size + j2);
        d_Hi_j1i2 = Hi.get_matrix()[j1][k * d_size + i2];//gsl_matrix_get(Hi, j1, k * d_size + i2);
        d_Hi_j1j2 = Hi.get_matrix()[j1][k * d_size + j2];//gsl_matrix_get(Hi, j1, k * d_size + j2);

        if (i1 == j1) {
            mat_dcdc.set_values_zero();
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] += d_Hi_j1i2*xHi.get_matrix()[l][k*d_size + i1]*
                                                                              xHi.get_matrix()[m][k*d_size + j2];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }

            if (i2 != j2) {
                mat_dcdc.set_values_zero();

                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] += d_Hi_j1j2*xHi.get_matrix()[l][k * d_size+i1]
                                                                                  *xHi.get_matrix()[m][k*d_size + i2];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
            }
        } else {
            mat_dcdc.set_values_zero();

            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] += d_Hi_j1i2*xHi.get_matrix()[l][k * d_size+i1]
                                                                              *xHi.get_matrix()[m][k*d_size + j2];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            mat_dcdc.set_values_zero();

            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] += d_Hi_i1i2*xHi.get_matrix()[l][k * d_size + j1]*
                                                                              xHi.get_matrix()[m][k * d_size + j2];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                }
            }
            for( l=0; l<dc_size; ++l ){
                for( m=0; m<dc_size; ++m ){
                    xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                }
            }

            if (i2 != j2) {

                mat_dcdc.set_values_zero();
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] += d_Hi_j1j2*xHi.get_matrix()[l][k * d_size + i1]*
                                                                                  xHi.get_matrix()[m][k * d_size + i2];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                mat_dcdc.set_values_zero();

                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] += d_Hi_i1j2*xHi.get_matrix()[l][k * d_size + j1]*
                                                                                  xHi.get_matrix()[m][k * d_size + i2];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ee.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_ge.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        mat_dcdc.get_matrix()[l][m] = mat_dcdc.get_matrix()[l][m]*delta;
                    }
                }
                for( l=0; l<dc_size; ++l ){
                    for( m=0; m<dc_size; ++m ){
                        xHiDHiDHix_gg.get_matrix()[l][m] += mat_dcdc.get_matrix()[l][m];
                    }
                }
            }
        }
    }
}

void Calc_traceHiD(const Eigen_result & eigen_r, const My_matrix<double> & Hi, const size_t & i, const size_t & j,
                   const size_t & n_size, const size_t & d_size,
                   double & tHiD_g, double & tHiD_e) {
    tHiD_g = 0.0;
    tHiD_e = 0.0;
    double delta, d;
    for (size_t k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];//gsl_vector_get(eval, k);
        d = Hi.get_matrix()[j][k * d_size + i];
        if (i == j) {
            tHiD_g += delta * d;
            tHiD_e += d;
        } else {
            tHiD_g += delta * d * 2.0;
            tHiD_e += d * 2.0;
        }
    }
}

void Calc_traceHiDHiD(const Eigen_result & eigen_r, const My_matrix<double> & Hi,
                      const size_t & i1, const size_t & j1, const size_t & i2,
                      const size_t & j2, const size_t & n_size, const size_t & d_size,
                      double & tHiDHiD_gg, double & tHiDHiD_ee,
                      double & tHiDHiD_ge) {
    tHiDHiD_gg = 0.0;
    tHiDHiD_ee = 0.0;
    tHiDHiD_ge = 0.0;

    double delta, d_Hi_i1i2, d_Hi_i1j2, d_Hi_j1i2, d_Hi_j1j2;

    for (size_t k = 0; k < n_size; k++) {
        delta = eigen_r.get_eigen_values().get_array()[k];

        d_Hi_i1i2 = Hi.get_matrix()[i1][k * d_size + i2];
        d_Hi_i1j2 = Hi.get_matrix()[i1][k * d_size + j2];
        d_Hi_j1i2 = Hi.get_matrix()[j1][k * d_size + i2];
        d_Hi_j1j2 = Hi.get_matrix()[j1][k * d_size + j2];

        if (i1 == j1) {
            tHiDHiD_gg += delta * delta * d_Hi_i1j2 * d_Hi_j1i2;
            tHiDHiD_ee += d_Hi_i1j2 * d_Hi_j1i2;
            tHiDHiD_ge += delta * d_Hi_i1j2 * d_Hi_j1i2;

            if (i2 != j2) {
                tHiDHiD_gg += delta * delta * d_Hi_i1i2 * d_Hi_j1j2;
                tHiDHiD_ee += d_Hi_i1i2 * d_Hi_j1j2;
                tHiDHiD_ge += delta * d_Hi_i1i2 * d_Hi_j1j2;
            }
        } else {
            tHiDHiD_gg +=
                    delta * delta * (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);
            tHiDHiD_ee += (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);
            tHiDHiD_ge += delta * (d_Hi_i1j2 * d_Hi_j1i2 + d_Hi_j1j2 * d_Hi_i1i2);

            if (i2 != j2) {
                tHiDHiD_gg +=
                        delta * delta * (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
                tHiDHiD_ee += (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
                tHiDHiD_ge += delta * (d_Hi_i1i2 * d_Hi_j1j2 + d_Hi_j1i2 * d_Hi_i1j2);
            }
        }
    }
}

// trace(PD) = trace((Hi-HixQixHi)D)=trace(HiD) - trace(HixQixHiD) // do it in a different way
void Calc_tracePD(const Eigen_result & eigen_r, const My_matrix<double> & Qi,
                  const My_matrix<double> & Hi, const My_matrix<double> & xHiDHix_all_g,
                  const My_matrix<double> & xHiDHix_all_e, const size_t & i,
                  const size_t & n_size, const size_t & d_size, const size_t & dc_size,
                  const size_t & j, double &tPD_g, double &tPD_e) {
    size_t v = GetIndex(i, j, d_size);

    double d = 0;
    Calc_traceHiD(eigen_r, Hi, i,  j, n_size, d_size, tPD_g, tPD_e);
    size_t k, l;
    for ( k = 0; k < dc_size; ++k) {
        d=0.0;
        for( l=0; l<Qi.get_num_column(); ++l ){
            d += Qi.get_matrix()[k][l]*xHiDHix_all_g.get_matrix()[l][v * dc_size + k];
        }
        tPD_g -= d;
        d=0.0;
        for( l=0; l<Qi.get_num_column(); ++l ){
            d += Qi.get_matrix()[k][l]*xHiDHix_all_e.get_matrix()[l][v * dc_size + k];
        }
        tPD_e -= d;
    }
}

// trace(PDPD) = trace((Hi-HixQixHi)D(Hi-HixQixHi)D)
//             = trace(HiDHiD) - trace(HixQixHiDHiD)
//               - trace(HiDHixQixHiD) + trace(HixQixHiDHixQixHiD)
void Calc_tracePDPD(const Eigen_result & eigen_r, const My_matrix<double> & Qi,
                    const My_matrix<double> & Hi, const My_matrix<double> & xHi,
                    const My_matrix<double> & QixHiDHix_all_g,
                    const My_matrix<double> & QixHiDHix_all_e,
                    const My_matrix<double> & xHiDHiDHix_all_gg,
                    const My_matrix<double> & xHiDHiDHix_all_ee,
                    const My_matrix<double> & xHiDHiDHix_all_ge, const size_t & i1,
                    const size_t & j1, const size_t & i2, const size_t & j2,
                    const size_t & n_size, const size_t & d_size,
                    double &tPDPD_gg, double & tPDPD_ee, double & tPDPD_ge) {
    const size_t dc_size = Qi.get_num_row();
    size_t v_size = d_size * (d_size + 1) / 2;
    size_t v1 = GetIndex(i1, j1, d_size), v2 = GetIndex(i2, j2, d_size);
    double d;
    // Calculate the first part: trace(HiDHiD).
    Calc_traceHiDHiD(eigen_r, Hi, i1, j1, i2, j2, n_size, d_size, tPDPD_gg, tPDPD_ee, tPDPD_ge);
//    std::cout << "line 1435 tPDPD_gg " << tPDPD_gg <<std::endl;
    size_t i, k;
    for (i = 0; i < dc_size; ++i) {
        d=0.0;
        for( k=0; k< Qi.get_num_column(); ++k ){
            d += Qi.get_matrix()[i][k]*xHiDHiDHix_all_gg.get_matrix()[k][(v1 * v_size + v2) * dc_size + i];
        }
        tPDPD_gg -= d * 2.0;
        d=0.0;
        for( k=0; k< Qi.get_num_column(); ++k ){
            d += Qi.get_matrix()[i][k]*xHiDHiDHix_all_ee.get_matrix()[k][(v1 * v_size + v2) * dc_size + i];
        }
        tPDPD_ee -= d * 2.0;
        d=0.0;
        for( k=0; k< Qi.get_num_column(); ++k ){
            d += Qi.get_matrix()[i][k]*xHiDHiDHix_all_ge.get_matrix()[k][(v1 * v_size + v2) * dc_size + i];
        }
        tPDPD_ge -= d * 2.0;
    }
//    std::cout << "line 1454 tPDPD_gg " << tPDPD_gg <<std::endl;
    // Calculate the fourth part: trace(HixQixHiDHixQixHiD).
    for (i = 0; i < dc_size; ++i) {
        d=0.0;
        for( k=0; k< dc_size; ++k ){
//            std::cout << "line 1459 d " << d <<std::endl;
            d += QixHiDHix_all_g.get_matrix()[i][v1 * dc_size+k]*QixHiDHix_all_g.get_matrix()[k][v2 * dc_size + i];
//            std::cout << "line 1461 d " << d <<std::endl;
        }
//        std::cout << "line 1463 tPDPD_gg " << tPDPD_gg <<std::endl;
//        std::cout << "line 1464 d " << d <<std::endl;
        tPDPD_gg += d;
        d=0.0;
        for( k=0; k<dc_size; ++k ){
            d += QixHiDHix_all_e.get_matrix()[i][v1 * dc_size+k]*QixHiDHix_all_e.get_matrix()[k][v2 * dc_size + i];
        }
        tPDPD_ee += d;
        d=0.0;
        for( k=0; k< dc_size; ++k ){
            d += QixHiDHix_all_g.get_matrix()[i][v1 * dc_size+k]*QixHiDHix_all_e.get_matrix()[k][v2 * dc_size + i];
        }
        tPDPD_ge += d;
    }
}

// Calculate (xHiDHiy) for every pair (i,j).
void Calc_xHiDHiy_all(const Eigen_result & eigen_r, const My_matrix<double> & xHi,
                      const My_matrix<double> & Hiy, const size_t & n_size, const size_t & d_size,
                      My_matrix<double> & xHiDHiy_all_g,
                      My_matrix<double> & xHiDHiy_all_e) {
    xHiDHiy_all_g.set_values_zero();
    xHiDHiy_all_e.set_values_zero();

    size_t v;
    size_t i, j, k;
    My_Vector<double> xHiDHiy_g(xHiDHiy_all_g.get_num_row());
    My_Vector<double> xHiDHiy_e(xHiDHiy_all_g.get_num_row());

    for (i = 0; i < d_size; ++i) {
        for (j = 0; j < d_size; ++j) {
            if (j < i) {
                continue;
            }
            v = GetIndex(i, j, d_size);
            for( k=0; k <xHiDHiy_all_g.get_num_row(); ++k ){
                xHiDHiy_g.get_array()[k] = xHiDHiy_all_g.get_matrix()[k][v];
                xHiDHiy_e.get_array()[k] = xHiDHiy_all_e.get_matrix()[k][v];
            }
            Calc_xHiDHiy(eigen_r, xHi, n_size, d_size, Hiy, i, j, xHiDHiy_g, xHiDHiy_e);
            for( k=0; k <xHiDHiy_all_g.get_num_row(); ++k ){
                xHiDHiy_all_g.get_matrix()[k][v]=xHiDHiy_g.get_array()[k];
                xHiDHiy_all_e.get_matrix()[k][v]=xHiDHiy_e.get_array()[k];
            }
        }
    }
}

// Calculate (xHiDHix) for every pair (i,j).
void Calc_xHiDHix_all(const Eigen_result & eigen_r, const My_matrix<double> & xHi,
                      const size_t & n_size, const size_t & d_size, const size_t & dc_size,
                      My_matrix<double> & xHiDHix_all_g, My_matrix<double> & xHiDHix_all_e) {
    xHiDHix_all_g.set_values_zero();
    xHiDHix_all_e.set_values_zero();

    size_t v;
    size_t k, l;
    My_matrix<double> xHiDHix_g(dc_size, dc_size);
    My_matrix<double> xHiDHix_e(dc_size, dc_size);
    size_t j;
    for (size_t i = 0; i < d_size; ++i) {
        for ( j = 0; j < d_size; ++j) {
            if (j < i) {
                continue;
            }
            v = GetIndex(i, j, d_size);

            for( k=0; k<dc_size; ++k ){
                for( l=0; l<dc_size; ++l ){
                    xHiDHix_g.get_matrix()[k][l] = xHiDHix_all_g.get_matrix()[k][v * dc_size+l];
                    xHiDHix_e.get_matrix()[k][l] = xHiDHix_all_e.get_matrix()[k][v * dc_size+l];
                }
            }
            Calc_xHiDHix( eigen_r, xHi, i, j, n_size, d_size, dc_size, xHiDHix_g, xHiDHix_e);

            for( k=0; k<dc_size; ++k ){
                for( l=0; l<dc_size; ++l ){
                    xHiDHix_all_g.get_matrix()[k][v * dc_size+l] = xHiDHix_g.get_matrix()[k][l];
                    xHiDHix_all_e.get_matrix()[k][v * dc_size+l] = xHiDHix_e.get_matrix()[k][l];
                }
            }
        }
    }
}

// Calculate (xHiDHiy) for every pair (i,j).
void Calc_xHiDHiDHiy_all(const size_t & v_size, const Eigen_result & eigen_r,
                         const My_matrix<double> & Hi, const My_matrix<double> & xHi,
                         const My_matrix<double> & Hiy, const size_t & n_size, const size_t & d_size,
                         My_matrix<double> & xHiDHiDHiy_all_gg,
                         My_matrix<double> & xHiDHiDHiy_all_ee,
                         My_matrix<double> & xHiDHiDHiy_all_ge) {
    xHiDHiDHiy_all_gg.set_values_zero();
    xHiDHiDHiy_all_ee.set_values_zero();
    xHiDHiDHiy_all_ge.set_values_zero();

//    size_t d_size = Hiy->size1;
    size_t v1, v2;
    size_t j1;
    size_t i2, j2;
    My_Vector<double> xHiDHiDHiy_gg(xHiDHiDHiy_all_gg.get_num_row());
    My_Vector<double> xHiDHiDHiy_ee(xHiDHiDHiy_all_gg.get_num_row());
    My_Vector<double> xHiDHiDHiy_ge(xHiDHiDHiy_all_gg.get_num_row());
    size_t i;
    for (size_t i1 = 0; i1 < d_size; ++i1) {
        for ( j1 = 0; j1 < d_size; ++j1) {
            if (j1 < i1) {
                continue;
            }
            v1 = GetIndex(i1, j1, d_size);

            for ( i2 = 0; i2 < d_size; ++i2) {
                for ( j2 = 0; j2 < d_size; ++j2) {
                    if (j2 < i2) {
                        continue;
                    }
                    v2 = GetIndex(i2, j2, d_size);
                    for( i=0; i<xHiDHiDHiy_all_gg.get_num_row(); ++i ){
                        xHiDHiDHiy_gg.get_array()[i] = xHiDHiDHiy_all_gg.get_matrix()[i][v1 * v_size + v2];
                        xHiDHiDHiy_ee.get_array()[i] = xHiDHiDHiy_all_ee.get_matrix()[i][v1 * v_size + v2];
                        xHiDHiDHiy_ge.get_array()[i] = xHiDHiDHiy_all_ge.get_matrix()[i][v1 * v_size + v2];
                    }
                    Calc_xHiDHiDHiy( eigen_r, Hi, xHi, Hiy, i1, j1, i2, j2, n_size, d_size, xHiDHiDHiy_gg,
                                     xHiDHiDHiy_ee, xHiDHiDHiy_ge);
                    for( i=0; i<xHiDHiDHiy_all_gg.get_num_row(); ++i ){
                        xHiDHiDHiy_all_gg.get_matrix()[i][v1 * v_size + v2] = xHiDHiDHiy_gg.get_array()[i];
                        xHiDHiDHiy_all_ee.get_matrix()[i][v1 * v_size + v2] = xHiDHiDHiy_ee.get_array()[i];
                        xHiDHiDHiy_all_ge.get_matrix()[i][v1 * v_size + v2] = xHiDHiDHiy_ge.get_array()[i];
                    }
                }
            }
        }
    }
}

// Calculate (xHiDHix) for every pair (i,j).
void Calc_xHiDHiDHix_all(const size_t & v_size, const Eigen_result & eigen_r,
                         const My_matrix<double> & Hi, const My_matrix<double> & xHi,
                         const size_t & n_size, const size_t & d_size, const size_t & dc_size,
                         My_matrix<double> & xHiDHiDHix_all_gg,
                         My_matrix<double> & xHiDHiDHix_all_ee,
                         My_matrix<double> & xHiDHiDHix_all_ge) {
    xHiDHiDHix_all_gg.set_values_zero();
    xHiDHiDHix_all_ee.set_values_zero();
    xHiDHiDHix_all_ge.set_values_zero();

//    size_t d_size = xHi->size2 / eval->size, dc_size = xHi->size1;
    My_matrix<double> xHiDHiDHix_gg1(dc_size, dc_size);
    My_matrix<double> xHiDHiDHix_ee1(dc_size, dc_size);
    My_matrix<double> xHiDHiDHix_ge1(dc_size, dc_size);

    size_t v1, v2;
    size_t j1, j2, i2, i, j;
    for (size_t i1 = 0; i1 < d_size; ++i1) {
        for ( j1 = 0; j1 < d_size; ++j1) {
            if (j1 < i1) {
                continue;
            }
            v1 = GetIndex(i1, j1, d_size);

            for ( i2 = 0; i2 < d_size; ++i2) {
                for ( j2 = 0; j2 < d_size; ++j2) {
                    if (j2 < i2) {
                        continue;
                    }
                    v2 = GetIndex(i2, j2, d_size);

                    if (v2 < v1) {
                        continue;
                    }

                    for( i=0; i<dc_size; ++i ){
                        for( j=0; j<dc_size; ++j ){
                            xHiDHiDHix_gg1.get_matrix()[i][j]=xHiDHiDHix_all_gg.get_matrix()[i][(v1*v_size+v2)*dc_size+j];
                            xHiDHiDHix_ee1.get_matrix()[i][j]=xHiDHiDHix_all_ee.get_matrix()[i][(v1*v_size+v2)*dc_size+j];
                            xHiDHiDHix_ge1.get_matrix()[i][j]=xHiDHiDHix_all_ge.get_matrix()[i][(v1*v_size+v2)*dc_size+j];
                        }
                    }

                    Calc_xHiDHiDHix(eigen_r, Hi, xHi, i1, j1, i2, j2, n_size, d_size, dc_size, xHiDHiDHix_gg1,
                                    xHiDHiDHix_ee1, xHiDHiDHix_ge1);

                    for( i=0; i<dc_size; ++i ){
                        for( j=0; j<dc_size; ++j ){
                            xHiDHiDHix_all_gg.get_matrix()[i][(v1*v_size+v2)*dc_size+j] = xHiDHiDHix_gg1.get_matrix()[i][j];
                            xHiDHiDHix_all_ee.get_matrix()[i][(v1*v_size+v2)*dc_size+j] = xHiDHiDHix_ee1.get_matrix()[i][j];
                            xHiDHiDHix_all_ge.get_matrix()[i][(v1*v_size+v2)*dc_size+j] = xHiDHiDHix_ge1.get_matrix()[i][j];
                        }
                    }

                    if (v2 != v1) {
                        for( i=0; i<dc_size; ++i ) {
                            for (j = 0; j < dc_size; ++j) {
                                xHiDHiDHix_all_gg.get_matrix()[i][(v2 * v_size + v1) * dc_size+j]
                                        = xHiDHiDHix_gg1.get_matrix()[i][j];
                                xHiDHiDHix_all_ee.get_matrix()[i][(v2 * v_size + v1) * dc_size+j]
                                        = xHiDHiDHix_ee1.get_matrix()[i][j];
                                xHiDHiDHix_all_ge.get_matrix()[i][(v2 * v_size + v1) * dc_size+j]
                                        = xHiDHiDHix_ge1.get_matrix()[i][j];
                            }
                        }
                    }
                }
            }
        }
    }
}

// Calculate (xHiDHix)Qi(xHiy) for every pair (i,j).
void Calc_xHiDHixQixHiy_all(const My_matrix<double> & xHiDHix_all_g,
                            const My_matrix<double> & xHiDHix_all_e,
                            const My_Vector<double> & QixHiy, const size_t & dc_size,
                            My_matrix<double> & xHiDHixQixHiy_all_g,
                            My_matrix<double> & xHiDHixQixHiy_all_e) {
    size_t v_size = xHiDHix_all_g.get_num_column() / dc_size;
    size_t j, k;
    /*
    My_matrix<double> xHiDHix_g(dc_size, dc_size);
    My_matrix<double> xHiDHix_e(dc_size, dc_size);

    My_Vector<double> xHiDHixQixHiy_g(xHiDHixQixHiy_all_g.get_num_row());
    My_Vector<double> xHiDHixQixHiy_e(xHiDHixQixHiy_all_e.get_num_row());
*/
    for (size_t i = 0; i < v_size; i++) {
        for( j=0; j<xHiDHixQixHiy_all_e.get_num_row(); ++j ) {
            xHiDHixQixHiy_all_g.get_matrix()[j][i] = 0;
            xHiDHixQixHiy_all_e.get_matrix()[j][i] = 0;

            for (k = 0; k < QixHiy.get_length(); ++k) {
                xHiDHixQixHiy_all_g.get_matrix()[j][i] += xHiDHix_all_g.get_matrix()[j][i * dc_size + k]
                                                          * QixHiy.get_array()[k];
                xHiDHixQixHiy_all_e.get_matrix()[j][i] += xHiDHix_all_e.get_matrix()[j][i * dc_size + k]
                                                          * QixHiy.get_array()[k];
            }
        }
    }
}

// Calculate Qi(xHiDHiy) and Qi(xHiDHix)Qi(xHiy) for each pair of i,j (i<=j).
void Calc_QiVec_all(const My_matrix<double> & Qi, const My_matrix<double> & vec_all_g,
                    const My_matrix<double> & vec_all_e, My_matrix<double> & Qivec_all_g,
                    My_matrix<double> & Qivec_all_e) {
//
//    std::cout << "Qi row " << Qi.get_num_row() <<std::endl;
//    std::cout << "Qi column " << Qi.get_num_column() <<std::endl;
//
//    std::cout << "vec_all_g row " << vec_all_g.get_num_row() <<std::endl;
//    std::cout << "vec_all_g column " << vec_all_g.get_num_column() <<std::endl;
//
//    std::cout << "vec_all_e row " << vec_all_e.get_num_row() <<std::endl;
//    std::cout << "vec_all_e column " << vec_all_e.get_num_column() <<std::endl;
//
//    std::cout << "Qivec_all_g row " << Qivec_all_g.get_num_row() <<std::endl;
//    std::cout << "Qivec_all_g column " << Qivec_all_g.get_num_column() <<std::endl;
//
//    std::cout << "Qivec_all_e row " << Qivec_all_e.get_num_row() <<std::endl;
//    std::cout << "Qivec_all_e column " << Qivec_all_e.get_num_column() <<std::endl;

    size_t j, k;
    for (size_t i = 0; i < vec_all_g.get_num_column(); ++i) {
        for( j=0; j<Qivec_all_g.get_num_row(); ++j ){
            Qivec_all_g.get_matrix()[j][i]=0.0;
            Qivec_all_e.get_matrix()[j][i]=0.0;
            for( k=0; k<vec_all_g.get_num_row(); ++k ){
                Qivec_all_g.get_matrix()[j][i] += Qi.get_matrix()[j][k]*vec_all_g.get_matrix()[k][i];
                Qivec_all_e.get_matrix()[j][i] += Qi.get_matrix()[j][k]*vec_all_e.get_matrix()[k][i];
            }
        }
    }
}

// Calculate Qi(xHiDHix) for each pair of i,j (i<=j).
void Calc_QiMat_all(const My_matrix<double> & Qi, const My_matrix<double> & mat_all_g,
                    const My_matrix<double> & mat_all_e, const size_t & dc_size,
                    My_matrix<double> & Qimat_all_g,
                    My_matrix<double> & Qimat_all_e) {
//    std::cout << "line 1801 dc_size: " << dc_size << std::endl;
    size_t v_size = mat_all_g.get_num_column() / mat_all_g.get_num_row();
    size_t j, k, l;
    for (size_t i = 0; i < v_size; ++i) {
        for( j=0; j<dc_size; ++j ){
            for( k=0; k<dc_size; ++k ){
                Qimat_all_g.get_matrix()[j][i* dc_size+k]=0.0;
                Qimat_all_e.get_matrix()[j][i* dc_size+k]=0.0;
                for( l=0; l<dc_size; ++l ){
                    Qimat_all_g.get_matrix()[j][i* dc_size+k] += Qi.get_matrix()[j][l] * mat_all_g.get_matrix()[l][i * dc_size+k];
                    Qimat_all_e.get_matrix()[j][i* dc_size+k] += Qi.get_matrix()[j][l] * mat_all_e.get_matrix()[l][i * dc_size+k];
                }
            }
        }
    }
}

// Calculate yPDPy
// yPDPy = y(Hi-HixQixHi)D(Hi-HixQixHi)y
//       = ytHiDHiy - (yHix)Qi(xHiDHiy) - (yHiDHix)Qi(xHiy)
//         + (yHix)Qi(xHiDHix)Qi(xtHiy)
void Calc_yPDPy(const Eigen_result &eigen_r, const My_matrix<double> & Hiy,
                const My_Vector<double> & QixHiy, const My_matrix<double> & xHiDHiy_all_g,
                const My_matrix<double> & xHiDHiy_all_e,
                const My_matrix<double> & xHiDHixQixHiy_all_g,
                const My_matrix<double> & xHiDHixQixHiy_all_e, const size_t & i,
                const size_t & j, const size_t & n_size, const size_t & d_size,
                double & yPDPy_g, double & yPDPy_e) {
    size_t v = GetIndex(i, j, d_size);

    double d=0.0;

    // First part: ytHiDHiy.
    Calc_yHiDHiy(eigen_r, Hiy, i, j, n_size, yPDPy_g, yPDPy_e);
    size_t k;
    for( k=0; k<xHiDHiy_all_g.get_num_row(); ++k ){
        d += QixHiy.get_array()[k]*xHiDHiy_all_g.get_matrix()[k][v];
    }
    yPDPy_g -= d * 2.0;
    d=0.0;
    for( k=0; k<xHiDHiy_all_e.get_num_row(); ++k ){
        d += QixHiy.get_array()[k]*xHiDHiy_all_e.get_matrix()[k][v];
    }
    yPDPy_e -= d * 2.0;

    d=0.0;
    for( k=0; k<xHiDHixQixHiy_all_g.get_num_row(); ++k ){
        d += QixHiy.get_array()[k]*xHiDHixQixHiy_all_g.get_matrix()[k][v];
    }
    yPDPy_g += d;
    d=0.0;
    for( k=0; k<xHiDHixQixHiy_all_e.get_num_row(); ++k ){
        d += QixHiy.get_array()[k]*xHiDHixQixHiy_all_e.get_matrix()[k][v];
    }
    yPDPy_e += d;
}


void Calc_yPDPDPy(
        const Eigen_result & eigen_r, const My_matrix<double> & Hi, const My_matrix<double> & xHi,
        const My_matrix<double> & Hiy, const  My_Vector<double> & QixHiy,
        const My_matrix<double> & xHiDHiy_all_g, const  My_matrix<double> & xHiDHiy_all_e,
        const My_matrix<double> & QixHiDHiy_all_g, const  My_matrix<double> & QixHiDHiy_all_e,
        const My_matrix<double> & xHiDHixQixHiy_all_g,
        const My_matrix<double> & xHiDHixQixHiy_all_e,
        const My_matrix<double> & QixHiDHixQixHiy_all_g,
        const My_matrix<double> & QixHiDHixQixHiy_all_e,
        const My_matrix<double> & xHiDHiDHiy_all_gg, const My_matrix<double> & xHiDHiDHiy_all_ee,
        const My_matrix<double> & xHiDHiDHiy_all_ge, const My_matrix<double> & xHiDHiDHix_all_gg,
        const My_matrix<double> & xHiDHiDHix_all_ee, const My_matrix<double> & xHiDHiDHix_all_ge,
        const size_t & i1, const size_t & j1, const size_t & i2, const size_t & j2,
        const size_t & n_size, /*const size_t & d_size, const size_t & dc_size,*/
        double &yPDPDPy_gg, double &yPDPDPy_ee, double &yPDPDPy_ge) {
    size_t d_size = Hi.get_num_row();
    size_t dc_size = xHi.get_num_row();
    size_t v1 = GetIndex(i1, j1, d_size), v2 = GetIndex(i2, j2, d_size);
    size_t v_size = d_size * (d_size + 1) / 2;

    My_Vector<double> xHiDHiDHixQixHiy(dc_size);

    // First part: yHiDHiDHiy.
    Calc_yHiDHiDHiy(eigen_r, Hi, Hiy, i1, j1, i2, j2, n_size, d_size,yPDPDPy_gg, yPDPDPy_ee,
                    yPDPDPy_ge);
//    std::cout << "line 1911 i1 " << i1 << " i2 " << i2 << " j1 " << j1 << " j2 " << j2 << " yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;
    // Second and third parts:
    // -(yHix)Qi(xHiDHiDHiy) - (yHiDHiDHix)Qi(xHiy).
    size_t i, j;
    double d=0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_gg.get_matrix()[i][v1 * v_size + v2];
    }
    yPDPDPy_gg -= d;
    d = 0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_ee.get_matrix()[i][v1 * v_size + v2];
    }
    yPDPDPy_ee -= d;
    d = 0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_ge.get_matrix()[i][v1 * v_size + v2];
    }
    yPDPDPy_ge -= d;


    d = 0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_gg.get_matrix()[i][v2 * v_size + v1];
    }
    yPDPDPy_gg -= d;
    d = 0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_ee.get_matrix()[i][v2 * v_size + v1];
    }
    yPDPDPy_ee -= d;
    d = 0.0;
    for( i=0; i<QixHiy.get_length(); ++i){
        d += QixHiy.get_array()[i]*xHiDHiDHiy_all_ge.get_matrix()[i][v2 * v_size + v1];
    }
    yPDPDPy_ge -= d;

//    std::cout << "line 1948 i1 " << i1 << " i2 " << i2 << " j1 " << j1 << " j2 " << j2 << " yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;
    // Fourth part: - (yHiDHix)Qi(xHiDHiy).
    d = 0.0;
    for( i=0; i<xHiDHiy_all_g.get_num_row(); ++i ){
        d += xHiDHiy_all_g.get_matrix()[i][v1]*QixHiDHiy_all_g.get_matrix()[i][v2];
    }
    yPDPDPy_gg -= d;
    d = 0.0;
    for( i=0; i<xHiDHiy_all_g.get_num_row(); ++i ){
        d += xHiDHiy_all_e.get_matrix()[i][v1]*QixHiDHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ee -= d;
    d = 0.0;
    for( i=0; i<xHiDHiy_all_g.get_num_row(); ++i ){
        d += xHiDHiy_all_g.get_matrix()[i][v1]*QixHiDHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ge -= d;

//    std::cout << "line 1966 i1 " << i1 << " i2 " << i2 << " j1 " << j1 << " j2 " << j2 << " yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;
    // Fifth and sixth parts:
    //   + (yHix)Qi(xHiDHix)Qi(xHiDHiy) +
    //   (yHiDHix)Qi(xHiDHix)Qi(xHiy)
    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row() ; ++i ){
        d += xHiDHixQixHiy_all_g.get_matrix()[i][v1] * QixHiDHiy_all_g.get_matrix()[i][v2];
    }
    yPDPDPy_gg += d;
    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row() ; ++i ){
        d += xHiDHixQixHiy_all_g.get_matrix()[i][v2] * QixHiDHiy_all_g.get_matrix()[i][v1];
    }
    yPDPDPy_gg += d;

    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row() ; ++i ){
        d += xHiDHixQixHiy_all_e.get_matrix()[i][v1] * QixHiDHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ee += d;
    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row() ; ++i ){
        d += xHiDHixQixHiy_all_e.get_matrix()[i][v2] * QixHiDHiy_all_e.get_matrix()[i][v1];
    }
    yPDPDPy_ee += d;

    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row() ; ++i ){
        d += xHiDHixQixHiy_all_g.get_matrix()[i][v1] * QixHiDHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ge += d;
    d = 0.0;
    for( i=0; i<xHiDHixQixHiy_all_g.get_num_row(); ++i ){
        d += xHiDHixQixHiy_all_e.get_matrix()[i][v2] * QixHiDHiy_all_g.get_matrix()[i][v1];
    }
    yPDPDPy_ge += d;
 //   std::cout << "line 2002 i1 " << i1 << " i2 " << i2 << " j1 " << j1 << " j2 " << j2 << " yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;
    // Seventh part: + (yHix)Qi(xHiDHiDHix)Qi(xHiy)
    xHiDHiDHixQixHiy.set_values_zero();
    for( i=0; i<dc_size; ++i ){
        for( j=0; j<dc_size; ++j ){
            xHiDHiDHixQixHiy.get_array()[i] += xHiDHiDHix_all_gg.get_matrix()[i][(v1 * v_size + v2) * dc_size+j]
                                                                              * QixHiy.get_array()[j];
        }
    }
    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += xHiDHiDHixQixHiy.get_array()[i]*QixHiy.get_array()[i];
    }
    yPDPDPy_gg += d;

    xHiDHiDHixQixHiy.set_values_zero();
    for( i=0; i<dc_size; ++i ){
        for( j=0; j<dc_size; ++j ){
            xHiDHiDHixQixHiy.get_array()[i] += xHiDHiDHix_all_ee.get_matrix()[i][(v1 * v_size + v2) * dc_size+j]
                                                                              * QixHiy.get_array()[j];
        }
    }
    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += xHiDHiDHixQixHiy.get_array()[i]*QixHiy.get_array()[i];
    }
    yPDPDPy_ee += d;

    xHiDHiDHixQixHiy.set_values_zero();
    for( i=0; i<dc_size; ++i ){
        for( j=0; j<dc_size; ++j ){
            xHiDHiDHixQixHiy.get_array()[i] += xHiDHiDHix_all_ge.get_matrix()[i][(v1 * v_size + v2) * dc_size+j]
                                                                              * QixHiy.get_array()[j];
        }
    }
    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += xHiDHiDHixQixHiy.get_array()[i]*QixHiy.get_array()[i];
    }
    yPDPDPy_ge += d;
//    std::cout << "line 2042 i1 " << i1 << " i2 " << i2 << " j1 " << j1 << " j2 " << j2 << " yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;
// Eighth part: - (yHix)Qi(xHiDHix)Qi(xHiDHix)Qi(xHiy).
    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += QixHiDHixQixHiy_all_g.get_matrix()[i][v1]*xHiDHixQixHiy_all_g.get_matrix()[i][v2];
    }
    yPDPDPy_gg -= d;
    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += QixHiDHixQixHiy_all_e.get_matrix()[i][v1]*xHiDHixQixHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ee -= d;

    d=0.0;
    for( i=0; i<xHiDHiDHixQixHiy.get_length(); ++i ){
        d += QixHiDHixQixHiy_all_g.get_matrix()[i][v1]*xHiDHixQixHiy_all_e.get_matrix()[i][v2];
    }
    yPDPDPy_ge -= d;
}

// Calculate Edgeworth correctation factors for small samples notation
// and method follows Thomas J. Rothenberg, Econometirca 1984; 52 (4)
// M=xHiDHix
void CalcCRT(const My_matrix<double> & Hessian_inv, const My_matrix<double> & Qi,
             const My_matrix<double> & QixHiDHix_all_g,
             const My_matrix<double> & QixHiDHix_all_e,
             const My_matrix<double> & xHiDHiDHix_all_gg,
             const My_matrix<double> & xHiDHiDHix_all_ee,
             const My_matrix<double> & xHiDHiDHix_all_ge,
             const size_t & d_size, const size_t & c_size, const size_t & dc_size,
             double &crt_a, double &crt_b, double &crt_c) {
    crt_a = 0.0;
    crt_b = 0.0;
    crt_c = 0.0;

    size_t v_size = Hessian_inv.get_num_row() / 2;
    double h_gg, h_ge, h_ee, d, B = 0.0, C = 0.0, D = 0.0;
    double trCg1, trCe1, trCg2, trCe2, trB_gg, trB_ge, trB_ee;
    double trCC_gg, trCC_ge, trCC_ee, trD_gg = 0.0, trD_ge = 0.0, trD_ee = 0.0;

    My_matrix<double> QiMQi_g1 (dc_size, dc_size);
    My_matrix<double> QiMQi_e1 (dc_size, dc_size);
    My_matrix<double> QiMQi_g2 (dc_size, dc_size);
    My_matrix<double> QiMQi_e2 (dc_size, dc_size);

    My_matrix<double> QiMQisQisi_g1 (d_size, d_size);
    My_matrix<double> QiMQisQisi_e1 (d_size, d_size);
    My_matrix<double> QiMQisQisi_g2 (d_size, d_size);
    My_matrix<double> QiMQisQisi_e2 (d_size, d_size);

    My_matrix<double> QiMQiMQi_gg (dc_size, dc_size);
    My_matrix<double> QiMQiMQi_ge (dc_size, dc_size);
    My_matrix<double> QiMQiMQi_ee (dc_size, dc_size);

    My_matrix<double> QiMMQi_gg (dc_size, dc_size);
    My_matrix<double> QiMMQi_ge (dc_size, dc_size);
    My_matrix<double> QiMMQi_ee (dc_size, dc_size);

    My_matrix<double> M_dd(d_size, d_size);
    My_matrix<double> M_dcdc(dc_size, dc_size);

    My_matrix<double> Qi_sub(d_size, d_size);

    size_t i, j;
    My_matrix<double> Qi_si(d_size, d_size);
    for( i=0; i<d_size; ++i ){
        for( j=0; j<d_size; ++j ){
            Qi_si.get_matrix()[i][j] = Qi.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1) * d_size+j];
        }
    }

    inverse_matrix(Qi_si);

    My_matrix<double> QiM_g1(dc_size, dc_size);
    My_matrix<double> QiM_e1(dc_size, dc_size);

    My_matrix<double> QiMQi_g1_s(dc_size, dc_size);
    My_matrix<double> QiMQi_e1_s(dc_size, dc_size);

    My_matrix<double> QiM_g2(dc_size, dc_size);
    My_matrix<double> QiM_e2(dc_size, dc_size);

    My_matrix<double> QiMQi_g2_s(dc_size, dc_size);
    My_matrix<double> QiMQi_e2_s(dc_size, dc_size);

    My_matrix<double> QiMQiMQi_gg_s(d_size, d_size);
    My_matrix<double> QiMQiMQi_ge_s(d_size, d_size);
    My_matrix<double> QiMQiMQi_ee_s(d_size, d_size);

    My_matrix<double> MM_gg(dc_size, dc_size);
    My_matrix<double> MM_ge(dc_size, dc_size);
    My_matrix<double> MM_ee(dc_size, dc_size);

    My_matrix<double> QiMMQi_gg_s(d_size, d_size);
    My_matrix<double> QiMMQi_ge_s(d_size, d_size);
    My_matrix<double> QiMMQi_ee_s(d_size, d_size);
    // Calculate correction factors.
    size_t k;
    for (size_t v1 = 0; v1 < v_size; ++v1) {
        for( i=0; i<dc_size; ++i ){
            for( j=0; j<dc_size; ++j ){
                QiM_g1.get_matrix()[i][j]=QixHiDHix_all_g.get_matrix()[i][v1 * dc_size+j];
                QiM_e1.get_matrix()[i][j]=QixHiDHix_all_e.get_matrix()[i][v1 * dc_size+j];
            }
        }
        trmul(QiM_g1, Qi, QiMQi_g1); //QiMQi_g1 is a new matrix here
        trmul(QiM_e1, Qi, QiMQi_e1);

        for( i=0; i<dc_size; ++i ){
            for( j=0; j<dc_size; ++j ){
                QiMQi_g1_s.get_matrix()[i][j]=QiMQi_g1.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1) * d_size+j];
                QiMQi_e1_s.get_matrix()[i][j]=QiMQi_e1.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1) * d_size+j];
            }
        }
        trmul(QiMQi_g1_s, Qi_si, QiMQisQisi_g1);
        // Calculate trCg1 and trCe1.
        trCg1 = 0.0;
        for ( k = 0; k < d_size; ++k) {
            trCg1 -= QiMQisQisi_g1.get_matrix()[k][k];//gsl_matrix_get(QiMQisQisi_g1, k, k);
        }
        trmul(QiMQi_e1_s, Qi_si, QiMQisQisi_e1);
        trCe1 = 0.0;
        for ( k = 0; k < d_size; ++k) {
            trCe1 -= QiMQisQisi_e1.get_matrix()[k][k];//gsl_matrix_get(QiMQisQisi_e1, k, k);
        }
        for (size_t v2 = 0; v2 < v_size; ++v2) {
            if (v2 < v1) {
                continue;
            }
            for( i=0; i<dc_size; ++i ){
                for( j=0; j<dc_size; ++j ){
                    QiM_g2.get_matrix()[i][j]=QixHiDHix_all_g.get_matrix()[i][v2 * dc_size+j];
                    QiM_e2.get_matrix()[i][j]=QixHiDHix_all_e.get_matrix()[i][v2 * dc_size+j];
                }
            }
            trmul(QiM_g2, Qi, QiMQi_g2);
            trmul(QiM_e2, Qi, QiMQi_e2);

            for( i=0; i<dc_size; ++i ){
                for( j=0; j<dc_size; ++j ){
                    QiMQi_g2_s.get_matrix()[i][j]=QiMQi_g2.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1) * d_size+j];
                    QiMQi_e2_s.get_matrix()[i][j]=QiMQi_e2.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1) * d_size+j];
                }
            }
            // Calculate trCg2 and trCe2.
            trmul(QiMQi_g2_s, Qi_si, QiMQisQisi_g2);

            trCg2 = 0.0;
            for ( k = 0; k < d_size; ++k) {
                trCg2 -= QiMQisQisi_g2.get_matrix()[k][k];//gsl_matrix_get(QiMQisQisi_g2, k, k);
            }
            trmul(QiMQi_e2_s, Qi_si, QiMQisQisi_e2);
            trCe2 = 0.0;
            for ( k = 0; k < d_size; ++k) {
                trCe2 -= QiMQisQisi_e2.get_matrix()[k][k]; //gsl_matrix_get(QiMQisQisi_e2, k, k);
            }
            trmul(QiMQisQisi_g1, QiMQisQisi_g2, M_dd);
            trCC_gg = 0.0;
            for ( k = 0; k < d_size; ++k) {
                trCC_gg += M_dd.get_matrix()[k][k]; //gsl_matrix_get(M_dd, k, k);
            }

            trmul(QiMQisQisi_g1, QiMQisQisi_e2, M_dd);
            for( i=0; i<d_size; ++i ){
                for( j=0; j<d_size; ++j ){
                    for( k=0; k<d_size; ++k ){
                        M_dd.get_matrix()[i][j] += QiMQisQisi_e1.get_matrix()[i][k]*QiMQisQisi_g2.get_matrix()[k][j];
                    }
                }
            }
            trCC_ge = 0.0;
            for ( k = 0; k < d_size; ++k) {
                trCC_ge += M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
            }
            trmul(QiMQisQisi_e1, QiMQisQisi_e2, M_dd);
            trCC_ee = 0.0;
            for ( k = 0; k < d_size; ++k) {
                trCC_ee += M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
            }

            // Calculate Qi(xHiDHix)Qi(xHiDHix)Qi, and subpart of it.
            trmul(QiM_g1, QiMQi_g2, QiMQiMQi_gg);
            trmul(QiM_g1, QiMQi_e2, QiMQiMQi_ge);
            for( i=0; i<d_size; ++i ){
                for( j=0; j<d_size; ++j ){
                    for( k=0; k<d_size; ++k ){
                        QiMQiMQi_ge.get_matrix()[i][j] += QiM_e1.get_matrix()[i][k]*QiMQi_g2.get_matrix()[k][j];
                    }
                }
            }
            trmul(QiM_e1, QiMQi_e2, QiMQiMQi_ee);

            for( i=0; i<d_size; ++i ) {
                for (j = 0; j < d_size; ++j) {
                    QiMQiMQi_gg_s.get_matrix()[i][j] = QiMQiMQi_gg.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                         * d_size+j];
                    QiMQiMQi_ge_s.get_matrix()[i][j] = QiMQiMQi_ge.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                         * d_size+j];
                    QiMQiMQi_ee_s.get_matrix()[i][j] = QiMQiMQi_ee.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                         * d_size+j];
                }
            }
            // and part of trB_gg, trB_ge, trB_ee.
            trmul(QiMQiMQi_gg_s, Qi_si, M_dd);
            trB_gg = 0.0;
            for ( k = 0; k < d_size; ++k) {
                d = M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
                trB_gg -= d;
            }
            trmul(QiMQiMQi_ge_s, Qi_si, M_dd);
            trB_ge = 0.0;
            for ( k = 0; k < d_size; ++k) {
                d = M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
                trB_ge -= d;
            }
            trmul(QiMQiMQi_ee_s, Qi_si, M_dd);
            trB_ee = 0.0;
            for ( k = 0; k < d_size; ++k) {
                d = M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
                trB_ee -= d;
            }

            for( i=0; i<dc_size; ++i ){
                for( j=0; j<dc_size; ++j ){
                    MM_gg.get_matrix()[i][j] = xHiDHiDHix_all_gg.get_matrix()[i][(v1 * v_size + v2) * dc_size+j];
                    MM_ge.get_matrix()[i][j] = xHiDHiDHix_all_ge.get_matrix()[i][(v1 * v_size + v2) * dc_size+j];
                    MM_ee.get_matrix()[i][j] = xHiDHiDHix_all_ee.get_matrix()[i][(v1 * v_size + v2) * dc_size+j];
                }
            }
            trmul(Qi, MM_gg, M_dcdc);

            trmul(M_dcdc, Qi, QiMMQi_gg);
            trmul(Qi, MM_ge, M_dcdc);
            trmul(M_dcdc, Qi, QiMMQi_ge);
            trmul(Qi, MM_ee, M_dcdc);
            trmul(M_dcdc, Qi, QiMMQi_ee);

            for( i=0; i<d_size; ++i ){
                for( j=0; j<d_size; ++j ){
                    QiMMQi_gg_s.get_matrix()[i][j] = QiMMQi_gg.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                     * d_size+j];
                    QiMMQi_ge_s.get_matrix()[i][j] = QiMMQi_ge.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                     * d_size+j];
                    QiMMQi_ee_s.get_matrix()[i][j] = QiMMQi_ee.get_matrix()[(c_size - 1) * d_size+i][(c_size - 1)
                                                                                                     * d_size+j];
                }
            }

            // Calculate the other part of trB_gg, trB_ge, trB_ee.
            trmul(QiMMQi_gg_s, Qi_si, M_dd);
            for ( k = 0; k < d_size; ++k) {
                trB_gg += M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
            }
            trmul(QiMMQi_ge_s, Qi_si, M_dd);
            for ( k = 0; k < d_size; ++k) {
                trB_ge += 2.0 * M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
            }
            trmul(QiMMQi_ee_s, Qi_si, M_dd);
            for ( k = 0; k < d_size; ++k) {
                trB_ee += M_dd.get_matrix()[k][k];//gsl_matrix_get(M_dd, k, k);
            }
            // Calculate trD_gg, trD_ge, trD_ee.
            trD_gg = 2.0 * trB_gg;
            trD_ge = 2.0 * trB_ge;
            trD_ee = 2.0 * trB_ee;

            // calculate B, C and D
            h_gg = -1.0 * Hessian_inv.get_matrix()[v1][v2];//gsl_matrix_get(Hessian_inv, v1, v2);
            h_ge = -1.0 * Hessian_inv.get_matrix()[v1][v2 + v_size];//gsl_matrix_get(Hessian_inv, v1, v2 + v_size);
            h_ee = -1.0 * Hessian_inv.get_matrix()[v1 + v_size][v2 + v_size];//gsl_matrix_get(Hessian_inv, v1 + v_size, v2 + v_size);

            B += h_gg * trB_gg + h_ge * trB_ge + h_ee * trB_ee;
            C += h_gg * (trCC_gg + 0.5 * trCg1 * trCg2) +
                 h_ge * (trCC_ge + 0.5 * trCg1 * trCe2 + 0.5 * trCe1 * trCg2) +
                 h_ee * (trCC_ee + 0.5 * trCe1 * trCe2);
            D += h_gg * (trCC_gg + 0.5 * trD_gg) + h_ge * (trCC_ge + 0.5 * trD_ge) +
                 h_ee * (trCC_ee + 0.5 * trD_ee);

            if (v1 != v2) {
                B += h_gg * trB_gg + h_ge * trB_ge + h_ee * trB_ee;
                C += h_gg * (trCC_gg + 0.5 * trCg1 * trCg2) +
                     h_ge * (trCC_ge + 0.5 * trCg1 * trCe2 + 0.5 * trCe1 * trCg2) +
                     h_ee * (trCC_ee + 0.5 * trCe1 * trCe2);
                D += h_gg * (trCC_gg + 0.5 * trD_gg) + h_ge * (trCC_ge + 0.5 * trD_ge) +
                     h_ee * (trCC_ee + 0.5 * trD_ee);
            }
        }
    }

    // Calculate a, b, c from B C D.
    crt_a = 2.0 * D - C;
    crt_b = 2.0 * B;
    crt_c = C;

}

//CalcDev(func_name, eigen_r, Qi, Hi_all, xHi_all, Hiy_all, QixHiy, n_size,
//        d_size, c_size, dc_size, gradient,
//        Hessian_inv, crt_a, crt_b, crt_c);

// Calculate first-order and second-order derivatives.
void CalcDev(const char func_name, const Eigen_result & eigen_r, const My_matrix<double> & Qi,
             const My_matrix<double> & Hi, const My_matrix<double> & xHi, const My_matrix<double> & Hiy,
             const My_Vector<double> & QixHiy, const size_t & n_size,
        /*const size_t & d_size, const size_t & c_size, const size_t & dc_size,*/
             My_Vector<double> & gradient,
             My_matrix<double> & Hessian_inv, double &crt_a, double &crt_b,
             double &crt_c) {
           // because se(Vg):  from the stand output is not correct

//    size_t n_size = ;

    size_t d_size = Hi.get_num_row();
    size_t dc_size = Qi.get_num_row();
    size_t c_size = dc_size/d_size;


    size_t v_size = d_size * (d_size + 1) / 2;
    size_t v1, v2;
    double dev1_g, dev1_e, dev2_gg, dev2_ee, dev2_ge;

    My_matrix<double> Hessian(v_size * 2, v_size * 2);

    My_matrix<double> xHiDHiy_all_g(dc_size, v_size);
    My_matrix<double> xHiDHiy_all_e(dc_size, v_size);
    My_matrix<double> xHiDHix_all_g(dc_size, v_size * dc_size);
    My_matrix<double> xHiDHix_all_e(dc_size, v_size * dc_size);
    My_matrix<double> xHiDHixQixHiy_all_g(dc_size, v_size);
    My_matrix<double> xHiDHixQixHiy_all_e(dc_size, v_size);

    My_matrix<double> QixHiDHiy_all_g(dc_size, v_size);
    My_matrix<double> QixHiDHiy_all_e(dc_size, v_size);
    My_matrix<double> QixHiDHix_all_g(dc_size, v_size * dc_size);
    My_matrix<double> QixHiDHix_all_e(dc_size, v_size * dc_size);

    My_matrix<double> QixHiDHixQixHiy_all_g(dc_size, v_size);
    My_matrix<double> QixHiDHixQixHiy_all_e(dc_size, v_size);

    My_matrix<double> xHiDHiDHiy_all_gg(dc_size, v_size * v_size);
    My_matrix<double> xHiDHiDHiy_all_ee(dc_size, v_size * v_size);
    My_matrix<double> xHiDHiDHiy_all_ge(dc_size, v_size * v_size);

    My_matrix<double> xHiDHiDHix_all_gg(dc_size, v_size * v_size * dc_size);
    My_matrix<double> xHiDHiDHix_all_ee(dc_size, v_size * v_size * dc_size);
    My_matrix<double> xHiDHiDHix_all_ge(dc_size, v_size * v_size * dc_size);

    Calc_xHiDHiy_all(eigen_r, xHi, Hiy, n_size, d_size, xHiDHiy_all_g, xHiDHiy_all_e);

    int i, j;
//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " xHiDHiy_all_g " << xHiDHiy_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHiy_all_e " << xHiDHiy_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_xHiDHix_all(eigen_r, xHi, n_size, d_size, dc_size, xHiDHix_all_g, xHiDHix_all_e);

//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " xHiDHix_all_g " << xHiDHix_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHix_all_e " << xHiDHix_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_xHiDHixQixHiy_all( xHiDHix_all_g, xHiDHix_all_e, QixHiy, dc_size, xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e);

//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " xHiDHixQixHiy_all_g " << xHiDHixQixHiy_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHixQixHiy_all_e " << xHiDHixQixHiy_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }


    Calc_xHiDHiDHiy_all(v_size, eigen_r, Hi, xHi, Hiy, n_size, d_size, xHiDHiDHiy_all_gg, xHiDHiDHiy_all_ee,
                        xHiDHiDHiy_all_ge);
//
//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " xHiDHiDHiy_all_gg " << xHiDHiDHiy_all_gg.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHiDHiy_all_ee " << xHiDHiDHiy_all_ee.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHiDHiy_all_ge " << xHiDHiDHiy_all_ge.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_xHiDHiDHix_all(v_size, eigen_r, Hi, xHi, n_size, d_size, dc_size, xHiDHiDHix_all_gg, xHiDHiDHix_all_ee,
                        xHiDHiDHix_all_ge);

//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " xHiDHiDHix_all_gg " << xHiDHiDHix_all_gg.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHiDHix_all_ee " << xHiDHiDHix_all_ee.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " xHiDHiDHix_all_ge " << xHiDHiDHix_all_ge.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_QiVec_all(Qi, xHiDHiy_all_g, xHiDHiy_all_e, QixHiDHiy_all_g,
                   QixHiDHiy_all_e);
//
//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " QixHiDHiy_all_g " << QixHiDHiy_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " QixHiDHiy_all_e " << QixHiDHiy_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_QiVec_all(Qi, xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e,
                   QixHiDHixQixHiy_all_g, QixHiDHixQixHiy_all_e);

//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " QixHiDHixQixHiy_all_g " << QixHiDHixQixHiy_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " QixHiDHixQixHiy_all_e " << QixHiDHixQixHiy_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }

    Calc_QiMat_all(Qi, xHiDHix_all_g, xHiDHix_all_e, dc_size, QixHiDHix_all_g,
                   QixHiDHix_all_e);
//
//    for( i=0; i<dc_size; ++i ){
//        for( j=0; j<v_size; ++j ){
//            std::cout << "i " << i << " j " << j << " QixHiDHix_all_g " << QixHiDHix_all_g.get_matrix()[i][j] << std::endl;
//            std::cout << "i " << i << " j " << j << " QixHiDHix_all_e " << QixHiDHix_all_e.get_matrix()[i][j] << std::endl;
//        }
//    }


    double tHiD_g, tHiD_e, tPD_g, tPD_e, tHiDHiD_gg, tHiDHiD_ee;
    double tHiDHiD_ge, tPDPD_gg, tPDPD_ee, tPDPD_ge;
    double yPDPy_g, yPDPy_e, yPDPDPy_gg, yPDPDPy_ee, yPDPDPy_ge;

    // Calculate gradient and Hessian for Vg.
    size_t j1, i2, j2;
    for (size_t i1 = 0; i1 < d_size; ++i1) {
        for ( j1 = 0; j1 < d_size; ++j1) {
            if (j1 < i1) {
                continue;
            }
            v1 = GetIndex(i1, j1, d_size);

            Calc_yPDPy(eigen_r, Hiy, QixHiy, xHiDHiy_all_g, xHiDHiy_all_e, xHiDHixQixHiy_all_g, xHiDHixQixHiy_all_e, i1,
                       j1, n_size, d_size, yPDPy_g, yPDPy_e);

            if (func_name == 'R' || func_name == 'r') {
                Calc_tracePD(eigen_r, Qi, Hi, xHiDHix_all_g, xHiDHix_all_e, i1, n_size, d_size, dc_size,
                                  j1, tPD_g, tPD_e);

                dev1_g = -0.5 * tPD_g + 0.5 * yPDPy_g;
                dev1_e = -0.5 * tPD_e + 0.5 * yPDPy_e;
                //std::cout << "line 2445 tPD_g: " << tPD_g << " yPDPy_g " << yPDPy_g << " tPD_e " << tPD_e << " yPDPy_e " << yPDPy_e << std::endl;
            } else {
                Calc_traceHiD(eigen_r, Hi, i1, j1, n_size, d_size, tHiD_g, tHiD_e);
                dev1_g = -0.5 * tHiD_g + 0.5 * yPDPy_g;
                dev1_e = -0.5 * tHiD_e + 0.5 * yPDPy_e;
            }

            gradient.get_array()[v1] = dev1_g;
            gradient.get_array()[v1 + v_size] = dev1_e;

            for ( i2 = 0; i2 < d_size; ++i2) {
                for ( j2 = 0; j2 < d_size; ++j2) {
                    if (j2 < i2) {
                        continue;
                    }
                    v2 = GetIndex(i2, j2, d_size);

                    if (v2 < v1) {
                        continue;
                    }

                    Calc_yPDPDPy(eigen_r, Hi, xHi, Hiy, QixHiy, xHiDHiy_all_g, xHiDHiy_all_e,
                                 QixHiDHiy_all_g, QixHiDHiy_all_e, xHiDHixQixHiy_all_g,
                                 xHiDHixQixHiy_all_e, QixHiDHixQixHiy_all_g,
                                 QixHiDHixQixHiy_all_e, xHiDHiDHiy_all_gg,
                                 xHiDHiDHiy_all_ee, xHiDHiDHiy_all_ge,
                                 xHiDHiDHix_all_gg, xHiDHiDHix_all_ee,
                                 xHiDHiDHix_all_ge, i1, j1, i2, j2,
                                 n_size, /*d_size, dc_size,*/ yPDPDPy_gg, yPDPDPy_ee, yPDPDPy_ge);
                    //std::cout << " line 2477 yPDPDPy_gg " << yPDPDPy_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;

                    // AI for REML.

                    if (func_name == 'R' || func_name == 'r') {

                        Calc_tracePDPD( eigen_r, Qi, Hi, xHi, QixHiDHix_all_g, QixHiDHix_all_e,
                                        xHiDHiDHix_all_gg, xHiDHiDHix_all_ee,
                                        xHiDHiDHix_all_ge, i1, j1, i2, j2, n_size, d_size,
                                        tPDPD_gg, tPDPD_ee, tPDPD_ge);


                        dev2_gg = 0.5 * tPDPD_gg - yPDPDPy_gg;
                        dev2_ee = 0.5 * tPDPD_ee - yPDPDPy_ee;
                        dev2_ge = 0.5 * tPDPD_ge - yPDPDPy_ge;

//                        std::cout << " line 2492 tPDPD_gg " << tPDPD_gg << " yPDPDPy_ee " << yPDPDPy_ee << " yPDPDPy_ge " << yPDPDPy_ge << std::endl;

                    } else {
                        Calc_traceHiDHiD(eigen_r, Hi, i1, j1, i2, j2, n_size, d_size,
                                         tHiDHiD_gg, tHiDHiD_ee, tHiDHiD_ge);
                        dev2_gg = 0.5 * tHiDHiD_gg - yPDPDPy_gg;
                        dev2_ee = 0.5 * tHiDHiD_ee - yPDPDPy_ee;
                        dev2_ge = 0.5 * tHiDHiD_ge - yPDPDPy_ge;
                    }

                    // Set up Hessian.
                    Hessian.get_matrix()[v1][v2] = dev2_gg;
                    Hessian.get_matrix()[v1 + v_size][v2 + v_size] = dev2_ee;
                    Hessian.get_matrix()[v1][v2 + v_size] = dev2_ge;
                    Hessian.get_matrix()[v2 + v_size][v1] = dev2_ge;

                    if (v1 != v2) {
                        Hessian.get_matrix()[v2][v1] = dev2_gg;
                        Hessian.get_matrix()[v2+v_size][v1+v_size] = dev2_ee;
                        Hessian.get_matrix()[v2][v1 + v_size] = dev2_ge;
                        Hessian.get_matrix()[v1 + v_size][v2] = dev2_ge;
                    }
                }
            }
        }
    }
//    for( i=0; i<Hessian.get_num_row(); ++i ){
//        for( j=0; j<Hessian.get_num_column(); ++j ){
//            std::cout << "line 2515 i " << i << " j " << j << " Hessian " << Hessian.get_matrix()[i][j] << std::endl;
//        }
//    }

    // Invert Hessian.
    Hessian_inv.value_copy(Hessian);
    inverse_matrix(Hessian_inv);
//
//    for( i=0; i<Hessian.get_num_row(); ++i ){
//        for( j=0; j<Hessian.get_num_column(); ++j ){
//            std::cout << "line 2515 i " << i << " j " << j << " Hessian_inv " << Hessian_inv.get_matrix()[i][j] << std::endl;
//        }
//    }

    // Calculate Edgeworth correction factors after inverting
    // Hessian.
    if (c_size > 1) {
        CalcCRT( Hessian_inv, Qi, QixHiDHix_all_g, QixHiDHix_all_e,
                 xHiDHiDHix_all_gg, xHiDHiDHix_all_ee, xHiDHiDHix_all_ge,
                 d_size, c_size, dc_size, crt_a, crt_b, crt_c);
    } else {
        crt_a = 0.0;
        crt_b = 0.0;
        crt_c = 0.0;
    }
}

// Update Vg, Ve.
void UpdateVgVe(const My_matrix<double> & Hessian_inv, const My_Vector<double> & gradient,
                const double & step_scale,
                My_matrix<double> & V_g, My_matrix<double> & V_e) {
    size_t v_size = gradient.get_length()/2, d_size = V_g.get_num_row();
    size_t v;
    My_Vector<double> vec_v(v_size * 2);

    double d;

    // Vectorize Vg and Ve.
    size_t i, j;
    for (i = 0; i < d_size; ++i) {
        for (j = 0; j < d_size; ++j) {
            if (j < i) {
                continue;
            }
            v = GetIndex(i, j, d_size);

            d = V_g.get_matrix()[i][j];//gsl_matrix_get(V_g, i, j);
            vec_v.get_array()[v]=d;//gsl_vector_set(vec_v, v, d);

            d = V_e.get_matrix()[i][j];//d = gsl_matrix_get(V_e, i, j);
            vec_v.get_array()[v+v_size]=d;//gsl_vector_set(vec_v, v + v_size, d);
        }
    }
    for (i = 0; i < Hessian_inv.get_num_row(); ++i) {
        for (j = 0; j < Hessian_inv.get_num_column(); ++j) {
            vec_v.get_array()[i] += -1.0 * step_scale*Hessian_inv.get_matrix()[i][j]*gradient.get_array()[j];
        }
    }

    // Save Vg and Ve.
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j < d_size; ++j) {
            if (j < i) {
                continue;
            }
            v = GetIndex(i, j, d_size);

            d = vec_v.get_array()[v];//gsl_vector_get(vec_v, v);
            V_g.get_matrix()[i][j]=d;
            V_g.get_matrix()[j][i]=d;

            d = vec_v.get_array()[v+v_size];//gsl_vector_get(vec_v, v + v_size);
            V_e.get_matrix()[i][j]=d;
            V_e.get_matrix()[j][i]=d;
        }
    }
}

double MphNR(const char & func_name, const size_t & max_iter, const double & max_prec,
             const Eigen_result & eigen_r, const My_matrix<double> & X, const My_matrix<double> & Y,
             const size_t & n_size,const size_t & d_size,const size_t & c_size,
             My_matrix<double> & Hi_all, My_matrix<double> & xHi_all, My_matrix<double> & Hiy_all,
             My_matrix<double> & V_g, My_matrix<double> & V_e, My_matrix<double> & Hessian_inv,
             double &crt_a, double &crt_b, double &crt_c) {
    size_t dc_size = d_size*c_size;
    size_t v_size = d_size * (d_size + 1) / 2;

    double logdet_H, logdet_Q, yPy, logl_const;
    double logl_old = 0.0, logl_new = 0.0, step_scale;
    size_t step_iter, flag_pd;

    My_matrix<double> Vg_save(d_size, d_size);
    My_matrix<double> Ve_save(d_size, d_size);
    My_matrix<double> V_temp(d_size, d_size);
    My_matrix<double> U_temp(d_size, d_size);
    My_Vector<double> D_temp(d_size);
    My_Vector<double> xHiy(dc_size);
    My_Vector<double> QixHiy(dc_size);
    My_matrix<double> Qi(dc_size, dc_size);
    My_matrix<double> Xt(n_size, c_size);
    My_matrix<double> XXt(c_size, c_size);
    My_Vector<double> gradient(v_size * 2);

    T_matrix(X, Xt);
    trmul(X, Xt, XXt);
    // Calculate |XXt| and (XXt)^{-1}.

    // Calculate the constant for logl.
    if (func_name == 'R' || func_name == 'r') { // here is correct
        logl_const =
                -0.5 * (double)(n_size - c_size) * (double)d_size * log(2.0 * M_PI) +
                0.5 * (double)d_size * log(fabs(determinant(XXt)));
        //std::cout << " line 2667  logl_const " << logl_const << std::endl;
    } else {
        logl_const = -0.5 * (double)n_size * (double)d_size * log(2.0 * M_PI);
        //std::cout << " line 2670  logl_const " << logl_const << std::endl;
    }
    // Optimization iterations.
    size_t i,j,t;

    for ( t = 0; t < max_iter; ++t) {
        Vg_save.value_copy(V_g);
        Ve_save.value_copy(V_e);

        step_scale = 1.0;
        step_iter = 0;
        do {
            V_g.value_copy(Vg_save);
            V_e.value_copy(Ve_save);

            // Update Vg, Ve, and invert Hessian.
            if (t != 0) {
                UpdateVgVe(Hessian_inv, gradient, step_scale, V_g, V_e);
//                std::cout << "updating Vg Ve" << std::endl;
            }
//            std::cout << "2485 Check if both Vg and Ve are positive definite." << std::endl;
            // Check if both Vg and Ve are positive definite.
            flag_pd = 1;
            V_temp.value_copy(V_e);
            Eigen_result eigen_temp = eigen(V_temp);// EigenDecomp(V_temp, U_temp, D_temp, 0);
            for ( i = 0; i < d_size; i++) {
                if (eigen_temp.get_eigen_values().get_array()[i] <= 0) {//gsl_vector_get(D_temp, i) <= 0
                    flag_pd = 0;
                }
            }

            V_temp.value_copy(V_g);
            eigen_temp = eigen(V_temp);
            for ( i = 0; i < d_size; i++) {
                if (eigen_temp.get_eigen_values().get_array()[i] <= 0) {
                    flag_pd = 0;
                }
            }
            // If flag_pd==1, continue to calculate quantities
            // and logl.
            //std::cout << "line 2509. t " << t << std::endl;
            if (flag_pd == 1) {
                CalcHiQi( eigen_r, X, n_size, d_size, c_size, dc_size, V_g, V_e, Hi_all, Qi, logdet_H, logdet_Q);
                Calc_Hiy_all(Y, Hi_all, n_size, d_size, Hiy_all);
/*
                for( int k=0; k<Hi_all.get_num_row(); ++k ){
                    for( int l=0; l<Hi_all.get_num_column(); ++l ){
                        std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " l " << l << " Hi_all: " << Hi_all.get_matrix()[k][l] << std::endl;
                    }
                }
                for( int k=0; k<Hiy_all.get_num_row(); ++k ){
                    for( int l=0; l<Hiy_all.get_num_column(); ++l ){
                        std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " l " << l << " Hiy_all: " << Hiy_all.get_matrix()[k][l] << std::endl;
                    }
                }
*/
// Calculate all xHi.
                Calc_xHi_all( X, Hi_all, n_size, d_size, c_size, xHi_all);
/*
                // Calculate QixHiy and yPy.
                for( int k=0; k<xHi_all.get_num_row(); ++k ){
                    for( int l=0; l<xHi_all.get_num_column(); ++l ){
                        std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " l " << l << " xHi_all: " << xHi_all.get_matrix()[k][l] << std::endl;
                    }
                }
                for( int k=0; k<Y.get_num_row(); ++k ){
                    for( int l=0; l<Y.get_num_column(); ++l ){
                        std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " l " << l << " Y: " << Y.get_matrix()[k][l] << std::endl;
                    }
                }
*/
                Calc_xHiy( Y, xHi_all, n_size, d_size, dc_size, xHiy);
                trmul(Qi, xHiy, QixHiy);
/*
                for( int k=0; k<Qi.get_num_row(); ++k ){
                    for( int l=0; l<Qi.get_num_column(); ++l ){
                        std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " l " << l << " Qi: " << Qi.get_matrix()[k][l] << std::endl;
                    }
                }
                for( int k=0; k<xHiy.get_length(); ++k ){
                    std::cout << "max_iter: " << max_iter << " t " << t << " k " << k << " xHiy: " << xHiy.get_array()[k] << std::endl;
                }
*/
                yPy=0.0;
                for( j=0; j<QixHiy.get_length(); ++j ){
                    yPy += QixHiy.get_array()[j]*xHiy.get_array()[j];
  //                  std::cout << "max_iter: " << max_iter << " t " << t << " line 2588 QixHiy: " << QixHiy.get_array()[j] << " xHiy " << xHiy.get_array()[j] << std::endl;
                }
 //               std::cout << "max_iter: " << max_iter << " t " << t << " line 2589 yPy: " << yPy << std::endl;
                yPy = Calc_yHiy(Y, Hiy_all, n_size, d_size) - yPy;
   //             std::cout << "max_iter: " << max_iter << " t " << t << " line 2591 yPy: " << yPy << std::endl;
                // Calculate log likelihood/restricted likelihood value.
//                std::cout << "line 2526" << std::endl;
                if (func_name == 'R' || func_name == 'r') {
                    logl_new = logl_const - 0.5 * logdet_H - 0.5 * logdet_Q - 0.5 * yPy;
//                    std::cout << "max_iter: " << max_iter << " t " << t << " line 2594 logl_new: " << logl_new <<
//                              " logl_const " << logl_const << " logdet_H " << logdet_H << " logdet_Q " << logdet_Q<< " yPy " << yPy << std::endl;
                } else {
                    logl_new = logl_const - 0.5 * logdet_H - 0.5 * yPy;
                }
            }
//            std::cout << "line 2533" << std::endl;
            step_scale /= 2.0;
            step_iter++;
        } while (
                (flag_pd == 0 || logl_new < logl_old || logl_new - logl_old > 10) &&
                step_iter < 10 && t != 0);
//        std::cout << "line 2644 Terminate if change is small. t " << t << std::endl;
        // Terminate if change is small.
        if (t != 0) {
            if (logl_new < logl_old || flag_pd == 0) {
                V_g.value_copy(Vg_save);
                V_e.value_copy(Ve_save);
                break;
            }

            if (logl_new - logl_old < max_prec) {
                break;
            }
        }

        logl_old = logl_new;
//        std::cout << "max_iter: " << max_iter << " t " << t << " line 2621 logl_new: " << logl_new << std::endl;
        CalcDev(func_name, eigen_r, Qi, Hi_all, xHi_all, Hiy_all, QixHiy, n_size,/*
                d_size, c_size, dc_size,*/ gradient,
                Hessian_inv, crt_a, crt_b, crt_c);

//        for( i=0; i<Hessian_inv.get_num_row(); ++i ){
//            for( j=0; j<Hessian_inv.get_num_column(); ++j ){
//                //Hessian_inv.get_matrix()[i][j] = -1 * Hessian_inv.get_matrix()[i][j];
//                std::cout << "line 2800 " << "i " << i << " j " << j << " Hessian_inv " << Hessian_inv.get_matrix()[i][j] << std::endl;
//            }
//        }
//        for( i=0; i<gradient.get_length(); ++i ){
//            std::cout << "line 2804 " << "i " << i << " gradient " << gradient.get_array()[i] << std::endl;
//        }
    }
    //std::cout << "Mutiply Hessian_inv with -1.0." << std::endl;
    // Mutiply Hessian_inv with -1.0.
    // Now Hessian_inv is the variance matrix.
    for( i=0; i<Hessian_inv.get_num_row(); ++i ){
        for( j=0; j<Hessian_inv.get_num_column(); ++j ){
            Hessian_inv.get_matrix()[i][j] = -1 * Hessian_inv.get_matrix()[i][j];
            //std::cout << "line 2813 i " << i << " j " << j << " Hessian_inv " << Hessian_inv.get_matrix()[i][j] << std::endl;
        }
    }
    return logl_new;
}

void MphInitial( const int & em_iter, const double & em_prec,
              const int & nr_iter, const double & nr_prec,
              const size_t & n_size, const size_t & d_size, const size_t & c_size,
              const Eigen_result & eigen_r, const My_matrix<double> & UtW,
              const My_matrix<double> & Y, const double & l_min, const double & l_max,
              const double & eps, const int & maxiter, const std::string & method,
              const int & ngrids, const Kinship_matrix_impl & k_i, My_matrix<double> & V_g, My_matrix<double> & V_e,
              My_matrix<double> & B, const My_matrix<double> & phenotype ){
    V_g.set_values_zero();
    V_e.set_values_zero();
    B.set_values_zero();

    //size_t n_s = X.get_num_row(), c_s = X.get_num_column(), d_s = Y.get_num_column();
    double lambda, vg, ve;

    // Initialize the diagonal elements of Vg and Ve using univariate
    // LMM and REML estimates.
    My_Vector<double> beta_temp(c_size);
    My_Vector<double> se_beta_temp(c_size);

    My_Vector<double> y(n_size);
    My_Vector<double> Uty(n_size);
    size_t i, j;

    for (i=0; i<d_size; ++i) {
        for(j=0; j<n_size; ++j){
            y.get_array()[j] = phenotype.get_matrix()[j][i];
        }

        trmul( eigen_r.get_eigen_vectors(), y,  Uty);
        My_matrix<double> Null_x(n_size, 0);
        My_matrix<double> tesp_x(n_size, 1);
        My_Vector<double> phe(n_size);
        for(j=0; j<n_size; ++j){
            tesp_x.get_matrix()[j][0]=1.0;
            phe.get_array()[j] = phenotype.get_matrix()[j][i];
        }

        gemma_estimates( y, Uty, Null_x, tesp_x, eigen_r, ngrids, l_min, l_max, eps, method,  maxiter, lambda);


        CalcLmmVgVe( y, eigen_r, UtW, Uty, lambda, vg, ve);
//        std::cout  << "lambda  "  << lambda << " vg: " << vg << " ve: " << ve << std::endl;
        V_g.get_matrix()[i][i] = vg;
        V_e.get_matrix()[i][i] = ve;
    }

    My_matrix<double> UltVehiY(d_size, n_size);
    My_Vector<double> D_l(d_size);
    My_matrix<double> UltVeh(d_size, d_size);
    My_matrix<double> UltVehi(d_size, d_size);

    My_matrix<double> Qi(d_size * c_size, d_size * c_size);
    My_Vector<double> XHiy(d_size * c_size);
    My_Vector<double> beta(d_size * c_size);

    XHiy.set_values_zero();

    double dl, d, delta, dx, dy;

    // Eigen decomposition

    EigenProc( V_g, V_e, d_size, D_l, UltVeh, UltVehi); // here h mean half and i means -1

    // Calculate Qi and log|Q|.
    // double logdet_Q = CalcQi(eval, D_l, X, Qi);
    CalcQi(eigen_r, D_l, n_size, d_size, c_size, UtW, Qi);
/*
    for( i=0; i<d_size*c_size; ++i ){
        for( j=0; j<d_size*c_size; ++j ){
            std::cout << "i: " <<i << " j: " << j << " QI: " << Qi.get_matrix()[i][j] << std::endl;
        }
    }
    for( i=0; i<d_size; ++i ) {
        for (j = 0; j < d_size; ++j) {
            std::cout << "i " << i << " j " << j << " UltVehi " << UltVehi.get_matrix()[i][j] << std::endl;
            std::cout << "i " << i << " j " << j << " UltVeh " << UltVeh.get_matrix()[i][j] << std::endl;
        }
    }
    */
    // Calculate UltVehiY.
    trmul(UltVehi, Y, UltVehiY);
    // calculate XHiy
    size_t k;
    for ( i = 0; i < d_size; ++i) {
        dl = D_l.get_array()[i];//gsl_vector_get(D_l, i);
        for ( j = 0; j < c_size; ++j) {
            d = 0.0;
            for ( k = 0; k < n_size; ++k) {
                delta = eigen_r.get_eigen_values().get_array()[k];//gsl_vector_get(eval, k);
                dx = UtW.get_matrix()[j][k];//gsl_matrix_get(X, j, k);
                dy = UltVehiY.get_matrix()[i][k];//gsl_matrix_get(UltVehiY, i, k);
                d += dy * dx / (delta * dl + 1.0);
                //std::cout << "line 2651 delta: " << delta << " dx: " << dx << " dy: " << dy << " another " << X.get_matrix()[j][k] << std::endl;
            }
//            std::cout << "line 2653 d: " << d << std::endl;
            XHiy.get_array()[j * d_size + i]=d;
        }
    }
    trmul(Qi, XHiy, beta);
//    for( i=0;i<d_size*c_size; ++i ){
//        std::cout << "line 2672 i: " << i << " beta.get_array()[i] " << beta.get_array()[i] << std::endl;
//    }
    // Multiply beta by UltVeh and save to B.
    B.set_values_zero();
    for ( i = 0; i < c_size; ++i) {
        for( j=0; j<B.get_num_row(); ++j){
            for( k=0; k< d_size; ++k){
                B.get_matrix()[j][i] += UltVeh.get_matrix()[k][j]*beta.get_array()[i*d_size+k];
            }
        }
    }
}

// p-value correction
// mode=1 Wald; mode=2 LRT; mode=3 SCORE;
double PCRT( const size_t & d_size, const double & p_value,
            const double & crt_a, const double & crt_b, const double & crt_c) {
    double p_crt = 0.0, chisq_crt = 0.0, q = (double)d_size;
    double chisq = critchi(p_value, (double)d_size);
    chisq_crt = chisq / (1.0 + crt_a / (2.0 * q));
    p_crt = 1- chii(chisq_crt, d_size);
    return p_crt;
}

void AnalyzePlink(const Eigen_result & eigen_r, const My_matrix<double> & UtW, const My_matrix<double> & UtY,
                  const char & method, Kinship_matrix_impl & k_i,
                  Genotype & genotype, const My_matrix<double> & phenotype, const double & man_l, const double & man_u ) {
    const double l_min = -5;
    const double l_max = 5;
    const double p_nr = 0.001;
    const size_t em_iter = 10000;
    const size_t nr_iter = 100;
    const size_t n_region = 100;

    const double em_prec = 0.0001;
    const double nr_prec = 0.0001;
    const double eps = 0.00001;
    const size_t maxiter = 1000;
    const int crt = 0;

    double logl_H0 = 0.0, logl_H1 = 0.0, p_lrt = 0;
    double crt_a, crt_b, crt_c;

    size_t c = 0;
    size_t n_size = UtY.get_num_row(), d_size = UtY.get_num_column(), c_size = UtW.get_num_column();
    size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

    My_matrix<double> Xlarge(n_size, 20000);
    Xlarge.set_values_zero();

    // Large matrices for EM.
    My_matrix<double> U_hat(d_size, n_size);
    My_matrix<double> E_hat(d_size, n_size);
    My_matrix<double> OmegaU(d_size, n_size);
    My_matrix<double> OmegaE(d_size, n_size);
    My_matrix<double> UltVehiY(d_size, n_size);
    My_matrix<double> UltVehiBX(d_size, n_size);
    My_matrix<double> UltVehiU(d_size, n_size);
    My_matrix<double> UltVehiE(d_size, n_size);

    // Large matrices for NR.
    // Each dxd block is H_k^{-1}.
    My_matrix<double> Hi_all(d_size, d_size * n_size);

    // Each column is H_k^{-1}y_k.
    My_matrix<double> Hiy_all(d_size, n_size);

    // Each dcxdc block is x_k\otimes H_k^{-1}.
    My_matrix<double> xHi_all(dc_size, d_size * n_size);
    My_matrix<double> Hessian(v_size * 2, v_size * 2);

    My_Vector<double> x(n_size);

    My_matrix<double> Y(d_size, n_size);
    My_matrix<double> X(c_size + 1, n_size);
    My_matrix<double> V_g(d_size, d_size);
    My_matrix<double> V_e(d_size, d_size);
    My_matrix<double> B(d_size, c_size + 1);
    My_Vector<double> beta(d_size);
    My_matrix<double> Vbeta(d_size, d_size);

    // Null estimates for initial values.
    My_matrix<double> V_g_null(d_size, d_size);
    My_matrix<double> V_e_null(d_size, d_size);
    My_matrix<double> B_null(d_size, c_size + 1);
    My_matrix<double> se_B_null(d_size, c_size);
    size_t i, j;
    size_t ti;
    size_t ii;
    My_matrix<double> X_sub(c_size, n_size);
    My_matrix<double> B_sub(d_size, c_size);
    My_matrix<double> xHi_all_sub(d_size * c_size, d_size * n_size);

    for( i=0; i<d_size; ++i ){
        for( j=0; j<c_size; ++j ) {
            B_sub.get_matrix()[i][j]=B.get_matrix()[i][j];
        }
    }
    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * n_size; ++j ){
            xHi_all_sub.get_matrix()[i][j] = xHi_all.get_matrix()[i][j];
        }
    }
    T_matrix(UtY, Y);
    T_matrix(UtW, X_sub);
/*
    for( i=0; i<c_size; ++i ){
        for( j=0; j<n_size; ++j ){
            X.get_matrix()[i][j] = UtW.get_matrix()[j][i];
        }
    }*/
    for( i=0; i<X.get_num_column(); ++i ){
        X.get_matrix()[c_size][i]=0.0;
    }
    for( i=0; i<B.get_num_row(); ++i ){
        B.get_matrix()[i][c_size]=0.0;
    }

    MphInitial( em_iter, em_prec, nr_iter, nr_prec, n_size, d_size, c_size, eigen_r, UtW,
                Y, l_min, l_max, eps, maxiter, "ML", n_region, k_i, V_g, V_e, B, phenotype);  // have been well tested

    logl_H0 = MphEM('R', em_iter, em_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, U_hat, E_hat,
                    OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                    V_e, B);

    logl_H0 = MphNR('R', nr_iter, nr_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, Hi_all,
                    xHi_all_sub, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                    crt_c);

    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * c_size; ++j ){
            xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
        }
    }
    MphCalcBeta(eigen_r, X_sub, n_size, d_size, c_size, Y, V_g, V_e, UltVehiY, B,
                se_B_null);
    c = 0;
    std::vector<double> Vg_remle_null, Ve_remle_null, Vg_mle_null, Ve_mle_null;
    std::vector<double> VVg_remle_null, VVe_remle_null, VVg_mle_null;
    std::vector<double> beta_remle_null, se_beta_remle_null, beta_mle_null;
    std::vector<double> VVe_mle_null;
    std::vector<double> se_beta_mle_null;

    Vg_remle_null.clear();
    Ve_remle_null.clear();
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            Vg_remle_null.push_back(V_g.get_matrix()[i][j]);
            Ve_remle_null.push_back(V_e.get_matrix()[i][j]);
            VVg_remle_null.push_back(Hessian.get_matrix()[c][c]);
            VVe_remle_null.push_back(Hessian.get_matrix()[c + v_size][c + v_size]);
            c++;
        }
    }
    beta_remle_null.clear();
    se_beta_remle_null.clear();
    for ( i = 0; i < se_B_null.get_num_row(); ++i) {
        for ( j = 0; j < se_B_null.get_num_column(); ++j) {
            beta_remle_null.push_back(B.get_matrix()[i][j]);
            se_beta_remle_null.push_back(se_B_null.get_matrix()[i][j]);
        }
    }
//    double logl_remle_H0 = logl_H0;

    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.precision(10);
    std::cout << "REMLE estimate for Vg in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_g.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Vg): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c][c]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "REMLE estimate for Ve in the null model: " << std::endl;
    for ( i = 0; i < d_size; i++) {
        for ( j = 0; j <= i; j++) {
            std::cout << V_e.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Ve): " << std::endl;
    for ( i = 0; i < d_size; i++) {
        for ( j = 0; j <= i; j++) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c + v_size][c + v_size]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "REMLE likelihood = " << logl_H0 <<std::endl;

    logl_H0 = MphEM('L', em_iter, em_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, U_hat, E_hat,
                    OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                    V_e, B);
    logl_H0 = MphNR('L', nr_iter, nr_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, Hi_all,
                    xHi_all_sub, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                    crt_c);
    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * c_size; ++j ){
            xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
        }
    }

    MphCalcBeta(eigen_r, X_sub, n_size, d_size, c_size, Y, V_g, V_e, UltVehiY, B,
                se_B_null);

    c = 0;
    Vg_mle_null.clear();
    Ve_mle_null.clear();
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            Vg_mle_null.push_back(V_g.get_matrix()[i][j]);
            Ve_mle_null.push_back(V_e.get_matrix()[i][j]);
            VVg_mle_null.push_back(Hessian.get_matrix()[c][c]);
            VVe_mle_null.push_back(Hessian.get_matrix()[c + v_size][c + v_size]);
            c++;
        }
    }
    beta_mle_null.clear();
    se_beta_mle_null.clear();
    for ( i = 0; i < se_B_null.get_num_row(); ++i) {
        for ( j = 0; j < se_B_null.get_num_column(); ++j) {
            beta_mle_null.push_back(B.get_matrix()[i][j]);
            se_beta_mle_null.push_back(se_B_null.get_matrix()[i][j]);
        }
    }
//    double logl_mle_H0 = logl_H0;

    std::cout << "MLE estimate for Vg in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_g.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Vg): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c][c]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "MLE estimate for Ve in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_e.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Ve): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c + v_size][c + v_size]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "MLE likelihood = " << logl_H0 << std::endl;

    std::vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
    for ( i = 0; i < d_size; i++) {
        v_beta.push_back(0.0);
    }
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            v_Vg.push_back(0.0);
            v_Ve.push_back(0.0);
            v_Vbeta.push_back(0.0);
        }
    }
    V_g_null.value_copy(V_g);
    V_e_null.value_copy(V_e);
    B_null.value_copy(B);
    //return;
    // Start reading genotypes and analyze.
    // Calculate n_bit and c, the number of bit for each snp.
    for( i=0; i<n_size; ++i ){
        for( j=0; j<c_size; ++j ){
            X.get_matrix()[j][i] = X_sub.get_matrix()[j][i];
        }
    }
    double sum;
    double count;
    double mean;
    bool has_missing;
    int32_t number_1;
    std::vector <int> indexs;
    size_t i_index;
    for (size_t t = 0; t < genotype.get_number_of_variant(); ++t) {
        indexs.clear();
        sum=0;
        count=0;
        number_1=0;
        has_missing=false;

        for( i=0; i< n_size; ++i ){
            if( genotype.get_genotype_matrix().get_matrix()[i][t] == missing_genotype ){
                indexs.push_back(i);
                has_missing=true;
            }else{
                sum += genotype.get_genotype_matrix().get_matrix()[i][t];
                count++;
                if( 1 == genotype.get_genotype_matrix().get_matrix()[i][t] ){
                    ++number_1;
                }
            }
            x.get_array()[i] = genotype.get_genotype_matrix().get_matrix()[i][t];
        }
        if( number_1>= man_l && number_1<=man_u) {
            if(has_missing){
                mean=sum/count;
                for( i_index=0;i_index<indexs.size();++i_index){
                    x.get_array()[indexs[i_index]] = mean;
                }
            }
            // Initial values.
            V_g.value_copy(V_g_null);
            V_e.value_copy(V_e_null);
            B.value_copy(B_null);

            My_Vector<double> X_row(X.get_num_column());
            trmul(eigen_r.get_eigen_vectors(), x, X_row);
            for (i = 0; i < n_size; ++i) {
                X.get_matrix()[c_size][i] = X_row.get_array()[i];
            }

            logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eigen_r, X, Y, n_size, d_size, c_size + 1,
                            U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                            UltVehiE, V_g, V_e, B);

//            p_lrt = MphCalcP(eigen_r, X_row, n_size, d_size, c_size, X_sub, Y, V_g, V_e,
//                             UltVehiY, beta, Vbeta);

            p_lrt = 1 - chii(2.0 * (logl_H1 - logl_H0), d_size);
            if (p_lrt < p_nr) {
                logl_H1 = MphNR('L', nr_iter / 10, nr_prec * 10, eigen_r, X, Y, n_size, d_size, c_size + 1,
                                Hi_all, xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

                for (i = 0; i < d_size * c_size; ++i) {
                    for (j = 0; j < d_size * c_size; ++j) {
                        xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
                    }
                }
                // Calculate beta and Vbeta.

//                p_lrt = MphCalcP(eigen_r, X_row, n_size, d_size, dc_size,
//                                 X, Y, V_g, V_e, UltVehiY, beta, Vbeta);
                p_lrt = 1 - chii(2.0 * (logl_H1 - logl_H0), d_size);
                //p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);
                if (crt == 1) {
                    p_lrt = PCRT(d_size, p_lrt, crt_a, crt_b, crt_c);
                }
            }

            // Store summary data.
            for (i = 0; i < d_size; ++i) {
                v_beta[i] = beta.get_array()[i];
            }

            c = 0;
            for (i = 0; i < d_size; ++i) {
                for (j = i; j < d_size; ++j) {
                    v_Vg[c] = V_g.get_matrix()[i][j];
                    v_Ve[c] = V_e.get_matrix()[i][j];
                    v_Vbeta[c] = Vbeta.get_matrix()[i][j];
                    c++;
                }
            }
            std::cout << genotype.get_variant_Vector()[t].getChromosome() << " " <<
                      genotype.get_variant_Vector()[t].getPosition() << " " << genotype.get_variant_Vector()[t].getId()
                      << "  " << p_lrt << std::endl;
        }
    }
}






















































































void AnalyzePlink_MultiAllic(const Eigen_result & eigen_r, const My_matrix<double> & UtW, const My_matrix<double> & UtY,
                             const char & method, Kinship_matrix_impl & k_i,
                             Genotype & genotype, const My_matrix<double> & phenotype, const double & man_l, const double & man_u ) {// todo the freedom to get p-value is not correct
    const double l_min = -5;
    const double l_max = 5;
    const double p_nr = 0.001;
    const size_t em_iter = 10000;
    const size_t nr_iter = 100;
    const size_t n_region = 100;

    const double em_prec = 0.0001;
    const double nr_prec = 0.0001;
    const double eps = 0.00001;
    const size_t maxiter=1000;
    const int crt = 0;

    double logl_H0 = 0.0, logl_H1 = 0.0, p_lrt = 0;
    double crt_a, crt_b, crt_c;

    size_t c = 0;
    size_t n_size = UtY.get_num_row(), d_size = UtY.get_num_column(), c_size = UtW.get_num_column();
    size_t dc_size = d_size * (c_size + 1), v_size = d_size * (d_size + 1) / 2;

    My_matrix<double> Xlarge(n_size, 20000);
    Xlarge.set_values_zero();

    // Large matrices for EM.
    My_matrix<double> U_hat(d_size, n_size);
    My_matrix<double> E_hat(d_size, n_size);
    My_matrix<double> OmegaU(d_size, n_size);
    My_matrix<double> OmegaE(d_size, n_size);
    My_matrix<double> UltVehiY(d_size, n_size);
    My_matrix<double> UltVehiBX(d_size, n_size);
    My_matrix<double> UltVehiU(d_size, n_size);
    My_matrix<double> UltVehiE(d_size, n_size);

    // Large matrices for NR.
    // Each dxd block is H_k^{-1}.
    My_matrix<double> Hi_all(d_size, d_size * n_size);

    // Each column is H_k^{-1}y_k.
    My_matrix<double> Hiy_all(d_size, n_size);

    // Each dcxdc block is x_k\otimes H_k^{-1}.
    My_matrix<double> xHi_all(dc_size, d_size * n_size);
    My_matrix<double> Hessian(v_size * 2, v_size * 2);

    My_matrix<double> Y(d_size, n_size);
//    My_matrix<double> X(c_size + 1, n_size);
    My_matrix<double> V_g(d_size, d_size);
    My_matrix<double> V_e(d_size, d_size);
    My_matrix<double> B(d_size, c_size + 1);
    My_Vector<double> beta(d_size);
    My_matrix<double> Vbeta(d_size, d_size);

    // Null estimates for initial values.
    My_matrix<double> V_g_null(d_size, d_size);
    My_matrix<double> V_e_null(d_size, d_size);
    My_matrix<double> B_null(d_size, c_size + 1);
    My_matrix<double> se_B_null(d_size, c_size);
    size_t i, j;
    size_t ti;
    size_t ii;
    My_matrix<double> X_sub(c_size, n_size);
    My_matrix<double> B_sub(d_size, c_size);
    My_matrix<double> xHi_all_sub(d_size * c_size, d_size * n_size);

    for( i=0; i<d_size; ++i ){
        for( j=0; j<c_size; ++j ) {
            B_sub.get_matrix()[i][j]=B.get_matrix()[i][j];
        }
    }
    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * n_size; ++j ){
            xHi_all_sub.get_matrix()[i][j] = xHi_all.get_matrix()[i][j];
        }
    }
    T_matrix(UtY, Y);
    T_matrix(UtW, X_sub);

    for( i=0; i<B.get_num_row(); ++i ){
        B.get_matrix()[i][c_size]=0.0;
    }

    MphInitial( em_iter, em_prec, nr_iter, nr_prec, n_size, d_size, c_size, eigen_r, UtW,
                Y, l_min, l_max, eps, maxiter, "ML", n_region, k_i, V_g, V_e, B, phenotype);  // have been well tested

    logl_H0 = MphEM('R', em_iter, em_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, U_hat, E_hat,
                    OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                    V_e, B);

    logl_H0 = MphNR('R', nr_iter, nr_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, Hi_all,
                    xHi_all_sub, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                    crt_c);

    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * c_size; ++j ){
            xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
        }
    }
    MphCalcBeta(eigen_r, X_sub, n_size, d_size, c_size, Y, V_g, V_e, UltVehiY, B,
                se_B_null);
    c = 0;
    std::vector<double> Vg_remle_null, Ve_remle_null, Vg_mle_null, Ve_mle_null;
    std::vector<double> VVg_remle_null, VVe_remle_null, VVg_mle_null;
    std::vector<double> beta_remle_null, se_beta_remle_null, beta_mle_null;
    std::vector<double> VVe_mle_null;
    std::vector<double> se_beta_mle_null;

    Vg_remle_null.clear();
    Ve_remle_null.clear();
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            Vg_remle_null.push_back(V_g.get_matrix()[i][j]);
            Ve_remle_null.push_back(V_e.get_matrix()[i][j]);
            VVg_remle_null.push_back(Hessian.get_matrix()[c][c]);
            VVe_remle_null.push_back(Hessian.get_matrix()[c + v_size][c + v_size]);
            c++;
        }
    }
    beta_remle_null.clear();
    se_beta_remle_null.clear();
    for ( i = 0; i < se_B_null.get_num_row(); ++i) {
        for ( j = 0; j < se_B_null.get_num_column(); ++j) {
            beta_remle_null.push_back(B.get_matrix()[i][j]);
            se_beta_remle_null.push_back(se_B_null.get_matrix()[i][j]);
        }
    }
//    double logl_remle_H0 = logl_H0;

    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout.precision(10);
    std::cout << "REMLE estimate for Vg in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_g.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Vg): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c][c]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "REMLE estimate for Ve in the null model: " << std::endl;
    for ( i = 0; i < d_size; i++) {
        for ( j = 0; j <= i; j++) {
            std::cout << V_e.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Ve): " << std::endl;
    for ( i = 0; i < d_size; i++) {
        for ( j = 0; j <= i; j++) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c + v_size][c + v_size]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "REMLE likelihood = " << logl_H0 <<std::endl;

    logl_H0 = MphEM('L', em_iter, em_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, U_hat, E_hat,
                    OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU, UltVehiE, V_g,
                    V_e, B);
    logl_H0 = MphNR('L', nr_iter, nr_prec, eigen_r, X_sub, Y, n_size, d_size, c_size, Hi_all,
                    xHi_all_sub, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b,
                    crt_c);
    for( i=0; i<d_size * c_size; ++i ){
        for( j=0; j<d_size * c_size; ++j ){
            xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
        }
    }

    MphCalcBeta(eigen_r, X_sub, n_size, d_size, c_size, Y, V_g, V_e, UltVehiY, B,
                se_B_null);

    c = 0;
    Vg_mle_null.clear();
    Ve_mle_null.clear();
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            Vg_mle_null.push_back(V_g.get_matrix()[i][j]);
            Ve_mle_null.push_back(V_e.get_matrix()[i][j]);
            VVg_mle_null.push_back(Hessian.get_matrix()[c][c]);
            VVe_mle_null.push_back(Hessian.get_matrix()[c + v_size][c + v_size]);
            c++;
        }
    }
    beta_mle_null.clear();
    se_beta_mle_null.clear();
    for ( i = 0; i < se_B_null.get_num_row(); ++i) {
        for ( j = 0; j < se_B_null.get_num_column(); ++j) {
            beta_mle_null.push_back(B.get_matrix()[i][j]);
            se_beta_mle_null.push_back(se_B_null.get_matrix()[i][j]);
        }
    }
//    double logl_mle_H0 = logl_H0;

    std::cout << "MLE estimate for Vg in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_g.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Vg): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c][c]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "MLE estimate for Ve in the null model: " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            std::cout << V_e.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "se(Ve): " << std::endl;
    for ( i = 0; i < d_size; ++i) {
        for ( j = 0; j <= i; ++j) {
            c = GetIndex(i, j, d_size);
            std::cout << sqrt(Hessian.get_matrix()[c + v_size][c + v_size]) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "MLE likelihood = " << logl_H0 << std::endl;

    std::vector<double> v_beta, v_Vg, v_Ve, v_Vbeta;
    for ( i = 0; i < d_size; i++) {
        v_beta.push_back(0.0);
    }
    for ( i = 0; i < d_size; ++i) {
        for ( j = i; j < d_size; ++j) {
            v_Vg.push_back(0.0);
            v_Ve.push_back(0.0);
            v_Vbeta.push_back(0.0);
        }
    }
    V_g_null.value_copy(V_g);
    V_e_null.value_copy(V_e);
    B_null.value_copy(B);
    //return;
    // Start reading genotypes and analyze.
    // Calculate n_bit and c, the number of bit for each snp.

    double sum;
    double count;
    double mean;
    bool has_missing;
    std::vector <int> indexs;
    size_t i_index;
    int size, p, q;
    int j_1;
    double this_this_state;
    for (size_t t = 0; t < genotype.get_number_of_variant(); ++t) {
        if (genotype.get_variant_Vector()[t].getId()[genotype.get_variant_Vector()[t].getId().size() - 1] == 'b') {

        } else {
            std::map<int, int> allThisStates;
            for (i = 0; i < n_size; ++i) {
                if (allThisStates.find(genotype.get_genotype_matrix().get_matrix()[i][t]) == allThisStates.end()) {
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][t]] = 1;
                } else {
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][t]] =
                            allThisStates[genotype.get_genotype_matrix().get_matrix()[i][t]] + 1;
                }
            }

            std::vector<int> allThisStates_vector;
            for (std::map<int, int>::iterator ittt = allThisStates.begin(); ittt != allThisStates.end(); ++ittt) {
                if (ittt->second >= man_l && ittt->second <= man_u && ittt->second != missing_genotype) {
                    allThisStates_vector.push_back(ittt->first);
                }
            }
            if (allThisStates_vector.size() > 0) {
                size = allThisStates_vector.size();
                q = c_size + size - 1; // here minute 1 is because we do not take all the allele as covariants
                p = n_size - q;

                My_matrix<double> full_x(n_size, size-1);
                for (i = 0; i < n_size; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (i = 0; i < n_size; ++i) {
                    for (j_1 = 0; j_1 < (size - 1); ++j_1) { // left one allele
                        if (allThisStates_vector[j_1] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                            this_this_state = 1.0;
                        } else {
                            this_this_state = 0.0;
                        }
                        full_x.get_matrix()[i][j_1 + 1] = this_this_state;
                    }
                }

                My_matrix<double> X(n_size, q);

                // Initial values.
                V_g.value_copy(V_g_null);
                V_e.value_copy(V_e_null);
                B.value_copy(B_null);
                X.set_values_zero();
//                My_Vector<double> X_row(X.get_num_column());
                for( int s1=0; s1<n_size; ++s1 ){
                    for( int s2=0; s2<full_x.get_num_column(); ++s2 ){
                        for( int s3=0; s3<n_size; ++s3 ){
                            X.get_matrix()[s1][s2] += eigen_r.get_eigen_vectors().get_matrix()[s1][s3] * full_x.get_matrix()[s2][s2];
                        }
                    }
                }
                for( int s1=0; s1<n_size; ++s1 ){
                    for( int s2=0; s2<UtW.get_num_column(); ++s2 ){
                        X.get_matrix()[s1][full_x.get_num_column()+s2]=UtW.get_matrix()[s1][s2];
                    }
                }

                logl_H1 = MphEM('L', em_iter / 10, em_prec * 10, eigen_r, X, Y, n_size, d_size, c_size + 1,
                                U_hat, E_hat, OmegaU, OmegaE, UltVehiY, UltVehiBX, UltVehiU,
                                UltVehiE, V_g, V_e, B);


                p_lrt = 1 - chii(2.0 * (logl_H1 - logl_H0), d_size);
                if (p_lrt < p_nr) {
                    logl_H1 = MphNR('L', nr_iter / 10, nr_prec * 10, eigen_r, X, Y, n_size, d_size, c_size + 1,
                                    Hi_all, xHi_all, Hiy_all, V_g, V_e, Hessian, crt_a, crt_b, crt_c);

                    for (i = 0; i < d_size * c_size; ++i) {
                        for (j = 0; j < d_size * c_size; ++j) {
                            xHi_all.get_matrix()[i][j] = xHi_all_sub.get_matrix()[i][j];
                        }
                    }
                    // Calculate beta and Vbeta.


                    p_lrt = 1 - chii(2.0 * (logl_H1 - logl_H0), d_size);
                    //p_lrt = gsl_cdf_chisq_Q(2.0 * (logl_H1 - logl_H0), (double)d_size);
                    if (crt == 1) {
                        p_lrt = PCRT(d_size, p_lrt, crt_a, crt_b, crt_c);
                    }
                }

                // Store summary data.
                for (i = 0; i < d_size; ++i) {
                    v_beta[i] = beta.get_array()[i];
                }

                c = 0;
                for (i = 0; i < d_size; ++i) {
                    for (j = i; j < d_size; ++j) {
                        v_Vg[c] = V_g.get_matrix()[i][j];
                        v_Ve[c] = V_e.get_matrix()[i][j];
                        v_Vbeta[c] = Vbeta.get_matrix()[i][j];
                        c++;
                    }
                }
                std::cout << genotype.get_variant_Vector()[t].getChromosome() << " " <<
                          genotype.get_variant_Vector()[t].getPosition() << " " << genotype.get_variant_Vector()[t].getId()
                          << "  " << p_lrt << std::endl;
            }
        }
    }
}
