//
// Created by baoxing on 2/21/18.
//

#include "My_linear_algebra_impl.h"

//matrix operation begin

double determinant( const My_matrix<double> & _x ){
    My_matrix<double> x(_x);
    int16_t i, j, k, is, js, n = x.get_num_column();
    double f=1.0, det=1.0, q, d;
    for( k=0; k<(n-1); ++k ){
        q=0.0;
        for( i=k; i<n;++i ){
            for( j=k; j<n;++j ){
                d=fabs(x.get_matrix()[i][j]);
                if( d>q ){
                    q=d;
                    is=i;
                    js=j;
                }
            }
        }
        if( q == 0.0 ){
            det=0.0;
            return det;
        }
        if( is != k ){
            f = -f;
            for( j=k; j<n; ++j ){
                d=x.get_matrix()[k][j];
                x.get_matrix()[k][j]=x.get_matrix()[is][j];
                x.get_matrix()[is][j]=d;
            }
        }
        if( js != k ){
            f = -f;
            for( i=k; i<n; ++i ){
                d=x.get_matrix()[i][js];
                x.get_matrix()[i][js]=x.get_matrix()[i][k];
                x.get_matrix()[i][k]=d;
            }
        }
        det = det * x.get_matrix()[k][k];
        for( i=k+1;i<n;++i ){
            d = x.get_matrix()[i][k]/x.get_matrix()[k][k];
            for( j=k+1;j<n;++j ){
                x.get_matrix()[i][j] -= d*x.get_matrix()[k][j];
            }
        }
    }
    det = f * det * x.get_matrix()[n-1][n-1];
    return det;
}

// inverse matrix
void inverse_matrix(My_matrix<double>& c){
//    assert( c.get_num_column() == c.get_num_row());
    int *is, *js, i, j, k, n=c.get_num_column();
    double d, p;
    is = new int[n];
    js = new int[n];
    for( k=0; k<n; ++k ){
        d=0.0;
        for( i=k; i<n; ++i ) {
            for (j = k; j < n; ++j) {
                p = fabs(c.get_matrix()[i][j]);
                if (p > d) {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }
        if( d == 0.0 && k != (n-1) ){
            delete [] is;
            delete [] js;
            std::cerr << ("everything is zero, error, matrix is not invertible") << std::endl;
            return;
            //exit(1);
        }
        if( is[k] != k ){
            for( j=0; j<n; ++j ){
                p = c.get_matrix()[k][j];
                c.get_matrix()[k][j] = c.get_matrix()[is[k]][j];
                c.get_matrix()[is[k]][j] = p;
            }
        }

        if( js[k] != k ){
            for( i=0; i<n; ++i ){
                p = c.get_matrix()[i][k];
                c.get_matrix()[i][k] = c.get_matrix()[i][js[k]];
                c.get_matrix()[i][js[k]] = p;
            }
        }

        c.get_matrix()[k][k]=1.0/c.get_matrix()[k][k];
        for( j=0; j<n; ++j ){
            if( j!=k ){
                c.get_matrix()[k][j]=c.get_matrix()[k][j]*c.get_matrix()[k][k];
            }
        }
        for( i=0; i<n; ++i){
            if( i != k ){
                for ( j=0; j<n; ++j ){
                    if( j != k ){
                        c.get_matrix()[i][j]=c.get_matrix()[i][j]-c.get_matrix()[i][k]*c.get_matrix()[k][j];
                    }
                }
            }
        }
        for( i=0; i<n; ++i ){
            if( i !=k ){
                c.get_matrix()[i][k]=-c.get_matrix()[i][k]*c.get_matrix()[k][k];
            }
        }
    }
    for( k=n-1; k>=0; --k ){
        if( js[k] != k ){
            for( j=0; j<=n-1; j++ ){
                p=c.get_matrix()[k][j];
                c.get_matrix()[k][j]=c.get_matrix()[js[k]][j];
                c.get_matrix()[js[k]][j]=p;
            }
        }
        if( is[k] != k ){
            for ( i=0; i<n; ++i ){
                p=c.get_matrix()[i][k];
                c.get_matrix()[i][k]=c.get_matrix()[i][is[k]];
                c.get_matrix()[i][is[k]]=p;
            }
        }
    }
    delete[] is;
    delete[] js;
}

// QR decomposition
// a[m][n] q[m][m] is the Q of QR decomposition
Qr_decomposition_result qr_decomposition( const My_matrix<double> & _a){

    My_matrix<double> a(_a);

    int i, j, k, nn, jj;
    double u, alpha, w, t;
    if( a.get_num_row() < a.get_num_column() ){
        std::cout << "QR decomposition failed" << std::endl;
        exit(1);
    }
    My_matrix<double> q(a.get_num_row(), a.get_num_row());

    for( i=0; i<a.get_num_row(); ++i ){
        for( j=0; j<a.get_num_row(); ++j ){
            q.get_matrix()[i][j] = 0.0;
            if( i == j ){
                q.get_matrix()[i][j] = 1.0;
            }
        }
    }
    nn=a.get_num_column();
    if( a.get_num_row() == a.get_num_column() ){
        nn = a.get_num_row() -1;
    }
    for ( k=0; k<nn; ++k){
        u=0.0;
        for ( i=k; i<a.get_num_row(); ++i ){
            w = fabs(a.get_matrix()[i][k]);
            if( w > u ){
                u = w;
            }
        }
        alpha = 0.0;
        for( i=k; i<a.get_num_row(); ++i ){
            t=a.get_matrix()[i][k]/u;
            alpha = alpha + t * t;
        }
        if( a.get_matrix()[k][k] > 0.0 ){
            u =-u;
        }
        alpha = u * sqrt(alpha);
        if( fabs(alpha) == 0.0 ){
            std::cout << "QR decomposition failed" << std::endl;
            exit(1);
        }
        u = sqrt(2.0*alpha*(alpha-a.get_matrix()[k][k]));
        if( u != 0.0 ){
            a.get_matrix()[k][k] = (a.get_matrix()[k][k]-alpha)/u;
            for( i=k+1; i<a.get_num_row(); ++i ){
                a.get_matrix()[i][k] = a.get_matrix()[i][k]/u;
            }
            for( j=0; j<a.get_num_row(); ++j ){
                t=0.0;
                for( jj=k; jj<a.get_num_row(); ++jj ){
                    t=t+a.get_matrix()[jj][k]*q.get_matrix()[jj][j];
                }
                for( i=k; i<a.get_num_row(); ++i ){
                    q.get_matrix()[i][j] = q.get_matrix()[i][j] - 2.0*t*a.get_matrix()[i][k];
                }
            }
            for( j=k+1; j<a.get_num_column();++j ){
                t=0.0;
                for( jj=k; jj<a.get_num_row(); ++jj ){
                    t = t + a.get_matrix()[jj][k]*a.get_matrix()[jj][j];
                }
                for( i=k; i<a.get_num_row(); ++i ){
                    a.get_matrix()[i][j]=a.get_matrix()[i][j]-2.0*t*a.get_matrix()[i][k];
                }
            }
            a.get_matrix()[k][k]=alpha;
            for( i=k+1; i<a.get_num_row(); ++i ){
                a.get_matrix()[i][k]=0.0;
            }
        }
    }
    for( i=0; i<(a.get_num_row()-1); ++i ){
        for( j=i+1; j<a.get_num_row(); ++j ){
            t=q.get_matrix()[i][j];
            q.get_matrix()[i][j]=q.get_matrix()[j][i];
            q.get_matrix()[j][i]=t;
        }
    }

    Qr_decomposition_result qr_decomposition_result(a, q);
    return qr_decomposition_result;
}
//
////eigen value from large to small
//// every line is a eigen vector
//Eigen_result eigen_2( const My_matrix<double> & _a){
//    assert(_a.get_num_row() == _a.get_num_column());
//    Eigen::MatrixXd A(_a.get_num_column(), _a.get_num_row());
//    int i, j;
//    for( i=0; i <_a.get_num_column(); ++i ){
//        for( j=0; j <_a.get_num_column(); ++j ) {
//            A(i,j)=_a.get_matrix()[i][j];
//        }
//    }
//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
//    My_Vector<double> eigen_values(_a.get_num_column());
//    My_matrix<double> v_o(_a.get_num_column(), _a.get_num_column());
//    for( i=0; i <_a.get_num_column(); ++i ){
//        eigen_values.get_array()[i] = es.eigenvalues()(i);
//        for( j=0; j <_a.get_num_column(); ++j ) {
//            v_o.get_matrix()[i][j]=es.eigenvectors().col(i)[j];
//        }
//    }
//    Eigen_result eigen_result(eigen_values, v_o);
//    return eigen_result;
//}



//
////Jacobi method
//Eigen_result eigen( const My_matrix & _a, const double & eps){ // it is too slow
//    My_matrix a(_a);
//    assert(a.get_num_row() == a.get_num_column());
//    My_matrix v(a.get_num_row(), a.get_num_row());
//    int i, j, p, q;
//    double ff, fm , cn, sn, omega, x, y, d;
//    v.get_matrix()[0][0]=1.0;
//    for ( i=1; i< a.get_num_row(); ++i ){
//        v.get_matrix()[i][i]=1.0;
//        for ( j=0; j<i; ++j ){
//            v.get_matrix()[i][j] = 0.0;
//            v.get_matrix()[j][i] = 0.0;
//        }
//    }
//    ff = 0.0;
//    for ( i=1; i< a.get_num_row(); ++i  ){
//        for( j=0; j<a.get_num_row(); ++j ){
//            d = a.get_matrix()[i][j];
//            ff += d*d;
//        }
//    }
//    ff = sqrt(2.0*ff);
//    loop0:
//    ff = ff/(1.0*a.get_num_row());
//    loop1:
//    for( i=1; i< a.get_num_row(); ++i ){
//        for( j=0; j<i; ++j ){
//            d = fabs(a.get_matrix()[i][j]);
//            if( d > ff ){
//                p =i;
//                q = j;
//                goto loop;
//            }
//        }
//    }
//    if(ff < eps){
//        My_Vector eigen_values(a.get_num_row());
//        int i, j;
//        for( i=0; i<a.get_num_row(); ++i ){
//            eigen_values.get_array()[i]= a.get_matrix()[i][i];
//        }
//        My_Vector eigen_values_o(eigen_values);
//        My_matrix v_o(v);
//
//        std::sort(eigen_values.get_array(), eigen_values.get_array() + eigen_values.get_length(), greater());
//        std::vector<int> order;
//        std::set<int> includedIndexs;
//        for( i=0; i<a.get_num_row(); ++i ){
//            for( j=0; j<a.get_num_row(); ++j ){
//                if( (eigen_values.get_array()[j] == eigen_values_o.get_array()[i]) && includedIndexs.find(j)==includedIndexs.end() ){
//                    order.push_back(j);
//                    includedIndexs.insert(j);
//                }
//            }
//        }
//        assert(order.size()== a.get_num_row());
//        for( i=0; i<a.get_num_row(); ++i ){
//            for( j=0; j<a.get_num_row(); ++j ) {
//                v.get_matrix()[order[j]][i] = v_o.get_matrix()[i][j]; // change order and transform
//            }
//        }
//        Eigen_result eigen_result(eigen_values, v);
//        return eigen_result;
//    }else{
//        goto loop0;
//    }
//    loop:
//    x = -a.get_matrix()[p][q];
//    y = (a.get_matrix()[q][q] - a.get_matrix()[p][p])/2.0;
//    omega = x/sqrt(x*x + y *y);
//    if ( y <0 ){
//        omega = - omega;
//    }
//    sn = 1.0+sqrt(1.0-omega*omega);
//    sn = omega/sqrt(2.0*sn);
//    cn = sqrt (1.0-sn*sn);
//    fm = a.get_matrix()[p][p];
//    a.get_matrix()[p][p] = fm*cn*cn+a.get_matrix()[q][q]*sn*sn+a.get_matrix()[p][q]*omega;
//    a.get_matrix()[q][q] = fm*sn*sn+a.get_matrix()[q][q]*cn*cn-a.get_matrix()[p][q]*omega;
//    a.get_matrix()[p][q] = 0.0;
//    a.get_matrix()[q][p] = 0.0;
//    for ( j=0; j<a.get_num_row(); ++j ){
//        if ( (j!=p) && (j!=q)){
//            fm = a.get_matrix()[p][j];
//            a.get_matrix()[p][j]=fm*cn+a.get_matrix()[q][j]*sn;
//            a.get_matrix()[q][j]=-fm*sn+a.get_matrix()[q][j]*cn;
//        }
//    }
//    for( i=0; i<a.get_num_row(); ++i ){
//        if( i!=p && i!=q ){
//            fm=a.get_matrix()[i][p];
//            a.get_matrix()[i][p]=fm*cn+a.get_matrix()[i][q]*sn;
//            a.get_matrix()[i][q]=-fm*sn+a.get_matrix()[i][q]*cn;
//        }
//
//
//    }
//    for ( i=0; i<=a.get_num_row()-1; ++i ){
//        fm=v.get_matrix()[i][p];
//        v.get_matrix()[i][p]=fm*cn+v.get_matrix()[i][q]*sn;
//        v.get_matrix()[i][q]=-fm*sn+v.get_matrix()[i][q]*cn;
//    }
//    goto loop1;
//}


// for a equation y=bx
//y is a vector with size n   x is a two dimension matrix with size n*m
// b is a vector with size m, when return it is the beta values of the function

void lsq(  const My_matrix<double> & x, const  My_Vector<double> & y, My_Vector<double> & b ){
    unsigned long n = y.get_length();
    unsigned long m = x.get_num_column();
    My_matrix<double> x_t(m, n);
    T_matrix(x, x_t);
    My_matrix<double> xt_x(m, m);
    trmul<double>(x_t, x, xt_x);
    inverse_matrix(xt_x);
    My_matrix<double> xt_x_inverse_xt(m, n);
    trmul<double>(xt_x, x_t, xt_x_inverse_xt);
    int i, j;
    for( i=0; i<m; ++i ){
        b.get_array()[i]=0;
        for( j=0; j<n; ++j ){
            b.get_array()[i] += xt_x_inverse_xt.get_matrix()[i][j]*y.get_array()[j];
        }
    }
}



// the result of lq2 is different with R, maybe there is something wrong
void lsq2(  const My_matrix<double> & x, const  My_Vector<double> & y, My_Vector<double> & b ){

    unsigned long n = y.get_length();
    unsigned long m = x.get_num_column();
//    assert(b.get_length() == m);
//    assert(y.get_length() == x.get_num_row());

    double * y_t = new double[y.get_length()];
    memcpy(y_t, y.get_array(), n * sizeof(double));
    int i, j;
    double d;

    double *c = new double[x.get_num_column()];
    Qr_decomposition_result qr_decomposition_result = qr_decomposition(x);
    My_matrix<double> q = qr_decomposition_result.get_q();
    My_matrix<double> r = qr_decomposition_result.get_r();

    for( i=0; i<x.get_num_column(); ++i ){
        d=0.0;
        for( j=0; j<n; ++j ){
            d += q.get_matrix()[i][j]*y_t[j];
        }
        c[i]=d;
    }
    y_t[m-1]=c[m-1]/r.get_matrix()[m-1][m-1];
    for ( i=m-2; i>=0; i-- ){
        d=0.0;
        for(j=i+1; j<=m-1; ++j){
            d=d+r.get_matrix()[i][j]*y_t[j];
        }
        y_t[i]=(c[i]-d)/r.get_matrix()[i][i];
    }
    for ( i=0; i<m; ++i ){
        b.get_array()[i]=y_t[i];
    }
    delete[] c;
    delete [] y_t;

}

void strq(const My_matrix<double> & a, double ** q, double *b, double * c) {
    int n = a.get_num_column();
    int i, j, k, v;
    double h, f, g, h2;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++) {
            q[i][j]= a.get_matrix()[i][j];
        }
    }
    for (i=n-1; i>0; i--) {
        h = 0.0;
        if (i > 1) {
            for (k = 0; k <= i - 1; k++) {
                h += q[i][k] * q[i][k];
            }
        }

        if (h + 1.0 == 1.0) {
            c[i] = 0.0;
            if (i == 1) {
                c[i] = q[i][i - 1];
            }
            b[i] = 0.0;
        } else {
            c[i] = sqrt(h);
            if (q[i][i - 1] > 0.0) {
                c[i] = -c[i];
            }
            h = h - q[i][i - 1] * c[i];
            q[i][i - 1] -= c[i];
            f = 0.0;
            for (j = 0; j <= i - 1; j++) {
                q[j][i] = q[i][j] / h;
                g = 0.0;
                for (k = 0; k <= j; k++) {
                    g = g + q[j][k] * q[i][k];
                }
                if (j + 1 <= i - 1) {
                    for (k = j + 1; k <= i - 1; k++) {
                        g = g + q[k][j] * q[i][k];
                    }
                }
                c[j] = g / h;
                f = f + g * q[j][i];
            }
            h2 = f / (h + h);
            for (j = 0; j <= i - 1; j++) {
                f = q[i][j];
                g = c[j] - h2 * f;
                c[j] = g;
                for (k = 0; k <= j; k++) { ;
                    q[j][k] = q[j][k] - f * c[k] - g * q[i][k];
                }
            }
            b[i] = h;
        }
    }
    for (i=0; i<=n-2; i++){
        c[i]=c[i+1];
    }
    c[n-1]=0.0;
    b[0]=0.0;
    for (i=0; i<=n-1; i++) {
        if ((b[i]!=0.0)&&(i-1>=0)){
            for (j=0; j<=i-1; j++) {
                g=0.0;

                for (k=0; k<=i-1; k++) {
                    g = g + q[i][ k] * q[k][j];
                }
                for (k=0; k<=i-1; k++){
                    q[k][j]=q[k][j]-g*q[k][i];
                }
            }
        }

        b[i]=q[i][i];
        q[i][i]=1.0;
        if (i-1>=0){
            for (j=0; j<=i-1; j++) {
                q[i][j]=0.0;
                q[j][i]=0.0;
            }
        }
    }
}

int sstq(int & n, double * b, double *c, double **q, double eps, int l) {
    int i,j,k,m,it;
    double d,f,h,g,p,r,e,s;
    c[n-1]=0.0;
    d=0.0;
    f=0.0;
    for (j=0; j<n; ++j) {
        it=0;
        h=eps*(fabs(b[j])+fabs(c[j]));
        if (h>d){
            d=h;
        }
        m=j;
        while ((m<n)&&(fabs(c[m])>d)){
            ++m;
        }
        if (m!=j){
            do {
                if (it==l) {
                    printf("fail\n");
                    return(-1);
                }
                it=it+1;
                g=b[j];
                p=(b[j+1]-g)/(2.0*c[j]);
                r=sqrt(p*p+1.0);
                if (p>=0.0){
                    b[j]=c[j]/(p+r);
                } else{
                    b[j]=c[j]/(p-r);
                }
                h=g-b[j];
                for (i=j+1; i<n; ++i){
                    b[i]=b[i]-h;
                }
                f+=h;
                p=b[m];
                e=1.0;
                s=0.0;
                for (i=m-1; i>=j; --i) {
                    g=e*c[i];
                    h=e*p;
                    if (fabs(p)>=fabs(c[i])) {
                        e=c[i]/p;
                        r=sqrt(e*e+1.0);
                        c[i+1]=s*p*r;
                        s=e/r;
                        e=1.0/r;
                    } else {
                        e=p/c[i];
                        r=sqrt(e*e+1.0);
                        c[i+1]=s*c[i]*r;
                        s=1.0/r;
                        e=e/r;
                    }
                    p=e*b[i]-s*g;
                    b[i+1]=h+s*(e*g+s*b[i]);
                    for (k=0; k<n; ++k) {
                        h=q[k][i+1];
                        q[k][i+1]=s*q[k][i]+e*h;
                        q[k][i]=e*q[k][i]-s*h;
                    }
                }
                c[j]=s*p;
                b[j]=e*p;
            } while (fabs(c[j])>d);
        }
        b[j]=b[j]+f;
    }
//    for (i=0; i<n; i++) {
//        k = i;
//        p = b[i];
//        if ((i + 1) < n) {
//            j = i + 1;
//            while ((j < n) && (b[j] <= p)) {
//                k = j;
//                p = b[j];
//                ++j;
//            }
//        }
//        if (k != i) {
//            b[k] = b[i];
//            b[i] = p;
//            for (j = 0; j < n; j++) {
//                p = q[j][i];
//                q[j][i] = q[j][k];
//                q[j][k] = p;
//            }
//        }
//    }
    return(1);
}

Eigen_result eigen( const My_matrix<double> & a, const double & eps, const int & l, const int & keep) {
    int t = a.get_num_column();
    int toDelete = t - keep;
    int i, j;
    double *b = new double[t];
    double *c = new double[t];
    double ** q = new double *[t];
    for( i=0; i<t; ++i ) {
        q[i] = new double[t];
    }

    strq(a,q,b,c);

    sstq(t,b,c,q,eps,l);

    My_Vector<double> eigen_values(t);
    for( i=0; i<t; ++i ) {
        eigen_values.get_array()[i] = fabs(b[i]);
    }

    std::sort(eigen_values.get_array(), eigen_values.get_array() + eigen_values.get_length());
    std::map<int, int> order;
    std::map<int, int> rev_order;
    std::set<int> includedIndexs;
    for( i=0; i<a.get_num_row(); ++i ){
        for( j=0; j<a.get_num_row(); ++j ){
            if( (eigen_values.get_array()[j] == fabs(b[i])) && includedIndexs.find(j)==includedIndexs.end() ){
                rev_order[j] = i;
                order[i] = j;
                includedIndexs.insert(j);
            }
        }
    }
    delete [] c;
    if( keep >0 && keep< t){
        My_Vector<double> eigen_values1(keep);
        My_matrix<double> v(keep, t);
        for( i=0; i<keep; ++i ){
            eigen_values1.get_array()[i] = b[rev_order[toDelete+i]];
            for( j=0; j<a.get_num_row(); ++j ) {
                v.get_matrix()[i][j] = q[j][rev_order[toDelete+i]];
            }
        }
        Eigen_result eigen_result(eigen_values1, v);

        for( i=0; i<t; ++i ) {
            delete [] q[i];
        }
        delete [] b;
        return eigen_result;
    }else{
        My_matrix<double> v(t, t);
        for( i=0; i<a.get_num_row(); ++i ){
            for( j=0; j<a.get_num_row(); ++j ) {
                v.get_matrix()[order[j]][i] = q[i][j]; // change order and transform
            }
        }
        My_Vector<double> eigen_values1(t);
        for( i=0; i<t; ++i ) {
            eigen_values1.get_array()[i] = b[rev_order[i]];
        }
        Eigen_result eigen_result(eigen_values1, v);

        for( i=0; i<t; ++i ) {
            delete [] q[i];
        }
        delete [] b;
        return eigen_result;
    }
}

Eigen_result eigen( const My_matrix<double> & _a, const int & keep){
    double eps = 0.000001;
    return eigen( _a, eps, 6000, keep);
}

Eigen_result eigen( const My_matrix<double> & _a, const double & eps){
    return eigen( _a, eps, 6000, _a.get_num_column());
}

Eigen_result eigen( const My_matrix<double> & _a) {
    double eps = 0.00000001;
    return eigen(_a, eps, 60000, _a.get_num_column());
}
