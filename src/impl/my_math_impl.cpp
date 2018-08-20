//
// Created by baoxing on 2/18/18.
//

#include "my_math_impl.h"
//
//double my_gamma( const double & x ){
//    int i;
//    double y, t, s, u;
//    static double a[11] = {0.0000677106, -0.0003442342,
//    0.0015397681, -0.0024467480, 0.0109736958,
//    -0.0002109075, 0.0742379071, 0.0815782188,
//    0.4118402518, 0.4227843370, 1.0};
//    if( x <= 0.0 ){
//        std::cout << "error, x<=0 for gamma test" << std::endl;
//        return -1.0;
//    }
//    y=x;
//    if( y<=1.0 ){
//        t = 1.0/(y*(y+1.0));
//        y += 2.0;
//    }else if(y<=2.0){
//        t=1.0/y;
//        y+=1.0;
//    }else if( y<= 3.0){
//        t=1.0;
//    }else{
//        t=1.0;
//        while( y > 3.0 ){
//            y -= 1.0;
//            t=t*y;
//        }
//    }
//    s = a[0];
//    u = y - 2.0;
//    for( i=0; i<11; ++i ){
//        s=s*u+a[i];
//    }
//    s=s*t;
//    return s;
//}
//
//static double bt(const double & a, const double & b, const double & x){
//    int k;
//    double d, p0=0.0, q0=1.0, p1=1.0, q1=1.0, s0, s1;
//    for( k=1; k <=1000000; ++k ){
//        d=(a+k)*(a+b+k)*x;
//        d=-d/((a+k+k)*(a+k+k+1.0));
//        p0=p1+d*p0;
//        q0=q1+d*q0;
//        s0=p0/q0;
//        d=k*(b-k)*x;
//        d=d/((a+k+k-1.0)*(a+k+k));
//        p1=p0+d*p1;
//        q1=q0+d*q1;
//        s1=p1/q1;
//        if( fabs(s1-s0) < fabs(s1)*1.0e-15 ){
//            return s1;
//        }
//    }
//    printf("a or b too big!");
//    return -1.0;
//}
//
//double beta( const double & a, const double & b, const double & x ){
//    double y;
//    double bt( const double &, const double &, const double &);
//    if( a<=0.0 ){
//        std::cout << "error, a<=0 for beta test" << std::endl;
//        return -1.0;
//    }
//    if( b<=0.0 ){
//        std::cout << "error, b<=0 for beta test" << std::endl;
//        return -1.0;
//    }
//    if( (x<0.0) || (x>1.0) ){
//        std::cout << "error, x<0 or x>1 for beta test" << std::endl;
//        return 1.0e+07;
//    }
//    if( (x==0.0) || (x==1.0) ){
//        y=0.0;
//    }else{
//        y=a*log(x)+b*log(1.0-x);
//        y=exp(y);
//        y=y*my_gamma(a+b)/(my_gamma(a)*my_gamma(b));
//    }
//    if( x < (a+1.0)/(a+b+2.0) ){
//        y = y * bt(a,b,x)/a;
//    }else{
//        y=1.0 - y * bt(b, a, 1.0-x)/b;
//    }
//    return y;
//}
//
//double sf(const double& f, const int& n1, const int &n2){
//    double y, f_t=f;
//    if( f_t<0.0 ){
//        f_t=-f_t;
//    }
//    y=beta(n2/2.0, n1/2.0, n2/(n2+n1*f_t));
//    return y;
//}


double c[11] = { 0.0000677106, -0.0003442342, 0.0015397681, -0.0024467480,
                 0.0109736958, -0.0002109075, 0.0742379071, 0.0815782188, 0.4118402518,
                 0.4227843370, 1.0000000000 };

double my_gamma(const double & xx) {
    double x, y, tmp, ser;
    static const double cof[6] = {
            76.18009172947146,
            -86.50532032941677,
            24.0140982408091,
            -1.231739572460155,
            0.1208650973866179e-2,
            -0.5395239384953e-5
    };
    y = x = xx;
    tmp = (x + 0.5) * log(x + 5.5) - (x + 5.5);
    ser = 1.000000000190015;
    for (int j = 0; j < 6; ++j) {
        ser += cof[j] / (y + 1);
        ++y;
    }
    return tmp + log(2.5066282746310005 * ser / x);
}

double beta(const double & x, const double & y) {
    if (x <= 0 || y <= 0) {
        return 0;
    }
    return exp(my_gamma(x)+my_gamma(y)-my_gamma(x + y));
}


double fi(const int & N, const double & x, const double & a, const double & b) {
    int n = N / 2;
    double f = 0.0, f1, s1, s2, tmpU, tmpV;
    int i;
    for (i = n; i >= 1; --i) {
        tmpU = (a + 2.0 * i - 1.0) * (a + 2.0 * i);
        s2 = i * (b - i) * x / tmpU;
        f1 = s2 / (1.0 + f);
        tmpV = (a + 2.0 * i - 2.0) * (a + 2.0 * i - 1.0);
        s1 = -(a + i - 1.0) * (b + a + i - 1.0) * x / tmpV;
        f = s1 / (1.0 + f1);
    }
    return 1.0 / (1.0 + f);
}

double incomBeta(const double & x, const double & a, const double & b) {
    double precise = 1.0e-30;// if change the value to 1.0e-20 the program does not obviously faster
    if (a <= 0.0 || b <= 0.0) {
        return 0.0;
    }
    if (fabs(x - 0.0) < precise || fabs(x - 1.0) < precise) {
        return 0.0;
    }

    double c1, c2, c3, f1, f2;
    int n;
    c1 = pow(x, a);
    c2 = pow(1.0 - x, b);
    c3 = beta(a, b);
    if (x < (a + 1.0) / (a + b + 2.0)) {
        n = 1;
        while (1) {
            f1 = fi(2 * n, x, a, b);
            f2 = fi(2 * n + 2, x, a, b);
            if (fabs(f2 - f1) < precise)
                return f2 * c1 * c2 / a / c3;
            else
                n++;
        }
    } else {
        if (fabs(x - 0.5) < precise && fabs(a - b) < precise)
            return 0.5;
        else {
            n = 1;
            while (1) {
                f1 = fi(2 * n, 1.0 - x, b, a);
                f2 = fi(2 * n + 2, 1.0 - x, b, a);
                if (fabs(f2 - f1) < precise)
                    return 1.0 - f2 * c1 * c2 / b / c3;
                else
                    n++;
            }
        }
    }
    return 0;
}

double sf(double & f, const int & n1, const int & n2) {
    if (f < 0.0)
        f = -f;
    return incomBeta(n2 / (n2 + n1 * f), n2 / 2.0, n1 / 2.0);
}
//---get P value from F value and df1 df2---------------------------end

double gam1(const double & x){
    int i;
    double y,t,s,u;
    static double a[11]={ 0.0000677106,-0.0003442342,
           0.0015397681,-0.0024467480,0.0109736958,
           -0.0002109075,0.0742379071,0.0815782188,
           0.4118402518,0.4227843370,1.0};
    if (x<=0.0){
        printf("err**x<=0!\n"); return(-1.0);
    }
    y=x;
    if (y<=1.0){
        t=1.0/(y*(y+1.0));
        y=y+2.0;
    }else if (y<=2.0){
        t=1.0/y;
        y=y+1.0;
    }else if (y<=3.0){
        t=1.0;
    }else{
        t=1.0;
        while (y>3.0){
            y=y-1.0;
            t=t*y;
        }
    }
    s=a[0];
    u=y-2.0;
    for (i=1; i<=10; i++){
        s=s*u+a[i];
    }
    s=s*t;
    return(s);
  }

double gam2(const double & a, const double & x){
    int n;
    double p,q,d,s,s1,p0,q0,p1,q1,qq;
    if ((a<=0.0)||(x<0.0)){
        if (a<=0.0) printf("err**a<=0!\n");
          if (x<0.0) printf("err**x<0!\n");
            return(-1.0);
    }
    if (x+1.0==1.0) return(0.0);
    if (x>1.0e+35) return(1.0);
    q=log(x); q=a*q; qq=exp(q);
    if (x<1.0+a){
        p=a;
        d=1.0/a;
        s=d;
        for (n=1; n<=100; n++){
            p=1.0+p;
            d=d*x/p;
            s=s+d;
            if (fabs(d)<fabs(s)*1.0e-07){
                s=s*exp(-x)*qq/gam1(a);
                return(s);
            }
        }
    }else{
        s=1.0/x;
        p0=0.0;
        p1=1.0;
        q0=1.0;
        q1=x;
        for (n=1; n<=100; n++){
            p0=p1+(n-a)*p0;
            q0=q1+(n-a)*q0;
            p=x*p0+n*p1;
            q=x*q0+n*q1;
            if (fabs(q)+1.0!=1.0){
                s1=p/q;
                p1=p;
                q1=q;
                if (fabs((s1-s)/s1)<1.0e-07){
                    s=s1*exp(-x)*qq/gam1(a);
                    return(1.0-s);
                }
                s=s1;
            }
            p1=p;
            q1=q;
        }
    }
    //printf("a too large !\n");
    s=1.0-s*exp(-x)*qq/gam1(a);
    return(s);
}

double chii(double x,const int &n){
    double y;
    if (x<0.0){
        x=-x;
    }
    y=gam2(n/2.0,x/2.0);
    return(y);
}

double critchi (const double & p, const double & df) {
    double  minchisq = 0.0;
    double  maxchisq = 99999.0;
    double  chisqval;

    if (p <= 0.0)
        return (maxchisq);
    else if (p >= 1.0)
        return (0.0);

    chisqval = df / sqrt (p);    /* fair first value */
    while (maxchisq - minchisq > 10e-7) {
        if ( (1-chii (chisqval, df)) < p) {
            maxchisq = chisqval;
        } else {
            minchisq = chisqval;
        }
        chisqval = (maxchisq + minchisq) * 0.5;
    }
    return chisqval;
}


//
//double chi(double x, int n) {
//    if (x < 0.0)
//        x = -x;
//    return incomBeta(n2 / (n2 + n1 * f), n2 / 2.0, n1 / 2.0);
//}
