//
// Created by baoxing on 2/24/18.
//

//
// Created by baoxing on 2/21/18.
//

#include "../../model/model.h"
#include "../../impl/impl.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

int add (int a, int b){
    return a+b;
}

TEST(sum_of_a_vector, c1){
    My_Vector<double> a(10);
    for (int i=0; i<a.get_length(); i++){
        a.get_array()[i]=i;
    }
    std::cout << sum_of_a_vector(a) << std::endl;
    ASSERT_EQ(45, sum_of_a_vector(a));
}


TEST(sum_of_powed_vector, c1){
    My_Vector<double> a(10);
//    EXPECT_EQ(0, sum_of_a_vector(a));
    for (int i=0; i<a.get_length(); i++){
        a.get_array()[i]=i;
    }
    int power = 5;
    ASSERT_EQ(120825, sum_of_powed_vector(a, 5));
}

TEST(sum_of_two_vectors, c1){
    My_Vector<double> a(10);
    My_Vector<double> b(10);
    for (int i=0; i<a.get_length(); i++){
        a.get_array()[i]=i;
        b.get_array()[9-i]=i;
    }
    My_Vector<double> sum(10);
    sum_of_two_vectors<double, double, double>(a, b, sum);
    for (int i=0; i<a.get_length(); i++) {
        ASSERT_EQ(9, sum.get_array()[i]);
    }
}

TEST(sum_of_a_vector_a_number, c1){
    My_Vector<double> a(10);
    for (int i=0; i<a.get_length(); i++){
        a.get_array()[i]=i;
    }
    My_Vector<double> result(10);
    for (int i=0; i<a.get_length(); i++) {
        sum_of_a_vector_a_number(a, 9-i, result);
        ASSERT_EQ(9, result.get_array()[i]);
    }
}

TEST(production_of_two_vectors, c1){
    My_Vector<double> a(10);
    My_Vector<double> b(10);
    for (int i=0; i<a.get_length(); i++){
        a.get_array()[i]=i;
        b.get_array()[9-i]=i;
    }
    My_Vector<double> c(10);
    production_of_two_vectors(a, b, c);
    for (int i=0; i<a.get_length(); i++) {
        ASSERT_EQ(i*(9-i), c.get_array()[i]);
    }
}

TEST(quotient_of_two_vectors, c1){
    My_Vector<double> a(10);
    My_Vector<double> b(10);
    for (int i=1; i<a.get_length()+1; i++){
        a.get_array()[i-1]=i;
        b.get_array()[i-1]=10*i;
    }
    My_Vector<double> c(10);
    quotient_of_two_vectors<double, double, double>(a, b, c);
    for (int i=0; i<a.get_length(); i++) {
        ASSERT_EQ(0.1, c.get_array()[i]);
    }
}


TEST (trmul, c1){
    My_matrix<double> a(4 ,5);
    a.get_matrix()[0][0] = 1.0;  a.get_matrix()[0][1] = 3.0;  a.get_matrix()[0][2] =  -2.0; a.get_matrix()[0][3] = 0.0; a.get_matrix()[0][4] = 4.0;
    a.get_matrix()[1][0] = -2.0; a.get_matrix()[1][1] = -1.0; a.get_matrix()[1][2] = 5.0;  a.get_matrix()[1][3] = -7.0; a.get_matrix()[1][4] = 2.0;
    a.get_matrix()[2][0] = 0.0;  a.get_matrix()[2][1] = 8.0;  a.get_matrix()[2][2] = 4.0;  a.get_matrix()[2][3] = 1.0; a.get_matrix()[2][4] = -5.0;
    a.get_matrix()[3][0] = 3.0;  a.get_matrix()[3][1] = -3.0; a.get_matrix()[3][2] = 2.0;  a.get_matrix()[3][3] = -4.0; a.get_matrix()[3][4] = 1.0;

    My_matrix<double> b(5 ,3);
    b.get_matrix()[0][0] = 4.0; b.get_matrix()[0][1]=5.0; b.get_matrix()[0][2]= -1.0;
    b.get_matrix()[1][0] = 2.0; b.get_matrix()[1][1]= -2.0; b.get_matrix()[1][2]= 6.0;
    b.get_matrix()[2][0] = 7.0; b.get_matrix()[2][1]= 8.0; b.get_matrix()[2][2]= 1.0;
    b.get_matrix()[3][0] = 0.0; b.get_matrix()[3][1]= 3.0; b.get_matrix()[3][2]= -5.0;
    b.get_matrix()[4][0] = 9.0; b.get_matrix()[4][1]= 8.0; b.get_matrix()[4][2]= -6.0;

    My_matrix<double> c(4 ,3);
    int m = 4;
    int n = 5;
    int k = 3;
    trmul(a, b, c);
    for( int i=0; i<4; ++i ){
        for( int j=0; j<3; ++j ){
            std::cout << c.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    ASSERT_EQ(c.get_matrix()[0][0], 32);
    ASSERT_EQ(c.get_matrix()[1][1], 27);
    ASSERT_EQ(c.get_matrix()[2][2], 77);
    ASSERT_EQ(0, 0);
}
TEST (inverse_matrix, c1){

    My_matrix<double> a(4,4);
    a.get_matrix()[0][0] = 0.2368; a.get_matrix()[0][1] = 0.2471; a.get_matrix()[0][2] = 0.2568; a.get_matrix()[0][3] = 1.2671;
    a.get_matrix()[1][0] = 1.1161; a.get_matrix()[1][1] = 0.1254; a.get_matrix()[1][2] = 0.1397; a.get_matrix()[1][3] = 0.1490;
    a.get_matrix()[2][0] = 0.1582; a.get_matrix()[2][1] = 1.1675; a.get_matrix()[2][2] = 0.1768; a.get_matrix()[2][3] = 0.1871;
    a.get_matrix()[3][0] = 0.1968; a.get_matrix()[3][1] = 0.2071; a.get_matrix()[3][2] = 1.2168; a.get_matrix()[3][3] = 0.2271;
    My_matrix<double> b(4,4);
    T_matrix(a,b);
    My_matrix<double> c(a);
    inverse_matrix(c);
    std::cout << "MAT A IS:" << std::endl;
    for ( int i=0; i<4; ++i ){
        for ( int j=0; j<4; ++j ) {
            std::cout << b.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    int i;
    std::cout << "MAT A- IS:" << std::endl;
    for ( i=0; i<4; ++i ){
        for ( int j=0; j<4; ++j ) {
            std::cout << a.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "MAT AA- IS:" << std::endl;
    trmul(b, a, c);
    for ( int i=0; i<4; ++i ){
        for ( int j=0; j<4; ++j ) {
            std::cout << c.get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    ASSERT_EQ(0, 0);
}

TEST (qr_decomposition, c1){
    My_matrix<double> a(4, 4);

    a.get_matrix()[0][0]=1.0;  a.get_matrix()[0][1]=1.0;  a.get_matrix()[0][2]=-1.0;
    a.get_matrix()[1][0]=2.0;  a.get_matrix()[1][1]=1.0;  a.get_matrix()[1][2]=0.0;
    a.get_matrix()[2][0]=1.0;  a.get_matrix()[2][1]=-1.0; a.get_matrix()[2][2]=0.0;
    a.get_matrix()[3][0]=-1.0; a.get_matrix()[3][1]=2.0;  a.get_matrix()[3][2]=1.0;
//    My_matrix _a(a);
    Qr_decomposition_result qr_decomposition_result = qr_decomposition(a);

    std::cout << "MAT Q IS:" << std::endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << qr_decomposition_result.get_q().get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "MAT R IS:" << std::endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << qr_decomposition_result.get_r().get_matrix()[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    ASSERT_TRUE(fabs(qr_decomposition_result.get_q().get_matrix()[0][0] - -0.377964) < 0.00001 );
    ASSERT_TRUE(fabs(qr_decomposition_result.get_q().get_matrix()[1][1] - -0.377964) < 0.00001 );
    ASSERT_TRUE(fabs(qr_decomposition_result.get_q().get_matrix()[2][2] - -0.377964) < 0.00001 );
    ASSERT_TRUE(fabs(qr_decomposition_result.get_q().get_matrix()[3][0] - 0.377964) < 0.00001 );

    ASSERT_TRUE(fabs(qr_decomposition_result.get_r().get_matrix()[0][0] - -2.64575) < 0.00001 );
    ASSERT_TRUE(fabs(qr_decomposition_result.get_r().get_matrix()[1][1] - -2.64575) < 0.00001 );
    ASSERT_TRUE(fabs(qr_decomposition_result.get_r().get_matrix()[2][2] - -1.13389) < 0.00001 );

    ASSERT_EQ(0, 0);
}

TEST (eigen, c1){
    My_matrix<double> a(5, 5);

    a.get_matrix()[0][0]=10.0; a.get_matrix()[0][1]=1.0;  a.get_matrix()[0][2]=2.0;  a.get_matrix()[0][3]=3.0;  a.get_matrix()[0][4]=4.0;
    a.get_matrix()[1][0]=1.0;  a.get_matrix()[1][1]=9.0;  a.get_matrix()[1][2]=-1.0; a.get_matrix()[1][3]=2.0;  a.get_matrix()[1][4]=-3.0;
    a.get_matrix()[2][0]=2.0;  a.get_matrix()[2][1]=-1.0; a.get_matrix()[2][2]=7.0;  a.get_matrix()[2][3]=3.0;  a.get_matrix()[2][4]=-5.0;
    a.get_matrix()[3][0]=3.0;  a.get_matrix()[3][1]=2.0;  a.get_matrix()[3][2]=3.0;  a.get_matrix()[3][3]=12.0; a.get_matrix()[3][4]=-1.0;
    a.get_matrix()[4][0]=4.0;  a.get_matrix()[4][1]=-3.0; a.get_matrix()[4][2]=-5.0; a.get_matrix()[4][3]=-1.0; a.get_matrix()[4][4]=15.0;
    My_matrix<double> c(a);
    for( int i=0; i<5; ++i ){
        c.get_matrix()[i][i]+=1;
    }
//    double eps=0.0000000000000000000001;
    Eigen_result r1 = eigen(a);
    Eigen_result r2 = eigen(c);
    std::cout << std::endl;
    int i, j;
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            printf("%12.25f ", r1.get_eigen_vectors().get_matrix()[i][j]);
            printf("%12.25f ", r2.get_eigen_vectors().get_matrix()[i][j]);
            std::cout << std::endl;
            ASSERT_TRUE(fabs(r1.get_eigen_vectors().get_matrix()[i][j]-r2.get_eigen_vectors().get_matrix()[i][j]) < 0.0000000001);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for( i=0; i<5; ++i ){
        printf("%12.25f\n", r1.get_eigen_values().get_array()[i]);
        ASSERT_TRUE(fabs(r1.get_eigen_values().get_array()[i] -r2.get_eigen_values().get_array()[i]+1) < 0.0000000001 );
    }

    std::cout << std::endl;

    My_matrix<double> this_eig_values(5, 5);
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            this_eig_values.get_matrix()[i][j]=0;
        }
        this_eig_values.get_matrix()[i][i]=r1.get_eigen_values().get_array()[i];
    }
    My_matrix<double> this_eig_vectors(5, 5);
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            this_eig_vectors.get_matrix()[i][j] = r1.get_eigen_vectors().get_matrix()[j][i];
        }
    }

    My_matrix<double> t1(5, 5);
    trmul(this_eig_vectors, this_eig_values, t1);
    My_matrix<double> t2(5, 5);
    trmul(t1, r1.get_eigen_vectors(), t2);
    for( i=0; i<5; ++i ) {
        for (j = 0; j < 5; ++j) {
            ASSERT_TRUE(fabs(t2.get_matrix()[i][j] - a.get_matrix()[i][j]) < 0.0000000001 );
        }
    }
}

TEST (eigen, cqr){
    My_matrix<double> a(5, 5);

    a.get_matrix()[0][0]=10.0; a.get_matrix()[0][1]=1.0;  a.get_matrix()[0][2]=2.0;  a.get_matrix()[0][3]=3.0;  a.get_matrix()[0][4]=4.0;
    a.get_matrix()[1][0]=1.0;  a.get_matrix()[1][1]=9.0;  a.get_matrix()[1][2]=-1.0; a.get_matrix()[1][3]=2.0;  a.get_matrix()[1][4]=-3.0;
    a.get_matrix()[2][0]=2.0;  a.get_matrix()[2][1]=-1.0; a.get_matrix()[2][2]=7.0;  a.get_matrix()[2][3]=3.0;  a.get_matrix()[2][4]=-5.0;
    a.get_matrix()[3][0]=3.0;  a.get_matrix()[3][1]=2.0;  a.get_matrix()[3][2]=3.0;  a.get_matrix()[3][3]=12.0; a.get_matrix()[3][4]=-1.0;
    a.get_matrix()[4][0]=4.0;  a.get_matrix()[4][1]=-3.0; a.get_matrix()[4][2]=-5.0; a.get_matrix()[4][3]=-1.0; a.get_matrix()[4][4]=15.0;
    My_matrix<double> c(a);
    for( int i=0; i<5; ++i ){
        c.get_matrix()[i][i]+=1;
    }
    double eps=0.0000000000001;
    Eigen_result r1 = eigen(a, eps, 1000000, a.get_num_column());
    Eigen_result r2 = eigen(c, eps, 1000000, c.get_num_column());
    std::cout << std::endl;
    int i, j;
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            printf("%12.25f ", r1.get_eigen_vectors().get_matrix()[i][j]);
            printf("%12.25f ", r2.get_eigen_vectors().get_matrix()[i][j]);
            std::cout << std::endl;
            ASSERT_TRUE(fabs(r1.get_eigen_vectors().get_matrix()[i][j]-r2.get_eigen_vectors().get_matrix()[i][j]) < 0.0000000001);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for( i=0; i<5; ++i ){
        printf("%12.25f\n", r1.get_eigen_values().get_array()[i]);
        ASSERT_TRUE(fabs(r1.get_eigen_values().get_array()[i] -r2.get_eigen_values().get_array()[i]+1) < 0.0000000001 );
    }

    std::cout << std::endl;
    My_matrix<double> this_eig_values(5, 5);
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            this_eig_values.get_matrix()[i][j]=0;
        }
        this_eig_values.get_matrix()[i][i]=r1.get_eigen_values().get_array()[i];
    }
    My_matrix<double> this_eig_vectors(5, 5);
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            this_eig_vectors.get_matrix()[i][j] = r1.get_eigen_vectors().get_matrix()[j][i];
        }
    }

    My_matrix<double> t1(5, 5);
    trmul(this_eig_vectors, this_eig_values, t1);
    My_matrix<double> t2(5, 5);
    trmul(t1, r1.get_eigen_vectors(), t2);
    for( i=0; i<5; ++i ) {
        for (j = 0; j < 5; ++j) {
            ASSERT_TRUE(fabs(t2.get_matrix()[i][j] - a.get_matrix()[i][j]) < 0.0000000001 );
        }
    }
}

TEST (eigen, c2){
    Kinship_matrix_impl k_i("/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf");
    My_matrix<double> c(k_i.getKin());
    int t = k_i.getKin().get_num_column();

    for( int i=0; i<t; ++i ){
        c.get_matrix()[i][i]+=1;
    }
    double eps=0.0000000000001;

    Eigen_result r1 = eigen(k_i.getKin(), eps, 1000000, k_i.getKin().get_num_column());
    Eigen_result r2 = eigen(c, eps, 1000000, c.get_num_column());
    int i, j;

    for( i=0; i<t; ++i ){
        for( j=0; j<t; ++j ){
            printf("%12.25f ", r1.get_eigen_vectors().get_matrix()[i][j]);
            printf("%12.25f ", r2.get_eigen_vectors().get_matrix()[i][j]);
            std::cout << std::endl;
            ASSERT_TRUE((fabs(r1.get_eigen_vectors().get_matrix()[i][j])-fabs(r2.get_eigen_vectors().get_matrix()[i][j])) < 0.0000001);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for( i=0; i<t; ++i ){
        printf("%12.25f\n", r1.get_eigen_values().get_array()[i]);
        ASSERT_TRUE(fabs(r1.get_eigen_values().get_array()[i] -r2.get_eigen_values().get_array()[i]+1.0) < 0.0000001 );
    }

    for( i=0; i<t; ++i ){
        printf("%12.25f\n", r1.get_eigen_values().get_array()[i]);
    }
    printf("\n\n\n\n\n\n");
    for( i=0; i<t; ++i ){
        printf("%12.25f\n", r1.get_eigen_vectors().get_matrix()[0][i]);
    }
    std::cout << std::endl;
    My_matrix<double> this_eig_values(t, t);
    for( i=0; i<t; ++i ){
        for( j=0; j<t; ++j ){
            this_eig_values.get_matrix()[i][j]=0;
        }
        this_eig_values.get_matrix()[i][i]=r1.get_eigen_values().get_array()[i];
    }
    My_matrix<double> this_eig_vectors(t, t);
    T_matrix(r1.get_eigen_vectors(), this_eig_vectors);

    My_matrix<double> t1(t, t);
    trmul(this_eig_vectors, this_eig_values, t1);
    My_matrix<double> t2(t, t);
    trmul(t1, r1.get_eigen_vectors(), t2);
    for( i=0; i<t; ++i ) {
        for (j = 0; j < t; ++j) {
            std::cout <<  t2.get_matrix()[i][j] << " " << k_i.getKin().get_matrix()[i][j] << std::endl;
            ASSERT_TRUE(fabs(t2.get_matrix()[i][j] - k_i.getKin().get_matrix()[i][j]) < 0.00001 );
        }
    }
}

TEST(T_matrix, c1){
    int n=5, m = 6, i, j;
    My_matrix<double> a(n, m);

    for( i=0; i<n; ++i ){
        for( j=0; j<m; ++j ){
            a.get_matrix()[i][j] = i*j;
        }
    }
    My_matrix<double> b(m, n);
    T_matrix(a, b);
    for( i=0; i<n; ++i ){
        for( j=0; j<m; ++j ){
            ASSERT_EQ(a.get_matrix()[i][j], b.get_matrix()[j][i]);
        }
    }
    ASSERT_EQ(0, 0);
}

TEST(lsq, c1){
    My_matrix<double> x (4, 3);

    x.get_matrix()[0][0]=1.0; x.get_matrix()[0][1]=1.0; x.get_matrix()[0][2]=-1.0;
    x.get_matrix()[1][0]=2.0; x.get_matrix()[1][1]=1.0; x.get_matrix()[1][2]=0.0;
    x.get_matrix()[2][0]=1.0; x.get_matrix()[2][1]=-1.0; x.get_matrix()[2][2]=0.0;
    x.get_matrix()[3][0]=-1.0; x.get_matrix()[3][1]=2.0; x.get_matrix()[3][2]=1.0;

    My_Vector<double> y(4);
    y.get_array()[0] = 2.0;
    y.get_array()[1]=-3.0;
    y.get_array()[2]=1.0;
    y.get_array()[3]=4.0;

    My_Vector<double> b(x.get_num_column());
    lsq(x, y, b);
    for( int i=0; i<b.get_length(); ++i ){
        printf("%10.20f\n", b.get_array()[i]);
    }
    ASSERT_TRUE(fabs(-1.19047619047619068766 - b.get_array()[0]) < 0.0000000001 );
    ASSERT_TRUE(fabs(0.95238095238095243911 - b.get_array()[1]) < 0.0000000001 );
    ASSERT_TRUE(fabs(-0.66666666666666696273 - b.get_array()[2]) < 0.0000000001 );

    My_matrix<double> x1 (4, 2);
    x1.get_matrix()[0][0]=1.0;
    x1.get_matrix()[1][0]=1.0;
    x1.get_matrix()[2][0]=1.0;
    x1.get_matrix()[3][0]=1.0;

    x1.get_matrix()[0][1]=1.0;
    x1.get_matrix()[1][1]=2.0;
    x1.get_matrix()[2][1]=2.0;
    x1.get_matrix()[3][1]=1.0;
    My_Vector<double> b1(x1.get_num_column());
    lsq(x1, y, b1);

    My_matrix<double> x2 (4, 2);
    x2.get_matrix()[0][0]=1.0;
    x2.get_matrix()[1][0]=1.0;
    x2.get_matrix()[2][0]=1.0;
    x2.get_matrix()[3][0]=1.0;

    x2.get_matrix()[0][1]=1.0;
    x2.get_matrix()[1][1]=0.0;
    x2.get_matrix()[2][1]=0.0;
    x2.get_matrix()[3][1]=1.0;

    My_Vector<double> b2(x2.get_num_column());
    lsq(x2, y, b2);

    for( int i=0; i<x1.get_num_row(); ++i ){
        printf("intercept: %10.50f\t%10.50f\n", b1.get_array()[0], b2.get_array()[0]);
        printf("slop: %10.50f\t%10.50f\n", b1.get_array()[1], b2.get_array()[1]);
        printf("%10.50f\t%10.50f\n", b1.get_array()[0]+b1.get_array()[1]*x1.get_matrix()[i][1], b2.get_array()[0]+b2.get_array()[1]*x2.get_matrix()[i][1]);
    }
    // TODO there is no approximation, why they are not identical with each other
    ASSERT_EQ(b1.get_array()[0]+b1.get_array()[1]*x1.get_matrix()[0][1], b2.get_array()[0]+b2.get_array()[1]*x2.get_matrix()[0][1]);
    ASSERT_EQ(b1.get_array()[0]+b1.get_array()[1]*x1.get_matrix()[1][1], b2.get_array()[0]+b2.get_array()[1]*x2.get_matrix()[1][1]);
    ASSERT_EQ(b1.get_array()[0]+b1.get_array()[1]*x1.get_matrix()[2][1], b2.get_array()[0]+b2.get_array()[1]*x2.get_matrix()[2][1]);
    ASSERT_EQ(b1.get_array()[0]+b1.get_array()[1]*x1.get_matrix()[3][1], b2.get_array()[0]+b2.get_array()[1]*x2.get_matrix()[3][1]);
}

TEST(strq, c1){
    int i, j;
    double *b = new double[5];
    double *c = new double[5];
    double ** q = new double *[5];
    for( i=0; i<5; ++i ) {
        q[i] = new double[5];
    }
    static double a[5][5]={ {10.0,1.0,2.0,3.0,4.0},
                            {1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},
                            {3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0}};
    My_matrix<double> a_(5,5);
    for( i=0; i<5; ++i ){
        for( j=0; j<5; ++j ){
            a_.get_matrix()[i][j]=a[i][j];
        }
    }
    strq(a_,q,b,c);
    printf("MAT A IS:\n");
    for (i=0; i<=4; i++) {
        for (j=0; j<=4; j++){
            printf("%13.7e ",a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("MAT Q IS:\n");
    for (i=0; i<=4; i++) {
        for (j=0; j<=4; j++){
            printf("%13.7e ",q[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("MAT B IS:\n");
    for (i=0; i<=4; i++){
        printf("%13.7e ",b[i]);
    }
    printf("\n\n");
    printf("MAT C IS:\n");
    for (i=0; i<=4; i++){
        printf("%13.7e ",c[i]);
    }
    printf("\n\n");

    ASSERT_EQ(0, 0);
}

TEST(sstq, c1){
    int i, j,k,l=100000;;
    double *b = new double[5];
    double *c = new double[5];
    double ** q = new double *[5];
    for( i=0; i<5; ++i ) {
        q[i] = new double[5];
    }
    My_matrix<double> a(5, 5);

    a.get_matrix()[0][0]=10.0; a.get_matrix()[0][1]=1.0;  a.get_matrix()[0][2]=2.0;  a.get_matrix()[0][3]=3.0;  a.get_matrix()[0][4]=4.0;
    a.get_matrix()[1][0]=1.0;  a.get_matrix()[1][1]=9.0;  a.get_matrix()[1][2]=-1.0; a.get_matrix()[1][3]=2.0;  a.get_matrix()[1][4]=-3.0;
    a.get_matrix()[2][0]=2.0;  a.get_matrix()[2][1]=-1.0; a.get_matrix()[2][2]=7.0;  a.get_matrix()[2][3]=3.0;  a.get_matrix()[2][4]=-5.0;
    a.get_matrix()[3][0]=3.0;  a.get_matrix()[3][1]=2.0;  a.get_matrix()[3][2]=3.0;  a.get_matrix()[3][3]=12.0; a.get_matrix()[3][4]=-1.0;
    a.get_matrix()[4][0]=4.0;  a.get_matrix()[4][1]=-3.0; a.get_matrix()[4][2]=-5.0; a.get_matrix()[4][3]=-1.0; a.get_matrix()[4][4]=15.0;

    double eps=0.0000000000001;

    strq(a,q,b,c);
    int t=5;
    k=sstq(t,b,c,q,eps,l);
    printf("MAT A IS:\n");
    for (i=0; i<=4; i++)
    { for (j=0; j<=4; j++)
            printf("%13.7e ",a.get_matrix()[i][j]);
        printf("\n");
    }
    printf("\n");
    if (k>0)
    { printf("MAT B IS:\n");
        for (i=0; i<=4; i++)
            printf("%13.7e ",b[i]);
        printf("\n\n");
        printf("MAT Q IS:\n");
        for (i=0; i<=4; i++)
        { for (j=0; j<=4; j++)
                printf("%13.7e ",q[i][j]);
            printf("\n");
        }
        printf("\n");
    }
}
TEST(determinant, c1){
    My_matrix<double> a(4, 4);
    a.get_matrix()[0][0]=1.0; a.get_matrix()[0][1]=2.0; a.get_matrix()[0][2]=3.0; a.get_matrix()[0][3]=4.0;
    a.get_matrix()[1][0]=5.0; a.get_matrix()[1][1]=6.0; a.get_matrix()[1][2]=7.0; a.get_matrix()[1][3]=8.0;
    a.get_matrix()[2][0]=9.0; a.get_matrix()[2][1]=10.0; a.get_matrix()[2][2]=11.0; a.get_matrix()[2][3]=12.0;
    a.get_matrix()[3][0]=13.0; a.get_matrix()[3][1]=14.0; a.get_matrix()[3][2]=15.0; a.get_matrix()[3][3]=16.0;

    My_matrix<double> b(4, 4);
    b.get_matrix()[0][0]=3.0; b.get_matrix()[0][1]=-3.0; b.get_matrix()[0][2]=-2.0; b.get_matrix()[0][3]=4.0;
    b.get_matrix()[1][0]=5.0; b.get_matrix()[1][1]=-5.0; b.get_matrix()[1][2]=1.0; b.get_matrix()[1][3]=8.0;
    b.get_matrix()[2][0]=11.0; b.get_matrix()[2][1]=8.0; b.get_matrix()[2][2]=5.0; b.get_matrix()[2][3]=-7.0;
    b.get_matrix()[3][0]=5.0;b.get_matrix()[3][1]=-1.0; b.get_matrix()[3][2]=-3.0; b.get_matrix()[3][3]=-1.0;
    std::cout << determinant(a) << std::endl;
    std::cout << determinant(b) << std::endl;
    ASSERT_EQ(determinant(a), 0.0);
    ASSERT_TRUE(determinant(b) - 595.0 < 0.0000000001 );
}
