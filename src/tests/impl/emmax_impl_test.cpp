//
// Created by baoxing on 2/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/impl.h"

TEST(get_eigen_L_test, c1){
    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    Kinship_matrix_impl k_i(kinship_file);
    std::cout << "line 25 " << currentDateTime() << std::endl;
    Eigen_result e_r = _get_eigen_L_(k_i.getKin());
    std::cout << "line 27 " << currentDateTime() << std::endl;
    for( int i=0; i< e_r.get_eigen_values().get_length(); ++i ){
        printf("%10.20f\n" , e_r.get_eigen_values().get_array()[i]);
    }
    ASSERT_EQ(0, 0);
}

TEST(get_eigen_R_test_function, c0){

    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    Kinship_matrix_impl k_i(kinship_file);

    My_matrix<double> x(k_i.getKin().get_num_row(), 1);
    int i, j;
    for( i=0; i< k_i.getKin().get_num_row(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    My_matrix<double> xs(k_i.getKin().get_num_row(), 0);
    My_matrix<double> full_x(x.get_num_row(), x.get_num_column()+xs.get_num_column());
    append_two_matrix( xs, x, full_x);
    Eigen_result e_r = _get_eigen_R_(full_x, k_i.getKin());
    for( i=0; i< e_r.get_eigen_values().get_length(); ++i ){
        printf("%10.20f\n" , e_r.get_eigen_values().get_array()[i]);
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for( i=0; i< e_r.get_eigen_vectors().get_num_column(); ++i ) {
        printf("%10.20f\n", e_r.get_eigen_vectors().get_matrix()[1][i]);
    }
    std::cout << e_r.get_eigen_values().get_array()[0] << std::endl;
    std::cout << e_r.get_eigen_values().get_array()[1] << std::endl;
    ASSERT_EQ(0, 0);
}

TEST(get_eigen_R_test_function, c1){
    std::string format = "tfam";
    std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype =  Genotype(number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);

    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    Kinship_matrix_impl k_i(kinship_file);


    My_matrix<double> x(genotype.get_number_of_individual(), 1);
    int i, j;
    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    My_matrix<double> xs(genotype.get_number_of_individual(), 0);
    My_matrix<double> full_x(x.get_num_row(), x.get_num_column()+xs.get_num_column());
    append_two_matrix( xs, x, full_x);
    Eigen_result e_r = _get_eigen_R_(full_x, k_i.getKin());
    for( int i=0; i< e_r.get_eigen_values().get_length(); ++i ){
        printf("%10.20f\n" , e_r.get_eigen_values().get_array()[i]);
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for( int i=0; i< e_r.get_eigen_vectors().get_num_column(); ++i ) {
        printf("%10.20f\n", e_r.get_eigen_vectors().get_matrix()[1][i]);
    }
    std::cout << e_r.get_eigen_values().get_array()[0] << std::endl;
    std::cout << e_r.get_eigen_values().get_array()[1] << std::endl;
    ASSERT_EQ(0, 0);
}

TEST(emma_estimates, c1){
    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype_2";
    std::string format = "tfam";
    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    double maf = 0.1;
    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;

    My_matrix<double> x(pi.getPhenotypes().get_length(), 1);
    int i;
    for( i=0; i<pi.getPhenotypes().get_length(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    int ngrids = 100;
    double llim = -10.0;
    double ulim = 10.0;
    double eps = 1e-10;
    std::string method="ML";
    int maxiter = 100;
    Eigen_result eigen_L = _get_eigen_L_( k_i.getKin());

    int j;
    int n = pi.getPhenotypes().get_length();
    int q = 1;
    int p = n - q;
    int m = ngrids+1;

    My_matrix<double> full_x(x.get_num_row(), 1);
    full_x.value_copy(x);
    std::cout << " line 243 " << std::endl;
    Eigen_result eigen_R = _get_eigen_R_(full_x, k_i.getKin());
    std::cout << " line 245 " << std::endl;
    My_Vector<double> etas(p);
    double t;
    for( i=0; i < p; ++i ){
        t = 0.0;
        std::cout << " eigen_R: " << i << " " << eigen_R.get_eigen_values().get_array()[i] << " " << std::endl;
        for( j=0; j < n; ++j ){
            t = t + eigen_R.get_eigen_vectors().get_matrix()[i][j]*pi.getPhenotypes().get_array()[j];
        }
        etas.get_array()[i]=t;
        std::cout << " etas i " << i << " " << t << std::endl;
    }

    My_Vector<double> sq_etas(p);
    for( i=0; i < p; ++i ){
        sq_etas.get_array()[i] = etas.get_array()[i] * etas.get_array()[i];
        std::cout << i << " sq_etas.get_array()[i] " << sq_etas.get_array()[i] << std::endl;
    }

}
TEST(emma_estimates, c2){
    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype_2";
    std::string format = "tfam";
    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    double maf = 0.1;
    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;

    My_matrix<double> x(pi.getPhenotypes().get_length(), 1);
    int i;
    for( i=0; i<pi.getPhenotypes().get_length(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    int ngrids = 100;
    double llim = -10.0;
    double ulim = 10.0;
    double eps = 1e-10;
    std::string method="ML";
    int maxiter = 100;
    Eigen_result eigen_L = _get_eigen_L_( k_i.getKin());

    int j;
    int n = pi.getPhenotypes().get_length();
    int q = 1;
    int p = n - q;
    int m = ngrids+1;

    My_matrix<double> full_x(x.get_num_row(), 1);
    full_x.value_copy(x);
    std::cout << " line 243 " << std::endl;
    Eigen_result eigen_R = _get_eigen_R_(full_x, k_i.getKin());
    std::cout << " line 245 " << std::endl;
    My_Vector<double> etas(p);
    double t;
    for( i=0; i < p; ++i ){
        t = 0.0;
        std::cout << " eigen_R: " << i << " " << eigen_R.get_eigen_values().get_array()[i] << " " << std::endl;
        for( j=0; j < n; ++j ){
            t = t + eigen_R.get_eigen_vectors().get_matrix()[i][j]*pi.getPhenotypes().get_array()[j];
        }
        etas.get_array()[i]=t;
        std::cout << " etas i " << i << " " << t << std::endl;
    }

    My_Vector<double> sq_etas(p);
    for( i=0; i < p; ++i ){
        sq_etas.get_array()[i] = etas.get_array()[i] * etas.get_array()[i];
        std::cout << i << " sq_etas.get_array()[i] " << sq_etas.get_array()[i] << std::endl;
    }

    My_Vector<double> log_deltas(m); // the space for deltas to search
    My_Vector<double> deltas(m);
    for ( i=0; i < m; ++i ){
        log_deltas.get_array()[i] = (double(i) / ngrids)*(ulim - llim) + llim;
        deltas.get_array()[i] = exp(log_deltas.get_array()[i]);
    }

    for( i=0; i<m; ++i ) {
        std::cout << " _ll_ i " << i << " deltas.get_array()[i] " << deltas.get_array()[i] << " " <<  _ll_(deltas.get_array()[i], eigen_R, eigen_L, sq_etas) << std::endl;
    }

    for( i=0; i<m; ++i ){
        std::cout << " _dll_ i " << i << " deltas.get_array()[i] " << deltas.get_array()[i] << " " <<  _dll_(deltas.get_array()[i], eigen_R, eigen_L, sq_etas) << std::endl;
    }
}
