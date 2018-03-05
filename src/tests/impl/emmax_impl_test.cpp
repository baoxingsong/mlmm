//
// Created by baoxing on 2/25/18.
//

#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/impl.h"
/**
import numpy as np
from scipy import linalg
K = np.loadtxt("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf") # this is the kinship matrix
K = np.mat(K)
evals, evecs = linalg.eigh(K)
evals = np.array(evals, dtype="double")
result= {'values': evals, 'vectors': np.mat(evecs, dtype="double").T}

data <- read.table("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf")
data <- as.matrix(data)
re = eigen(data)
 */

TEST(get_eigen_L_test, c1){
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    Kinship_matrix_impl k_i(kinship_file);
    std::cout << "line 25 " << currentDateTime() << std::endl;
    Eigen_result e_r = _get_eigen_L_(k_i.getKin());
    std::cout << "line 27 " << currentDateTime() << std::endl;
    for( int i=0; i< e_r.get_eigen_values().get_length(); ++i ){
        printf("%10.20f\n" , e_r.get_eigen_values().get_array()[i]);
    }
    ASSERT_EQ(0, 0);
}


/**
import numpy as np
from scipy import linalg
K = np.loadtxt("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf") # this is the kinship matrix
K = np.mat(K)

def _get_eigen_R_( X=None, K=None, dtype='double'):
	q = X.shape[1]
    X_squared_inverse = linalg.pinv(X.T * X)  # (X.T*X)^{-1}
	# linalg.pinv: Calculate the generalized inverse of a matrix using
	# its singular-value decomposition (SVD) and including all large singular values.
	hat_matrix = X * X_squared_inverse * X.T
	S = np.mat(np.identity(a.shape[0])) - hat_matrix  # S=I-X(X'X)^{-1}X'
	M = np.mat(S * (K + np.matrix(np.identity(a.shape[0]))) * S, dtype='double')
	evals, evecs = linalg.eigh(M)  # eigen of S(K+I)S
	eig_values = np.array(evals[q:], dtype=dtype) - 1  # Don't know what minute one, the emma R code did this
	return {'values': eig_values, 'vectors': (np.mat(evecs, dtype=dtype).T[q:])}


a = np.matrix([[2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 2, 1, 1],[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]).T

e_r = _get_eigen_R_(a, K)
e_r['values']

*/

TEST(get_eigen_R_test_function, c1){
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    Kinship_matrix_impl k_i(kinship_file);

    Genotype genotype= Read_tped_file("/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf");

    My_matrix<double> x(genotype.get_number_of_individual(), 1);
    int i, j;
    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    My_matrix<double> xs(genotype.get_number_of_individual(), 1);
    j=0;
    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        xs.get_matrix()[i][0]=genotype.get_genotype_matrix().get_matrix()[i][j];
        //std::cout << xs.get_matrix()[i][0] << std::endl;
    }
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
    ASSERT_TRUE(fabs(e_r.get_eigen_values().get_array()[0] - 5.0497597) < 0.00001 );
    ASSERT_TRUE(fabs(e_r.get_eigen_values().get_array()[1] - 3.68871906) < 0.00001 );
    //ASSERT_TRUE(fabs(e_r.get_eigen_vectors().get_matrix()[1][1] - -0.116992663) < 0.00001 );
    ASSERT_EQ(0, 0);
}

//
//TEST(EigenSolver, c1){
//    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
//    Kinship_matrix_impl k_i(kinship_file);
//    using Eigen::MatrixXd;
//    MatrixXd A(k_i.getKin().get_num_column(), k_i.getKin().get_num_column());
//    for( int i=0; i<k_i.getKin().get_num_column(); ++i ){
//        for( int j=0; j<k_i.getKin().get_num_column(); ++j ){
//            A(i,j)=k_i.getKin().get_matrix()[i][j];
//        }
//    }
//    //std::cout << "Here is a random 6x6 matrix, A:" << std::endl << A << std::endl << std::endl;
//
//    Eigen::EigenSolver<MatrixXd> es(A);
//    std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
//    //std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
//
//    ASSERT_EQ(0, 0);
//}
