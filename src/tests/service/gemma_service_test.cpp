//
// Created by baoxing on 2/25/18.
//

#include "../../service/gemma_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(gemma_test_test, c1){
    double maf=0.1;
//    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
//    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
//    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    std::string phenotype_path = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/phenotype.tfam";
    std::string genotype_path = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/snp.356.scan_test";
    std::string kinship_file = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/snp.356.cXX.txt";
    gemma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
