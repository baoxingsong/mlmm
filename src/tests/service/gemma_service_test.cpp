//
// Created by baoxing on 2/25/18.
//

#include "../../service/gemma_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(gemma_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    gemma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
