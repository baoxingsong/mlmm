//
// Created by baoxing on 2/25/18.
//

#include "../../service/gemma_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(gemma_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype_2";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    gemma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
