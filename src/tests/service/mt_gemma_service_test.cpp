//
// Created by Baoxing song on 14.06.18.
//
#include "../../service/mt_gemma_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(mt_gemma_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    mt_gemma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
