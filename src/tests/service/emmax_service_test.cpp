//
// Created by baoxing on 2/25/18.
//

#include "../../service/emmax_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(emma_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    emma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    emmax_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_multi_allic_single_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/indel_snp_from_msa";
    //std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    emmax_test_multi_allic_single_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_multi_allic_multi_test_null_model, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/indel_snp_from_msa";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    emmax_test_multi_allic_multi_test_null_model ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_multi_allic_multi_test_full_model, c1){
    double maf=0.1;
    std::string phenotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/indel_snp_from_msa";
    std::string kinship_file = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    emmax_test_multi_allic_multi_test_full_model ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
