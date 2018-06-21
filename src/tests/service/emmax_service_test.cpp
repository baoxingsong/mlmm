//
// Created by baoxing on 2/25/18.
//

#include "../../service/emmax_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(emma_test_test, c1){
    double maf=0.1;
//    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
//    std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
//    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
    std::string phenotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype_2";
    std::string genotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940";
    std::string kinship_file = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";
    emma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_test, c1){
    double maf=0.1;
//    std::string phenotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
//    std::string genotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/orf";
//    std::string kinship_file = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
//    std::string phenotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mt_phenotype_2";
//    std::string genotype_path = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940";
//    std::string kinship_file = "/home/who/Dropbox/mlmm_cpp/src/tests/testData/mouse_hs1940.BN.kinf";

    std::string phenotype_path = "/media/who/8t2/22092016/107/34/phenotype.tfam";
    std::string genotype_path = "/media/who/8t2/22092016/107/34/snp";
    std::string kinship_file = "/media/who/8t2/22092016/107/34/snp.aBN.kinf";


    emmax_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}

TEST(emmax_test_multi_allic_single_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/indel_snp_from_msa";
    //std::string genotype_path = "/Users/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
    std::string kinship_file = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/snp.aBN.kinf";
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
