//
// Created by Baoxing song on 14.06.18.
//
#include "../../service/mt_gemma_service.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(mt_gemma_test_test, c1){
    double maf=0.1;
    std::string phenotype_path = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/phenotype.tfam";
    std::string genotype_path = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/snp.356.scan_test";
    std::string kinship_file = "/netscratch/dep_tsiantis/grp_gan/song/multipleTritGwasTestData/snp.356.cXX.txt";
    mt_gemma_test ( phenotype_path, genotype_path, kinship_file, maf );
    ASSERT_EQ(0, 0);
}
