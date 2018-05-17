//
// Created by Baoxing song on 16.05.18.
//

#include "../../impl/phenotype_impl.h"
#include "../../model/model.h"
#include "../../impl/impl.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"
#include "../../impl/Read_plink_file_impl.h"


TEST(Read_tped_file, c1){
    for( int32_t i =0 ; i<1000; ++i) {
        std::string format = "tfam";
        std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
        std::string tfamFile = genotype_path + ".tfam";
        std::string tpedFile = genotype_path + ".tped";
        uint64_t number_of_individuals = getFileLineNumber(tfamFile);
        uint64_t number_of_variants = getFileLineNumber(tpedFile);
        Genotype genotype = Genotype(number_of_individuals, number_of_variants);
        Read_tped_file(tfamFile, tpedFile, genotype);
        ASSERT_EQ(0, 0);
    }
}

TEST(Read_ped_file, c1){
    for( int32_t i =0 ; i<1000; ++i){
        std::string genotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/orf";
        std::string mapFile = genotype_path + ".map";
        std::string pedFile = genotype_path + ".ped";
        uint64_t number_of_individuals = getFileLineNumber ( pedFile );
        uint64_t number_of_variants  = getFileLineNumber ( mapFile );
        Genotype genotype =  Genotype(number_of_individuals, number_of_variants);
        Read_ped_file(mapFile, pedFile, genotype);
        ASSERT_EQ(0, 0);
    }
}
