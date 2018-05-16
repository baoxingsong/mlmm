//
// Created by song on 5/15/18.
//

#include "../../service/PerMutationTest.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"

TEST(PerMutationTest, c1){
    std::string format = "tfam";
    std::string phenotype_path = "/home/song/Dropbox/mlmm_cpp/src/tests/testData/phenotype.tfam";
    phenotype_impl pi(phenotype_path, format);
    size_t time=1000;
    size_t seed=100;
    std::string program = "emmax";
    PerMutationTest perMutationTest(time, seed, pi, program);
    ASSERT_EQ(0, 0);
}
