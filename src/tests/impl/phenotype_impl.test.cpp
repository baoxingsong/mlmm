//
// Created by baoxing on 2/25/18.
//

#include "../../impl/phenotype_impl.h"
#include "../../model/model.h"
#include "../../impl/impl.h"
#include "../../../googletest/googletest/include/gtest/gtest.h"


TEST(phenotype_impl, c1){
    std::string file_path = "./src/tests/testData/phenotype.tfam";
    std::string format = "tfam";
    phenotype_impl(file_path, format);
    ASSERT_EQ(0, 0);
}
