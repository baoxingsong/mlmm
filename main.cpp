#include "./googletest/googletest/include/gtest/gtest.h"

#include "./src/Control.h"

int main(int argc, char ** argv) {
//    testing::InitGoogleTest(&argc, argv);
//    RUN_ALL_TESTS();
    if( argc<=1 ){
        usage();
        return 1;
    }
    control(argc, argv);
    return 0;
}
