cmake_minimum_required(VERSION 3.8)
project(mlmm_cpp)

set(CMAKE_CXX_STANDARD 11 )
add_definitions(-fno-elide-constructors) # it is important for the My_matrix class

if (APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
    add_definitions(-D__GLIBCXX__)
endif (APPLE)


#for google test begin
add_subdirectory(./googletest/googletest)


include_directories(googletest/googletest/include googletest/googletest)

add_executable(mlmm main.cpp
        googletest
        googletest/googletest
        src/InputParser.cpp
        src/InputParser.h
        src/model/Genotype.cpp
        src/model/Genotype.h
        src/model/Variant.cpp
        src/model/Variant.h
        src/model/Individual.cpp
        src/model/Individual.h
        src/impl/Read_plink_file_impl.cpp
        src/impl/Read_plink_file_impl.h
        src/util/util.h
        src/impl/Kinship_matrix_impl.cpp
        src/impl/Kinship_matrix_impl.h
        src/util/myutil.cpp
        src/util/myutil.h
        src/impl/emmax_impl.cpp
        src/impl/emmax_impl.h
        src/impl/my_math_impl.cpp
        src/impl/my_math_impl.h
        src/model/My_matrix.h
        src/model/My_Vector.h
        src/impl/My_linear_algebra_impl.cpp
        src/impl/My_linear_algebra_impl.h
        src/model/Eigen_result.cpp
        src/model/Eigen_result.h

        src/tests/impl/My_linear_algebra_impl_test.cpp

        src/model/Qr_decomposition_result.cpp
        src/model/Qr_decomposition_result.h

        src/impl/impl.h
        src/model/model.h
        src/service/emmax_service.cpp src/service/emmax_service.h src/impl/phenotype_impl.cpp src/impl/phenotype_impl.h src/tests/impl/phenotype_impl.test.cpp src/tests/service/emmax_service_test.cpp src/tests/impl/emmax_impl_test.cpp src/tests/impl/my_math_impl_test.cpp src/model/Emma_result.cpp src/model/Emma_result.h src/service/service.h)


target_link_libraries( mlmm gtest gtest_main )
