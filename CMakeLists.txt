cmake_minimum_required(VERSION 3.0)
project(mlmm_cpp CXX)
# here four variables is defined ${mlmm_cpp_BINARY_DIR} and ${mlmm_cpp_SOURCE_DIR}
# and ${mlmm_cpp_BINARY_DIR} = ${PROJECT_BINARY_DIR}
# ${mlmm_cpp_SOURCE_DIR} = ${PROJECT_SOURCE_DIR}
# so it is better to use the latter two variables

# set grammer SET(VAR [VALUE] [CACHE TYPE DOCSTRING [FORCE]])
SET(CMAKE_CXX_STANDARD 11 )
ADD_DEFINITIONS(-fno-elide-constructors) # it is important for the My_matrix class

set (CMAKE_CXX_FLAGS "-O3 -msse2 -mavx2 -Wall -std=c++0x -pthread -I ./ -std=c++11" CACHE STRING "compile flags" FORCE)
if (APPLE)
    add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
    add_definitions(-D__GLIBCXX__)
endif (APPLE)

if (WIN32)
    MESSAGE( FATAL_ERROR "we do not support windows platform" )
endif (WIN32)

if (WIN64)
    MESSAGE( FATAL_ERROR "we do not support windows platform" )
endif (WIN64)

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
    message("debug mode")
endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")


#for google test begin
add_subdirectory(./googletest/googletest)
include_directories(googletest/googletest/include googletest/googletest)

# https://cmake.org/cmake/help/v3.0/command/message.html
MESSAGE( STATUS "This is BINARY DIR " ${mlmm_cpp_BINARY_DIR} )
MESSAGE( STATUS "This is SOURCE DIR " ${mlmm_cpp_SOURCE_DIR} )

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
        #       src/impl/gemma_impl.cpp
        #      src/impl/gemma_impl.h
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
        src/service/emmax_service.cpp
        src/service/emmax_service.h
        src/impl/phenotype_impl.cpp
        src/impl/phenotype_impl.h
        src/tests/impl/phenotype_impl.test.cpp
        src/tests/service/emmax_service_test.cpp
        src/tests/impl/emmax_impl_test.cpp
        src/tests/impl/my_math_impl_test.cpp
        src/model/Emma_result.cpp
        src/model/Emma_result.h
        src/service/service.h
        src/Control.cpp src/Control.h
        src/impl/Average_information_impl.cpp
        src/impl/Average_information_impl.h
        src/service/PerMutationTest.cpp
        src/service/PerMutationTest.h
        src/tests/service/PermutationTest.cpp)


target_link_libraries( mlmm gtest gtest_main )
