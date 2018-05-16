//
// Created by song on 5/15/18.
//

#ifndef MLMM_CPP_PERMUTATIONTEST_H
#define MLMM_CPP_PERMUTATIONTEST_H

#include "../service/service.h"
class PerMutationTest{
    private:
        size_t time;
        size_t seed;
        phenotype_impl pi;
        std::string program;
    public:
        PerMutationTest(const size_t & _time, const size_t & _seed, const phenotype_impl & _pi, const std::string & _program);
};
#endif //MLMM_CPP_PERMUTATIONTEST_H
