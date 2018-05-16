//
// Created by song on 5/15/18.
//

#include "PerMutationTest.h"
PerMutationTest::PerMutationTest(const size_t & _time, const size_t & _seed, const phenotype_impl & _pi, const std::string & _program){
    this->time = _time;
    this->seed = _seed;
    this->pi = _pi;
    this->program = _program;
    uint32_t arrayIndex[this->pi.getPhenotypes().get_length()];
    size_t i;
    for( i=0; i<this->pi.getPhenotypes().get_length(); ++i ){
        arrayIndex[i]=i;
    }
    std::random_shuffle(arrayIndex, arrayIndex+this->pi.getPhenotypes().get_length());
    for( i=0; i<this->pi.getPhenotypes().get_length(); ++i ){
        std::cout << i << ": " << arrayIndex[i] << std::endl;
    }
}
