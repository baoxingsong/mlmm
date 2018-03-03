//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_READ_PLINK_FILE_IMPL_H
#define MLMM_CPP_READ_PLINK_FILE_IMPL_H
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <assert.h>

#include "../model/model.h"
#include "../util/util.h"
Genotype Read_ped_file(const std::string& filePath);
Genotype Read_tped_file(const std::string& filePath);
#endif //MLMM_CPP_READ_PLINK_FILE_IMPL_H
