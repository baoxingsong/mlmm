/*
 * =====================================================================================
 *
 *       Filename:  InputParser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string InputParser::getCmdOption( std::string &option) {
    std::vector<std::string>::iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return "NA";
}

std::string InputParser::getCmdOption( const char* o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}

bool InputParser::cmdOptionExists( const char* o){
    std::string option = o;
    return cmdOptionExists(option);
}

void usage( ) {
    std::string progName = "mlmm_cpp";
    std::cout << "Program " << progName << std::endl <<
              "Usage:  " << progName << " <command> [options]" << std::endl <<
              "Commands:    emma/emmax/emmax_test_multi_allic_single_test/emmax_test_multi_allic_multi_test_null_model/emmax_test_multi_allic_multi_test_full_model" << std::endl <<
              "Options" << std::endl <<
              "   -p        phenotype file" << std::endl <<
              "   -g        genotype file" << std::endl <<
              "   -f        genotype file format ped/tped" << std::endl<<
              "   -m        minor allele frequency for filtering" << std::endl<<
              "   -k        kinship matrix file" << std::endl;

}
