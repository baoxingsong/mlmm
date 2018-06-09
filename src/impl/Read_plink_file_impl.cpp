//
// Created by baoxing on 2/17/18.
//
/*
 * Initially, I use the function split implemented in myutil to parse the genotype file while that involes some unnecessary RAM copies
 * The current implementation use more code and less friendly for code reading, but it is more efficient.
 * Parsing the initial version of the orf genotype data under the test folder for 1000 times takes about 1m 13s
 * The current version takes about 46s
 * **/
#include "Read_plink_file_impl.h"

void Read_ped_file(const std::string & mapFile, const std::string & pedFile, Genotype & genotype){ // todo it is not well test yet
    // check file begin
    int int_temp;
    uint64_t uint64_temp;
    int64_t int64_temp;
    double double_temp;

    std::string subs;
    std::istringstream iss;

    std::ifstream infile_map(mapFile);
    // since the file has been opened before, so do not check the success of opening file anymore
    uint64_t index=0;
    std::string line_map;
    while (std::getline(infile_map, line_map)){
        iss.clear();
        iss.str(line_map);
//        std::cout << line_map << std::endl;
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_variant_Vector()[index].setChromosome(subs);
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_variant_Vector()[index].setId(subs);
        iss >> double_temp;
//        std::cout << double_temp << std::endl;
        genotype.get_variant_Vector()[index].setGenetic_distance(double_temp);
        iss >> int64_temp;
//        std::cout << int64_temp << std::endl;
        genotype.get_variant_Vector()[index].setPosition(int64_temp);
        ++index;
    }
    infile_map.close();

    std::ifstream infile_ped(pedFile);
    index=0;
    std::string line_ped;
//    int8_t nt1;
    int i;
    std::istringstream iss2;
    while (std::getline(infile_ped, line_ped)){
//        std::cout << line_ped << std::endl;
        iss2.str(line_ped);
        iss2 >> subs;
        genotype.get_individual_Vector()[index].setFamily_id(subs);
        iss2 >> subs;
        genotype.get_individual_Vector()[index].setIndividual_id(subs);
        iss2 >> subs;
        genotype.get_individual_Vector()[index].setPaternal_id(subs);
        iss2 >> subs;
        genotype.get_individual_Vector()[index].setMaternal_id(subs);
        iss2 >> int_temp;
        genotype.get_individual_Vector()[index].setSex(int_temp);
        for (i=0; i< genotype.get_number_of_variant()-1; ++i){
            iss2 >> int64_temp;
            iss2 >> int64_temp;
            if (int64_temp == 0 || int64_temp>INT8_MAX || int64_temp<INT8_MIN){
                genotype.get_genotype_matrix().get_matrix()[index][i]=missing_genotype;
            }else{
                genotype.get_genotype_matrix().get_matrix()[index][i]=int64_temp;
            }
//            std::cout << int64_temp << "\t";
        }
        iss2 >> int64_temp;
        if (int64_temp == 0 || int64_temp>INT8_MAX || int64_temp<INT8_MIN){
            genotype.get_genotype_matrix().get_matrix()[index][i] = missing_genotype;
        }else{
            genotype.get_genotype_matrix().get_matrix()[index][i] = int64_temp;
        }
//        std::cout << int64_temp << std::endl;
        index++;
    }
    infile_ped.close();
}

void Read_tped_file(const std::string & tfamFile, const std::string & tpedFile, Genotype & genotype){
    std::string subs;
    std::istringstream iss;

    std::ifstream infile_tfam(tfamFile);
    uint64_t index=0;
    std::string line_tfam;
    int int_temp;
    uint64_t uint64_temp;
    int64_t int64_temp;
    double double_temp;
    while (std::getline(infile_tfam, line_tfam)){
        iss.str(line_tfam);
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_individual_Vector()[index].setFamily_id(subs);
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_individual_Vector()[index].setIndividual_id(subs);
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_individual_Vector()[index].setPaternal_id(subs);
        iss >> subs;
//        std::cout << subs << std::endl;
        genotype.get_individual_Vector()[index].setMaternal_id(subs);
        iss >> int_temp;
//        std::cout << int_temp << std::endl;
        genotype.get_individual_Vector()[index].setSex(int_temp);
        ++index;
    }
    infile_tfam.close();
    //    for( int_temp=0;int_temp<index;++int_temp ){
//        std::cout << genotype.get_individual_Vector()[int_temp].get_family_id() << std::endl;
//    }
//    std::cout << subs << std::endl;
    std::istringstream iss2;
    std::ifstream infile_tped(tpedFile);
    index=0;
    int i;
    std::string line_tped;
    while (std::getline(infile_tped, line_tped)){
//        std::cout << line_tped << std::endl;
        iss2.str(line_tped);
        iss2 >> subs;
//        std::cout << subs << std::endl;
        genotype.get_variant_Vector()[index].setChromosome(subs);
        iss2 >> subs;
//        std::cout << subs << std::endl;
        genotype.get_variant_Vector()[index].setId(subs);
        iss2 >> uint64_temp;
        genotype.get_variant_Vector()[index].setGenetic_distance(double_temp);
//        std::cout << double_temp << std::endl;
        iss2 >> uint64_temp;
        genotype.get_variant_Vector()[index].setPosition(uint64_temp);
//        std::cout << uint64_temp << std::endl;
//        std::cout << genotype.get_number_of_individual() << std::endl;
        for (i=0; i< genotype.get_number_of_individual()-1; ++i){
            iss2 >> int64_temp;
            iss2 >> int64_temp;
//            std::cout << int64_temp << "\t";
            if (int64_temp == 0 || int64_temp>INT8_MAX || int64_temp<INT8_MIN){
                genotype.get_genotype_matrix().get_matrix()[i][index] = missing_genotype;
            }else{
                genotype.get_genotype_matrix().get_matrix()[i][index] = int64_temp;
            }
        }
        iss2 >> int64_temp;
//        std::cout << int64_temp << "\t";
//        std::cout << std::endl;
        if (int64_temp == 0 || int64_temp>INT8_MAX || int64_temp<INT8_MIN){
            genotype.get_genotype_matrix().get_matrix()[i][index] = missing_genotype;
        }else{
            genotype.get_genotype_matrix().get_matrix()[i][index] = int64_temp;
        }
        ++index;
    }
    infile_tped.close();
}
