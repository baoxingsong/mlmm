//
// Created by baoxing on 2/17/18.
//

#include "Read_plink_file_impl.h"

void Read_ped_file(const std::string & mapFile, const std::string & pedFile, Genotype & genotype){
    // check file begin
    std::ifstream infile_map(mapFile);
    uint64_t index=0;
    std::string line_map;
    while (std::getline(infile_map, line_map)){
        std::vector<std::string> elements;
        split(line_map, elements);
        Variant variant( elements[0], elements[1], stof(elements[2]), stoul(elements[3]));
        genotype.get_variant_Vector()[index]=variant;
        ++index;
    }
    infile_map.close();

    std::ifstream infile_ped(pedFile);
    index=0;
    std::string line_ped;
    while (std::getline(infile_map, line_ped)){
        std::vector<std::string> elements;
        split(line_ped, elements);
        Individual individual(elements[0], elements[1], elements[2], elements[3], stoi(elements[4]) );
        genotype.get_individual_Vector()[index] = individual;
        for (int i=6; i< 2 * genotype.get_number_of_variant() + 6; i+=2){
            int nt1 = stoi(elements[i]);
            if (nt1 == 0 || nt1>INT8_MAX || nt1<INT8_MIN){
                genotype.get_genotype_matrix().get_matrix()[index][(i-6)/2]=missing_genotype;
            }else{
                genotype.get_genotype_matrix().get_matrix()[index][(i-6)/2]=nt1;
            }
        }
        index++;
    }
    infile_ped.close();
}

void Read_tped_file(const std::string & tfamFile, const std::string & tpedFile, Genotype & genotype){
    std::ifstream infile_tfam(tfamFile);
    uint64_t index=0;
    std::vector<std::string> elements;
    std::string line_tfam;
    while (std::getline(infile_tfam, line_tfam)){
        elements.clear();
        split(line_tfam, elements);
        Individual individual( elements[0], elements[1], elements[2], elements[3], stoi(elements[4]) );
        genotype.get_individual_Vector()[index]=individual;
        ++index;
    }
    infile_tfam.close();

    std::ifstream infile_tped(tpedFile);
    index=0;
    std::string line_tped;
    while (std::getline(infile_tped, line_tped)){
        elements.clear();
        split(line_tped, elements);
        Variant variant( elements[0], elements[1], stof(elements[2]), stoul(elements[3]));
        genotype.get_variant_Vector()[index] = variant;
        for (int i=4; i< 2 * genotype.get_number_of_individual() + 4; i+=2){
            int nt1 = stoi(elements[i]);
            if (nt1 == 0 || nt1>INT8_MAX || nt1<INT8_MIN){
                genotype.get_genotype_matrix().get_matrix()[(i-4)/2][index] = missing_genotype;
            }else{
                genotype.get_genotype_matrix().get_matrix()[(i-4)/2][index] = nt1;
            }
        }
        ++index;
    }
    infile_tped.close();
}
