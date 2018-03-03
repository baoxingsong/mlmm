//
// Created by baoxing on 2/17/18.
//

#include "Read_plink_file_impl.h"

Genotype Read_ped_file(const std::string& filePath){
    std::string mapFile = filePath + ".map";
    std::string pedFile = filePath + ".ped";

    std::ifstream infile_map(mapFile);
    if( ! infile_map.good()){
        std::cerr << "error in opening map file " << filePath << std::endl;
        exit (1);
    }
    std::vector<Variant> variants;
    std::string line_map;
    while (std::getline(infile_map, line_map)){
        std::vector<std::string> elements;
        split(line_map, elements);
        Variant variant( elements[0], elements[1], stof(elements[2]), stoul(elements[3]));
        variants.push_back(variant);
    }
    infile_map.close();


    std::ifstream infile_ped(pedFile);
    if( ! infile_ped.good()){
        std::cerr << "error in opening ped file " << filePath << std::endl;
        exit (1);
    }

    unsigned long number_of_variant = variants.size();

    std::vector<int*> genotypes;
    std::vector<Individual> individuals;
    std::string line_ped;
    while (std::getline(infile_map, line_ped)){
        std::vector<std::string> elements;
        split(line_ped, elements);
        Individual individual(elements[0], elements[1], elements[2], elements[3], stoi(elements[4]) );
        individuals.push_back(individual);
        int* snps = new int[number_of_variant];
        for (int i=6; i< 2 * number_of_variant + 6; i+=2){
            int nt1 = stoi(elements[i]);
            if (nt1 == 0){
                snps[(i-6)/2] = missing_genotype;
            }else{
                snps[(i-6)/2]=nt1;
            }
        }
        genotypes.push_back(snps);
    }
    std::string format = "ped";
    assert(genotypes.size() == individuals.size());
    Genotype genotype(genotypes.size(), number_of_variant, variants, individuals, genotypes, format);
    infile_ped.close();

    for( int i=0; i<genotypes.size(); ++i ){
        delete [] genotypes[i];
    }

    return genotype;
}


Genotype Read_tped_file(const std::string& filePath){
    std::string tfamFile = filePath + ".tfam";
    std::string tpedFile = filePath + ".tped";

    std::ifstream infile_tfam(tfamFile);
    if( ! infile_tfam.good()){
        std::cerr << "error in opening map file " << filePath << std::endl;
        exit (1);
    }
    std::vector<Individual> individuals;
    std::string line_tfam;
    while (std::getline(infile_tfam, line_tfam)){
        std::vector<std::string> elements;
        split(line_tfam, elements);
        Individual individual( elements[0], elements[1], elements[2], elements[3], stoi(elements[4]) );
        individuals.push_back(individual);
    }
    infile_tfam.close();

    unsigned long number_of_individual = individuals.size();

    std::ifstream infile_tped(tpedFile);
    if( ! infile_tped.good()){
        std::cerr << "error in opening ped file " << filePath << std::endl;
        exit (1);
    }
    std::vector<int*> genotypes;
    std::vector<Variant> variants;
    std::string line_tped;
    while (std::getline(infile_tped, line_tped)){
        std::vector<std::string> elements;
        split(line_tped, elements);
        Variant variant( elements[0], elements[1], stof(elements[2]), stoul(elements[3]));
        variants.push_back(variant);

        int* snps = new int[individuals.size()];
        for (int i=4; i< 2 * number_of_individual + 4; i+=2){
            int nt1 = stoi(elements[i]);
            if (nt1 == 0){
                snps[(i-4)/2] = missing_genotype;
            }else{
                snps[(i-4)/2]=nt1;
            }
        }
        genotypes.push_back(snps);
    }
    std::string format = "tped";
    assert(genotypes.size() == variants.size());
    Genotype genotype(number_of_individual, genotypes.size(), variants, individuals, genotypes, format);
    infile_tped.close();
    for( int i=0; i<genotypes.size(); ++i ){
        delete [] genotypes[i];
    }
    return genotype;
}
