//
// Created by baoxing on 3/5/18.
//

#include "Control.h"

int control(int argc, char** argv){

    std::string program = argv[1];
    if (program.compare("-h") == 0 || program.compare("--help") == 0) {
        usage();
        exit(1);
    }

    InputParser inputParser(argc, argv);
    std::string phenotype_file = inputParser.getCmdOption("-p");
    std::string genotype_file = inputParser.getCmdOption("-g");
    std::string genotype_file_format = inputParser.getCmdOption("-f");
    std::string kinship_matrix_file = inputParser.getCmdOption("-k");
    double maf = stod(inputParser.getCmdOption("-m"));
    std::string format = "tfam";

    Genotype genotype;

    phenotype_impl pi(phenotype_file, format);
    Kinship_matrix_impl k_i(kinship_matrix_file);

    if (genotype_file_format.compare("ped") == 0) {
        genotype = Read_ped_file(genotype_file);
    } else if (genotype_file_format.compare("tped") == 0) {
        genotype = Read_tped_file(genotype_file);
    } else {
        std::cerr << "unknown genotype file format" << std::endl;
        usage();
        exit(1);
    }
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    std::set<std::string> validateMehtods;
    validateMehtods.insert("emma");
    validateMehtods.insert("emmax");
    validateMehtods.insert("emmax_test_multi_allic_single_test");
    validateMehtods.insert("emmax_test_multi_allic_multi_test_null_model");
    validateMehtods.insert("emmax_test_multi_allic_multi_test_full_model");

    if( validateMehtods.find(program) == validateMehtods.end() ){
        std::cerr << "unknown regression model" << std::endl;
        usage();
        exit(1);
    }
    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;

    if (program.compare("emma") == 0) {
        emma_test(pi, k_i, genotype, man_l, man_u);
    } else if (program.compare("emmax") == 0) {
        emmax_test(pi, k_i, genotype, man_l, man_u);
    } else if (program.compare("emmax_test_multi_allic_single_test") == 0) {
        emmax_test_multi_allic_single_test(pi, k_i, genotype, man_l, man_u);
    } else if (program.compare("emmax_test_multi_allic_multi_test_null_model") == 0) {
        emmax_test_multi_allic_multi_test_null_model(pi, k_i, genotype, man_l, man_u);
    } else if (program.compare("emmax_test_multi_allic_multi_test_full_model") == 0) {
        emmax_test_multi_allic_multi_test_full_model(pi, k_i, genotype, man_l, man_u);
    } else {
        std::cerr << "unknown regression model" << std::endl;
        usage();
        exit(1);
    }
    return 0;
}
