//
// Created by baoxing on 3/5/18.
//

#include "Control.h"

int control(int argc, char** argv){
    //check parameter begin
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


    std::string mapFile = genotype_file + ".map";
    std::string pedFile = genotype_file + ".ped";
    std::string tfamFile = genotype_file + ".tfam";
    std::string tpedFile = genotype_file + ".tped";


    std::set<std::string> validateMehtods;
    validateMehtods.insert("emma");
    validateMehtods.insert("emmax");
    validateMehtods.insert("emmax_test_multi_allic_single_test");
    validateMehtods.insert("emmax_test_multi_allic_multi_test_null_model");
    validateMehtods.insert("emmax_test_multi_allic_multi_test_full_model");
    validateMehtods.insert("gemma");
    validateMehtods.insert("gemma_mt");
    validateMehtods.insert("gemma_mt_ma");

    if( validateMehtods.find(program) == validateMehtods.end() ){
        std::cerr << "unknown regression model" << std::endl;
        usage();
        exit(1);
    }
    //check parameter end


    if (program.compare("gemma_mt") == 0) {
        mt_gemma_test(phenotype_file, genotype_file, kinship_matrix_file, maf);
    } else if (program.compare("gemma_mt_ma") == 0) {
        mt_gemma_test_MultiAllic(phenotype_file, genotype_file, kinship_matrix_file, maf);
    } else {

        std::string format = "tfam";
        phenotype_impl pi(phenotype_file, format);
        Kinship_matrix_impl k_i(kinship_matrix_file);

        size_t number_of_individuals;
        size_t number_of_variants;
        if (genotype_file_format.compare("ped") == 0) {
            number_of_individuals = getFileLineNumber ( pedFile );
            number_of_variants  = getFileLineNumber ( mapFile );
        } else if (genotype_file_format.compare("tped") == 0) {
            number_of_individuals = getFileLineNumber ( tfamFile );
            number_of_variants  = getFileLineNumber ( tpedFile );
        } else {
            std::cerr << "unknown genotype file format" << std::endl;
            usage();
            exit(1);
        }

        // prepare data begin
        Genotype genotype = Genotype(number_of_individuals, number_of_variants);
        if (genotype_file_format.compare("ped") == 0) {
            Read_ped_file(mapFile, pedFile, genotype);
        } else if (genotype_file_format.compare("tped") == 0) {
            Read_tped_file(tfamFile, tpedFile, genotype);
        }
        genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());
        //prepare data end
        double man_l = maf * pi.getIndividual_ids().size();
        double man_u = (double)pi.getIndividual_ids().size() - man_l;
        if (program.compare("emma") == 0) {
            emma_test(pi, k_i, genotype, man_l, man_u);
        } else if (program.compare("gemma") == 0) {
            gemma_test(pi, k_i, genotype, man_l, man_u);
        }  else if (program.compare("emmax") == 0) {
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
    }

    //permutation test begin todo
    if( inputParser.cmdOptionExists("-m") ){
        int time = 100;
        if( inputParser.cmdOptionExists("-t") ) {
            time = std::stoi(inputParser.getCmdOption("-t"));
        }
        int seed = 1000;
        if( inputParser.cmdOptionExists("-s") ) {
            seed = std::stoi(inputParser.getCmdOption("-s"));
        }
    }
    //permutation test end
    return 0;
}
