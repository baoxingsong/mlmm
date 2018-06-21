//
// Created by Baoxing song on 14.06.18.
//

#include "mt_gemma_service.h"

void mt_gemma_test ( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf){

    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    Kinship_matrix_impl k_i(kinship_file);

    Eigen_result eigen_r = _get_eigen_L_( k_i.getKin());
    const char method = 'R';
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype =  Genotype (number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);

    My_matrix<double> w(genotype.get_number_of_individual(), 1);
    int i, j;
    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        w.get_matrix()[i][0]=1;
    }

    My_matrix<double> UtW(genotype.get_number_of_individual(), 1);
    trmul(eigen_r.get_eigen_vectors(), w, UtW);

    std::ifstream infile(phenotype_path);
    if( ! infile.good()){
        std::cerr << "error in opening phenotype file " << phenotype_path << std::endl;
        exit(1);
    }

    std::vector<double> phenotypes_t;
    std::string line;
    std::getline(infile, line);
    std::vector<std::string> elements_temp;
    split(line, elements_temp);
    My_matrix<double> Y(number_of_individuals, elements_temp.size());
    infile.close();
    infile=std::ifstream(phenotype_path);
    i=0;

    while (std::getline(infile, line)){
        std::vector<std::string> elements;
        split(line, elements);
        for (j=0; j<elements_temp.size();++j) {
            Y.get_matrix()[i][j]=stod(elements[j]);
        }
        ++i;
    }
    My_matrix<double> UtY(number_of_individuals, elements_temp.size());
    trmul(eigen_r.get_eigen_vectors(), Y, UtY);
    AnalyzePlink( eigen_r, UtW, UtY, method, k_i, genotype, Y);
}