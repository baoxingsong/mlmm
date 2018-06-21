//
// Created by baoxing on 2/25/18.
//

#include "phenotype_impl.h"
phenotype_impl::phenotype_impl(){

}
phenotype_impl::phenotype_impl(const std::string & file_path, const std::string & format){
    std::ifstream infile(file_path);
    if( ! infile.good()){
        std::cerr << "error in opening phenotype file " << file_path << std::endl;
        exit(1);
    }
    std::vector<double> phenotypes_t;
    std::string line;
    while (std::getline(infile, line)){
        std::vector<std::string> elements;
        split(line, elements);
        this->individual_ids.push_back(elements[1]);
        this->family_ids.push_back(elements[0]);
        phenotypes_t.push_back(stod(elements[2]));
    }
    this->phenotypes = My_Vector<double>(phenotypes_t);
}

const std::vector<std::string> &phenotype_impl::getIndividual_ids() const {
    return individual_ids;
}

void phenotype_impl::setIndividual_ids(const std::vector<std::string> &individual_ids) {
    this->individual_ids = individual_ids;
}

const std::vector<std::string> &phenotype_impl::getFamily_ids() const {
    return family_ids;
}

void phenotype_impl::setFamily_ids(const std::vector<std::string> &family_ids) {
    this->family_ids = family_ids;
}

const My_Vector<double> &phenotype_impl::getPhenotypes() const {
    return phenotypes;
}

void phenotype_impl::setPhenotypes(const My_Vector<double> &phenotypes) {
    this->phenotypes = phenotypes;
}
