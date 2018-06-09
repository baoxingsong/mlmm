//
// Created by Baoxing song on 06.06.18.
//

#include "phenotype_mul_impl.h"
phenotype_mul_impl::phenotype_mul_impl(){

}
phenotype_mul_impl::phenotype_mul_impl(const std::string & file_path, const std::string & format){
    std::ifstream infile(file_path);
    if( ! infile.good()){
        std::cerr << "error in opening kinship matrix file " << file_path << std::endl;
        exit(1);
    }
    std::vector<double> phenotypes_t;
    std::string line;
    while (std::getline(infile, line)){
        std::vector<std::string> elements;
        split(line, elements);
        this->individual_ids.push_back(elements[0]);
        this->family_ids.push_back(elements[1]);
        phenotypes_t.push_back(stod(elements[2]));
    }
    //this->phenotypes = My_Vector<double>(phenotypes_t); todo
}
const std::vector<std::string> &phenotype_mul_impl::getIndividual_ids() const{
    return individual_ids;
}
void phenotype_mul_impl::setIndividual_ids(const std::vector<std::string> &individual_ids){
    this->individual_ids = individual_ids;
}
const std::vector<std::string> &phenotype_mul_impl::getFamily_ids() const{
    return family_ids;
}
void phenotype_mul_impl::setFamily_ids(const std::vector<std::string> &family_ids){
    this->family_ids = family_ids;
}
const My_matrix<double> &phenotype_mul_impl::getPhenotypes() const{
    return phenotypes;
}
void phenotype_mul_impl::setPhenotypes(const My_matrix<double> &phenotypes){
    this->phenotypes = phenotypes;
}