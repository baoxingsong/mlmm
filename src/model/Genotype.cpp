//
// Created by baoxing on 2/17/18.
//

#include "Genotype.h"
Genotype::Genotype(){
}

Genotype::Genotype(const size_t & number_of_individuals, const size_t & number_of_variants){
    this->genotype_matrix = My_matrix<int8_t> (number_of_individuals, number_of_variants);
    this->number_of_individual=number_of_individuals;
    this->number_of_variant=number_of_variants;
    variant_Vector = new Variant[number_of_variants];
    individual_Vector = new Individual[number_of_individuals];
}

My_matrix<int8_t > & Genotype::get_genotype_matrix(){
    return this->genotype_matrix;
}

Variant * Genotype::get_variant_Vector(){
    return this->variant_Vector;
}
Individual * Genotype::get_individual_Vector(){
    return this->individual_Vector;
}
void Genotype::onlyKeepThoseIndividuls(const std::vector <std::string> & individual_ids) {
    size_t number_of_matched_individual=0;
    int i, j, z;
    std::string individual_id;
    bool ifFound;
    for (i = 0; i < individual_ids.size(); ++i) {
        ifFound = false;
        individual_id = individual_ids[i];
        for (j = 0; j < this->number_of_individual; ++j) {
            Individual individual = this->individual_Vector[j];
            if (individual.get_individual_id().compare(individual_id) == 0) {
                ++number_of_matched_individual;
                ifFound = true;
            }
        }
        if (!ifFound) {
            std::cerr << individual_id << " had phenotype available, but the genotype is absent " << std::endl;
            exit(1);
        }
    }

    Individual * _new_individuals = new Individual[number_of_matched_individual];
    My_matrix<int8_t> _new_genotype_matrix(number_of_matched_individual, this->number_of_variant);
    int new_individual_id=0;
    for (i = 0; i < individual_ids.size(); ++i) {
        individual_id = individual_ids[i];
        for (j = 0; j < this->number_of_individual; ++j) {
            Individual individual = this->individual_Vector[j];
            if (individual.get_individual_id().compare(individual_id) == 0) {
                _new_individuals[new_individual_id] = individual;
                for( z=0; z < this->number_of_variant; ++z ) {
                    _new_genotype_matrix.get_matrix()[new_individual_id][z] = this->genotype_matrix.get_matrix()[j][z];
                }
                ++new_individual_id;
            }
        }
    }
    delete [] this->individual_Vector;
    this->number_of_individual = number_of_matched_individual;
    this->genotype_matrix = _new_genotype_matrix;
    this->individual_Vector = _new_individuals;
}

size_t Genotype::get_number_of_individual(){
    return this->number_of_individual;
}
size_t Genotype::get_number_of_variant(){
    return this->number_of_variant;
}
