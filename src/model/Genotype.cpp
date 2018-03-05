//
// Created by baoxing on 2/17/18.
//

#include "Genotype.h"


Genotype::Genotype(){
    this->number_of_individual = 0;
    this->number_of_variant = 0;
}
Genotype::Genotype(const unsigned long& _number_of_individual,
                   const unsigned long& _number_of_variant, const std::vector<Variant>& _variant_Vector,
                   const std::vector<Individual>& _individual_Vector, const std::vector<int*>& genotypes, const std::string& format){
    this->number_of_individual = _number_of_individual;
    this->number_of_variant = _number_of_variant;
    this->variant_Vector = _variant_Vector;
    this->individual_Vector = _individual_Vector;

    this->genotype_matrix = My_matrix<unsigned int > (_number_of_individual, _number_of_variant);
    if( format == "ped" ){
        for (int i = 0; i < _number_of_individual; ++i ){
            for( int j=0; j <number_of_variant; ++j ){
                this->genotype_matrix.get_matrix()[i][j] = genotypes[i][j];
            }
        }
    }else if( format == "tped"){
        for (int i = 0; i < _number_of_individual; ++i){
            for ( int j=0; j<_number_of_variant; ++j ){
                this->genotype_matrix.get_matrix()[i][j] = genotypes[j][i];
            }
        }
    }
}

My_matrix<unsigned int > & Genotype::get_genotype_matrix(){
    return this->genotype_matrix;
}

const unsigned long& Genotype::get_number_of_individual() const{
    return this->number_of_individual;
}

const unsigned long& Genotype::get_number_of_variant() const{
    return this->number_of_variant;
}

const std::vector<Variant>& Genotype::get_variant_Vector() const{
    return this->variant_Vector;
}

const std::vector<Individual>& Genotype::get_individual_Vector() const{
    return this->individual_Vector;
}
void Genotype::onlyKeepThoseIndividuls(const std::vector <std::string> & individual_ids){
    this->number_of_individual = individual_ids.size();
    std::vector<Individual> _new_individuals;
    My_matrix<unsigned int > _new_genotype_matrix(this->number_of_individual, this->number_of_variant);
    int i, j;
    std::string individual_id;
    bool ifFound;
    for( i=0; i<individual_ids.size() ; ++i){
        ifFound = false;
        individual_id = individual_ids[i];
        for( j=0; j<this->individual_Vector.size(); ++j ){
            Individual individual = this->individual_Vector[j];
            if( individual.get_individual_id().compare(individual_id) == 0 ){
                ifFound = true;
                _new_individuals.push_back(individual);
                _new_genotype_matrix.get_matrix()[i] = this->genotype_matrix.get_matrix()[i];
            }
        }
        if( ! ifFound ){
            std::cerr << "could not find the genotype of individual." << std::endl;
            exit(1);
        }
    }

    this->genotype_matrix = _new_genotype_matrix;
    this->individual_Vector = _new_individuals;
}
