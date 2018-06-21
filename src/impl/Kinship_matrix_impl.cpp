//
// Created by baoxing on 2/17/18.
//

#include "Kinship_matrix_impl.h"

double get_ibs_score( const int& i, const int& j, Genotype& genotype){
    unsigned int total_number = 0;
    unsigned int identical_number = 0;
    for( int ij=0; ij<genotype.get_number_of_variant(); ++ij ){
        if( genotype.get_genotype_matrix().get_matrix()[i][ij] != missing_genotype && genotype.get_genotype_matrix().get_matrix()[j][ij] != missing_genotype){
            ++total_number;
            if( genotype.get_genotype_matrix().get_matrix()[i][ij] == genotype.get_genotype_matrix().get_matrix()[j][ij] ){
                ++identical_number;
            }
        }
    }
    return double(identical_number)/double(total_number);
}

Kinship_matrix_impl::Kinship_matrix_impl( Genotype& genotype ){
    My_matrix<double> k(genotype.get_number_of_individual(), genotype.get_number_of_individual());
    this->kin=k;
    for( int i=1; i< genotype.get_number_of_individual(); ++i){
        kin.get_matrix()[i][i]=1;
        for( int j=0; j<i; ++j ){
            kin.get_matrix()[i][j] = get_ibs_score(i, j, genotype);
            kin.get_matrix()[j][i] = kin.get_matrix()[i][j];
        }
    }
}

void Kinship_matrix_impl::scale_k(){
    assert(this->kin.get_num_row() == this->kin.get_num_column());
    My_matrix<double> first_half(this->kin.get_num_row(), this->kin.get_num_row());
    int i,j;
    double element = -1.0/this->kin.get_num_column();
    for ( i=0; i<this->kin.get_num_column(); ++i ){
        first_half.get_matrix()[i][i] = 1+element;
        for( j=0; j<i; ++j ){
            first_half.get_matrix()[i][j] = element;
            first_half.get_matrix()[j][i] = element;
        }
    }
    My_matrix<double> c(first_half.get_num_row(), this->kin.get_num_column());
    double c_d = 0.0;
    for ( i=0; i<this->kin.get_num_column(); ++i ){
        for( j=0; j<this->kin.get_num_column(); ++j ){
            c_d += c.get_matrix()[i][j];
        }
    }
    double scalar = (this->kin.get_num_column() - 1.0) / c_d;
    for ( i=0; i<this->kin.get_num_column(); ++i ){
        for( j=0; j<this->kin.get_num_column(); ++j ){
            this->kin.get_matrix()[i][j]=this->kin.get_matrix()[i][j]*scalar;
        }
    }
}

Kinship_matrix_impl::Kinship_matrix_impl(const std::string & kin_file_path){
    std::ifstream infile(kin_file_path);
    if( ! infile.good()){
        std::cerr << "error in opening kinship matrix file " << kin_file_path << std::endl;
        exit (1);
    }
    std::vector<Variant> variants;
    std::string line;

    std::getline(infile, line);
    std::vector<std::string> elements;
    split(line, elements);
    My_matrix<double> k(elements.size(), elements.size());
    this->kin=k;

    infile.clear();
    infile.seekg(0, std::ios::beg);

    int line_index=0, i;
    while (std::getline(infile, line)){
        std::vector<std::string> elements_t;
        split(line, elements_t);
//        assert(elements_t.size() == k.get_num_column());
        for( i=0; i<this->kin.get_num_column(); ++i  ){
            this->kin.get_matrix()[line_index][i] = stod(elements_t[i]);
        }
        ++line_index;
    }
}

const My_matrix<double> &Kinship_matrix_impl::getKin() const {
    return kin;
}

void Kinship_matrix_impl::setKin(const My_matrix<double> &kin) {
    Kinship_matrix_impl::kin = kin;
}
