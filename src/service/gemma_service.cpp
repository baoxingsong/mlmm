//
// Created by song on 5/29/18.
//

#include "gemma_service.h"

void gemma_test ( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){

    My_matrix<double> w(genotype.get_number_of_individual(), 1);
    int i, j;

    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        w.get_matrix()[i][0]=1;
    }
    int ngrids = 100;
    double llim = -5.0;
    double ulim = 5.0;
    double eps = 0.00001;
    std::string method="ML";
    int maxiter = 1000;

    Eigen_result eigen_L = _get_eigen_L_( k_i.getKin());

    My_matrix<double> x(genotype.get_number_of_individual(), 1);
    double sum;
    double count;
    double mean;
    double p_val;
    std::vector <int> indexs;
    bool has_missing;
    int i_index;

    My_Vector<double> Uty(genotype.get_number_of_individual());
    My_matrix<double> xNUll(genotype.get_number_of_individual(), 0);
    _get_Uty_( eigen_L, pi.getPhenotypes(), pi.getPhenotypes().get_length(), Uty);

    double l0 = gemma_estimates (pi.getPhenotypes(), Uty, xNUll, w, eigen_L, ngrids, llim, ulim, eps, method, maxiter);

    double l1;
    double dlrt;
    for( j=0; j<genotype.get_number_of_variant(); ++j){
        indexs.clear();
        sum=0;
        count=0;
        has_missing=false;
        for( i=0; i< genotype.get_number_of_individual(); ++i ){
            if( genotype.get_genotype_matrix().get_matrix()[i][j] == missing_genotype ){
                indexs.push_back(i);
                has_missing=true;
            }else{
                sum += genotype.get_genotype_matrix().get_matrix()[i][j];
                count++;
            }
            x.get_matrix()[i][0]=genotype.get_genotype_matrix().get_matrix()[i][j];
        }
        if( has_missing ){
            mean=sum/count;
            for( i_index=0; i_index<indexs.size(); ++i_index ){
                x.get_matrix()[indexs[i_index]][0]=mean;
            }
        }
        std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t" << genotype.get_variant_Vector()[j].getPosition() << "\t" << genotype.get_variant_Vector()[j].getId();


        l1 = gemma_estimates (pi.getPhenotypes(), Uty, x, w, eigen_L, ngrids, llim, ulim, eps, method, maxiter);
//        printf("\t%10.20f\t%10.20f", l0, l1);
        //dlrt = 2 * log(l1/l0);
        dlrt = 2 * (l1-l0); // this function on the original publication is confusing
        p_val = 1- chii(dlrt, 1);

        if( p_val == 0.0 ){
            p_val = 0.00000000000000000001;
        }
        printf("\t%10.20f\t%10.20f\n", dlrt, p_val);
        //return;
    }
}

void gemma_test ( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf){
    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype =  Genotype (number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;

    gemma_test ( pi, k_i, genotype, man_l, man_u );
}