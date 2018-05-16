// =====================================================================================
//
//       Filename:  emmax_service.cpp
//
//    Description:
//
//        Version:  1.0
//        Created:  02/27/2018 05:53:23 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
//        Company:  unemployed
//
// =====================================================================================
/*************************************************************************
 *
 *
 * compine several functions and run the linear mixed model under the
 * statistical framework of EMMA/EMMAX
 *
 *
 *
************************************************************************/

#include "emmax_service.h"


void emma_test ( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){
    My_matrix<double> x(genotype.get_number_of_individual(), 1);
    int i, j;

    for( i=0; i< genotype.get_number_of_individual(); ++i ){
        x.get_matrix()[i][0]=1;
    }

    int ngrids = 50;
    double llim = -10.0;
    double ulim = 10.0;
    double eps = 0.0000001;
    std::string method="REML";
    int maxiter = 100;

    Eigen_result eigen_L = _get_eigen_L_( k_i.getKin());

    My_matrix<double> xs(genotype.get_number_of_individual(), 1);
    double sum;
    double count;
    double mean;
    double p_val;
    std::vector <int> indexs;
    bool has_missing;


    for( j=0; j<genotype.get_number_of_variant(); ++j ){
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
            xs.get_matrix()[i][0]=genotype.get_genotype_matrix().get_matrix()[i][j];
        }
        if( has_missing ){
            mean=sum/count;
            for( int i_index : indexs ){
                xs.get_matrix()[i_index][0]=mean;
            }
        }
        std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t" << genotype.get_variant_Vector()[j].getPosition() << "\t" << genotype.get_variant_Vector()[j].getId();
        p_val =get_estimates ( pi.getPhenotypes(), k_i.getKin(), x, xs, eigen_L, ngrids, llim, ulim, eps, method, maxiter).getPvalue();

        if( p_val == 0.0 ){
            p_val = 0.00000000000000000001;
        }
        printf("\t%10.20f\n", p_val);
    }
}

Emma_result estimate( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype ){
    My_matrix<double> x(genotype.get_number_of_individual(), 1);
    int i;
    for( i=0; i<genotype.get_number_of_individual(); ++i ){
        x.get_matrix()[i][0]=1;
    }
    int ngrids = 50;
    double llim = -10.0;
    double ulim = 10.0;
    double eps = 0.0000001;
    std::string method="REML";
    int maxiter = 100;
    Eigen_result eigen_L = _get_eigen_L_( k_i.getKin());
    My_matrix<double> x0(genotype.get_number_of_individual(), 0);
    return get_estimates ( pi.getPhenotypes(), k_i.getKin(), x, x0, eigen_L, ngrids, llim, ulim, eps, method, maxiter); // do some thing here and get several parameters

}

void emmax_test( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){
    int i, j, j_1;
    Emma_result e_er = estimate( pi, k_i, genotype );
    int covariantsNumber = 1;
    int n = genotype.get_number_of_individual();
    int q = covariantsNumber + 1; //test each genotypic variants one by one
    int p = n - q;

    // claim those variables firstly could save RAM operation and save computational time begin
    My_matrix<double> full_x(n, 2);
    for( i=0; i< n; ++i ){
        full_x.get_matrix()[i][0]=1;
    }
    int freedome1 = 1;
    double mahalanobis_rss, f_stat, p_val;
    double h0_rss = e_er.getMahalanobis_rss();
    My_matrix<double> X_t(n, 2);
    My_Vector<double> beta_est(2);
    My_Vector<double> x_beta(n);
    // claim those variables firstly could save RAM operation and save computational time end
    printf("Chr\tposition\tid\tp_val\n");
    double sum;
    double count;
    double mean;
    std::vector <int> indexs;
    bool has_missing;
    for( j=0; j<genotype.get_number_of_variant(); ++j ) {
        indexs.clear();
        sum=0.0;
        count=0.0;
        has_missing=false;
        for( i=0; i< n; ++i ){
            if( genotype.get_genotype_matrix().get_matrix()[i][j] == missing_genotype ){
                indexs.push_back(i);
                has_missing=true;
            }else{
                sum += genotype.get_genotype_matrix().get_matrix()[i][j];
                count++;
            }
            full_x.get_matrix()[i][1] = genotype.get_genotype_matrix().get_matrix()[i][j];
        }
        if(has_missing){
            mean=sum/count;
            for( int i_index : indexs ){
                full_x.get_matrix()[i_index][1] = mean;
            }
        }
        trmul(e_er.getH_sqrt_inv(), full_x, X_t);
        lsq(X_t, e_er.getY_t(), beta_est);

        for( i=0; i<n; ++i ){
            x_beta.get_array()[i] = 0;
            for( j_1=0; j_1<q; ++j_1 ){
                x_beta.get_array()[i] += X_t.get_matrix()[i][j_1]*beta_est.get_array()[j_1];
            }
        }

        mahalanobis_rss = 0.0;
        for ( i=0; i<n; ++i ){
            mahalanobis_rss += pow( e_er.getY_t().get_array()[i]-x_beta.get_array()[i], 2);
        }
        f_stat = ( h0_rss / mahalanobis_rss - 1) * p / freedome1;
        p_val = sf(f_stat, freedome1, p);
        std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t" << genotype.get_variant_Vector()[j].getPosition() << "\t" << genotype.get_variant_Vector()[j].getId();
        if( p_val == 0.0 ){
            p_val = 0.00000000000000000001;
        }
        printf("\t%10.20f\n", p_val);
        //printf("%s\t%d\t%s\t%10.20f\n", genotype.get_variant_Vector()[j].getChromosome(), genotype.get_variant_Vector()[j].getPosition(), genotype.get_variant_Vector()[j].getId(),  p_val);
    }
}


void emmax_test_multi_allic_single_test(phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u  ){

    int i, j, j_1;
    Emma_result e_er = estimate( pi, k_i, genotype );
    int covariantsNumber = 1;
    int n = genotype.get_number_of_individual();

    // claim those variables firstly could save RAM operation and save computational time begin

    double mahalanobis_rss, f_stat, p_val;
    double h0_rss = e_er.getMahalanobis_rss();


    My_Vector<double> x_beta(n);
    // claim those variables firstly could save RAM operation and save computational time end

    printf("Chr\tposition\tid\tp_val\n");
    int this_this_state, q, p, size, freedome1;
    for( j=0; j<genotype.get_number_of_variant(); ++j ) {
        if (genotype.get_variant_Vector()[j].getId()[genotype.get_variant_Vector()[j].getId().size() - 1] == 'b') {

        } else {
            std::map<int, int> allThisStates;
            for (i = 0; i < n; ++i) {
                if( allThisStates.find(genotype.get_genotype_matrix().get_matrix()[i][j]) == allThisStates.end() ){
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = 1;
                }else{
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] + 1;
                }
            }

            std::vector<int> allThisStates_vector;
            for (std::map<int, int>::iterator ittt = allThisStates.begin(); ittt != allThisStates.end(); ++ittt) {
                if( ittt->second >= man_l && ittt->second <= man_u && ittt->second!=missing_genotype){
                    allThisStates_vector.push_back(ittt->first);
                }
            }
            if( allThisStates_vector.size() > 1 ){

                size = allThisStates_vector.size();
                q = covariantsNumber + size - 1;
                p = n - q;

                My_matrix<double> full_x(n, q);
                My_matrix<double> X_t(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (i = 0; i < n; ++i) {
                    for (j_1 = 0; j_1 < (size - 1); ++j_1) { // left one allele
                        if (allThisStates_vector[j_1] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                            this_this_state = 1;
                        } else {
                            this_this_state = 0;
                        }
                        full_x.get_matrix()[i][j_1 + 1] = this_this_state;
                    }
                }
                My_Vector<double> beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, X_t);
                lsq(X_t, e_er.getY_t(), beta_est);

                for (i = 0; i < n; ++i) {
                    x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        x_beta.get_array()[i] += X_t.get_matrix()[i][j_1] * beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0.0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                }
                freedome1 = q - covariantsNumber;
                f_stat = (h0_rss / mahalanobis_rss - 1) * p / (freedome1);
                p_val = sf(f_stat, freedome1, p);
                std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                          << genotype.get_variant_Vector()[j].getPosition() << "\t"
                          << genotype.get_variant_Vector()[j].getId();
                if (p_val == 0.0) {
                    p_val = 0.00000000000000000001;
                }
                printf("\t%10.20f\n", p_val); // it is good to use \n than std::endl. Since std::endl output the result immediately and less efficient
                //printf("%s\t%d\t%s\t%10.20f\n", genotype.get_variant_Vector()[j].getChromosome(), genotype.get_variant_Vector()[j].getPosition(), genotype.get_variant_Vector()[j].getId(),  p_val);
            }else if ( allThisStates_vector.size() == 1 && allThisStates_vector[0]!=missing_genotype ){
                size = 1;
                q = covariantsNumber+ size;
                p = n - q;

                My_matrix<double> full_x(n, q);
                My_matrix<double> X_t(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (i = 0; i < n; ++i) { // here we take missing value as a state
                    if (allThisStates_vector[0] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                        this_this_state = 1;
                    } else {
                        this_this_state = 0;
                    }
                    full_x.get_matrix()[i][1] = this_this_state;
                }
                My_Vector<double> beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, X_t);
                lsq(X_t, e_er.getY_t(), beta_est);

                for (i = 0; i < n; ++i) {
                    x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        x_beta.get_array()[i] += X_t.get_matrix()[i][j_1] * beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0.0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                }
                freedome1 = 1;
                f_stat = (h0_rss / mahalanobis_rss - 1) * p / (freedome1);
                p_val = sf(f_stat, freedome1, p);
                std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                          << genotype.get_variant_Vector()[j].getPosition() << "\t"
                          << genotype.get_variant_Vector()[j].getId();
                if (p_val == 0.0) {
                    p_val = 0.00000000000000000001;
                }
                printf("\t%10.20f\n", p_val);
            }
        }
    }
}


void summaryLocus( Genotype & genotype, const int & j, std::map<int, int> & allThisStates,
                   std::vector<int> & allThisStates_vector, const int & n, const double & man_l, const double & man_u){
    for (int i = 0; i < n; ++i) {
        if( allThisStates.find(genotype.get_genotype_matrix().get_matrix()[i][j]) == allThisStates.end() ){
            allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = 1;
        }else{
            allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] + 1;
        }
    }
    for (std::map<int, int>::iterator ittt = allThisStates.begin(); ittt != allThisStates.end(); ++ittt) {
        if( ittt->second >= man_l && ittt->second <= man_u && ittt->first!=missing_genotype){
            allThisStates_vector.push_back(ittt->first);
        }
    }
}


void emmax_test_multi_allic_multi_test_null_model( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){

    int i, j, j_1, j_2;
    Emma_result e_er = estimate( pi, k_i, genotype );
    int covariantsNumber = 1;

    // claim those variables firstly could save RAM operation and save computational time begin
    int n = genotype.get_number_of_individual();
    int q = 2;
    int p = n - q;
    int freedome1 = 1;
    double mahalanobis_rss, f_stat, p_val;
    double h0_rss = e_er.getMahalanobis_rss();
    My_matrix<double> X_t(n, 2);
    My_matrix<double> full_x(n, 2);
    for( i=0; i< n; ++i ){
        full_x.get_matrix()[i][0]=1;
    }
    My_Vector<double> beta_est(2);
    My_Vector<double> x_beta(n);
    // claim those variables firstly could save RAM operation and save computational time end

    printf("Chr\tposition\tid\tp_val\n");
    int this_this_state, size;
    for( j=0; j<genotype.get_number_of_variant(); ++j ) {
        if( genotype.get_variant_Vector()[j].getId()[genotype.get_variant_Vector()[j].getId().size()-1] == 'b' ) {

        }else {
            std::map<int, int> allThisStates;
            std::vector<int> allThisStates_vector;
            summaryLocus( genotype, j, allThisStates, allThisStates_vector, n, man_l, man_u);

            if( allThisStates_vector.size() > 0 ){
                size = allThisStates_vector.size();
                if( allThisStates_vector.size() == allThisStates.size() && allThisStates_vector.size()==2 ){ // there are two and only two states, it is enough to test only one of them
                    --size;
                }
                for (j_1 = 0; j_1 < size; ++j_1) {
                    if (allThisStates_vector[j_1] != missing_genotype) {
                        for (i = 0; i < n; ++i) {
                            if (allThisStates_vector[j_1] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                                this_this_state = 0;
                            } else {
                                this_this_state = 1;
                            }
                            full_x.get_matrix()[i][1] = this_this_state;
                        }
                        trmul(e_er.getH_sqrt_inv(), full_x, X_t);
                        lsq(X_t, e_er.getY_t(), beta_est);

                        for (i = 0; i < n; ++i) {
                            x_beta.get_array()[i] = 0;
                            for (j_2 = 0; j_2 < q; ++j_2) {
                                x_beta.get_array()[i] += X_t.get_matrix()[i][j_2] * beta_est.get_array()[j_2];
                            }
                        }
                        mahalanobis_rss = 0.0;
                        for (i = 0; i < n; ++i) {
                            mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                        }
                        f_stat = (h0_rss / mahalanobis_rss - 1) * p / freedome1;
                        p_val = sf(f_stat, freedome1, p);
                        std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                                  << genotype.get_variant_Vector()[j].getPosition() << "\t"
                                  << genotype.get_variant_Vector()[j].getId();
                        if (p_val == 0.0) {
                            p_val = 0.00000000000000000001;
                        }
                        printf("\t%10.20f\n", p_val);
                    }
                }
            }
        }
    }
}


void emmax_test_multi_allic_multi_test_full_model( phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){

    int i, j, j_1, j_2;
    Emma_result e_er = estimate( pi, k_i, genotype );
    int covariantsNumber = 1;

    // claim those variables firstly could save RAM operation and save computational time begin
    int n = genotype.get_number_of_individual();
    int freedome1 = 1, x_index, q, p;
    double mahalanobis_rss, h0_rss, f_stat, p_val;
    My_Vector<double> full_x_beta(n);
    My_Vector<double> x_beta(n);
    // claim those variables firstly could save RAM operation and save computational time end
    printf("Chr\tposition\tid\tp_val\n");
    for( j=0; j < genotype.get_number_of_variant(); ++j ) {
        if( genotype.get_variant_Vector()[j].getId()[genotype.get_variant_Vector()[j].getId().size()-1] == 'b' ) {

        }else{
            std::map<int, int> allThisStates;
            std::vector<int> allThisStates_vector;
            summaryLocus( genotype, j, allThisStates, allThisStates_vector, n, man_l, man_u);

            if (allThisStates_vector.size() > 1) {

                q = covariantsNumber + allThisStates_vector.size() - 1;
                p = n - q;
                My_matrix<double> full_X_t(n, q); // there is and intercept, the last state is not included in the full model
                My_matrix<double> full_x(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (j_1 = 0; j_1 < (allThisStates_vector.size() - 1); ++j_1) {
                    for (i = 0; i < n; ++i) {
                        if (allThisStates_vector[j_1] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                            full_x.get_matrix()[i][j_1 + covariantsNumber] = 1;
                        } else {
                            full_x.get_matrix()[i][j_1 + covariantsNumber] = 0;
                        }
                    }
                }
                My_Vector<double> full_beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, full_X_t);
                lsq(full_X_t, e_er.getY_t(), full_beta_est);
                for (i = 0; i < n; ++i) {
                    full_x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        full_x_beta.get_array()[i] += full_X_t.get_matrix()[i][j_1] * full_beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - full_x_beta.get_array()[i], 2);
                } // this is the full model

                My_Vector<double> beta_est(q - 1);
                My_matrix<double> X_t(n, q - 1);
                My_matrix<double> x_1(n, q - 1);
                for (i = 0; i < n; ++i) {
                    x_1.get_matrix()[i][0] = 1;
                }
                for (j_1 = 0; j_1 < (allThisStates_vector.size() - 1); ++j_1) {
                    x_index = covariantsNumber;
                    for (j_2 = 0; j_2 < (allThisStates_vector.size() - 1); ++j_2) {
                        if (j_1 != j_2) { // exclude j_1 each time
                            for (i = 0; i < n; ++i) {
                                if (allThisStates_vector[j_2] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                                    x_1.get_matrix()[i][x_index] = 1;
                                } else {
                                    x_1.get_matrix()[i][x_index] = 0;
                                }
                            }
                            ++x_index;
                        }
                    }
                    trmul(e_er.getH_sqrt_inv(), x_1, X_t);
                    lsq(X_t, e_er.getY_t(), beta_est);
                    for (i = 0; i < n; ++i) {
                        x_beta.get_array()[i] = 0;
                        for (j_2 = 0; j_2 < (q - 1); ++j_2) {
                            x_beta.get_array()[i] += X_t.get_matrix()[i][j_2] * beta_est.get_array()[j_2];
                        }
                    }
                    h0_rss = 0.0;
                    for (i = 0; i < n; ++i) {
                        h0_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                    }
                    f_stat = (h0_rss / mahalanobis_rss - 1) * p / freedome1;
                    p_val = sf(f_stat, freedome1, p);
                    std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                              << genotype.get_variant_Vector()[j].getPosition() << "\t"
                              << genotype.get_variant_Vector()[j].getId();
                    if (p_val == 0.0) {
                        p_val = 0.00000000000000000001;
                    }
                    printf("\t%10.20f\n", p_val);
                }
            }else if ( allThisStates_vector.size() == 1 ){
                q = covariantsNumber + 1;
                p = n - q;
                My_matrix<double> full_X_t(n, q);
                My_matrix<double> full_x(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }

                for (i = 0; i < n; ++i) {
                    if (allThisStates_vector[0] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                        full_x.get_matrix()[i][1] = 0;
                    } else {
                        full_x.get_matrix()[i][1] = 1;
                    }
                }
                My_Vector<double> full_beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, full_X_t);
                lsq(full_X_t, e_er.getY_t(), full_beta_est);
                for (i = 0; i < n; ++i) {
                    full_x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        full_x_beta.get_array()[i] += full_X_t.get_matrix()[i][j_1] * full_beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - full_x_beta.get_array()[i], 2);
                }

                h0_rss = e_er.getMahalanobis_rss();

                f_stat = (h0_rss / mahalanobis_rss - 1) * p / freedome1;
                p_val = sf(f_stat, freedome1, p);
                std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                          << genotype.get_variant_Vector()[j].getPosition() << "\t"
                          << genotype.get_variant_Vector()[j].getId();
                if (p_val == 0.0) {
                    p_val = 0.00000000000000000001;
                }
                printf("\t%10.20f\n", p_val);
            }
        }
    }
}

void emma_test ( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf){
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

    emma_test ( pi, k_i, genotype, man_l, man_u );
}
void emmax_test( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf ){
    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype  genotype =  Genotype(number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;
    emmax_test ( pi, k_i, genotype, man_l, man_u );
}
void emmax_test_multi_allic_single_test( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf ){
    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype = Genotype(number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;
    emmax_test_multi_allic_single_test ( pi, k_i, genotype, man_l, man_u  );
}
void emmax_test_multi_allic_multi_test_null_model( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf ){
    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype =  Genotype(number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;
    emmax_test_multi_allic_multi_test_null_model ( pi, k_i, genotype, man_l, man_u  );
}
void emmax_test_multi_allic_multi_test_full_model( const std::string & phenotype_path, const std::string & genotype_path, const std::string & kinship_file, const double & maf ){
    std::string format = "tfam";
    std::string tfamFile = genotype_path + ".tfam";
    std::string tpedFile = genotype_path + ".tped";

    phenotype_impl pi(phenotype_path, format);
    Kinship_matrix_impl k_i(kinship_file);
    uint64_t number_of_individuals = getFileLineNumber ( tfamFile );
    uint64_t number_of_variants  = getFileLineNumber ( tpedFile );
    Genotype genotype = Genotype(number_of_individuals, number_of_variants);
    Read_tped_file(tfamFile, tpedFile, genotype);
    genotype.onlyKeepThoseIndividuls(pi.getIndividual_ids());

    double man_l = maf * pi.getIndividual_ids().size();
    double man_u = (double)pi.getIndividual_ids().size() - man_l;
    emmax_test_multi_allic_multi_test_full_model ( pi, k_i, genotype,man_l, man_u );
}



// if I could develop a model that take fixed model as a special case of mixed model
// then we could use multi allic variants as random variables and get the p-values
void emmax_test_multi_allic_single_random_test(phenotype_impl & pi, Kinship_matrix_impl & k_i, Genotype & genotype, const double & man_l, const double & man_u ){

    int i, j, j_1;
    Emma_result e_er = estimate( pi, k_i, genotype );
    int covariantsNumber = 1;
    int n = genotype.get_number_of_individual();

    // claim those variables firstly could save RAM operation and save computational time begin

    double mahalanobis_rss, f_stat, p_val;
    double h0_rss = e_er.getMahalanobis_rss();


    My_Vector<double> x_beta(n);
    // claim those variables firstly could save RAM operation and save computational time end

    printf("Chr\tposition\tid\tp_val\n");
    int this_this_state, q, p, size, freedome1;
    for( j=0; j<genotype.get_number_of_variant(); ++j ) {
        if (genotype.get_variant_Vector()[j].getId()[genotype.get_variant_Vector()[j].getId().size() - 1] == 'b') {

        } else {
            std::map<int, int> allThisStates;
            for (i = 0; i < n; ++i) {
                if( allThisStates.find(genotype.get_genotype_matrix().get_matrix()[i][j]) == allThisStates.end() ){
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = 1;
                }else{
                    allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] = allThisStates[genotype.get_genotype_matrix().get_matrix()[i][j]] + 1;
                }
            }

            std::vector<int> allThisStates_vector;
            for (std::map<int, int>::iterator ittt = allThisStates.begin(); ittt != allThisStates.end(); ++ittt) {
                if( ittt->second >= man_l && ittt->second <= man_u && ittt->second!=missing_genotype){
                    allThisStates_vector.push_back(ittt->first);
                }
            }
            if( allThisStates_vector.size() > 1 ){

                size = allThisStates_vector.size();
                q = covariantsNumber + size - 1;
                p = n - q;

                My_matrix<double> full_x(n, q);
                My_matrix<double> X_t(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (i = 0; i < n; ++i) {
                    for (j_1 = 0; j_1 < (size - 1); ++j_1) { // left one allele
                        if (allThisStates_vector[j_1] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                            this_this_state = 1;
                        } else {
                            this_this_state = 0;
                        }
                        full_x.get_matrix()[i][j_1 + 1] = this_this_state;
                    }
                }
                My_Vector<double> beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, X_t);
                lsq(X_t, e_er.getY_t(), beta_est);

                for (i = 0; i < n; ++i) {
                    x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        x_beta.get_array()[i] += X_t.get_matrix()[i][j_1] * beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0.0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                }
                freedome1 = q - covariantsNumber;
                f_stat = (h0_rss / mahalanobis_rss - 1) * p / (freedome1);
                p_val = sf(f_stat, freedome1, p);
                std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                          << genotype.get_variant_Vector()[j].getPosition() << "\t"
                          << genotype.get_variant_Vector()[j].getId();
                if (p_val == 0.0) {
                    p_val = 0.00000000000000000001;
                }
                printf("\t%10.20f\n", p_val);
                //printf("%s\t%d\t%s\t%10.20f\n", genotype.get_variant_Vector()[j].getChromosome(), genotype.get_variant_Vector()[j].getPosition(), genotype.get_variant_Vector()[j].getId(),  p_val);
            }else if ( allThisStates_vector.size() == 1 && allThisStates_vector[0]!=missing_genotype ){
                size = 1;
                q = covariantsNumber+ size;
                p = n - q;

                My_matrix<double> full_x(n, q);
                My_matrix<double> X_t(n, q);
                for (i = 0; i < n; ++i) {
                    full_x.get_matrix()[i][0] = 1;
                }
                for (i = 0; i < n; ++i) { // here we take missing value as a state
                    if (allThisStates_vector[0] == genotype.get_genotype_matrix().get_matrix()[i][j]) {
                        this_this_state = 1;
                    } else {
                        this_this_state = 0;
                    }
                    full_x.get_matrix()[i][1] = this_this_state;
                }
                My_Vector<double> beta_est(q);
                trmul(e_er.getH_sqrt_inv(), full_x, X_t);
                lsq(X_t, e_er.getY_t(), beta_est);

                for (i = 0; i < n; ++i) {
                    x_beta.get_array()[i] = 0;
                    for (j_1 = 0; j_1 < q; ++j_1) {
                        x_beta.get_array()[i] += X_t.get_matrix()[i][j_1] * beta_est.get_array()[j_1];
                    }
                }
                mahalanobis_rss = 0.0;
                for (i = 0; i < n; ++i) {
                    mahalanobis_rss += pow(e_er.getY_t().get_array()[i] - x_beta.get_array()[i], 2);
                }
                freedome1 = 1;
                f_stat = (h0_rss / mahalanobis_rss - 1) * p / (freedome1);
                p_val = sf(f_stat, freedome1, p);
                std::cout << genotype.get_variant_Vector()[j].getChromosome() << "\t"
                          << genotype.get_variant_Vector()[j].getPosition() << "\t"
                          << genotype.get_variant_Vector()[j].getId();
                if (p_val == 0.0) {
                    p_val = 0.00000000000000000001;
                }
                printf("\t%10.20f\n", p_val);
            }
        }
    }
}
