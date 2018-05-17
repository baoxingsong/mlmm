//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_INDIVIDUAL_H
#define MLMM_CPP_INDIVIDUAL_H
#include <string>

class Individual {
    private:
        std::string family_id;
        std::string individual_id;
        std::string paternal_id;
        std::string maternal_id;
        int sex;
    public:
        Individual();
        Individual(const std::string& _family_id, const std::string& _individual_id, const std::string _paternal_id,
            const std::string & _maternal_id, const int& _sex);
        const std::string& get_family_id() const;
        const std::string& get_individual_id() const;
        const std::string& get_paternal_id() const;
        const std::string& get_maternal_id() const;
        const int& get_sex() const;
        void setFamily_id(const std::string &family_id);
        void setIndividual_id(const std::string &individual_id);
        void setPaternal_id(const std::string &paternal_id);
        void setMaternal_id(const std::string &maternal_id);
        void setSex(int sex);
};
#endif //MLMM_CPP_INDIVIDUAL_H
