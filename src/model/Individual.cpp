//
// Created by baoxing on 2/17/18.
//

#include "Individual.h"
Individual::Individual(){
    this->family_id="NA";
    this->individual_id="NA";
    this->paternal_id="NA";
    this->maternal_id="NA";
    this->sex = 999;
}
Individual::Individual(const std::string& _family_id, const std::string& _individual_id, const std::string _paternal_id,
           const std::string & _maternal_id, const int & _sex){
    this->family_id=_family_id;
    this->individual_id=_individual_id;
    this->paternal_id=_paternal_id;
    this->maternal_id=_maternal_id;
    this->sex = _sex;
}
const std::string& Individual::get_family_id() const{
    return this->family_id;
}
const std::string& Individual::get_individual_id() const{
    return this->individual_id;
}
const std::string& Individual::get_paternal_id() const{
    return this->paternal_id;
}
const std::string& Individual::get_maternal_id() const{
    return this->maternal_id;
}
const int& Individual::get_sex() const{
    return this->sex;
}
