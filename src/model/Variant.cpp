//
// Created by baoxing on 2/17/18.
//

#include "Variant.h"
Variant::Variant(const std::string& _chromosome, const std::string& _id,
                 const float & _geneticDistance, const unsigned long & _position){
    this->id = _id;
    this->chromosome=_chromosome;
    this->genetic_distance=_geneticDistance;
    this->position=_position;
}

const std::string &Variant::getId() const{
    return this->id;
}

void Variant::setId(const std::string &id) {
    this->id = id;
}

const std::string & Variant::getChromosome() const{
    return this->chromosome;
}

void Variant::setChromosome(const std::string &chromosome) {
    this->chromosome = chromosome;
}

const float Variant::getGenetic_distance() const {
    return this->genetic_distance;
}

void Variant::setGenetic_distance(float genetic_distance) {
    this->genetic_distance = genetic_distance;
}

const unsigned long Variant::getPosition() const{
    return this->position;
}

void Variant::setPosition(unsigned long position) {
    this->position = position;
}
