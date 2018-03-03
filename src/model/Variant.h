//
// Created by baoxing on 2/17/18.
//

#ifndef MLMM_CPP_Variant_H
#define MLMM_CPP_Variant_H

#include <string>
#include <iostream>

class Variant {
    private:
        std::string id;
        std::string chromosome;
        float genetic_distance;
        unsigned long position;
    public:
        Variant(const std::string& _chromosome, const std::string& _id,
                const float & _geneticDistance, const unsigned long & _position);

        const std::string &getId() const;
        void setId(const std::string &id);

        const std::string &getChromosome() const;
        void setChromosome(const std::string &chromosome);

        const float getGenetic_distance() const;
        void setGenetic_distance(float genetic_distance);

        const unsigned long getPosition() const;
        void setPosition(unsigned long position);

};

#endif //MLMM_CPP_Variant_H
