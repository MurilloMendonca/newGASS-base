//---> STRUCTURE ATOMO - GENE OF THE INDIVIDUAL <---//
#pragma once
#include <string>
#include "AtomoCompat.hpp"

struct Atomo
{
	std::string amino;
	std::string atomo;
    int amino_ID;
	char cadeia;
	int atomo_ID;
	float x, y, z;
    int cluster;

    Atomo()
    {
        amino = {};
        atomo = {};
        cadeia = '\0';
        atomo_ID = 0;
        x = 0.0;
        y = 0.0;
        z = 0.0;
        cluster=0;

    }

    Atomo(AtomoCompat aux)
    {
        amino = aux.amino;
        atomo = aux.atomo;
        cadeia = aux.cadeia;
        atomo_ID = aux.atomo_ID;
        x = aux.x;
        y = aux.y;
        z = aux.z;
        cluster=0;
    }


    bool operator==(const Atomo& other) const
    {
        return (amino == other.amino) &&
               (atomo == other.atomo) &&
               (cadeia == other.cadeia) &&
               (atomo_ID == other.atomo_ID) &&
               (x == other.x) &&
               (y == other.y) &&
               (z == other.z);
    }

    bool operator!=(const Atomo& other) const
    {
        return !(*this == other);
    }
};