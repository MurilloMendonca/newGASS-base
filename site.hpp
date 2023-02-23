#pragma once
#include "Atomo.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct site
{
    std::vector<Atomo> residuos;
    std::vector<float> distances;
    std::string pdbId;
    std::string uniprotId;
    std::string resolution;
    std::string ecNumber;
    std::vector<std::vector<int>> substitutions;
    int size;
    site(){
        this->size=0;
        substitutions = std::vector<std::vector<int>>(20);
    };
    site(std::vector<Atomo> residuos){
        this->residuos=residuos;
        this->size=residuos.size();
        calculateDistances();
    }
    void calculateDistances(){
        distances = std::vector<float>(size*(size-1)/2);
        int pos=0;
        for(int i=0;i<size;i++){
            for(int j=i+1;j<size;j++){
                distances[pos++]=sqrt((residuos[i].x - residuos[j].x)*(residuos[i].x - residuos[j].x) +
                    (residuos[i].y - residuos[j].y)*(residuos[i].y - residuos[j].y)  +  
                    (residuos[i].z - residuos[j].z)*(residuos[i].z - residuos[j].z));
            }
        }
        std::sort(distances.begin(), distances.end());

    }
};

