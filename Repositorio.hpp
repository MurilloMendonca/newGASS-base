#pragma once
#include "Atomo.hpp"
#include "AtomoCompat.hpp"
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <functional>
#include <fstream>
#include <iostream>

struct Repositorio{
    std::vector<Atomo> atoms;
    std::vector<std::vector<std::vector<Atomo*>>> rep; // [cluster][amino][atomo]
    std::string pdbId;
    float centroide[3];

    std::vector<std::vector<Atomo*>> operator [](int key){
        return rep[key];
    }
    
    auto getRepository(){
        const std::map<std::string, int> AMINO_CODES_MAP = {
            {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4}, {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9}, {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
        };
        std::vector<std::vector<std::vector<Atomo*>>> aux (5);
        for(int i = 0; i < 5; i++){
            aux[i].resize(20);
        }
        for(auto& atomo : this->atoms){
            if(AMINO_CODES_MAP.find(atomo.amino) != AMINO_CODES_MAP.end()){
            aux[atomo.cluster][AMINO_CODES_MAP.at(atomo.amino)].push_back(&atomo);
            aux[0][AMINO_CODES_MAP.at(atomo.amino)].push_back(&atomo);
            }
        }
        return aux;
    }
    void readRepository(std::string path){
        
        AtomoCompat auxCompat;
        float centroideMax[3] = {-1000,-1000,-1000};
        float centroideMin[3] = {1000,1000,1000};
        std::ifstream inFile(path,std::ios::binary);
        if(!inFile){
            throw std::runtime_error("Protein file can not be open: " + path);
        }
        int contador = 0;
        char linha[100];
        while(inFile){
            if(contador<=3){
                inFile.read((char*)&linha, sizeof(linha));
                contador++;
            }
            else{
                inFile.read((char *)&auxCompat, sizeof(AtomoCompat));
                Atomo aux(auxCompat);
                centroideMax[0] = std::max(centroideMax[0], aux.x);
                centroideMax[1] = std::max(centroideMax[1], aux.y);
                centroideMax[2] = std::max(centroideMax[2], aux.z);
                centroideMin[0] = std::min(centroideMin[0], aux.x);
                centroideMin[1] = std::min(centroideMin[1], aux.y);
                centroideMin[2] = std::min(centroideMin[2], aux.z);
                this->atoms.push_back(aux);
            }
        }
        this->centroide[0] = (centroideMax[0] + centroideMin[0])/2;
        this->centroide[1] = (centroideMax[1] + centroideMin[1])/2;
        this->centroide[2] = (centroideMax[2] + centroideMin[2])/2;
        for(auto& amino : this->atoms){
            amino.cluster = getPartition(this->centroide, amino);
        }
        inFile.close();
        rep = this->getRepository();
        
    }

    void readRepository(std::string path, std::function<int(std::vector<Atomo>, Atomo)> clusterFunction){
        AtomoCompat auxCompat;
        float centroideMax[3] = {-1000,-1000,-1000};
        float centroideMin[3] = {1000,1000,1000};
        std::ifstream inFile(path,std::ios::binary);
        if(!inFile){
            throw std::runtime_error("Protein file can not be open: " + path);
        }
        int contador = 0;
        char linha[100];
        while(inFile){
            if(contador<=3){
                inFile.read((char*)&linha, sizeof(linha));
                contador++;
            }
            else{
                inFile.read((char *)&auxCompat, sizeof(AtomoCompat));
                Atomo aux(auxCompat);
                centroideMax[0] = std::max(centroideMax[0], aux.x);
                centroideMax[1] = std::max(centroideMax[1], aux.y);
                centroideMax[2] = std::max(centroideMax[2], aux.z);
                centroideMin[0] = std::min(centroideMin[0], aux.x);
                centroideMin[1] = std::min(centroideMin[1], aux.y);
                centroideMin[2] = std::min(centroideMin[2], aux.z);
                this->atoms.push_back(aux);
            }
        }
        this->centroide[0] = (centroideMax[0] + centroideMin[0])/2;
        this->centroide[1] = (centroideMax[1] + centroideMin[1])/2;
        this->centroide[2] = (centroideMax[2] + centroideMin[2])/2;
        for(auto& amino : this->atoms){
            amino.cluster = clusterFunction(this->atoms, amino);
        }
        inFile.close();
        rep = this->getRepository();

    }
    

    int getPartition(float centroide[3], Atomo res){
        if(res.x > centroide[0] && res.y > centroide[1]){ //primeiro quadrante
            return 1;
        }
        else if(res.x > centroide[0] && res.y < centroide[1]){
            return 2;
        }
        else if(res.x < centroide[0] && res.y < centroide[1]){
            return 3;
        }
        else{
            return 4;
        }   
    }

    ~Repositorio(){
        this->rep.clear();
        this->atoms.clear();
    }
};