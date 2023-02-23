#pragma once

#include "Atomo.hpp"
#include "site.hpp"
#include "Repositorio.hpp"
#include "GA.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <thread>
#include <chrono>
#include <queue>
#include <mutex>


namespace GASS{
    struct Parameters{
        int AG_POPULATION_SIZE;
        int AG_NUMBER_OF_GENERATIONS;
        int AG_NUMBER_OF_ELITE_INDIVIDUALS;
        int AG_TOURNAMENT_SIZE;
        float AG_MUTATION_RATE;
        float AG_CROSSOVER_RATE;
        int GASS_NUMBER_OF_TEMPLATES;
    };

// GASS FUNCTIONS
void readConfigFile(std::string configFileName);


void readTemplateFile(std::string templateFileName, site &ref_site);
void readSetupTemplateFile(std::string setupTemplateFileName, std::vector<site>& templates);

auto getSitesByPdbId(std::string pdbId, std::vector<site>& sites);
void readSubstitutuionMatrix(std::string fileName, std::vector<site>& sites);

void run(Repositorio& repositorio, std::vector<site>& templates, std::string outputName);
void setup(std::vector<site>& templates, std::vector<Repositorio>& repositorios, std::string configFileName, std::string templatesFolderName,std::string runListFileName);

void run(std::vector<site>& templates, Repositorio& repositorio);
void run(std::vector<site>& templates, Repositorio& repositorio, std::string outputName);
void run(std::vector<site>& templates, Repositorio& repositorio, std::string outputName, Parameters parameters);
void readRunFile(std::string fileName, std::vector<std::string>& proteinNames);
void runOneToOne(site* temp, Repositorio* repositorio, std::set<Individuo>* results);
void runSanityTest(Parameters param, std::set<Individuo>* results);
}