#include "newGASS.hpp"
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

//---> GLOBAL VARIABLES <---//
namespace GASS{

    
    


    
    // AG PARAMETERS
    int AG_POPULATION_SIZE;
    int AG_NUMBER_OF_GENERATIONS;
    int AG_NUMBER_OF_ELITE_INDIVIDUALS;
    int AG_TOURNAMENT_SIZE;
    float AG_MUTATION_RATE;
    float AG_CROSSOVER_RATE;

    // GASS PARAMETERS
    int GASS_NUMBER_OF_TEMPLATES;


    // GASS FUNCTIONS
    void readConfigFile(std::string configFileName)
    {
        std::ifstream configFile(configFileName);
        char buff[255];
        int nlinha = 0;
        if (configFile.is_open())
        {
            while (!configFile.eof())
            {
                configFile.getline(buff, 255);
                switch (nlinha)
                {
                case 1:
                    AG_NUMBER_OF_GENERATIONS = atoi(buff);
                    break;
                case 3:
                    AG_POPULATION_SIZE = atoi(buff);
                    break;
                case 5:
                    AG_TOURNAMENT_SIZE = atoi(buff);
                    break;
                case 7:
                    AG_CROSSOVER_RATE = atoi(buff)/100.0;
                    break;
                case 9:
                    AG_MUTATION_RATE = atoi(buff)/100.0;
                    break;
                case 11:
                    AG_NUMBER_OF_ELITE_INDIVIDUALS = atoi(buff);
                    break;
                case 13:
                    GASS_NUMBER_OF_TEMPLATES = atoi(buff);
                    break;
                }
                nlinha++;
            }
            configFile.close();
        }
        else
        {
            std::cout << "Error opening config file" << std::endl;
            throw std::runtime_error("Configuration file not found");
        }
    }
    void readConfigFile(std::string configFileName, Parameters &params)
    {
        std::ifstream configFile(configFileName);
        char buff[255];
        int nlinha = 0;
        if (configFile.is_open())
        {
            while (!configFile.eof())
            {
                configFile.getline(buff, 255);
                switch (nlinha)
                {
                case 1:
                    params.AG_NUMBER_OF_GENERATIONS = atoi(buff);
                    break;
                case 3:
                    params.AG_POPULATION_SIZE = atoi(buff);
                    break;
                case 5:
                    params.AG_TOURNAMENT_SIZE = atoi(buff);
                    break;
                case 7:
                    params.AG_CROSSOVER_RATE = atoi(buff)/100.0;
                    break;
                case 9:
                    params.AG_MUTATION_RATE = atoi(buff)/100.0;
                    break;
                case 11:
                    params.AG_NUMBER_OF_ELITE_INDIVIDUALS = atoi(buff);
                    break;
                case 13:
                    params.GASS_NUMBER_OF_TEMPLATES = atoi(buff);
                    break;
                }
                nlinha++;
            }
            configFile.close();
        }
        else
        {
            std::cout << "Error opening config file" << std::endl;
            throw std::runtime_error("Configuration file not found");
        }
    }
    void readTemplateFile(std::string templateFileName, site &ref_site)
    {
        std::ifstream templateFile;
        templateFile.open(templateFileName, std::ios::binary);
        char buff[100];
        int nlinha = 0;
        AtomoCompat aux;
        std::cout.flush();
        if (templateFile.is_open())
        {
            int i = 0;
            int contador = 0;
            while ((!templateFile.eof()) && (i < ref_site.size))
            {
                if (contador <= 3)
                {
    
                    templateFile.read((char *)&buff, sizeof(buff));
                    if (contador == 0)
                        ref_site.pdbId = buff;
                    else if (contador == 1)
                        ref_site.ecNumber = buff;
                    else if (contador == 2)
                        ref_site.uniprotId = buff;
                    else if (contador == 3)
                        ref_site.resolution = buff;
                    contador++;
                }
                else
                {
                    templateFile.read((char *)&aux, sizeof(AtomoCompat));
                    for (int j = 0; j < ref_site.size; j++)
                        if ((aux.atomo_ID == ref_site.residuos[j].atomo_ID) && (aux.cadeia == ref_site.residuos[j].cadeia))
                        {
                            ref_site.residuos[j].x = aux.x;
                            ref_site.residuos[j].y = aux.y;
                            ref_site.residuos[j].z = aux.z;
                            i++;
                        }
                }
                
            }
            templateFile.close();
            ref_site.calculateDistances();
        }
        else
        {
            std::cout << "Error opening template file" << std::endl;
            throw std::runtime_error("Template file not found");
        }
    }

    void readSetupTemplateFile(std::string setupTemplateFileName, std::vector<site>& templates){
        std::ifstream fin;
        char buff[255];
        int nlinha = 0;
        int linhaarquivoref = 0;

        fin.open(setupTemplateFileName);
        if(fin.is_open()){
            for(int i = 0; i < GASS_NUMBER_OF_TEMPLATES; i++){
                nlinha = 0;
                templates.emplace_back();
                site &ref_site = templates.back();
                while(nlinha <= 4){
                    fin.getline(buff, 25);
                    switch (nlinha) {
                        case 1: ref_site.pdbId=buff; break;
                        case 3: ref_site.size=std::stoi(buff);break;
                    }
                    nlinha++;
                }
                ref_site.residuos.resize(ref_site.size);
                for(Atomo& atomo : ref_site.residuos){
                    fin.getline(buff, 25);
                    atomo.amino=buff;
                    fin.getline(buff, 25);
                    atomo.atomo_ID=std::stoi(buff);
                    fin.getline(buff, 25);
                    atomo.cadeia=buff[0];
                    nlinha+=3;
                }
            }
            fin.close();
        }
        
        else{
            std::cout << "File templates.txt can not be open.\n";
            throw std::runtime_error("File templates.txt can not be open.");
        }
    }


    auto getSitesByPdbId(std::string pdbId, std::vector<site>& sites){
        auto it = std::find_if(sites.begin(), sites.end(), [&pdbId](const site& s){
            return s.pdbId == pdbId;
        });
        if(it == sites.end()){
            return sites.end();
        }
        return it;
    }
    void readSubstitutuionMatrix(std::string fileName, std::vector<site>& sites){
        const std::map<std::string, int> AMINO_CODES_MAP = {
            {"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3}, {"CYS", 4}, {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9}, {"LEU", 10}, {"LYS", 11}, {"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}
        };
        std::ifstream arquivotexto;  // variável para controlar o fluxo de entrada

    // fin.open("Templates.txt");

        arquivotexto.open(fileName);
        int l1 = 0;
        if(!arquivotexto.is_open()){
            std::cout<<"File can not be open: SubstitutionMatrix.txt.\n";
            throw std::runtime_error("File can not be open: "+fileName);
        }
        else{
            char linha[30];
            while(!arquivotexto.eof()){
                arquivotexto.getline(linha,30);
                std::string temp="";
                int c2 = 0;
                while(linha[c2]!='\0'){
                    while(linha[c2]!=','){
                        temp += toupper(linha[c2]);
                        c2++;
                    }
                    std::string templateName=temp;
                    c2++;
                    temp="";
                    while(linha[c2]!=','){
                        temp+=linha[c2];
                        c2++;
                    }
                    std::string amino1=temp;
                    c2++;
                    temp="";
                    while(linha[c2]!='\0'){
                        temp+= linha[c2];
                        c2++;
                    }
                    std::string amino2=temp;
                    auto sitesToAddSubstitution=getSitesByPdbId(templateName, sites);
                    while(sitesToAddSubstitution != sites.end()){
                        sitesToAddSubstitution->substitutions[AMINO_CODES_MAP.at(amino1)].push_back(AMINO_CODES_MAP.at(amino2));
                        sitesToAddSubstitution++;
                    }
                }
                l1++;
            }
        }
    }

    void setup(std::vector<site>& templates, std::vector<Repositorio>& repositorios, std::string configFileName, std::string templatesFolderName,std::string runListFileName){
        readConfigFile(configFileName);
        readSetupTemplateFile(templatesFolderName+"Templates.txt", templates);
        for(site& s : templates){
            readTemplateFile(templatesFolderName+s.pdbId, s);
            std::cout <<"PDB: "<< s.pdbId << std::endl;
            std::cout <<"Template Size: "<< s.size << std::endl;
            std::cout <<"EC Number: "<< s.ecNumber << std::endl;
            std::cout <<"UniprotID: "<< s.uniprotId << std::endl;
            for(Atomo& a : s.residuos){
                std::cout << a.amino << std::endl;
                std::cout << a.atomo_ID << std::endl;
                std::cout << a.cadeia << std::endl;
                std::cout << a.x << std::endl;
                std::cout << a.y << std::endl;
                std::cout << a.z << std::endl;
            }
        }

        
        std::vector<std::string> proteinNames;
        readRunFile(runListFileName, proteinNames);
        std::cout<<"Proteinas: "<<std::endl;
        for(std::string& proteinName : proteinNames){
            std::cout<<proteinName<<std::endl;
        }
        for(std::string& proteinName : proteinNames){
            Repositorio repositorio;
            //getRepository(repositorio, "../cache/" + proteinName + "/targ_.dat");
            repositorio.readRepository("../cache/" + proteinName + "/targ_.dat");
            repositorios.push_back(repositorio);
        }

        readSubstitutuionMatrix(templatesFolderName+"SubstitutionMatrix.txt", templates);
        
    }

    void run(std::vector<site>& templates, Repositorio& repositorio){
        run(std::ref(templates), std::ref(repositorio), "output.txt");
    }
    void run(std::vector<site>& templates, Repositorio& repositorio, std::string outputName){
        std::vector<std::set<Individuo>> populacao(templates.size());
        std::vector<std::thread> threads;
        std::mutex mutex;
        for(int quad=1;quad<=4;quad++){
            threads.push_back(std::thread([quad, &populacao, &repositorio, &templates, &mutex]{
                //std::cout<<"Criando GA para o quadrante "<<quad<<std::endl;
                GA ga(AG_POPULATION_SIZE, 
                    AG_TOURNAMENT_SIZE,
                    AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    AG_MUTATION_RATE, 
                    AG_CROSSOVER_RATE, 
                    std::ref(repositorio), 
                    std::ref(templates));
                ga.run(AG_NUMBER_OF_GENERATIONS,quad);

                mutex.lock();
                for(int i=0;i<templates.size();i++){
                    populacao[i].insert(ga.melhoresIndividuos[i].begin(),ga.melhoresIndividuos[i].end());
                }
                mutex.unlock();
            }));
        }
        for(int i=0;i<threads.size();i++){
            threads[i].join();
        }
        std::vector<std::vector<Individuo>> populacao_final (templates.size());
        for(int i=0;i<templates.size();i++){
            if(populacao[i].size()>AG_POPULATION_SIZE){
                populacao[i].erase(std::next(populacao[i].begin(),AG_POPULATION_SIZE),populacao[i].end());
            }
            populacao_final[i].insert(populacao_final[i].begin(),populacao[i].begin(),populacao[i].end());
        }
        GA ga(AG_POPULATION_SIZE, 
                AG_TOURNAMENT_SIZE,
                AG_NUMBER_OF_ELITE_INDIVIDUALS,
                AG_MUTATION_RATE, 
                AG_CROSSOVER_RATE, 
                std::ref(repositorio), 
                std::ref(templates),
                outputName);
        ga.run(AG_NUMBER_OF_GENERATIONS, populacao_final);
        //std::cout << "Elapsed time running final GA: " << elapsed.count() << " s";
    }

    void run(std::vector<site>& templates, Repositorio& repositorio, std::string outputName, Parameters parameters){
        std::vector<std::set<Individuo>> populacao(templates.size());
        std::vector<std::thread> threads;
        std::mutex mutex;
        for(int quad=1;quad<=4;quad++){
            threads.push_back(std::thread([parameters, quad, &populacao, &repositorio, &templates, &mutex]{
                //std::cout<<"Criando GA para o quadrante "<<quad<<std::endl;
                GA ga(parameters.AG_POPULATION_SIZE, 
                    parameters.AG_TOURNAMENT_SIZE,
                    parameters.AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    parameters.AG_MUTATION_RATE, 
                    parameters.AG_CROSSOVER_RATE, 
                    std::ref(repositorio), 
                    std::ref(templates));
                ga.run(parameters.AG_NUMBER_OF_GENERATIONS,quad);

                mutex.lock();
                for(int i=0;i<templates.size();i++){
                    populacao[i].insert(ga.melhoresIndividuos[i].begin(),ga.melhoresIndividuos[i].end());
                }
                mutex.unlock();
            }));
        }
        for(int i=0;i<threads.size();i++){
            threads[i].join();
        }
        std::vector<std::vector<Individuo>> populacao_final (templates.size());
        for(int i=0;i<templates.size();i++){
            if(populacao[i].size()>parameters.AG_POPULATION_SIZE){
                populacao[i].erase(std::next(populacao[i].begin(),parameters.AG_POPULATION_SIZE),populacao[i].end());
            }
            populacao_final[i].insert(populacao_final[i].begin(),populacao[i].begin(),populacao[i].end());
        }
        GA ga(parameters.AG_POPULATION_SIZE, 
                parameters.AG_TOURNAMENT_SIZE,
                parameters.AG_NUMBER_OF_ELITE_INDIVIDUALS,
                parameters.AG_MUTATION_RATE, 
                parameters.AG_CROSSOVER_RATE, 
                std::ref(repositorio), 
                std::ref(templates),
                outputName);
        ga.run(parameters.AG_NUMBER_OF_GENERATIONS, populacao_final);
        //std::cout << "Elapsed time running final GA: " << elapsed.count() << " s";
    }
    
    void runOneToOne(site* temp, Repositorio* repositorio, std::set<Individuo>* results){
        std::vector<std::thread> threads;
        std::mutex mutex;
        for(int i=1;i<=4;i++){
            threads.push_back(std::thread([i, &temp, &repositorio, &results, &mutex]{
            std::set <Individuo> melhoresIndividuos;
                GA ga(AG_POPULATION_SIZE, 
                    AG_TOURNAMENT_SIZE,
                    AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    AG_MUTATION_RATE,
                    AG_CROSSOVER_RATE);
                ga.run(repositorio, temp, AG_NUMBER_OF_GENERATIONS, i, &melhoresIndividuos);
                mutex.lock();
                results->insert(melhoresIndividuos.begin(),melhoresIndividuos.end());
                mutex.unlock();
            }));
        }
        for(int i=0;i<threads.size();i++){
            threads[i].join();
        }
        std::vector<Individuo> populacao_final(results->begin(),results->end());
        if (populacao_final.size() > AG_POPULATION_SIZE){
            populacao_final.resize(AG_POPULATION_SIZE);
        }
        GA ga(AG_POPULATION_SIZE, 
                    AG_TOURNAMENT_SIZE,
                    AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    AG_MUTATION_RATE,
                    AG_CROSSOVER_RATE);
        ga.setPopulation(populacao_final);
        ga.run(repositorio, temp, AG_NUMBER_OF_GENERATIONS, 0, results);

    }

    void runOneToOne(site* temp, Repositorio* repositorio, std::set<Individuo>* results, Parameters param){
        std::vector<std::thread> threads;
        std::mutex mutex;
        for(int i=1;i<=4;i++){
            threads.push_back(std::thread([&param,i, &temp, &repositorio, &results, &mutex]{
            std::set <Individuo> melhoresIndividuos;
                GA ga(param.AG_POPULATION_SIZE, 
                    param.AG_TOURNAMENT_SIZE,
                    param.AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    param.AG_MUTATION_RATE,
                    param.AG_CROSSOVER_RATE);
                ga.run(repositorio, temp, param.AG_NUMBER_OF_GENERATIONS, i, &melhoresIndividuos);
                mutex.lock();
                results->insert(melhoresIndividuos.begin(),melhoresIndividuos.end());
                mutex.unlock();
            }));
        }
        for(int i=0;i<threads.size();i++){
            threads[i].join();
        }
        std::vector<Individuo> populacao_final(results->begin(),results->end());
        if (populacao_final.size() > param.AG_POPULATION_SIZE){
            populacao_final.resize(param.AG_POPULATION_SIZE);
        }
        GA ga(param.AG_POPULATION_SIZE, 
                    param.AG_TOURNAMENT_SIZE,
                    param.AG_NUMBER_OF_ELITE_INDIVIDUALS,
                    param.AG_MUTATION_RATE,
                    param.AG_CROSSOVER_RATE);
        ga.setPopulation(populacao_final);
        ga.run(repositorio, temp, param.AG_NUMBER_OF_GENERATIONS, 0, results);

    }

    void readRunFile(std::string fileName, std::vector<std::string>& proteinNames){
        std::ifstream file(fileName);
        std::string line;
        while(std::getline(file, line)){
            proteinNames.push_back(line);
        }
    }

    void runSanityTest(Parameters param, std::set<Individuo>* results){
        Repositorio repositorio;
        std::vector<site> temp;
        GASS_NUMBER_OF_TEMPLATES =1;
        GASS::readSetupTemplateFile("../templates/Templates_Zn/Template_3nos.txt", temp);
        GASS::readTemplateFile("../templates/Templates_Zn/3nos_.dat", temp[0]);
        repositorio.readRepository("../cache/3NOS/targ_.dat");
        runOneToOne(&temp[0], &repositorio, results, param);
    }
}