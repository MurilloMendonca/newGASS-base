#include "GA.hpp"
#include "Individuo.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <set>
template <typename T>
auto GA::melhor(T vec){
    return *(std::min_element(vec.begin(), vec.end()));
}
GA::GA(int tamanhoDaPopulacao, 
        int torneioTam,
        int elitismo,
        double mut, 
        double cruzamento, 
        Repositorio& repositorio, 
        std::vector<site>& templates){
    rd = new std::random_device;
    gen = std::mt19937((*rd)());
    this->taxaMutacao=mut;
    this->taxaCruzamento=cruzamento;
    this->tamanhoPopulacao=tamanhoDaPopulacao;
    this->tamanhoTorneio=torneioTam;
    this->repositorio=&repositorio;
    this->templates=&templates;
    this->quantidadeElitismo=elitismo;
    this->melhoresIndividuos=std::vector<std::set<Individuo>>(templates.size());
    this->outputFilePath = "";
}

GA::GA(int tamanhoDaPopulacao, 
        int torneioTam,
        int elitismo,
        double mut, 
        double cruzamento){
    rd = new std::random_device;
    gen = std::mt19937((*rd)());
    this->taxaMutacao=mut;
    this->taxaCruzamento=cruzamento;
    this->tamanhoPopulacao=tamanhoDaPopulacao;
    this->tamanhoTorneio=torneioTam;
    this->quantidadeElitismo=elitismo;
}



GA::GA(int tamanhoDaPopulacao, 
        int torneioTam,
        int elitismo,
        double mut, 
        double cruzamento, 
        Repositorio& repositorio, 
        std::vector<site>& templates,
        std::string outputFilePath){
    rd = new std::random_device;
    gen = std::mt19937((*rd)());
    this->taxaMutacao=mut;
    this->taxaCruzamento=cruzamento;
    this->tamanhoPopulacao=tamanhoDaPopulacao;
    this->tamanhoTorneio=torneioTam;
    this->repositorio=&repositorio;
    this->templates=&templates;
    this->quantidadeElitismo=elitismo;
    this->melhoresIndividuos=std::vector<std::set<Individuo>>(templates.size());
    this->outputFilePath=outputFilePath;
}

void GA::mostraPopulacao(){
    int x=0;
    for(Individuo ind:populacao){
        std::cout<<"\nIndividuo "<<x++<<": ";
        for(Atomo atom:ind.cromossomo)
            std::cout<<" "<<atom.amino<<" "<<atom.atomo_ID<<" "<<atom.cadeia<<";";
        std::cout<<" fitness="<<ind.getFitness();
    }
}

bool GA::valido(Individuo i){
    std::vector<Atomo> cromossomo = i.cromossomo;
    for(int i=0;i<cromossomo.size();i++){
        if(cromossomo[i].amino!=temp.residuos[i].amino)
            return false;
    }
    return true;
}

Individuo GA::geraIndividuoValido(int quadrante){
    std::vector<Atomo> cromossomo (temp.residuos.size()) ;
    int i=0;
    for(Atomo atom:temp.residuos){
        std::uniform_int_distribution<int> distribution(0,(*repositorio)[quadrante][AMINO_CODES_MAP.at(atom.amino)].size()-1);
        Atomo a = *(*repositorio)[quadrante][AMINO_CODES_MAP.at(atom.amino)][distribution(gen)];
        cromossomo[i++]=a;
    }
    //std::random_shuffle(cromossomo.begin(),cromossomo.end());
    return Individuo(cromossomo);
}

void GA::run(int geracoes){
    int quadrante=0;
    int contTemp=0;
    double setupTime=0;
    double eliteTime=0;
    double crossTime=0;
    double mutTime=0;
    double trocaTime=0;
    double fitTime=0;
    double topTime=0;
    double saveTime=0;
    double totalTime=0;
    double repTime=0;
    for(site temp:*templates){
        
        auto start = std::chrono::high_resolution_clock::now();
        this->temp=temp;
        if(!possivel(quadrante)){
            contTemp++;
            continue;
        }
        geraPopulacao(quadrante);
        calculaFitPopulacao();
        std::set<Individuo> top;
        melhorPorGeracao.clear();

        auto end = std::chrono::high_resolution_clock::now();
        setupTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        int contMelhor=0;
        for(int i=0;i<geracoes;i++){
            //mostraPopulacao();
            
            auto start1 = std::chrono::high_resolution_clock::now();
            auto elite = melhores(populacao, quantidadeElitismo);
            end = std::chrono::high_resolution_clock::now();
            eliteTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            //std::cout<<"\nChamando Cruzamento\n";
            start1 = std::chrono::high_resolution_clock::now();
            cruzaPopulacao();
            end = std::chrono::high_resolution_clock::now();
            crossTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            //std::cout<<"\nChamando Mutacao\n";
            start1 = std::chrono::high_resolution_clock::now();
            mutaPopulacao(quadrante);
            end = std::chrono::high_resolution_clock::now();
            mutTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();

            start1 = std::chrono::high_resolution_clock::now();
            trocaAminosEquivPopulacao(quadrante);
            end = std::chrono::high_resolution_clock::now();
            trocaTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();

            //std::cout<<"\nChamando CalculaFitPopulacao\n";
            start1 = std::chrono::high_resolution_clock::now();
            if(i>1){
                if(melhorPorGeracao[i-1].getFitness()==melhorPorGeracao[i-2].getFitness())
                    contMelhor++;
                else
                    contMelhor=0;
            }
            if(contMelhor>5){
                //std::cout<<"\n\nRemovendo o melhor individuo da populacao\nMelhor:"<<melhorPorGeracao[i].getFitness()<<"\nAnterior:"<<melhorPorGeracao[i-5].getFitness();
                auto it = std::find(populacao.begin(),populacao.end(),melhorPorGeracao[i-1]);
                while(it!=populacao.end()){
                    //std::cout<<"Removendo o: "<<(*it).getFitness()<<std::endl;
                    Individuo ind = geraIndividuoValido(quadrante);
                    populacao[it-populacao.begin()]=ind;
                    it = std::find(it,populacao.end(),melhorPorGeracao[i-1]);
                }
                contMelhor=0;
            }
            end = std::chrono::high_resolution_clock::now();
            repTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            start1 = std::chrono::high_resolution_clock::now();
            calculaFitPopulacao();
            end = std::chrono::high_resolution_clock::now();
            fitTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            //std::cout<<"\nChamando Elite\n";
            start1 = std::chrono::high_resolution_clock::now();
            for(auto ind:elite)
                populacao.push_back(ind);
            end = std::chrono::high_resolution_clock::now();
            eliteTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            //std::cout<<"\nChamando melhores\n";
            melhorPorGeracao.push_back(melhordaPopulacao());
            //piorPorGeracao.push_back(piordaPopulacao());
            //mediaPorGeracao.push_back(fitnessMedia());
            //desvioPadraoPorGeracao.push_back(desvioPadraoDaPopulacao());

            //std::cout<<"\nChamando top\n";
            start1 = std::chrono::high_resolution_clock::now();
            top.insert(populacao.begin(),populacao.end());
            //std::cout<<"\nInseri\n";
            
            //std::cout<<"\nPeguei os meelhores\n";
            // std::vector<Individuo> newTop;
            // for(auto ind:top){
            //     if(std::find(newTop.begin(),newTop.end(),ind)==newTop.end())
            //         newTop.push_back(ind);
            // }
            //top=newTop;
            
            if(top.size()>tamanhoPopulacao){
                top.erase(std::next(top.begin(),tamanhoPopulacao),top.end());
            }

            end = std::chrono::high_resolution_clock::now();
            topTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
            //std::cout<<"\nRemovi repetidos\n";
           

            //std::cout<<"\nVendo repeticao\n";
            // std::cout<<"\n\nPopulacao antes:\n";
            // mostraPopulacao();
            //Remove o melhor indivíduo da população se ele não mudar por 5 gerações
            
            //std::cout<<"\n\nPopulacao depois:\n";
            //mostraPopulacao();

            //std::cout<<"\n\nGeracao "<<i<<" fitness="<<melhorPorGeracao[i].getFitness()<<std::endl;
            //std::cout.flush();
        }
        auto start1 = std::chrono::high_resolution_clock::now();
        salvaPopulacao(top, contTemp>0);
        end = std::chrono::high_resolution_clock::now();
        saveTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start1).count();
        totalTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        // std::cout<<"Tempo total: "<<totalTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de setup: "<<setupTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de fitness: "<<fitTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de elite: "<<eliteTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de top: "<<topTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de repeticao: "<<repTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de save: "<<saveTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de crossover: "<<crossTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de mutacao: "<<mutTime/1000000000.0<<"s"<<std::endl;
        // std::cout<<"Tempo de troca: "<<trocaTime/1000000000.0<<"s"<<std::endl;
        // //Mostrar tempos em porcentagem do total
        // std::cout<<"Tempo de setup: "<<(setupTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de fitness: "<<(fitTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de elite: "<<(eliteTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de top: "<<(topTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de repeticao: "<<(repTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de save: "<<(saveTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de crossover: "<<(crossTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de mutacao: "<<(mutTime/totalTime)*100<<"%"<<std::endl;
        // std::cout<<"Tempo de troca: "<<(trocaTime/totalTime)*100<<"%"<<std::endl;
        //mostraPopulacao();
        populacao.clear();
        contTemp++;
    }
}

void GA::run(Repositorio *repo, site* temp, int geracoes, int quadrante, std::set<Individuo>* top){
    outputFilePath = "output/"+temp->pdbId+"_"+std::to_string(quadrante)+".txt";
    this->repositorio=repo;
    this->temp=*temp;
    if(!possivel(quadrante)){
        return;
    }
    geraPopulacao(quadrante);
    calculaFitPopulacao();
    int contMelhor=0;
    for(int i=0;i<geracoes;i++){

        auto elite = melhores(populacao, quantidadeElitismo);

        cruzaPopulacao();

        mutaPopulacao(quadrante);

        trocaAminosEquivPopulacao(quadrante);

        calculaFitPopulacao();

        for(auto ind:elite)
            populacao.push_back(ind);
        

        melhorPorGeracao.push_back(melhordaPopulacao());

        //std::cout<<"\nChamando CalculaFitPopulacao\n";
        if(i>1){
            if(melhorPorGeracao[i-1].getFitness()==melhorPorGeracao[i-2].getFitness())
                contMelhor++;
            else
                contMelhor=0;
        }
        if(contMelhor>5){
            //std::cout<<"\n\nRemovendo o melhor individuo da populacao\nMelhor:"<<melhorPorGeracao[i].getFitness()<<"\nAnterior:"<<melhorPorGeracao[i-5].getFitness();
            auto it = std::find(populacao.begin(),populacao.end(),melhorPorGeracao[i-1]);
            while(it!=populacao.end()){
                //std::cout<<"Removendo o: "<<(*it).getFitness()<<std::endl;
                Individuo ind = geraIndividuoValido(quadrante);
                populacao[it-populacao.begin()]=ind;
                it = std::find(it,populacao.end(),melhorPorGeracao[i-1]);
            }
            contMelhor=0;
        }
        


        top->insert(elite.begin(),elite.end());

        
        if(top->size()>quantidadeElitismo){
            top->erase(std::next(top->begin(),quantidadeElitismo),top->end());
        }
    }
    salvaPopulacao(*top, false);

}


void GA::run(int geracoes, int quadrante){
    int contTemp=0;
    for(site temp:*templates){
        //std::cout<<"\n\n#######Rodando para o template "<<contTemp<<std::endl;
        this->temp=temp;
        if(!possivel(quadrante)){
            contTemp++;
            continue;
        }
        geraPopulacao(quadrante);
        calculaFitPopulacao();
        std::set<Individuo> top;
        melhorPorGeracao.clear();


        int contMelhor=0;
        for(int i=0;i<geracoes;i++){
            //mostraPopulacao();
            
            auto start1 = std::chrono::high_resolution_clock::now();
            auto elite = melhores(populacao, quantidadeElitismo);
            //std::cout<<"\nChamando Cruzamento\n";
            cruzaPopulacao();
            //std::cout<<"\nChamando Mutacao\n";
            mutaPopulacao(quadrante);

            trocaAminosEquivPopulacao(quadrante);

            //std::cout<<"\nChamando CalculaFitPopulacao\n";
            // if(i>1){
            //     if(melhorPorGeracao[i-1].getFitness()==melhorPorGeracao[i-2].getFitness())
            //         contMelhor++;
            //     else
            //         contMelhor=0;
            // }
            // if(contMelhor>5){
            //     //std::cout<<"\n\nRemovendo o melhor individuo da populacao\nMelhor:"<<melhorPorGeracao[i].getFitness()<<"\nAnterior:"<<melhorPorGeracao[i-5].getFitness();
            //     auto it = std::find(populacao.begin(),populacao.end(),melhorPorGeracao[i-1]);
            //     while(it!=populacao.end()){
            //         //std::cout<<"Removendo o: "<<(*it).getFitness()<<std::endl;
            //         Individuo ind = geraIndividuoValido(quadrante);
            //         populacao[it-populacao.begin()]=ind;
            //         it = std::find(it,populacao.end(),melhorPorGeracao[i-1]);
            //     }
            //     contMelhor=0;
            // }
            calculaFitPopulacao();
            //std::cout<<"\nChamando Elite\n";
            for(auto ind:elite)
                populacao.push_back(ind);
            //std::cout<<"\nChamando melhores\n";
            melhorPorGeracao.push_back(melhordaPopulacao());
            //piorPorGeracao.push_back(piordaPopulacao());
            //mediaPorGeracao.push_back(fitnessMedia());
            //desvioPadraoPorGeracao.push_back(desvioPadraoDaPopulacao());

            //std::cout<<"\nChamando top\n";
            if(i%10==0){
                top.insert(populacao.begin(),populacao.end());
                //std::cout<<"\nInseri\n";
                
                if(top.size()>tamanhoPopulacao){
                    top.erase(std::next(top.begin(),tamanhoPopulacao),top.end());
                }

                //std::cout<<"\nRemovi repetidos\n";
            }
           

            //std::cout<<"\nVendo repeticao\n";
            // std::cout<<"\n\nPopulacao antes:\n";
            // mostraPopulacao();
            //Remove o melhor indivíduo da população se ele não mudar por 5 gerações
            
            //std::cout<<"\n\nPopulacao depois:\n";
            //mostraPopulacao();

            //std::cout<<"\n\nGeracao "<<i<<" fitness="<<melhorPorGeracao[i].getFitness()<<std::endl;
            //std::cout.flush();
        }
        //salvaPopulacao(top, contTemp>0);
        melhoresIndividuos[contTemp].insert(top.begin(),top.end());
        if(melhoresIndividuos[contTemp].size()>tamanhoPopulacao){
            melhoresIndividuos[contTemp].erase(std::next(melhoresIndividuos[contTemp].begin(),tamanhoPopulacao),melhoresIndividuos[contTemp].end());
        }
        //mostraPopulacao();
        populacao.clear();
        contTemp++;
    }
}

void GA::run(int geracoes, std::vector<std::vector<Individuo>> pop){
    int contTemp=0;

    int quadrante=0;
    for(site temp:*templates){
        //std::cout<<"\n\n#######Rodando para o template "<<contTemp<<std::endl;
        this->temp=temp;
        if(pop[contTemp].size()<tamanhoPopulacao)
            geraPopulacao(quadrante);
        else
            populacao=pop[contTemp];
        geraPopulacao(quadrante);
        calculaFitPopulacao();
        std::set<Individuo> top;
        melhorPorGeracao.clear();


        int contMelhor=0;
        for(int i=0;i<geracoes;i++){
            //mostraPopulacao();
            
            auto elite = melhores(populacao, quantidadeElitismo);
            //std::cout<<"\nChamando Cruzamento\n";
            cruzaPopulacao();
            //std::cout<<"\nChamando Mutacao\n";
            mutaPopulacao(quadrante);

            trocaAminosEquivPopulacao(quadrante);

            //std::cout<<"\nChamando CalculaFitPopulacao\n";
            if(i>1){
                if(melhorPorGeracao[i-1].getFitness()==melhorPorGeracao[i-2].getFitness())
                    contMelhor++;
                else
                    contMelhor=0;
            }
            if(contMelhor>5){
                //std::cout<<"\n\nRemovendo o melhor individuo da populacao\nMelhor:"<<melhorPorGeracao[i].getFitness()<<"\nAnterior:"<<melhorPorGeracao[i-5].getFitness();
                auto it = std::find(populacao.begin(),populacao.end(),melhorPorGeracao[i-1]);
                while(it!=populacao.end()){
                    //std::cout<<"Removendo o: "<<(*it).getFitness()<<std::endl;
                    Individuo ind = geraIndividuoValido(quadrante);
                    populacao[it-populacao.begin()]=ind;
                    it = std::find(it,populacao.end(),melhorPorGeracao[i-1]);
                }
                contMelhor=0;
            }
            calculaFitPopulacao();
            //std::cout<<"\nChamando Elite\n";
            for(auto ind:elite)
                populacao.push_back(ind);
            //std::cout<<"\nChamando melhores\n";
            melhorPorGeracao.push_back(melhordaPopulacao());
            //piorPorGeracao.push_back(piordaPopulacao());
            //mediaPorGeracao.push_back(fitnessMedia());
            //desvioPadraoPorGeracao.push_back(desvioPadraoDaPopulacao());

            //std::cout<<"\nChamando top\n";
            if(i%10==0){
                top.insert(populacao.begin(),populacao.end());
                //std::cout<<"\nInseri\n";
                
                if(top.size()>quantidadeElitismo){
                    top.erase(std::next(top.begin(),quantidadeElitismo),top.end());
                }

                //std::cout<<"\nRemovi repetidos\n";
            }
           

            //std::cout<<"\nVendo repeticao\n";
            // std::cout<<"\n\nPopulacao antes:\n";
            // mostraPopulacao();
            //Remove o melhor indivíduo da população se ele não mudar por 5 gerações
            
            //std::cout<<"\n\nPopulacao depois:\n";
            //mostraPopulacao();

            //std::cout<<"\n\nGeracao "<<i<<" fitness="<<melhorPorGeracao[i].getFitness()<<std::endl;
            //std::cout.flush();
        }
        salvaPopulacao(top, contTemp>0);
        // melhoresIndividuos[contTemp].insert(top.begin(),top.end());
        // if(melhoresIndividuos[contTemp].size()>quantidadeElitismo){
        //     melhoresIndividuos[contTemp].erase(std::next(melhoresIndividuos[contTemp].begin(),quantidadeElitismo),melhoresIndividuos[contTemp].end());
        // }
        
        //mostraPopulacao();
        populacao.clear();
        contTemp++;
    }
}

void GA::salvaPopulacao( std::set<Individuo> populacao, bool append){
    std::ofstream ofs(outputFilePath, append ? std::ios::app : std::ios::out);
    int contador=1;
    int tam_sitiof = temp.residuos.size();
    int i=0;
    for (Individuo ind:populacao)
    {
        ofs << contador << "\t" << i << "\t";
        ofs << std::setprecision(3) << ind.getFitness() << "\t";
        // ofs << nome_proteina << " " << fittop[i] << " ";
        if (ind.getFitness() != 1000)
        {
            ofs << ind.cromossomo[0].amino << " " << ind.cromossomo[0].atomo_ID << " " << ind.cromossomo[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << ind.cromossomo[j].amino << " " << ind.cromossomo[j].atomo_ID << " " << ind.cromossomo[j].cadeia;
            ofs << "\t" << temp.pdbId << "\t";
            ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
        }
        else
        {
            ofs << " - "
                << " "
                << " - "
                << " "
                << " - ";
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";"
                    << " - "
                    << " "
                    << " - "
                    << " "
                    << " - ";
            ofs << "\t" << temp.pdbId << "\t";
            ofs << temp.residuos[0].amino << " " << temp.residuos[0].atomo_ID << " " << temp.residuos[0].cadeia;
            for (int j = 1; j < tam_sitiof; j++)
                ofs << ";" << temp.residuos[j].amino << " " << temp.residuos[j].atomo_ID << " " << temp.residuos[j].cadeia;
        
        }
        //ofs << "\t" << ecnumber << "\t" << uniprot << "\t" << resolution << "\t" << tam_sitiof;
        ofs << std::endl;
        i++;
    }
}

void GA::calculaFitPopulacao(){
    for(Individuo& ind:populacao)
        ind.setFitness(calculaFit(ind));
}

void GA::geraPopulacao(int quadrante){
    for(int i=0;i<this->tamanhoPopulacao;i++)
        populacao.push_back(geraIndividuoValido(quadrante));
}

double GA::calculaFit(Individuo ind){
    int tam_sitiof = ind.cromossomo.size();
    float fit = 0.0f, aux;
    std::vector<float> distanciasitio(tam_sitiof*(tam_sitiof-1)/2);
    int pos=0;
    for(int i=0; i<tam_sitiof; i++)
        for(int j=i+1; j<tam_sitiof; j++){
            aux = sqrt(pow(ind.cromossomo[i].x - ind.cromossomo[j].x, 2)
                    + pow(ind.cromossomo[i].y - ind.cromossomo[j].y, 2)  
                    + pow(ind.cromossomo[i].z - ind.cromossomo[j].z, 2));
            distanciasitio[pos++] =aux;
        }
    std::sort(distanciasitio.begin(), distanciasitio.end());
    //std::sort(distanciasTemplate.begin(), distanciasTemplate.end());
	for(int i=0;i<distanciasitio.size();i++)
        fit += fabs(distanciasitio[i]-temp.distances[i]);


	for(int j=0; j<tam_sitiof; j++)
	  for(int i=j+1; i < tam_sitiof; i++)
		if ((ind.cromossomo[j].atomo_ID == ind.cromossomo[i].atomo_ID) && 
            (ind.cromossomo[j].cadeia == ind.cromossomo[i].cadeia)){
			fit += 100.0;
			i = tam_sitiof;
			j = tam_sitiof;
		}
    return fit;
}

void GA::cruzaPopulacao(){
    std::vector<Individuo> novapop(tamanhoPopulacao-quantidadeElitismo);
    std::vector<Individuo> filhos;
    int i=0;
    while (i<(tamanhoPopulacao-quantidadeElitismo)){
        Individuo pai1 = torneio();
        Individuo pai2 = torneio();

            filhos = cruza(pai1, pai2);
            novapop[i++] = filhos[0];
            if(i<(tamanhoPopulacao-quantidadeElitismo))
                novapop[i++] = filhos[1];
    }
    populacao=std::move(novapop);
}

template <typename T>
std::vector<Individuo> GA::melhores( T vec, int n){
    std::vector<Individuo> v(n);
    std::partial_sort_copy(
        std::begin(vec), std::end(vec), //.begin/.end in C++98/C++03
        std::begin(v), std::end(v),
        std::less<Individuo>()
    );
    return v;
}

std::vector<Individuo> GA::cruza( Individuo pai1,  Individuo pai2){
    
    std::uniform_real_distribution<double> dist(0,1);
    int size = pai1.cromossomo.size();
    if(dist(gen)<taxaCruzamento){
        std::uniform_int_distribution<int> distribution(0,size-1);
        std::vector<std::vector<bool>> genFilhos(2);
        std::vector<Atomo> genesFilho1 (size);
        std::vector<Atomo> genesFilho2 (size);
        bool nval;



        int corte = distribution(gen);

        for(int i=0;i<corte;i++){
            genesFilho1[i]=(pai1.cromossomo[i]);
            genesFilho2[i]=(pai2.cromossomo[i]);
        }
        for(int i=corte;i<size;i++){
            genesFilho1[i]=pai2.cromossomo[i];
            genesFilho2[i]=(pai1.cromossomo[i]);
        }
        Individuo f1(genesFilho1);
        Individuo f2(genesFilho2);

        std::vector<Individuo>res (2);
        res[0]=f1;
        res[1]=f2;
        return res;
    }
    std::vector<Individuo>res(2);
    res[0]=pai1;
    res[1]=pai2;
    return res;
    
}

void GA::mutaPopulacao(int quadrante){
    std::uniform_real_distribution<double> distribution(0,1);
    int cont=0;
    for(Individuo& ind:populacao)
        if(distribution(gen)<taxaMutacao){
            muta(ind, quadrante);
            cont++;
        }
    //std::cout<<"Mutacoes: "<<cont<<std::endl;
    //std::cout<<"Taxa: "<<100.0*cont/populacao.size()<<"%"<<std::endl;
}

bool GA::possivel(int quadrante){
    //verifica se o repositorio tem a quantidade de aminoacidos necessarios para o template
    const char AMINO_CODES[20][4] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    
    for(int i=0;i<20;i++){
        int cont=0;
        for(Atomo& res:temp.residuos)
            if(res.amino==AMINO_CODES[i])
                cont++;
        if(cont>(*repositorio)[quadrante][AMINO_CODES_MAP.at(AMINO_CODES[i])].size())
            return false;
    }
    return true;
}

void GA::trocaAminosEquivPopulacao(int quadrante){
    std::uniform_real_distribution<double> distribution(0,1);
    for(Individuo& ind:populacao)
        if(distribution(gen)<taxaMutacao)
        trocaAminosEquivalentes(ind, quadrante);
}

void GA::muta(Individuo& ind, int quadrante){
    std::uniform_int_distribution<int> sorteiaPos(0,ind.cromossomo.size()-1);
    int pos = sorteiaPos(gen);

    std::uniform_int_distribution<int> dist(0,(*repositorio)[quadrante][AMINO_CODES_MAP.at(ind.cromossomo[pos].amino)].size()-1);
    ind.cromossomo[pos] = *(*repositorio).rep.at(quadrante).at(AMINO_CODES_MAP.at(ind.cromossomo[pos].amino)).at(dist(gen));

    ind.setFitness(calculaFit(ind));

}

void GA::trocaAminosEquivalentes(Individuo& ind, int quadrante){

    std::uniform_real_distribution<double> distribution(0,1);
    auto vec = ind.cromossomo;
    Individuo mutado(vec);

    vec = ind.cromossomo;
    int tam = vec.size();
    for(int i =0;i<tam;i++){
        if(distribution(gen)<taxaMutacao && temp.substitutions[AMINO_CODES_MAP.at(vec[i].amino)].size()!=0){
            std::uniform_int_distribution<int> sub(0,temp.substitutions[AMINO_CODES_MAP.at(vec[i].amino)].size()-1);
            int subAmino = temp.substitutions[AMINO_CODES_MAP.at(vec[i].amino)][sub(gen)];
            std::uniform_int_distribution<int> dist(0,(*repositorio)[quadrante][subAmino].size()-1);
            vec[i] = *(*repositorio)[quadrante][subAmino][dist(gen)];
        }
    }
    mutado.setCromossomo(vec);
 
    mutado.setFitness(calculaFit(mutado));
    ind = mutado;
}

Individuo GA::torneio(){
    std::uniform_int_distribution<int> distribution(0,tamanhoPopulacao-1 -tamanhoTorneio);
    int posInicial = distribution(gen);
    Individuo aux=populacao[posInicial];
    for(int i=0;i<tamanhoTorneio;i++){
        int pos = posInicial+i;
        if(populacao[pos].getFitness()<aux.getFitness())
            aux=populacao[pos];
    }
    return aux;
}

Individuo GA::pior (std::vector<Individuo>vec){
    return *(std::max_element(vec.begin(), vec.end()));;
}

Individuo GA::piordaPopulacao (){
    return pior(populacao);
}

Individuo GA::melhordaPopulacao(){
    return melhor(populacao);
}

double GA::fitnessMedia(){
    double media = 0;
    for(auto ind:populacao)
        media+=ind.getFitness()/populacao.size();
    return media;
}

void GA::setPopulation(std::vector<Individuo> nova){
    this->populacao=nova;
}

std::vector<Individuo> GA::getmelhorPorGeracao(){
    return this->melhorPorGeracao;
}

std::vector<Individuo> GA::getpiorPorGeracao(){
    return this->piorPorGeracao;
}

std::vector<double> GA::getMediaPorGeracao(){
    return this->mediaPorGeracao;
}

std::vector<double> GA::getDesvioPadraoPorGeracao(){
    return this->desvioPadraoPorGeracao;
}

double GA::desvioPadraoDaPopulacao(){
    double media = fitnessMedia();
    double dp=0;
    for(auto ind:populacao){
        dp+=pow((ind.getFitness() - media),2)/populacao.size();
    }
    return dp;
}

GA::~GA(){
    this->populacao.clear();
    this->melhorPorGeracao.clear();
    this->piorPorGeracao.clear();
    this->mediaPorGeracao.clear();
    this->desvioPadraoPorGeracao.clear();
    this->melhoresIndividuos.clear();

    this->populacao.shrink_to_fit();
    this->melhorPorGeracao.shrink_to_fit();
    this->piorPorGeracao.shrink_to_fit();
    this->mediaPorGeracao.shrink_to_fit();
    this->desvioPadraoPorGeracao.shrink_to_fit();
    this->melhoresIndividuos.shrink_to_fit();
}
