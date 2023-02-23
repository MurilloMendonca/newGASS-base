
#include "Individuo.hpp"
#include <algorithm>
Individuo::Individuo(){

}
Individuo::Individuo(std::vector<Atomo> cromossomo ){
    this->cromossomo=cromossomo;
}
std::vector<Atomo> Individuo::getCromossomo() const{
    return this->cromossomo;
}

double Individuo::getFitness(){
    return this->fitness;
}
void Individuo::setFitness(double x){
    this->fitness=x;
}
void Individuo::setCromossomo(std::vector<Atomo> x){
    this->cromossomo=x;
}

bool operator<( Individuo &ind1,  Individuo &ind2){
    return ind1.getFitness()<ind2.getFitness();
}
bool operator==( const Individuo ind1,  const Individuo ind2){
    int s = ind1.getCromossomo().size();
    for(int i=0;i<ind1.getCromossomo().size();i++)
        if(ind1.getCromossomo()[i]!=ind2.getCromossomo()[i]) return false;
    return true;
}
// bool operator==( const Individuo ind1,  const Individuo ind2){
//     if(ind1.cromossomo.size()!=ind2.cromossomo.size()) return false;
//     for(Atomo a:ind1.cromossomo){
//         if(std::find(ind2.cromossomo.begin(),ind2.cromossomo.end(),a) == ind2.cromossomo.end()) return false;
//     }
//     return true;
// }
bool operator<( const Individuo &ind1,  const Individuo &ind2){
    return ind1.fitness<ind2.fitness;
}