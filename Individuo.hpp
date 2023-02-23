#pragma once

#include <cmath>
#include <vector>
#include <random>
#include "Atomo.hpp"


class Individuo
{
protected:
    

    /**
     * @brief Aptidao deste individuo
     */
    double fitness;


public:
    /**
     * @brief Vetor de genes pertencentes a este individuo que formam seu genotipo
     */
    std::vector<Atomo> cromossomo;
    /**
     * @brief Constroi um novo objeto da classe Individuo
     * 
     */
    Individuo();

    /**
     * @brief Constroi um novo objeto da classe Individuo passando seu genotipo
     * 
     * @param cromossomo genotipo deste individuo
     */
    Individuo(std::vector<Atomo> cromossomo );

    /**
     * @brief Retorna o cromossomo(genotipo) deste individuo
     * 
     * @return std::vector<Atomo> genotipo
     */
    virtual std::vector<Atomo> getCromossomo() const;

    /**
     * @brief Retorna a aptidao deste individuo
     * 
     * @return double aptidao
     */
    virtual double getFitness();
    
    /**
     * @brief Define a aptidao deste individuo
     * 
     * @param x valor a ser definido como a aptidao
     */
    virtual void setFitness(double x);

    /**
     * @brief Define o genotipo deste individuo
     * 
     * @param x valor a ser definido como o genotipo
     */
    virtual void setCromossomo(std::vector<Atomo> x);

    /**
     * @brief Operador de comparacao entre individuos
     * 
     * @param ind1 individuo a esquerda do sinal de comparacao
     * @param ind2 individuo a direita
     * @return true caso a aptidao do individuo a esquerda for menor do que a do individuo a direita
     * @return false caso contrario
     */
    friend bool operator<( Individuo &ind1,  Individuo &ind2);

    friend bool operator<( const Individuo &ind1,  const Individuo &ind2);

    /**
     * @brief Operador de igualdade entre individuos
     * 
     * @param ind1 individuo a esquerda do sinal de comparacao
     * @param ind2 individuo a direita
     * @return true caso a aptidao de ambos os individuos for a mesma
     * @return false caso contrario
     */
    friend bool operator==(const  Individuo ind1,  const Individuo ind2);


};

