/*! \mainpage Trabalho de Introdução à Computação Científica 2019/02

 \section introducao Introdução
*		O objetivo deste trabalho é implementar um programa computacional
*		para calcular a solução discreta para uma Equação Diferencial Parcial
*		com duas variáveis indepententes utilizando Diferenças Finitas centrais
*		de primeira ordem e o método de Gauss-Seidel.
*
*\subsection autores Autores
* 		GRR20163049 Bruno Henrique Labres\n
*		GRR20182981 Giovani Gurkevicz Marciniak
*
*\subsection professor Professor
* 		Armando Luiz Nicolini Delgado
*/

/**
*	@file pdeSolver.c
*	@author GRR20163049 Bruno Henrique Labres
*	@author GRR20182981 Giovani Gurkevicz Marciniak
*	@date 6 de Outubro de 2019
*	@brief Arquivo de programa principal.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "edp_lib.h"

/**
	@brief Função principal do programa. Cria as estruturas de dados e chama funções
		para calcular a solução discreta para uma Equação Diferencial Parcial
 		com duas variáveis indepententes utilizando Diferenças Finitas centrais
 		de primeira ordem e o método de Gauss-Seidel.
	@param argc Numero de argumentos da linha de comando
	@param argv Argumentos da linha de comando
	@return Retorna 0 em caso de sucesso ou EXIT_FAILURE caso ocorra algum erro.
*/


int main(int argc, char *argv[]){
	int nx, ny, maxIter;
	FILE * arquivo_saida;
	real_t lx = PI, ly = PI;
	real_t mediaTempo;
	EDP_t *e;

	arquivo_saida = trataArgumentos(argc, argv, &nx, &ny, &maxIter); // trata argumentos da linha de comando
	real_t *vetorResiduos = (real_t*) malloc(maxIter*(sizeof(real_t)));

	e = criaTipoEDP(nx, ny, lx, ly, maxIter); // cria estruturas de dados necessarias para os calculos envolvendo EDPs

	mediaTempo = calculaGaussSeidel(e, vetorResiduos);
	escreveSolucao(arquivo_saida, e, vetorResiduos, mediaTempo) ;
}
