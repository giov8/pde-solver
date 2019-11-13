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
#include <sys/time.h>
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
	int i, j, iter, pos;
	int nx, ny, maxIter, k;
	FILE * arquivo_saida;
	real_t *b, *x;
	real_t lx = PI, ly = PI;
	real_t hx, hy, xi, yi;
	real_t dp, ds, di, dsa, dia;
	real_t norma;
	// trata argumentos da linha de comando
	arquivo_saida = trataArgumentos(argc, argv, &nx, &ny, &maxIter);

	hx = lx/(nx+1);
	hy = ly/(ny+1);

	x = (real_t *) calloc (nx * ny, sizeof (real_t)); // vetor solucao
	b = (real_t *) calloc (nx * ny, sizeof (real_t)); // vetor de termos independentes


	// diagonais da matriz A
	dp = 4*hy*hy + 4*hx*hx + 8*PI*PI*hx*hx*hy*hy; // diagonal principal
	ds =  -2*hy*hy + hx*hy*hy; // diagonal superior
	di = -2*hy*hy - hx*hy*hy; // diagonal inferior
	dsa = -2*hx*hx + hx*hx*hy;  // diagonal superior afastada
	dia = -2*hx*hx - hx*hx*hy; // diagonal inferior afastada

	// gera os termos independentes
	pos = 0; // posicao no vetor de bias
	yi = 0.0;
	for(j = 0 ; j < ny; ++j) {
		yi += hy;
		xi = 0.0;
		for(i = 0 ; i < nx ; ++i) {
			xi += hx;	
			// verifica se faz operacoes com regiao de fronteira
			if(j == 0){
				// aplica condicao de contorno
				b[pos] = (2*hx*hx*hy*hy)*F(xi,yi) - dia*sin(2*PI*(PI - xi))*sinh(PI*PI);
			}
			else if(j == (ny-1)){
				// aplica condicao de contorno
				b[pos] = (2*hx*hx*hy*hy)*F(xi,yi) - dsa*sin(2*PI*xi)*sinh(PI*PI);
			}
			else
				b[pos] = (2*hx*hx*hy*hy)*F(xi,yi); // calcula vetor de termos independentes
	 		pos++; // prox posicao
	 	}
	}

	double *x_ant = calloc(nx*ny,sizeof(double));
	double tempo_inicio, tempo_fim;
	double media_tempo = 0.0;

	for(iter= 0; iter< maxIter; iter++)
	{			
		tempo_inicio = 0.0;// timestamp();

		k = 0;
		for(j = 0; j < ny; j++){
			for(i = 0; i < nx; i++){
				x[k] = b[k];
				
				// analisa os casos de fronteira
				if(i != 0)    x[k] -= di*x[k-1];
				if(i != nx-1) x[k] -= ds*x_ant[k+1];
				if(j != 0)    x[k] -= dia*x[k-nx];
				if(j != ny-1) x[k] -= dsa*x_ant[k+nx];
				
				x[k] = x[k]/dp; // guarda valor de z na malha
				
				x_ant[k] = x[k]; // guarda it anterior
				
				k++;
			}
		}				

		tempo_fim = 0.0; //timestamp();
		media_tempo += (tempo_fim - tempo_inicio);
		
		//residuo[it] = normaL2Residuo(SL,nx,ny,x);

	}
	
	media_tempo = media_tempo/maxIter;

	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"# Tempo Método GS: %g ms\n#\n", media_tempo);
	fprintf(arquivo_saida,"# Norma L2 do Residuo\n");
	
	for(int i = 0; i < maxIter; i++)
		//fprintf( arquivo_saida, "#i = %d: %g\n" , i+1 , residuo[i]);
	fprintf(arquivo_saida,"###########\n");

	
	
	k = 0;

	for(j = 0; j <= ny+1; j++)
	{
		yi = j*hy;	
		for(i = 0; i <= nx+1 ; i++)
		{
			xi = i*hx;	
	
			fprintf(arquivo_saida,"%g %g ", xi, yi);
			
			if(i == 0 || i == nx+1)
				fprintf(arquivo_saida,"0\n");
			else
				if(j == 0)
				{
					double cont_x_0 = sin(2*PI*(PI - xi))*sinh(PI*PI);//U(x,0) = sin(2pi(pi-xi))*sinH(pi²)
					fprintf(arquivo_saida,"%g\n", cont_x_0);
				}	
				else
					if(j == nx+1)
					{
						double cont_x_pi = sin(2*PI*xi)*sinh(PI*PI); //U(x,pi) = sin(2piX)*sinH(pi²)
						fprintf(arquivo_saida,"%g\n", cont_x_pi);
					}
					else
					{	
						fprintf(arquivo_saida,"%g\n", x[k]);
						k++;
					}
	 	}
	 	fprintf(arquivo_saida, "\n");
	}

	fclose(arquivo_saida);

}