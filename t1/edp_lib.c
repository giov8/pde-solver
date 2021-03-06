/**
*	@file edp_lib.c
*	@author GRR20163049 Bruno Henrique Labres
*	@author GRR20182981 Giovani Gurkevicz Marciniak 
*	@date 6 de Outubro de 2019
*	@brief Biblioteca com as funções envolvendo equações diferenciais parciais.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "edp_lib.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

/*!
  @brief Essa função retorna o tempo atual em milisegundos
  @return o tempo em milisegundos no formato real_t
*/
real_t timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((real_t)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}


/**
*  @brief Essa função cria um estrutura de dados necessarias para o calculo de EDPs (Equação Diferencial Parcial) 
*  @param nx Número de pontos no eixo x da malha
*  @param ny Número de pontos no eixo y da malha
*  @param lx Limite do eixo x na malha
*  @param ly Limite do eixo y na malha
*  @return e Ponteiro para uma estrutura EDP_t com suas variaveis inicializadas
*/

EDP_t *criaTipoEDP(int nx, int ny, real_t lx, real_t ly, int maxIter){
	EDP_t *e = (EDP_t *) malloc (sizeof(EDP_t));
	real_t hx, hy;

	// recebe a qtd desejada de iteracoes
	e->maxIter = maxIter;

	// define numero de equacoes em x e y
	e->nx = nx;
	e->ny = ny;

	// define limites em x e y
	e->lx = lx;
	e->ly = ly;

	e->x = (real_t *) calloc (nx * ny, sizeof (real_t)); // vetor solucao
	e->x_prev = (real_t *) calloc (nx * ny, sizeof (real_t)); // vetor solucao
	e->b = (real_t *) calloc (nx * ny, sizeof (real_t)); // vetor de termos independentes


	// passo dado na malha
	hx = lx/(nx+1);
	hy = ly/(ny+1);

	// diagonais da matriz A
	e->dp = 4*hy*hy + 4*hx*hx + 8*PI*PI*hx*hx*hy*hy; // diagonal principal
	e->ds = (hy*hy)*(-2 + hx); // diagonal superior
	e->di = (hy*hy)*(-2 - hx); // diagonal inferior
	e->dsa = (hx*hx)*(-2 + hy);  // diagonal superior afastada
	e->dia = (hx*hx)*(-2 - hy); // diagonal inferior afastada

	e->hx = hx;
	e->hy = hy;
	return e;
}

/**
*	@brief Função que trata os argumentos dados na linha de comando
*	@param argc Numero de argumentos da linha de comando
*	@param argv Argumentos da linha de comandos
*	@param nx Numero de pontos a serem calculados na dimensão X
*	@param ny Numero de pontos a serem calculados na dimensão Y
*	@param maxIter Numero maximo de iteracoes a serem executadas
*	@return Retorna um FILE * que representa o arquivo de saida
*/
FILE * trataArgumentos(int argc, char *argv[], int *nx, int *ny, int *maxIter){
	int i;
	int contador_arg = 0;
	FILE * arquivo_saida = stdout;
	// busca por nx, ny, maxIter, arquivo_saida
	for (i = 0 ; i < argc ; i++){
		if (!strcmp(argv[i], "-nx")){
			*nx = atoi(argv[i+1]);
			contador_arg++;
		}
		
		if (!strcmp(argv[i], "-ny")){
			*ny = atoi(argv[i+1]);
			contador_arg++;
		}
	
		if (!strcmp(argv[i], "-i")){
			*maxIter = atoi(argv[i+1]);
			contador_arg++;
		}

		if (!strcmp(argv[i], "-o")){
			arquivo_saida = fopen(argv[i+1], "w+");
		}
	}

	if (contador_arg < 3) {
			fprintf(stderr, "Argumentos faltando no programa. O programa não será executado. Verifique se está no padrão: pdeSolver -nx <Nx> -ny <Ny> -i <maxIter> -o arquivo_saida\n");
			exit(EXIT_FAILURE);
		}

	return arquivo_saida;
}

/**
*	@brief Função que calcula os termos independentes, aplicando condicoes de fronteira
*	@param e Estrutura de dados com as variaveis para os calculos envolvendo EDPs
*/
void calculaMatrizCoef(EDP_t *e){
	int i, j;
	real_t hx, hy, xi, yi;

	// passos da malha
	hx = e->lx/(e->nx+1);
	hy = e->ly/(e->ny+1);

	// gera coeficientes
	yi = 0.0;
	for(j = 0; j < e->ny; j++) {
		yi += hy;
		xi = 0.0;
		for(i = 0; i < e->nx ; i++) {
			xi += hx;	
			e->b[j*e->nx+i] = (2*hx*hx*hy*hy)*F(xi,yi); // calcula vetor de termos independentes
	 	
			if(j == 0){
				// aplica condicao de fronteira
				e->b[j*e->nx+i] -= e->dia*sin(2*PI*(PI - xi))*sinh(PI*PI);
			}
			if(j == e->ny-1){
				// / aplica condicao de fronteira
				e->b[j*e->nx+i] -= e->dsa*sin(2*PI*xi)*sinh(PI*PI);
			}
	 	}
	}
}

/**
*  \brief Essa função calcula o método de Gauss Seidel para a EDP e retorna o tempo médio
*  \param e Estrutura de dados com as variaveis para os calculos envolvendo EDPs
*  \param r Vetor de resíduos a cada iteração
*  \return mediaTempo Tempo medio do calculo de Gauss Seidel
*/
real_t calculaGaussSeidel(EDP_t *e, real_t *r){
	calculaMatrizCoef(e);
	
	int i,j,pos,iter;

	if (e->nx <= 0) {
		fprintf(stderr, "Parametro nx inválido.\n");
		exit(EXIT_FAILURE);
	}

	if (e->ny <= 0) {
		fprintf(stderr, "Parametro ny inválido.\n");
		exit(EXIT_FAILURE);
	}

	if (e->maxIter <= 0) {
		fprintf(stderr, "Parametro i inválido.\n");
		exit(EXIT_FAILURE);
	}


	LIKWID_MARKER_INIT;
	for(iter = 0 ; iter < e->maxIter ; iter++)
	{
		LIKWID_MARKER_START("Gauss");
		for(j = 0 ; j < e->ny ; j++){
			for(i = 0; i < e->nx; i++){
				pos = j*e->nx+i; // iterar pelo vetor solucao

				e->x[pos] = e->b[pos];
				
				// analisa os casos de fronteira
				if(i != 0)    e->x[pos] -= e->di*e->x[pos-1];
				if(i != (e->nx-1)) e->x[pos] -= e->ds*e->x_prev[pos+1];
				if(j != 0)    e->x[pos] -= e->dia*e->x[pos-e->nx];
				if(j != (e->ny-1)) e->x[pos] -= e->dsa*e->x_prev[pos+e->nx];
				
				e->x[pos] = e->x[pos]/e->dp; // guarda valor de z na malha
				
				e->x_prev[pos] = e->x[pos]; // guarda it anterior
			}
		}
		LIKWID_MARKER_STOP("Gauss");

		//LIKWID_MARKER_START("Residuo");
		r[iter] = calculaResiduo(e);
		//LIKWID_MARKER_STOP("Residuo");
	}

	LIKWID_MARKER_CLOSE;
	return 1;
}


/**
*  \brief Essa função calcula do resíduo de uma iteração do GaussSeidel 
*  \param e Estrutura de dados com as variaveis para os calculos envolvendo EDPs
*  \return sqrt(normaquadrada) a norma L2 da iteração
*/
real_t calculaResiduo(EDP_t *e) {
	int i, j, pos;
	real_t residuo;
	real_t normaquadrada = 0.0;

	for(j = 0 ; j < e->ny ; j++){
		for(i = 0; i < e->nx; i++){
			pos = j*e->nx+i; // iterar pelo vetor 
			
			residuo = (e->b[pos]) - (e->dp*e->x[pos]);
		
			//analisa os casos de fronteira
			if(i != 0)    residuo -= e->di * e->x[pos-1];
			if(i != (e->nx-1)) residuo -= e->ds * e->x[pos+1];
			if(j != 0)    residuo -= e->dia * e->x[pos-e->nx];
			if(j != (e->ny-1)) residuo -= e->dsa * e->x[pos+e->nx];
			
			normaquadrada += residuo * residuo;
		}
	}				

	return sqrt(normaquadrada);

}

/*!
  \brief Essa função escreve a solucao da EDP em um arquivo

  \param arquivo_saida Arquivo a ter a solucao escrita
  \param nx Tamanho da grade
  \param ny Tamanho da grade
  \param s Vetor de solucao
*/
void escreveSolucao(FILE *arquivo_saida, EDP_t *e, real_t *r, real_t mediaTempo){
	int i,j,k;
	real_t xi, yi;

	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"# Tempo Método GS: %g ms\n#\n", mediaTempo);
	fprintf(arquivo_saida,"# Norma L2 do Residuo\n");

	for (i = 0; i < e->maxIter; i++) {
		fprintf(arquivo_saida,"# i=%d: %g\n", i,r[i]);
	}
	fprintf(arquivo_saida,"###########\n\n");

	xi = 0.0;
	k = 0;
	for (i = 0 ; i < e->nx ; i++){
		xi+= e->hx;
		yi = 0.0;
		for (j = 0 ; j < e->ny ; j++){
			yi += e->hy;
			fprintf(arquivo_saida, "%g %g %g \n",xi,yi,e->x[k]);
			k++;
		}
		fprintf(arquivo_saida, "\n");
	}
	fclose(arquivo_saida);
}
