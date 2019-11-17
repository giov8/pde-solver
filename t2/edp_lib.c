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
	unsigned int i;
	unsigned int contador_arg = 0;
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
void calculaMatrizCoef(EDP_t *e, XB_t *xb){
	unsigned int i, j;
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
			xb[j*e->nx+i].b = (2*hx*hx*hy*hy)*F(xi,yi); // calcula vetor de termos independentes

			if(j == 0){
				// aplica condicao de fronteira
				xb[j*e->nx+i].b -= e->dia*sin(2*PI*(PI - xi))*sinh(PI*PI);
			}
			if(j == e->ny-1){
				// / aplica condicao de fronteira
				xb[j*e->nx+i].b -= e->dsa*sin(2*PI*xi)*sinh(PI*PI);
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

real_t calculaGaussSeidel(EDP_t* restrict e, real_t* restrict r, XB_t* restrict xb){
	/*register*/
  unsigned int nx = e->nx;
	unsigned int ny = e->ny;
	calculaMatrizCoef(e, xb);

	real_t diag[5];
	diag[DP] = e->dp;
	diag[DS] = e->ds;
	diag[DI] = e->di;
	diag[DSA] = e->dsa;
	diag[DIA] = e->dia;


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
	for (unsigned int iter = 0 ; iter < e->maxIter ; iter++) {
		LIKWID_MARKER_START("Gauss");

		// calcula por fora o valor da "matriz" pos[0,0]
		xb[0].x = xb[0].b - diag[DS]*xb[1].x - diag[DSA]*xb[nx].x;
	  	xb[0].x = xb[0].x/diag[DP];

		// calcula a primeira linha da "matriz"
		for (unsigned int i = 1 ; i <= nx-3  ; i+=2) {
			xb[i].x = xb[i].b - diag[DI]*xb[i-1].x;
			real_t aux = diag[DS]*xb[i+1].x + diag[DSA]*xb[i+nx].x;
			xb[i].x -= aux;
	    	xb[i].x = xb[i].x/diag[DP];

	    	xb[i+1].x = xb[i+1].b - diag[DI]*xb[i].x;
			real_t aux2 = diag[DS]*xb[i+2].x + diag[DSA]*xb[i+1+nx].x;
			xb[i+1].x -= aux2;
	    	xb[i+1].x = xb[i+1].x/diag[DP];

		}

		// calcula por fora o ultimo valor da "matriz" pos[0,nx]
		xb[nx-1].x = xb[nx-1].b - diag[DI]*xb[nx-2].x - diag[DSA]*xb[2*nx-1].x;
		xb[nx-1].x = xb[nx-1].x/diag[DP];

		// calcula a "matriz" interna
		for (unsigned int j = 1 ; j <= ny -2 ; j++) {
			// calculando o primeiro valor
	  		unsigned int k = j*nx; // só pra não calcular toda vez
	  		xb[k].x = xb[k].b - diag[DS]*xb[k+1].x;
	  		real_t aux = diag[DIA]*xb[k-nx].x + diag[DSA]*xb[k+e->nx].x;
	 		xb[k].x -= aux;
	 		xb[k].x = xb[k].x/diag[DP];

	        // calculando valores centrais
	      	for (unsigned int i = 1 ; i <= nx -3 ; i +=2) {
				xb[k+i].x = xb[k+i].b - diag[DI]*xb[k+i-1].x;
				aux = diag[DS]*xb[k+i+1].x +  diag[DIA]*xb[k+i-nx].x;
				xb[k+i].x = xb[k+i].x - aux - diag[DSA]*xb[k+i+nx].x;
				xb[k+i].x = xb[k+i].x/diag[DP];

				xb[k+i+1].x = xb[k+i+1].b - diag[DI]*xb[k+i].x;
				aux = diag[DS]*xb[k+i+2].x +  diag[DIA]*xb[k+i-nx+1].x;
				xb[k+i+1].x = xb[k+i+1].x - aux - diag[DSA]*xb[k+i+nx+1].x;
				xb[k+i+1].x = xb[k+i+1].x/diag[DP];
			}

			// calculando o ultimo valor
			xb[k+nx-1].x = xb[k+nx-1].b - diag[DI]*xb[k+nx-2].x;
			aux = diag[DIA]*xb[k-1].x + diag[DSA]*xb[k+2*nx-1].x;
			xb[k+nx-1].x = xb[k+nx-1].x - aux;
			xb[k+nx-1].x = xb[k+nx-1].x/diag[DP];
		}

	  	// calcula pos[0,ny-1]
	  	unsigned int k = nx*(ny-1); // só pra calcular menos
	  	xb[k].x = xb[k].b - diag[DS]*xb[k+1].x - diag[DIA]*xb[k-nx].x;
	  	xb[k].x = xb[k].x/diag[DP];

	  	// calcula última linha da matriz "fronteira"
		for (unsigned int i = 1 ; i <= nx-2 ; i+=2) {
			xb[k+i].x = xb[k+i].b -  diag[DI]*xb[k+i-1].x;
			real_t aux = diag[DS]*xb[k+i+1].x + diag[DIA]*xb[k+i-nx].x;
			xb[k+i].x -= aux;
		  	xb[k+i].x = xb[k+i].x/diag[DP];

			xb[k+i+1].x = xb[k+i+1].b -  diag[DI]*xb[k+i].x;
			real_t aux2 = diag[DS]*xb[k+i+2].x + diag[DIA]*xb[k+i+1-nx].x;
			xb[k+i+1].x -= aux2;
		  	xb[k+i+1].x = xb[k+i+1].x/diag[DP];
		}

	  	// calcula pos[nx-1,ny-1]
	  	xb[(nx)*(ny)-1].x = xb[(nx)*(ny)-1].b - diag[DI]*xb[(nx)*(ny)-2].x - diag[DIA]*xb[(nx)*(ny)-1-nx].x;
	  	xb[(nx)*(ny)-1].x = xb[(nx)*(ny)-1].x/diag[DP];
		LIKWID_MARKER_STOP("Gauss");

		LIKWID_MARKER_START("Residuo");
		r[iter] = calculaResiduo(e, xb);
		LIKWID_MARKER_STOP("Residuo");

	}

	LIKWID_MARKER_CLOSE;
	return 0;
}


/**
*  \brief Essa função calcula do resíduo de uma iteração do GaussSeidel
*  \param e Estrutura de dados com as variaveis para os calculos envolvendo EDPs
*  \return sqrt(normaquadrada) a norma L2 da iteração
*/
real_t calculaResiduo(EDP_t* restrict e, XB_t* restrict xb) {
	real_t residuo, residuo2;
	real_t normaquadrada = 0.0;

	real_t diag[5];
	diag[DP] = e->dp;
	diag[DS] = e->ds;
	diag[DI] = e->di;
	diag[DSA] = e->dsa;
	diag[DIA] = e->dia;

	int nx = e->nx;
	int ny = e->ny;

	// calcula por fora o valor da "matriz" pos[0,0]
	residuo = xb[0].b - (diag[DP]*xb[0].x) - diag[DSA]*xb[nx].x;
	normaquadrada = residuo * residuo;

	// calcula a primeira linha da "matriz"
	for (unsigned int i = 1 ; i <= nx-3  ; i+=2) {
		residuo =  xb[i].b - xb[i].x*diag[DP] - diag[DI]*xb[i-1].x;
		real_t aux = diag[DS]*xb[i+1].x + diag[DSA]*xb[i+nx].x;
    	residuo = residuo - aux;
    	normaquadrada += residuo * residuo;

		residuo2 = xb[i+1].b - xb[i+1].x*diag[DP] - diag[DI]*xb[i].x;
		real_t aux2 = diag[DS]*xb[i+2].x + diag[DSA]*xb[i+1+nx].x;
		residuo2 = residuo2 - aux2;
		normaquadrada += residuo2 * residuo2;
	}

	// calcula por fora o ultimo valor da "matriz" pos[0,nx]
	residuo =  xb[nx-1].b - xb[nx-1].x*diag[DP] - diag[DI]*xb[nx-2].x - diag[DSA]*xb[2*nx-1].x;
	normaquadrada =  residuo * residuo;

	// calcula a "matriz" interna
	for (unsigned int j = 1 ; j <= ny -2 ; j++) {
		// calculando o primeiro valor
  		unsigned int k = j*nx; // só pra não calcular toda vez
  		residuo = xb[k].b - xb[k].x*diag[DP] - diag[DS]*xb[k+1].x;
  		real_t aux = diag[DIA]*xb[k-nx].x + diag[DSA]*xb[k+e->nx].x;
 		residuo = residuo - aux;
    	normaquadrada += residuo * residuo;

        // calculando valores centrais
      	for (unsigned int i = 1 ; i <= nx -3 ; i +=2) {
			residuo = xb[k+i].b - xb[k+i].x*diag[DP] - diag[DI]*xb[k+i-1].x;
			aux = diag[DS]*xb[k+i+1].x +  diag[DIA]*xb[k+i-nx].x;
			residuo = residuo - aux - diag[DSA]*xb[k+i+nx].x;
    		normaquadrada += residuo * residuo;

			residuo2 = xb[k+i+1].b - xb[k+i+1].x*diag[DP] - diag[DI]*xb[k+i].x;
			real_t aux2 = diag[DS]*xb[k+i+2].x +  diag[DIA]*xb[k+i-nx+1].x;
			residuo2 = residuo2 - aux2 - diag[DSA]*xb[k+i+nx+1].x;
			normaquadrada += residuo2 * residuo2;
		}

		// calculando o ultimo valor
		residuo = xb[k+nx-1].b - xb[k+nx-1].x*diag[DP] - diag[DI]*xb[k+nx-2].x;
		aux = diag[DIA]*xb[k-1].x + diag[DSA]*xb[k+2*nx-1].x;
		residuo =  residuo - aux;
		normaquadrada += residuo * residuo;
	}

	// calcula pos[0,ny-1]
	unsigned int k = nx*(ny-1); // só pra calcular menos
	residuo = xb[k].b - xb[k].x*diag[DP] - diag[DS]*xb[k+1].x - diag[DIA]*xb[k-nx].x;
	normaquadrada += residuo * residuo;

	// calcula última linha da matriz "fronteira"
	for (unsigned int i = 1 ; i <= nx-2 ; i+=2) {
		residuo = xb[k+i].b - xb[k+i].x*diag[DP] - diag[DI]*xb[k+i-1].x;
		real_t aux = diag[DS]*xb[k+i+1].x + diag[DIA]*xb[k+i-nx].x;
		residuo = residuo - aux;
	  	normaquadrada += residuo * residuo;

		residuo2 = xb[k+i+1].b - xb[k+i+1].x*diag[DP] - diag[DI]*xb[k+i].x;
		real_t aux2 = diag[DS]*xb[k+i+2].x + diag[DIA]*xb[k+i+1-nx].x;
		residuo2 = residuo2 - aux2;
	  	normaquadrada += residuo2 * residuo2;
	}
	// calcula pos[nx-1,ny-1]
	residuo = xb[(nx)*(ny)-1].b -  xb[(nx)*(ny)-1].x*diag[DP]  - diag[DI]*xb[(nx)*(ny)-2].x - diag[DIA]*xb[(nx)*(ny)-1-nx].x;
	normaquadrada += residuo * residuo;

	return sqrt(normaquadrada);
}

/*!
  \brief Essa função escreve a solucao da EDP em um arquivo

  \param arquivo_saida Arquivo a ter a solucao escrita
  \param nx Tamanho da grade
  \param ny Tamanho da grade
  \param s Vetor de solucao
*/
void escreveSolucao(FILE *arquivo_saida, EDP_t *e, real_t *r, real_t mediaTempo, XB_t *xb){
	unsigned int i,j,k;
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
			fprintf(arquivo_saida, "%g %g %g \n",xi,yi,xb[k].x);
			k++;
		}
		fprintf(arquivo_saida, "\n");
	}
	fclose(arquivo_saida);
}
