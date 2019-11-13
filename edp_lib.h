/**
*	@file edp_lib.h
*	@author GRR20163049 Bruno Henrique Labres
*	@author GRR20182981 Giovani Gurkevicz Marciniak 
*	@date 6 de Outubro de 2019
*	@brief Arquivo de header da biblioteca de equações diferenciais parciais.
*/

#ifndef EDP_LIB
#define EDP_LIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>


/**
 * @brief Funcao a ser calculada
 */
#define F(x,y) (4*PI*PI*((sin(2*PI*x)*sinh(PI*y)) + (sin(2*PI*(PI-x))*sinh(PI*(PI-y)))))

/**
 * @brief Tipo generico para numeros reais
 */
typedef double real_t;

/**
 * @brief Struct que guarda variaveis para o calculo de EDPs
 */
typedef struct EDP_t {
	//real_t *x;
	//real_t *b;
	real_t dp, ds, di, dsa, dia; // diagonais
	real_t lx, ly, hx, hy;
	int nx, ny;
	int maxIter; // numero de iteracoes que serao usadas para o calculo iterativo
} EDP_t; 

typedef struct XB_t {
	real_t x;
	real_t b;
} XB_t;

#define PI 3.14159265358979323846

#define DP 0
#define DS 1
#define DI 2
#define DSA 3
#define DIA 4 

real_t timestamp(void);
EDP_t *criaTipoEDP(int nx, int ny, real_t lx, real_t ly, int maxIter);
FILE * trataArgumentos(int argc, char *argv[], int *nx, int *ny, int *maxIter);
void calculaMatrizCoef(EDP_t *e, XB_t *xb);
real_t calculaGaussSeidel(EDP_t *e, real_t *r, XB_t *xb);
real_t calculaResiduo(EDP_t *e);
void escreveSolucao(FILE *arquivo_saida, EDP_t *e, real_t *r, real_t mediaTempo, XB_t *xb);
double timestamp(void);
#endif