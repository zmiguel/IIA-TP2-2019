#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pLocal.h"
#include "util.h"

int main(int argc, char *argv[]){
    char ficheiro[100];
	int vert, itr, k, runs, custo, bCusto, opt=0;
	int *matrix, *sol, *best;
	float mbf = 0.0;

    if(argc != 3){
        printf("ERRO: numero de argumentos inválido!\n");
        exit(1);
    }

    strcpy(ficheiro,argv[1]);
    runs = atoi(argv[2]);

    prep_random();

    matrix = init_dados_local(ficheiro, &vert, &itr);
	sol = malloc(sizeof(int) * vert);
	best = malloc(sizeof(int) * vert);

    if (sol == NULL || best == NULL)
	{
		printf("download more ram");
		exit(1);
	}
	fflush(stdin);
    do{
		printf("vizinhança 1 ou 2?: ");
		scanf("%d", &opt);
	}while(opt != 1 && opt != 2);

    for (k = 0; k < runs; k++)
	{
		// Gerar solini
		gSolIni_local(sol, vert);
		pSol_local(sol, vert);

		// Trepa colinas
		if(opt == 1){
			custo = tCol1(sol, matrix, vert, itr);
		} else if (opt == 2){
			custo = tCol2(sol, matrix, vert, itr);
		} else{
			printf("vizinhança errada????\n");
        }

        printf("\nRep %d:\n", k);
        pSol_local(sol, vert);
        printf("Custo final: %2d\n", custo);

		mbf += custo;
		if (k == 0 || bCusto > custo)
		{
			bCusto = custo;
			sub(best, sol, vert);
		}
	}

	printf("\n");
	printf("\nMBF: %f\n", mbf / k);
	
	printf("Iterações: %d\n File: %s \n", atoi(argv[2]), argv[1]);
	printf("\nMelhor solucao encontrada:\n");
	pSol_local(best, vert);
	printf("Custo final: %2d\n", bCusto);
	free(matrix);
	free(sol);
	free(best);
	return 0;

}

int tCol1(int sol[], int *mat, int vert, int num_iter){
    int *nova_sol, custo, custo_viz, i;
	nova_sol = malloc(sizeof(int)*vert);
    if(nova_sol == NULL)
    {
        printf("Erro na alocacao de memoria");
        exit(1);
    }
	// Avalia solucao inicial
    custo = cFit_local(sol, mat, vert);
    for(i=0; i<num_iter; i++)
    {
		// Gera vizinho
		viz(sol, nova_sol, vert);
		// Avalia vizinho
		custo_viz = cFit_local(nova_sol, mat, vert);
		// Aceita vizinho se o custo aumentar
        if(custo_viz >= custo)
        {
			sub(sol, nova_sol, vert);
			custo = custo_viz;
        }
    }
    printf("\nCusto: %d\n\n", custo);
    free(nova_sol);
    return custo;
}

void viz(int a[], int b[], int n){
    int i, p1, p2;

    for(i=0; i<n; i++)
        b[i]=a[i];
	// Encontra posicao com valor 0
    do
        p1=random_l_h(0, n-1);
    while(b[p1] != 0);
	// Encontra posicao com valor 0
    do
        p2=random_l_h(0, n-1);
    while(b[p2] != 1);
	// Troca
    b[p1] = 1;
    b[p2] = 0;
}

int tCol2(int sol[], int *mat, int vert, int num_iter){
    int *nova_sol, custo, custo_viz, i;

	nova_sol = malloc(sizeof(int)*vert);
    if(nova_sol == NULL)
    {
        printf("download more ram");
        exit(1);
    }
	// Avalia solucao inicial
    custo = cFit_local(sol, mat, vert);
    for(i=0; i<num_iter; i++)
    {
		// Gera vizinho
		viz2(sol, nova_sol, vert);
		// Avalia vizinho
		custo_viz = cFit_local(nova_sol, mat, vert);
		// Aceita vizinho se o custo diminuir (problema de minimizacao)
        if(custo_viz < custo)
        {
			sub(sol, nova_sol, vert);
			custo = custo_viz;
        }
    }
    free(nova_sol);
    return custo;
}

void viz2(int a[], int b[], int n){
    int i, p1, p2, p3, p4;

    for(i=0; i<n; i++)
        b[i]=a[i];

	// Encontra posicao com valor 0
    do
        p1=random_l_h(0, n-1);
    while(b[p1] != 0);
	// Encontra posicao com valor 0
    do
        p2=random_l_h(0, n-1);
    while(b[p2] != 1);

    // Encontra posicao com valor 0
    do
        p3=random_l_h(0, n-1);
    while(b[p3] != 0);
    // Encontra posicao com valor 0
    do
        p4=random_l_h(0, n-1);
    while(b[p4] != 1);
	// Troca
    b[p1] = 1;
    b[p2] = 0;

    b[p3] = 1;
    b[p4] = 0;
}