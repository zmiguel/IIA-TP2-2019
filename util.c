#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"

void prep_random(){
	srand((unsigned)time(NULL));
}

int* init_dados_local(char *ficheiro, int *vert, int *iter){
	FILE *f;
	int *p, *q;
	int i, j, col, lin;
	char rip[50];

	f=fopen(ficheiro, "r");
	if(!f)
	{
		printf("ficheiro nao encontrado\n");
		exit(1);
	}
	fscanf(f,"%99[^\n]", &rip);

	fscanf(f, "%d %d %d", vert, vert, iter);

	// Alocacao da matriz
	p = malloc(sizeof(int)*(*vert)*(*vert));
	if(!p)
	{
	    printf("nao tens memoria\n");
	    exit(1);
	}
	// matriz
    for(i=0; i<*vert; i++)
        for(j=0; j<*vert; j++)
            *(p+(*vert)*i+j)=0; //0x0
    for(i=0; i<*iter; i++){
        fscanf(f, "%d %d", &lin, &col);
        *(p+(*vert)*(lin-1)+col-1)=1;
        *(p+(*vert)*(col-1)+lin-1)=1;
    }
	fclose(f);
	return p;
}

void gSolIni_local(int *sol, int v){
	int i, x;
	int flag[v-1];
	int aux;

	for(i = 0; i<v; i++){
		flag[i] = 0;
	}

	for(i=0; i<v; i++){
		do{
			aux = random_l_h(0, v-1); 
			sol[i] = aux;
		} while (flag[aux] != 0);
		flag[aux] = 1;
	}
}

void pSol_local(int *sol, int vert){
	int i;

	printf("Vert: ");
	for(i=0; i<vert; i++)
		printf("%2d  ", i+1);
	printf("\nVal:  ");
	for(i=0; i<vert; i++)
		printf("%2d  ", sol[i]);
	printf("\n");
}

int random_l_h(int min, int max){
	return min + rand() % (max-min+1);
}

void sub(int a[], int b[], int n){
    int i;
    for(i=0; i<n; i++)
        a[i]=b[i];
}

float rand_01(){
	return ((float)rand())/RAND_MAX;
}

int cFit_local(int a[], int *mat, int vert){
	int total=0;
	int max;
	int i, j;

	for(i=0; i<vert; i++)
		for(j=0; j<vert;j++)
			if(*(mat+i*vert+j)==1 && i != j){
				total = abs(a[i]-a[j]);
				if(total > max)
					max = total;
			}
	return max;
}

int flip(){
	if ((((float)rand()) / RAND_MAX) < 0.5)
		return 0;
	else
		return 1;
}


