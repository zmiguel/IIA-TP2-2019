#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pEvol.h"
#include "util.h"

int main(int argc, char *argv[]){

    char        ficheiro[100];
	struct info_evol EA_param;
	pchrom      pop = NULL, parents = NULL;
	chrom       best_run, best_ever;
	int         gen_actual, r, runs, i, inv, mat[MAX_OBJ][MAX_OBJ];
	float       mbf = 0.0;
    int         opt1=0, opt2=0;

    if(argc != 3){
        printf("ERRO: numero de argumentos inválido\n%s <nome_ficheiro> <iteracoes>\n",argv[0]);
        exit(1);
    }

    strcpy(ficheiro,argv[1]);
    runs = atoi(argv[2]);

    prep_random();

    EA_param = init_data_evol(ficheiro, mat);
    fflush(stdin);
    while(opt1<=0 || opt1 >3){
        printf("Escolha forma de avaliacao:\n");
        printf("1 - Penalizacao\n");
        printf("2 - Reparacao 1\n");
        printf("3 - Reparacao 2\n");
        scanf("%d", &opt1);
    }
    fflush(stdin);
    while(opt2 <= 0 || opt2 >3){
        printf("\nEscolha o operador genetico:\n");
        printf("1 - Recombinacao 1 ponto\n");
        printf("2 - Recombinacao 2 pontos\n");
        printf("3 - Recombinacao Uniforme\n");
        scanf("%d", &opt2);
    }

    for (r=0; r<runs; r++){
		printf("Repeticao %d\n",r+1);
        // Geração da população inicial
		pop = init_pop_evol(EA_param);
        // Avalia a população inicial
		if (opt1 == 1){
			eval(pop, EA_param, mat);
		}else if (opt1 == 2){
			eval2(pop, EA_param, mat);
		}else if (opt1 == 3){
			eval3(pop, EA_param, mat);
		}else{
			printf("Erro na selecao da avaliacao\n");
		}
        // Aplicação do algoritmo trepa colinas para refinar a população inicial
        // Exercício 4.6(i)
		// trepa_colinas(pop, EA_param, mat);
		// Como ainda não existe, escolhe-se como melhor solução a primeira da população (poderia ser outra qualquer)
		best_run = pop[0];
        // Encontra-se a melhor solução dentro de toda a população
		best_run = get_best_evol(pop, EA_param, best_run);
        // Reserva espaço para os pais da população seguinte
		parents = malloc(sizeof(chrom)*EA_param.popsize);
        // Caso não consiga fazer a alocação, envia aviso e termina o programa
		if (parents==NULL)
		{
			printf("Erro na alocacao de memoria\n");
			exit(1);
		}
		// Ciclo de optimização
		gen_actual = 1;
		while (gen_actual <= EA_param.numGenerations)
		{
            // Torneio binário para encontrar os progenitores (ficam armazenados no vector parents)
			tournament(pop, EA_param, parents);
            // Torneio de k elementos, com k >= 2, para encontrar os progenitores (ficam armazenados no vector parents)
            // Exercício 4.5
            // tournament_geral(pop, EA_param, parents);
            // Aplica os operadores genéticos aos pais (os descendentes ficam armazenados na estrutura pop)
			if (opt2 == 1){
			genetic_operators(parents, EA_param, pop);
			}else if (opt2 == 2){
			genetic_operators2(parents, EA_param, pop);
			}else if (opt2 == 3){
			genetic_operators3(parents, EA_param, pop);
			}else{
				printf("Erro na selecao de operador genetico\n");
			}

            // Avalia a nova população (a dos filhos)
			if (opt1 == 1){
				eval(pop, EA_param, mat);
			}else if (opt1 == 2){
				eval2(pop, EA_param, mat);
			}else if (opt1 == 3){
				eval3(pop, EA_param, mat);
			}else{
				printf("Erro na selecao da avaliacao\n");
			}
            // Aplicação do algoritmo trepa colinas para refinar a população final
            // Exercício 4.6(iii)
            // trepa_colinas(pop, EA_param, mat);
            // Actualiza a melhor solução encontrada
			best_run = get_best_evol(pop, EA_param, best_run);
			gen_actual++;
		}
        // Aplicação do algoritmo trepa colinas para refinar a população final
        // Exercício 4.6(ii)
        // trepa_colinas(pop, EA_param, mat);
		// Contagem das soluções inválidas
		for (inv=0, i=0; i<EA_param.popsize; i++)
			if (pop[i].valido == 0)
				inv++;
		// Escreve resultados da repetição que terminou
		printf("\nRepeticao %d:", r+1);
		write_best_evol(best_run, EA_param);
		printf("\nPercentagem Invalidos: %f\n\n", 100*(float)inv/EA_param.popsize);
		mbf += best_run.fitness;
		if (r==0 || best_run.fitness < best_ever.fitness)
			best_ever = best_run;
        // Liberta a memória usada
		free(parents);
		free(pop);
	}
	// Escreve resultados globais
	printf("\n\nMBF: %f\n", mbf/r);
	printf("\nMelhor solucao encontrada");
	write_best_evol(best_ever, EA_param);
	return 0;
}

float eval_individual_penalizado(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v){
	int     i, j;
	float   sum_weight, sum_profit;
    int total = 0;
    int max = 0;

    for(i=0; i<d.numGenes; i++){
        for(j=0; j<d.numGenes;j++){
            if(mat[i-1][j-1]==1 && i != j){
                total = abs(sol[i]-sol[j]);
                if(total > max)
                    max = total;
                if (d.ro < (float)mat[i-1][j-1]/mat[i-1][j-1])
                    d.ro = (float)mat[i-1][j-1]/mat[i-1][j-1];
            }            
        }
    }

	if (max <= 0)
	{
        // Solução inválida
		*v = 0;
        //return -total; // Solucao com penalização
		return 0;
        //return max-d.ar*d.ro; // Solucao com penalização
	}
	else
	{
        // Solução válida
		*v = 1;
		return max;
	}
}

float eval_individual_reparado1(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v){
	int     i;
	float   sum_weight, sum_profit;

	sum_weight = sum_profit = 0;
	// Percorre todos os objectos
	for (i=0; i < d.numGenes; i++)
	{
        // Verifica se o objecto i esta na mochila
		if (sol[i] == 1)
		{
            // Actualiza o peso total
			sum_weight += mat[i][0];
            // Actualiza o lucro total
			sum_profit += mat[i][1];
		}
	}
	// Processo de reparacao
    while (sum_weight > d.ar)
    {
        // escolhe um objeto aleatoriamente
        i = random_l_h(0, d.numGenes-1);
        // Se esse objeto estiver na mochila, retira-o e ajusta os somat�rios do peso e lucro
        if (sol[i] == 1)
        {
            sol[i] = 0;
            sum_weight -= mat[i][0];
            sum_profit -= mat[i][1];
        }
    }
    *v = 1;
	return sum_profit;
}

float eval_individual_reparado2(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v){
	int     i, mv, pos;
	float   sum_weight, sum_profit;

	sum_weight = sum_profit = 0;
	// Percorre todos os objectos
	for (i=0; i < d.numGenes; i++)
	{
        // Verifica se o objecto i esta na mochila
		if (sol[i] == 1)
		{
            // Actualiza o peso total
			sum_weight += mat[i][0];
            // Actualiza o lucro total
			sum_profit += mat[i][1];
		}
	}
	// Processo de reparacao 2
    while (sum_weight > d.ar)
    {
        pos = -1;
        for (i=0; i < d.numGenes; i++)
        {
            if (sol[i] == 1)
            {
                if  (pos == -1 || mv > mat[i][1])
                {
                    mv = mat[i][1];
                    pos = i;
                }
            }
        }
        sol[pos] = 0;
        sum_weight -= mat[pos][0];
        sum_profit -= mat[pos][1];
    }
    *v = 1;
	return sum_profit;
}

void eval(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]){
	int i;

	for (i=0; i<d.popsize; i++){
		//pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(a)
		pop[i].fitness = eval_individual_penalizado(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(b)
		//pop[i].fitness = eval_individual_reparado1(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(c)
        //pop[i].fitness = eval_individual_reparado2(pop[i].p, d, mat, &pop[i].valido);
    }
}

void eval2(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]){
	int i;

	for (i=0; i<d.popsize; i++)
		//pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(a)
		//pop[i].fitness = eval_individual_penalizado(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(b)
		pop[i].fitness = eval_individual_reparado1(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(c)
        //pop[i].fitness = eval_individual_reparado2(pop[i].p, d, mat, &pop[i].valido);
}

void eval3(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]){
	int i;

	for (i=0; i<d.popsize; i++)
		//pop[i].fitness = eval_individual(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(a)
		//pop[i].fitness = eval_individual_penalizado(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(b)
		//pop[i].fitness = eval_individual_reparado1(pop[i].p, d, mat, &pop[i].valido);
        // Exercício 4.2(c)
        pop[i].fitness = eval_individual_reparado2(pop[i].p, d, mat, &pop[i].valido);
}

void tournament(pchrom pop, struct info_evol d, pchrom parents){
	int i, x1, x2;

	// Realiza popsize torneios
	for (i = 0; i < d.popsize; i++)
	{
		x1 = random_l_h(0, d.popsize - 1);
		do
			x2 = random_l_h(0, d.popsize - 1);
		while (x1 == x2);
		if (pop[x1].fitness < pop[x2].fitness) // Problema de maximizacao
			parents[i] = pop[x1];
		else
			parents[i] = pop[x2];
	}
}

void tournament_geral(pchrom pop, struct info_evol d, pchrom parents){
	int i, j, k, sair, best, *pos;

	pos = malloc(d.tsize * sizeof(int));
	// Realiza popsize torneios
	for (i = 0; i < d.popsize; i++)
	{
		// Seleciona tsize soluções diferentes para entrarem em torneio de seleção
		for (j = 0; j < d.tsize; j++)
		{
			do
			{
				pos[j] = random_l_h(0, d.popsize - 1);
				// Verifica se a nova posição escolhida é igual a alguma das outras posições escolhidas
				sair = 0;
				for (k = 0; k < j; k++)
				{
					if (pos[k] == pos[j])
						sair = 1;
				}
			} while (sair);
			// Guarda a posição da melhor solução de todas as que entraram em torneio
			if (j == 0 || pop[pos[j]].fitness > pop[pos[best]].fitness) // Problema de maximizacao
				best = j;
		}
		parents[i] = pop[pos[best]];
	}
	free(pos);
}

void genetic_operators(pchrom parents, struct info_evol d, pchrom offspring) {
	// Recombinação com um ponto de corte
	crossover(parents, d, offspring);
	// Recombinação com dois pontos de corte
	// Exercício 4.4(a)
	//recombinacao_dois_pontos_corte(parents, d, offspring);
	// Recombinação uniforme
	// Exercício 4.4(b)
	//recombinacao_uniforme(parents, d, offspring);
	// Mutação binária
	//mutation(offspring, d);
	// Mutação por troca
	// Exercício 4.3
	//mutacao_por_troca(offspring, d);
}

void genetic_operators2(pchrom parents, struct info_evol d, pchrom offspring) {
	// Recombinação com um ponto de corte
	//crossover(parents, d, offspring);
	// Recombinação com dois pontos de corte
	// Exercício 4.4(a)
	recombinacao_dois_pontos_corte(parents, d, offspring);
	// Recombinação uniforme
	// Exercício 4.4(b)
	//recombinacao_uniforme(parents, d, offspring);
	// Mutação binária
	mutation(offspring, d);
	// Mutação por troca
	// Exercício 4.3
	//mutacao_por_troca(offspring, d);
}

void genetic_operators3(pchrom parents, struct info_evol d, pchrom offspring) {
	// Recombinação com um ponto de corte
	//crossover(parents, d, offspring);
	// Recombinação com dois pontos de corte
	// Exercício 4.4(a)
	//recombinacao_dois_pontos_corte(parents, d, offspring);
	// Recombinação uniforme
	// Exercício 4.4(b)
	recombinacao_uniforme(parents, d, offspring);
	// Mutação binária
	mutation(offspring, d);
	// Mutação por troca
	// Exercício 4.3
	//mutacao_por_troca(offspring, d);
}

void crossover(pchrom parents, struct info_evol d, pchrom offspring){
	int i, j, point;

	for (i=0; i<d.popsize; i+=2)
	{
		if (rand_01() < d.pr)
		{
			point = random_l_h(0, d.numGenes-1);
			for (j=0; j<point; j++)
			{
				offspring[i].p[j] = parents[i].p[j];
				offspring[i+1].p[j] = parents[i+1].p[j];
			}
			for (j=point; j<d.numGenes; j++)
			{
				offspring[i].p[j]= parents[i+1].p[j];
				offspring[i+1].p[j] = parents[i].p[j];
			}
		}
		else
		{
			offspring[i] = parents[i];
			offspring[i+1] = parents[i+1];
		}
	}
}

void recombinacao_dois_pontos_corte(pchrom parents, struct info_evol d, pchrom offspring) {
	int i, j, point1, point2;

	for (i = 0; i < d.popsize; i += 2)
	{
		if (rand_01() < d.pr)
		{
			point1 = random_l_h(0, d.numGenes - 2);
			point2 = random_l_h(point1 + 1, d.numGenes - 1);
			for (j = 0; j < point1; j++)
			{
				offspring[i].p[j] = parents[i].p[j];
				offspring[i + 1].p[j] = parents[i + 1].p[j];
			}
			for (j = point1; j < point2; j++)
			{
				offspring[i].p[j] = parents[i + 1].p[j];
				offspring[i + 1].p[j] = parents[i].p[j];
			}
			for (j = point2; j < d.numGenes; j++)
			{
				offspring[i].p[j] = parents[i].p[j];
				offspring[i + 1].p[j] = parents[i + 1].p[j];
			}
		}
		else
		{
			offspring[i] = parents[i];
			offspring[i + 1] = parents[i + 1];
		}
	}
}

void recombinacao_uniforme(pchrom parents, struct info_evol d, pchrom offspring){
	int i, j;

	for (i = 0; i < d.popsize; i += 2)
	{
		if (rand_01() < d.pr)
		{
			for (j = 0; j < d.numGenes; j++)
			{
				if (flip() == 1)
				{
					offspring[i].p[j] = parents[i].p[j];
					offspring[i + 1].p[j] = parents[i + 1].p[j];
				}
				else
				{
					offspring[i].p[j] = parents[i + 1].p[j];
					offspring[i + 1].p[j] = parents[i].p[j];
				}
			}
		}
		else
		{
			offspring[i] = parents[i];
			offspring[i + 1] = parents[i + 1];
		}
	}
}

void mutation(pchrom offspring, struct info_evol d) {
	int i, j;

	for (i = 0; i < d.popsize; i++)
		for (j = 0; j < d.numGenes; j++)
			if (rand_01() < d.pm)
				offspring[i].p[j] = !(offspring[i].p[j]);
}

void mutacao_por_troca(pchrom offspring, struct info_evol d) {
	int i, pos1, pos2, aux;

	for (i = 0; i < d.popsize; i++)
		if (rand_01() < d.pm)
		{
			do
				pos1 = random_l_h(0, d.numGenes - 1);
			while (offspring[i].p[pos1] == 1);
			do
				pos2 = random_l_h(0, d.numGenes - 1);
			while (offspring[i].p[pos2] == 0);
			aux = offspring[i].p[pos1];
			offspring[i].p[pos1] = offspring[i].p[pos2];
			offspring[i].p[pos2] = aux;
		}
}

int my_score(int p[], int **mat, int l){
	int max = 0;
	for(int i=0; i<l; i++){
		for(int j = 0; j<l; j++){
			int cand = abs( p[i] - p[j] * mat[i][j]);
			if(cand>max) max =cand; 
		} 
	}
	
	return max; 
}

struct info_evol init_data_evol(char *filename, int mat[][MAX_OBJ]){
	struct  info_evol x;
	FILE    *f;
	int     i, j, lin, col;
	char lixo[100]; //variável que vai "tirar" o comentário inicial do ficheiro

	f = fopen(filename, "rt");
	if (!f)
	{
		printf("File not found\n");
		exit(1);
	}
	fscanf(f,"%99[^\n]", lixo); //Ler o comentário e adicionar ao lixo.

	// Atribuição dos parâmetros do problema
	x.popsize = 100; 
	x.pm = 0.01;
	x.pr = 1;
	x.tsize = 3;
	x.numGenerations = 2000;
	
	//x.numGenes = 100;
	//x.capacity = 250;

	fscanf(f, "%d %d %d", &x.numGenes , &x.numGenes, &x.ar);

	if (x.numGenes > MAX_OBJ)
	{
		printf("Number of itens is superior to MAX_OBJ\n");
		exit(1);
	}
	x.ro = 0.0;
	// Leitura dos dados do KSP (peso e lucro)
	for (i=0; i<x.numGenes; i++)
		for(j=0; j<x.numGenes; j++)
			mat[i][j] = 0;
	for(i=0; i<x.ar; i++){
		fscanf(f, " %d %d", &lin, &col);

		mat[lin-1][col-1]=1;
		mat[col-1][lin-1]=1;
	}
	fclose(f);
	// Devolve a estrutura com os parâmetros
	return x;
}

pchrom init_pop_evol(struct info_evol d){
	int     i, j;
	int flag[d.numGenes - 1];
	int aux;
	pchrom  indiv;
	indiv = malloc(sizeof(chrom)*d.popsize);
	if (indiv==NULL)
	{
		printf("Erro na alocacao de memoria\n");
		exit(1);
	}

	for(i=0;i<d.numGenes; i++)
		flag[i] = 0;

	for (i=0; i<d.popsize; i++)
	{
		for (j=0; j<d.numGenes; j++){
			do{
				aux = random_l_h(0, d.numGenes-1);
				indiv[i].p[j] = aux;
			} while(flag[aux] != 0);
			flag[aux] = 1;
		}
		for(j=0; j<d.numGenes; j++)
			flag[j] = 0;
	}
	return indiv;
}

chrom get_best_evol(pchrom pop, struct info_evol d, chrom best){
	int i;
	best.fitness = 9999;
	for (i=0; i<d.popsize; i++)
	{
		if (best.fitness > pop[i].fitness && pop[i].fitness > 0)
			best=pop[i];
	}
	return best;
}

void write_best_evol(chrom x, struct info_evol d){
	int i;

	printf("\nBest individual: %4.1f\n", x.fitness);
	for (i=0; i<d.numGenes; i++)
		printf("%d ", x.p[i]);
	putchar('\n');
}