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
    int         opt1=0, tofile=0;

    if(argc <= 3){
        printf("ERRO: numero de argumentos inválido\n%s <nome_ficheiro> <iteracoes> [opt] [tofile (1)]\n",argv[0]);
        exit(1);
    }

	if(argc == 4){
        opt1 = atoi(argv[3]);
    }
    if(argc == 5){
        opt1 = atoi(argv[3]);
        tofile=1;
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
			genetic_operators(parents, EA_param, pop);

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
	if(tofile==1){
		char fname[100];
        sprintf(fname,"./out/pEvol_%sitt_%s_opt%d",argv[2],argv[1],opt1);
        FILE *f;
        f = fopen(fname,"wt");
        if(f==NULL){
            printf("nao abri o ficheiro\n");
        }
		fprintf(f,"\n\nMBF: %f\n", mbf/r);
		fprintf(f,"\nBest individual: %4.1f\n", best_ever.fitness);
		fclose(f);
	}else{
		printf("\n\nMBF: %f\n", mbf/r);
		printf("\nMelhor solucao encontrada");
		write_best_evol(best_ever, EA_param);
	}
	
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

	if (max >= d.numGenes)
	{
        // Solução inválida
		*v = 0;
        //return -total; // Solucao com penalização
		//return 0;
        return -(max-d.ar*d.ro); // Solucao com penalização
	}
	else
	{
        // Solução válida
		*v = 1;
		return max;
	}
}

float eval_individual_reparado1(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v){
	int     i, j;
	float   sum_weight, sum_profit;
    int total = 0;
    int max = 0;
    int max_temp;
    max_temp = max;

    for(i=0; i<d.numGenes; i++){
        for(j=0; j<d.numGenes;j++){
            if(mat[i-1][j-1]==1 && i != j){
                total = abs(sol[i]-sol[j]);
                if(total > max)
                    max = total;
            }            
        }
    }

	// Processo de reparacao
    while (max >= d.numGenes){
        // escolhe um objeto aleatoriamente
        i = random_l_h(0, d.numGenes-1);
        // Se esse objeto estiver na mochila, retira-o e ajusta os somatórios do peso e lucro
        for(j = 0;j<d.numGenes;j++){
            if(mat[i-1][j-1] == 1 && i != j){
                total = abs(sol[i]-sol[j]);
                if(total > max && total < max_temp)
                    max = total;
            }
        }
    }
    *v = 1;
	return max;
}

float eval_individual_reparado2(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v){
	int     i, j, j_max, i_max, mv, pos;
	float   sum_weight, sum_profit;
    int total = 0;
    int max = 0;
    int max_temp = max;
    int rdm, tmp;

    for(i=0; i<d.numGenes; i++){
        for(j=0; j<d.numGenes;j++){
            if(mat[i-1][j-1]==1 && i != j){
                total = abs(sol[i]-sol[j]);
                if(total > max)
                    max = total;
                    i_max = i;
                    j_max = j;
            }            
        }
    }

	// Processo de reparacao 2
    while (max >= d.numGenes)
    {
        // Se esse objeto estiver na mochila, retira-o e ajusta os somatórios do peso e lucro
        for(j = 0;j<d.numGenes;j++){
            //Gera um número diferente que não esteja ocupado
            do{
                rdm = random_l_h(0, d.numGenes-1);
            }while(rdm == i_max || rdm == j_max);

            tmp = sol[j_max];
            sol[j_max] = sol[rdm];
            sol[rdm] = tmp;

            if(mat[i_max-1][j-1] == 1 && i != j){
                total = abs(sol[i_max]-sol[j_max]);
                if(total > max && total < max_temp)
                    max = total;
            }
        }
    }
    *v = 1;
	return max;
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
	for (i=0; i<d.popsize;i++)
	{
		x1 = random_l_h(0, d.popsize-1);
		do
			x2 = random_l_h(0, d.popsize-1);
		while (x1==x2);
		if (pop[x1].fitness < pop[x2].fitness)	// Problema de minimização
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

void crossover(pchrom parents, struct info_evol d, pchrom offspring){
	int i, j, point;
	int flag[d.numGenes];
	pchrom pai, mae;
	int mask_pai[d.numGenes], mask_mae[d.numGenes];

	for (i = 0; i < d.popsize; i += 2){
		if (rand_01() < d.pr){
			pai = &parents[i];
			mae = &parents[i + 1];

			//Reiniciar a flag (definir todos os números como disponíveis)
			for (j = 0; j < d.numGenes; j++){
				flag[j] = 0;
			}
			//Gerar aleatóriamente uma mask para o pai e a inversa para a mãe (ex 11010 | 00101) -> Testado e a funcionar
			for (j = 0; j < d.numGenes; j++){
				mask_pai[j] = flip();
				mask_mae[j] = (int)!mask_pai[j];
			}
			//Se a mae tiver uma fitness mais baixa, trocar (pai tem sempre a fitness mais baixa)
			if (pai->fitness < mae->fitness){
				const pchrom temp = pai;
				pai = mae;
				mae = temp;
			}
			//Offspring herda valores do pai, conforme a mask_pai
			//Flag = 1 para todos os valores usados, impedindo assim repetidos.
			for (j = 0; j < d.numGenes; j++){
				if (mask_pai[j]){
					offspring[i].p[j] = pai->p[j];
					flag[pai->p[j]] = 1;
				}
			}
			//Offspring herda valores restantes da mãe, conforme a mask_mae
			//Caso o valor herdade da mãe seja repetido, continua a procurar por um valor disponível.
			for (j = 0; j < d.numGenes; j++){
				if (mask_mae[j]){
					int k = mae->p[j];
					while (flag[k]){
						k++;
						k = k % d.numGenes;
					}
					flag[k] = 1;
					offspring[i].p[j] = k; 
				}
			}
		}else{
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
	x.popsize = 150; 
	x.pm = 0.01;
	x.pr = 0.7;
	x.tsize = 3;
	x.numGenerations = 2500;
	
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