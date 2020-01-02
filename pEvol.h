#define MAX_OBJ 1000

struct info_evol {
    // Tamanho da população
    int     popsize;
    // Probabilidade de mutação
    float   pm;
    // Probabilidade de recombinação
    float   pr;
    // Tamanho do torneio para seleção do pai da próxima geração
	int     tsize;
	// Constante para avaliação com penalização
	float   ro;
	// Número de objetos que se podem colocar na mochila
    int     numGenes;
	// Número de arestas
	int     ar;
	// Número de gerações
    int     numGenerations;
};

typedef struct individual chrom, *pchrom;

struct individual{
    // Solução (objetos que estão dentro da mochila)
    int     p[MAX_OBJ];
    // Valor da qualidade da solução
	float   fitness;
    // 1 se for uma solução válida e 0 se não for
	int     valido;
};

float eval_individual_penalizado(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v);
float eval_individual_reparado1(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v);
float eval_individual_reparado2(int sol[], struct info_evol d, int mat[][MAX_OBJ], int *v);
void eval(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]);
void eval2(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]);
void eval3(pchrom pop, struct info_evol d, int mat[][MAX_OBJ]);
void tournament(pchrom pop, struct info_evol d, pchrom parents);
void tournament_geral(pchrom pop, struct info_evol d, pchrom parents);
void genetic_operators(pchrom parents, struct info_evol d, pchrom offspring);
void genetic_operators2(pchrom parents, struct info_evol d, pchrom offspring);
void genetic_operators3(pchrom parents, struct info_evol d, pchrom offspring);
void crossover(pchrom parents, struct info_evol d, pchrom offspring);
void recombinacao_dois_pontos_corte(pchrom parents, struct info_evol d, pchrom offspring);
void recombinacao_uniforme(pchrom parents, struct info_evol d, pchrom offspring);
void mutation(pchrom offspring, struct info_evol d);
void mutacao_por_troca(pchrom offspring, struct info_evol d);
int my_score(int p[], int **mat, int l);

struct info_evol init_data_evol(char *filename, int mat[][MAX_OBJ]);
pchrom init_pop_evol(struct info_evol d);
chrom get_best_evol(pchrom pop, struct info_evol d, chrom best);
void write_best_evol(chrom x, struct info_evol d);