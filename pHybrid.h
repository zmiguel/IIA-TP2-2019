#define MAX_OBJ 1000

struct info_hybrid{
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
	// Capacidade da mochila
	int     capacity;
	// Número de gerações
    int     numGenerations;
};

// Individuo (solução)
typedef struct individual chrom, *pchrom;

struct individual{
    // Solução (objetos que estão dentro da mochila)
    int     p[MAX_OBJ];
    // Valor da qualidade da solução
	float   fitness;
    // 1 se for uma solução válida e 0 se não for
	int     valido;
};

void tournament(pchrom, struct info_hybrid, pchrom);
void tournament_geral(pchrom, struct info_hybrid, pchrom);
void genetic_operators(pchrom, struct info_hybrid, pchrom);
void genetic_operators2(pchrom, struct info_hybrid, pchrom);
void genetic_operators3(pchrom, struct info_hybrid, pchrom);
void crossover(pchrom, struct info_hybrid, pchrom);
void mutation(pchrom, struct info_hybrid);
void recombinacao_dois_pontos_corte(pchrom parents, struct info_hybrid d, pchrom offspring);
void recombinacao_uniforme(pchrom parents, struct info_hybrid d, pchrom offspring);
void mutacao_por_troca(pchrom, struct info_hybrid);

void eval(pchrom, struct info_hybrid, int mat[][1000]);
void eval2(pchrom, struct info_hybrid, int mat[][1000]);
void eval3(pchrom, struct info_hybrid, int mat[][1000]);
void eval4(pchrom, struct info_hybrid, int mat[][1000]);

void trepa_colinas(pchrom, struct info_hybrid, int mat[][1000]);

struct info_hybrid init_data_hybrid(char *filename, int mat[][1000]);
pchrom init_pop_hybrid(struct info_hybrid d);
chrom get_best_hybrid(pchrom pop, struct info_hybrid d, chrom best);
void write_best_hybrid(chrom x, struct info_hybrid d);