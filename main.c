#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "dec.h"

int main(int argc, char *argv[]){
    int opt1=0, sair=0;
    char ficheiro[100],nitr[9];

    printf("IIA TP2\n\n");

    while(sair!=1){
        showMenu();
        scanf("%d",&opt1);
        switch (opt1)
        {
        case 1://pesquisa local
            printf("Indique o nome do ficheiro: ");
            scanf("%s",&ficheiro);
            printf("Indique o numero de iterações: ");
            scanf("%s",&nitr);
            execl("pLocal","./pLocal",ficheiro, nitr, (char *) NULL);
            break;
        case 2://pesquisa evolutiva
            printf("Indique o nome do ficheiro: ");
            scanf("%s",&ficheiro);
            printf("Indique o numero de iterações: ");
            scanf("%s",&nitr);
            execl("pEvol","./pEvol",ficheiro, nitr, (char *) NULL);
            break;
        case 3://pesquisa hybrida
            printf("Indique o nome do ficheiro: ");
            scanf("%s",&ficheiro);
            printf("Indique o numero de iterações: ");
            scanf("%s",&nitr);
            execl("pHybrid","./pHybrid",ficheiro, nitr, (char *) NULL);
            break;
        case 4://correr testes
            doTestes();
            break;
        case 9://sair
            sair = 1;
            break;
        default://sair por defeito
            sair = 1;
            break;
        } 
    }

    return 0;
}

void showMenu(){
    printf("1 - Pesquisa Local\n");
    printf("2 - Pesquisa Evolutiva\n");
    printf("3 - Pesquisa Hybrida\n");
    printf("4 - Correr testes\n");
    printf("9 - Sair!\n\n");
    printf("O que queres fazer? ");
}

void doTestes(){
    char ficheiro[100], line[128];
    printf("Ficheiro com info testes: ");
    scanf("%s", &ficheiro);

    char instname[100];
    char itt[9];
    char opt1[9];
    int mode;
    int val=1;

    FILE *f;
    f = fopen(ficheiro, "r+");
    if(f==NULL){
        printf("ficheiro nao encontrado!\n");
        exit(2);
    }

    while(!feof(f)){
        fscanf(f,"%d %s %s %s\n",&mode, &itt, &opt1, &instname);

        switch (fork())
        {
        case 0:             //no filho
            switch (mode)
            {
            case 1:
                execl("pLocal","./pLocal", instname, itt, opt1, "1", (char *)NULL);
                break;
            case 2:
                execl("pEvol","./pEvol", instname, itt, opt1, "1", (char *)NULL);
                break;
            case 3:
                execl("pHybrid","./pHybrid", instname, itt, opt1, "1", (char *)NULL);
                break;
            default:
                printf("ficheiro com dados errados!\n");
                break;
            }
            break;
        case 1:             //no pai
            continue;
        default:
            break;
        }
        
    }
    fclose(f);
    exit(10);
}