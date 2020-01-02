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
            execl("pEvol","./pLocal",ficheiro, nitr, (char *) NULL);
            break;
        case 3://pesquisa hybrida
            printf("Indique o nome do ficheiro: ");
            scanf("%s",&ficheiro);
            printf("Indique o numero de iterações: ");
            scanf("%s",&nitr);
            execl("pHybrid","./pLocal",ficheiro, nitr, (char *) NULL);
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
    printf("9 - Sair!\n\n");
    printf("O que queres fazer? ");
}