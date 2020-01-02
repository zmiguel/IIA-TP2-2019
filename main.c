#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "dec.h"

int main(int argc, char *argv[]){
    int opt1=0, sair=0;

    printf("IIA TP2\n\n");

    while(sair!=1){
        showMenu();
        scanf("%d",&opt1);
        switch (opt1)
        {
        case 1://pesquisa local
            execlp("pLocal","./pLocal","filename_goes_here","number os runs goes here", (char *) NULL);
            break;
        case 2://pesquisa evolutiva

            break;
        case 3://pesquisa hybrida
            
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