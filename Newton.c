#include <stdio.h>
#include <stdlib.h>
#include "fonction.c"


int main(){
    int taille, i, j;           // taille = nombre de points dans le tableau associant x et y
    double valeur;              // variable permettant de remplir x et y
    printf("Combien de valeurs connues ?:\n");
    scanf("%d",&taille);

    double *x = malloc(taille*sizeof(double));      //tableau des abscisses
    double *y = malloc(taille*sizeof(double));      //tableau des ordonnées


    for (i=0; i<taille; i++){                       //remplissage du tableau des x
        printf("Quelle valeur pour x[%d] : \n", i);
        scanf("%lf",&valeur);
        x[i] = valeur;
    }
    for (i=0; i<taille; i++){                       //remplissage du tableau des y
        printf("Quelle valeur pour y[%d] : \n", i);
        scanf("%lf",&valeur);
        y[i] = valeur;
    }

    double *result, *result2, *resfinal = malloc(taille*sizeof(double));  //definition de plusieurs pointeurs qui serviront à stocker certains résultats pour les réutiliser plus tard dans les calculs
                                                                          //result, result 2 existent pour des calculs intermédiares / resfinal pour le resultat final a savoir le polynome de Newton
    result = newton(x,y,taille);                                          //newton() calcule les coefficiants b0, ..., bn
    printf("\nOn obtient les coefficients b0, ..., bn : \n");             //le resultat est un tableau pointé par *result
    for (i=0; i<taille; i++){
        printf("%.15lf ", result[i]);
        resfinal[i] = 0;                                                  //Etant donne que l'on va ajoute des valeurs a resfinal plusieurs fois on doit l'initialiser à 0
    }


    for (j=1; j<taille; j++){
        result2 = polynomeDeNewton(result, x,taille,j);                   //polynomedeNewton() sert a organiser la formule du polynome dans un tableau : b0 (x-x0) + ... + bn (x-x0)...(x-xn)
        if  (j>1){                                                        //il calcule d'abord b0(x-x0) dans result2, ensuite b1(x-x0)(x-x1)...
            result2 = multi(result2,j+1);
        }
        result2 = reverse(result2, j+1);                                  //on inverse le tableau juste pour que a0 soit a la case d'indice 0, a1 a 1...
        for (i=0; i<j+1; i++){
            resfinal[i] += result2[i];                                    //resfinal accumule a chaque fois les valeurs de result2 c'est à dire : b0(x-x0) + ...
        }
    }

    printf("\nVoici le polynome obtenu apres l'interpolation de Newton :\n");
    resfinal[0] +=y[0];                                                  //on rajoute a la fin y[0]
    resfinal = reverse(resfinal, taille);                                //on inverse le tableau pour qu'il soit conforme a notre sens d'écriture
    for (i=0; i<taille; i++){
        printf("%.15lf ", resfinal[i]);                                     //affichage du polynome dans la console
    }
    printf("\n");

    free(x);                                                             //Partie libération de la mémoire
    free(y);
    free(resfinal);

    x = NULL;
    y = NULL;
    result = NULL;
    result2 = NULL;
    resfinal = NULL;

    return 0;
}
