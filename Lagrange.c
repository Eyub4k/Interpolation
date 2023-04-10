#include <stdio.h>
#include <stdlib.h>
#include "fonction.c"

int main(){

    int i, j, n;                                    //n le nombre de points dans le tableau
    double valeur;                                  //variable permettant de remplir x et y
    PolyBox poly;                                   //poly de type PolyBox : structure détaille dans fonction.h
    printf("Combien de points ?\n");
    scanf("%d", &n);

    double *result, *result2, *resfinal;   //idem que dans newton.c
    double *x = malloc(n*sizeof(int));              //tableau des x
    double *y = malloc(n*sizeof(int));              //tableau des x
    printf("\n");

    for (i=0; i<n; i++){
        printf("Quelle valeur pour x[%d] : \n", i);     //idem que dans newton.c
        scanf("%lf", &valeur);
        x[i] = valeur;
    }
    for (i=0; i<n; i++){
        printf("Quelle valeur pour y[%d] : \n", i);
        scanf("%lf", &valeur);
        y[i] = valeur;
    }

    printf("\nx[ ] = ");                             //affichage des x
    for (i=0; i<n; i++){
        printf("%.15lf ", x[i]);
    }
    printf("\ny[ ] = ");                             //affichage des y
    for (i=0; i<n; i++){
        printf("%.15lf ", y[i]);
    }

    for (i=0; i<n; i++){
        poly = lagrange(x,y,n,i);                   //la PolyBox poly stocke une PolyBox renvoye par lagrange() qui calcule L0, ,L1, L2, ...
        result2 = range(poly);                      //amenage le polynome pour l'utiliser dans la fonction multi()
        if (n<=2){
            result = result2;
        }
        else{
            result = multi(result2,n);
        }
        for (j=0; j<n; j++){
            result[j] = result[j] * poly.coeff;     //TRES IMPORTANT : poly correspond à (x-x1)(x-x2)...(x-xj) et coeff ) yi/(x-x1)(x-x2)...(x-xj)
        }
        if (i==0){
            resfinal = result;                      //initialisation du resultat final
        }
        else{
            for (j=0; j<n; j++){
                resfinal[j] = resfinal[j] + result[j];  //resfinal accumule toutes les valeurs, il effectue les sommes de y0L0, y1L1
            }
        }
    }
    printf("\n\nVoici le polynome obtenu apres l'interpolation de Lagrange :\n");
    for (i=0; i<n; i++){
        printf("%.15lf ", resfinal[i]);                //affichage du polynome de Lagrange
    }
    printf("\n");
    free(result);                                   //Partie libération de la mémoire
    free(result2);
    free(resfinal);
    free(x);
    free(y);
    result = NULL;
    result2 = NULL;
    resfinal = NULL;
    x= NULL;
    y = NULL;



}



















