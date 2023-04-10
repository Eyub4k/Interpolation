#include <stdio.h>
#include <stdlib.h>
#include "fonction.h"

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

PolyBox lagrange (double *x, double *y, int n, int ligne){      //Entrée : tableau des x et des y, le nombre de points n et la ligne i correspondant au i de la somme des Yi*Li dans la formule du cours
    PolyBox retour;                                             //Renvoie un PolyBox une structure détaillé dans fonction.h
    double coeff = 1;                                           //Coeff est l'élément neutre de la multiplication soit 1, on le multipliera plusieurs fois avec lui même
    double *poly = malloc (n*sizeof(double));
    int i, j, booleen = 0;                                      //booleen nous servira a sauter le terme xi a ne pas prendre dans Li
    for (j=0; j<n; j++){
        if (ligne!=j){                                          //si nous ne sommes pas sur une diagonale, si i!=j dans la formule du cours
            if (booleen == 1){                                  //si booleen a ete active, on saute une case de xj, dans ce cas la on doit donc remplir la case j-1 de poly sinon elle demeurera vide
                poly[j-1] = x[j];
            }
            else{
                poly[j] = x[j];                                 //si booleen n'a pas ete active alors on procède comme d'habitude
            }
            coeff = coeff * ( x[ligne] - x[j] );                //calcul du coeff qui correspond à Yi/(x-x0)...(x-xj)
        }
        else if((j+1)<n){                                       //si i==j et que j+1<n c'est à dire que l'on a pas dépassé le tableau
            j++;                                                //on saute la case xj donc j prend +1
            poly[j-1] = x[j];
            booleen = 1;                                        //booleen s'active
            coeff = coeff * ( x[ligne] - x[j] );
        }
    }
    coeff = y[ligne]/coeff;                                     //calcul final de coeff
    retour.coeff = coeff;                                       //rangement des variable dans la structure
    retour.polynome = poly;
    retour.taille = n;
    return retour;                                              //retourne une PolyBox
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double *range (PolyBox retour){                                 //amenage le polynome obtenu par lagrange pour les multiplier
    int i, n = retour.taille;                                   //n prend retour.taille c'est à dire la taille du polynôme obtenu par lagrange()
                                                                // C'est sur ce modèle que tout notre programme se base pour la gestion des polynomes :
    double *resultat = malloc (2*(n-1)*sizeof(double));         // on double la taille-1 du tableau : (x-1)(x+2) sera implemente ainsi p[X] = (1, -1) et  q[X] = (1, 2) dans un tableau (1, -1, 1, 2)
    for (i=0; i<(2*n-1); i++){
        if (i%2==0){                                            //si nous sommes sur une case pair alors elle correspond à x, donc la case prend 1
            resultat[i] = 1;
        }
        else{
            resultat[i] = -retour.polynome[i/2];
        }
    }
    return resultat;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

PolyBox multiP(double *p, double*q, int n, int m){              //Largement inspirée de la fonction multiplication de deux polynome dans un autre cours Mathématiques à l'usage des informaticiens
    int i, j, l;
    PolyBox z;
    if ( n == 1 ){                                              //Soient 2 polynome P et Q de degré n et m
        l = m;                                                  //Si n=1 alors degre(p.q)=m vice-versa
    }                                                           // Si n et m sont supérieurs a 1 alors degre(p.q)=m+n;
    else if (m == 1){
        l = n;
    }
    else if (m > 1 && n >1){
        l = m+n-1;
    }
    double *result = malloc( (l*sizeof(double)));               //Algorithme du cours de Mathématiques à l'suage des Informaticiens
    for (i=0; i<l; i++){
        result[i] = 0;}

    for (i=0; i<n; i++){
        for (j=0; j<m; j++){
            result[i+j] = result[i+j] + (p[i]*q[j]);
        }
    }
    z.polynome = result;                                        //polynome stocke dans PolyBox
    z.taille = l;                                               //idem pour le degre du polynome
    return z;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double *multi (double *poly, int n){                            //fonction qui choisit les polynomes a multiplier dans L avec L0, L1, ...
    PolyBox z;
    int nprim = n;                                              //taille de tableau modulable
    double *result = malloc (nprim*sizeof(double));
    double premier[2] = {0};                                    //on va multiplier les polynomes 2 a 2, c'est à dire que (x-1)(x-2)(x-5)(x-6), on va d'abord faire P1*P2 ensuite le rezsultat de P1*P2 *P3 ...
    double deuxieme[2] = {0};                                   //Ils sont de taille 2 car les polynomes sont toujours de taille 2 (x- ?)(x-?)...
    int i,j,compteur = 0;
    for (i=0; i<n; i++){                                        //on initialise result à n'importe quelle valeur pour pas qu'il soit vide
        result[i] = 2;
    }
    if (n>3){
        n = n + (n-3);
    }
    for (i=0; i<n; i+=2){                                       //i avance par pas de 2 car les polynomes sont de tailles 2

            if (compteur%2==0 && i<3){
                premier[0] = poly[i];
                premier[1] = poly[i+1];
                compteur++;
            }
            else if (compteur !=0 && i<3){
                deuxieme[0] = poly[i];
                deuxieme[1] = poly[i+1];
                compteur++;
            }
            else if (i>=4){

                deuxieme[0] = poly[i];
                deuxieme[1] = poly[i+1];
                compteur++;
            }

            if (i==2){
                z = multiP(premier,deuxieme, 2, 2);
                nprim = z.taille;
                result = z.polynome;

            }
            else if(i%2==0 && i!=0 && i!=n-2){
                z = multiP(deuxieme,result, 2, nprim);
                nprim = z.taille;
                result = z.polynome;
            }
            else if (i==n-2 && n%4!=0){
                z = multiP(deuxieme,result, 2, nprim);
                nprim = z.taille;
                result = z.polynome;

            }
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double *newton(double *x, double *y, int n){                //newton renvoie les coefficiants b0, b1, ...
    int i, j, k;
    int nprim = n-1;                                        //taille de tableau modulable
    double *inter = malloc(nprim*sizeof(double));
    double *inter2 = malloc(nprim*sizeof(double));
    double *result = malloc(n*sizeof(double));

    printf("\nx[ ] = ");                        //affichage de x
    for (i=0; i<n; i++){
        printf("%.15lf ", x[i]);
    }

    printf("\ny[ ] = ");                        //affichage de y
    for (i=0; i<n; i++){
        printf("%.15lf ", y[i]);
    }

    printf("\n\nVoici le tableau des differences divisees :\n");
    result[0] = y[0];                                       //algorithme du cours pour le calcul des differences divisées
    for (j=0; j<n; j++){
        nprim = n-j-1;
        for (k=0; k<nprim; k++){
            if (j==0){
                inter[k] = ( y[k+1] - y[k] ) / ( x[k+j+1] - x[k] );
                printf("%.15lf ",inter[k]);
                result[j+1] = inter[0];
            }
            else if (j%2!=0){
                inter2[k] = ( inter[k+1] - inter[k] ) / ( x[k+j+1] - x[k] );
                printf("%.15lf ",inter2[k]);
                result[j+1] = inter2[0];
            }
            else{
                inter[k] = ( inter2[k+1] - inter2[k] ) / ( x[k+j+1] - x[k] );
                printf("%.15lf ",inter[k]);
                result[j+1] = inter[0];
            }
        }
        printf("\n");
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double *polynomeDeNewton(double *coeffB, double *x, int n, int degre){      //arrange la formule du polynome de Newton b0 + b1(x-x1) +... dans un tableau selon le  même modele que Lagrange
    int i, compteur = 1;                                                    //(x-x1) sera implmente ainsi (1, -x1), le coefficient b sera directement distribué dans une paranthese ainsi on aura, b(x-x1) -> (b*1, b*x1)
    double *result = malloc( (2*degre) * sizeof(double) );

    for (i=0; i<(2*degre); i++){                                            //la taille du tableau s'en retrouve donc doublée
        if (i%2==0){                                                        //si on se trouve sur une case d'indice paire alors elle prend soit coeffB si on se situe avant la 1ere paranthese comme dans l'exemple au dessus
            if (i<2){
                result[i] = coeffB[degre];                                  //sinon elle prend 1
            }
            else{
                result[i] = 1;
            }
        }
        else{
            if (i<2){                                                       //si la case est d'indice impair alors si on se situe avant la 1ere paranthese on multiplie toujours le terme xi par b
                result[i] = coeffB[degre] * -x[i-compteur];
                compteur++;
            }
            else{
                result[i] = -x[i-compteur];                                 //sinon elle prend juste -xi
                compteur++;
            }
        }
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

double *reverse(double *poly, int n){                                       //fonction qui inverse l'ordre des cases dans un tableau
    int i;
    double *result = malloc(n*sizeof(double));
    for (i=0; i<n; i++){
        result[i] = poly[n-i-1];
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

