typedef struct {            //Structure permettant de stocker : la taille d'un polyome, le polynome, et le coefficient dans Lagrange Yi/(x-x1)(x-x2)...(x-xj)
    int taille;
    double *polynome;
    double coeff;
}PolyBox;

PolyBox lagrange (double *x, double *y, int n, int ligne);

double *range (PolyBox retour);

PolyBox multiP(double *p, double*q, int n, int m);

double *multi (double *poly, int n);

double *multi (double *poly, int n);

double *polynomeDeNewton(double *coeffB, double *x, int n, int degre);

double *reverse(double *poly, int n);
