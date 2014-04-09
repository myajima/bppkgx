/* file mymod.c */
#include <R.h>
static double parms[16];
#define p01 parms[0]
#define p02 parms[1]
#define p03 parms[2]
#define p04 parms[3]
#define p05 parms[4]
#define p06 parms[5]
#define p07 parms[6]
#define p08 parms[7]
#define p09 parms[8]
#define p10 parms[9]
#define p11 parms[10]
#define p12 parms[11]
#define p13 parms[12]
#define p14 parms[13]
#define p15 parms[14]
#define p16 parms[15]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=16;
odeparms(&N, parms);
}
/* Derivatives and 1 output variable */
// void derivs (int *neq, double *t, double *y, double *ydot,
//                 double *yout, int *ip)
// {
// if( ip[0] <1 ) error("nout should be at least 1");
// ydot[0] = -k1*y[0] + k2*y[1]*y[2];
// ydot[2] = k3 * y[1]*y[1];
// ydot[1] = -ydot[0]-ydot[2];
// yout[0] = y[0]+y[1]+y[2];
// }


//EHRpk  <- function  (  t, y, p ) {
void ehrpk (int *neq, double *t, double *y, double *ydot,
                double *yout, int *ip)
{
    //if( ip[0] <1 ) error("nout should be at least 1");
    if(t[0] <= 1.5) ydot[0] = p03*y[1] - (p02+p04+p08+p01)*y[0] + p16;
    if(t[0] >  1.5) ydot[0] = p03*y[1] - (p02+p04+p08+p01)*y[0];
                 ydot[1] = p02*y[0] -  p03*y[1] ;
                 ydot[2] = p04*y[0] - (p05+p06+p10)*y[2] + p14*y[6]; 
                 ydot[3] = p06*y[2] -           p07*y[3];
                 ydot[4] = p08*y[0] -           p09*y[4];
    if( t[0] > p11 && t[0] < ( p11 + 1.0 ) ){ 
                 ydot[5] = p10*y[2] - ( p12 + p13 )*y[5];
                 ydot[6] =            ( p12 + p13 )*y[5] - p14*y[6];
    } else {
                 ydot[5] = p10*y[2] -  p12*y[5];
                 ydot[6] =             p12*y[5]          - p14*y[6];
    }
    //yout = 0;
}







// /* The Jacobian matrix */
// void jac(int *neq, double *t, double *y, int *ml, int *mu,
// double *pd, int *nrowpd, double *yout, int *ip)
// {
// pd[0] = -k1;
// pd[1] = k1;
// pd[2] = 0.0;
// pd[(*nrowpd)] = k2*y[2];
// pd[(*nrowpd) + 1] = -k2*y[2] - 2*k3*y[1];
// pd[(*nrowpd) + 2] = 2*k3*y[1];
// pd[(*nrowpd)*2] = k2*y[1];
// pd[2*(*nrowpd) + 1] = -k2 * y[1];
// pd[2*(*nrowpd) + 2] = 0.0;
// }
// /* END file mymod.c */



//# Define Jacobian (optional) -------------------------------------------------------- #
//jac1 <- function(t, y, p){
void jac(int *neq, double *t, double *y, int *ml, int *mu,
double *pd, int *nrowpd, double *yout, int *ip)
{
    pd[              0] =  -(p02+p04+p08+p01);
    pd[              1] =  p02;
    pd[              2] =  p04;
    pd[              3] =  0;
    pd[              4] =  p08;
    pd[              5] =  0;
    pd[              6] =  0;
    pd[  (*nrowpd)    ] =   p03;
    pd[  (*nrowpd) + 1] =  -p03;
    pd[  (*nrowpd) + 2] =  0;
    pd[  (*nrowpd) + 3] =  0;
    pd[  (*nrowpd) + 4] =  0;
    pd[  (*nrowpd) + 5] =  0;
    pd[  (*nrowpd) + 6] =  0;
    pd[2*(*nrowpd)    ] =  0;
    pd[2*(*nrowpd) + 1] =  0;
    pd[2*(*nrowpd) + 2] =  -(p05+p06+p10);
    pd[2*(*nrowpd) + 3] =  p06;
    pd[2*(*nrowpd) + 4] =  0;
    pd[2*(*nrowpd) + 5] =  p10;
    pd[2*(*nrowpd) + 6] =  0;
    pd[3*(*nrowpd)    ] =  0;
    pd[3*(*nrowpd) + 1] =  0;
    pd[3*(*nrowpd) + 2] =  0;
    pd[3*(*nrowpd) + 3] =  -p07;
    pd[3*(*nrowpd) + 4] =  0;
    pd[3*(*nrowpd) + 5] =  0;
    pd[3*(*nrowpd) + 6] =  0;
    pd[4*(*nrowpd)    ] =  0;
    pd[4*(*nrowpd) + 1] =  0;
    pd[4*(*nrowpd) + 2] =  0;
    pd[4*(*nrowpd) + 3] =  0;
    pd[4*(*nrowpd) + 4] =  -p09;
    pd[4*(*nrowpd) + 5] =  0;
    pd[4*(*nrowpd) + 6] =  0;
    pd[5*(*nrowpd)    ] =  0;
    pd[5*(*nrowpd) + 1] =  0;
    pd[5*(*nrowpd) + 2] =  0;
    pd[5*(*nrowpd) + 3] =  0;
    pd[5*(*nrowpd) + 4] =  0;
    if(t[0] > p11 && t[0] < (p11+1)){ 
        pd[5*(*nrowpd) + 5] =  -(p12 + p13);
        pd[5*(*nrowpd) + 6] =   (p12 + p13);
    }else{
        pd[5*(*nrowpd) + 5] =  -p12;
        pd[5*(*nrowpd) + 6] =   p12;
    }
    pd[6*(*nrowpd)    ] =  0;
    pd[6*(*nrowpd) + 1] =  0;
    pd[6*(*nrowpd) + 2] =  p14;
    pd[6*(*nrowpd) + 3] =  0;
    pd[6*(*nrowpd) + 4] =  0;
    pd[6*(*nrowpd) + 5] =  0;
    pd[6*(*nrowpd) + 6] = -p14;

}