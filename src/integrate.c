#include "fem.h"


double femJacDPhi(femDiscrete* space, femIntegration* rule, double* x, double* y, double* dphidx, double* dphidy){
    double* xsi = rule->xsi;
    double* eta = rule->eta;
    double* phi = calloc(space->n,sizeof(double));
    double* dphidxsi = calloc(space->n,sizeof(double));
    double* dphideta = calloc(space->n,sizeof(double));
    double dxdxsi = 0.0;
    double dxdeta = 0.0;
    double dydxsi = 0.0; 
    double dydeta = 0.0;
    for(int i=0; i<rule->n; i++){
        femDiscretePhi2(space,xsi[i],eta[i],phi);
        femDiscreteDphi2dx(space,xsi[i],eta[i],dphidxsi,dphideta);
    }
    for (int i = 0; i < space->n; i++) {  
        dxdxsi += x[i]*dphidxsi[i];       
        dxdeta += x[i]*dphideta[i];   
        dydxsi += y[i]*dphidxsi[i];   
        dydeta += y[i]*dphideta[i]; }
    double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

    for (int i = 0; i < space->n; i++) {    
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
    }
    return jac; 
}

double femIntegrate(femIntegration* rule, double* points, double jac){
    double I = 0;
    for(int i=0; i<rule->n; i++){
        I += rule->weight[i]*points[i];
    }
    return I*jac;
};

femLocal(femDiscrete* space, double **A, double *dx, double *dy, double weight, double jac, double* coeff){
    int n = space->n;
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            A[i][j] += jac * weight * (dx[i]*coeff[0]*dx[j] + dy[i]*coeff[2]*dy[j]);
        }
        for(int j=(n-2); j<n; j++){
            A[i][j] += jac * weight * (dx[i]*coeff[1]*dy[j] + dy[i]*coeff[2]*dx[j]);
        }
        for(int j=0; j<2; j++){
            A[i+(n-2)][j] += jac * weight * (dy[i+(n-2)]*coeff[1]*dx[j] + dx[i+(n-2)]*coeff[2]*dy[j]);
        }
        for(int j=(n-2); j<n; j++){
            A[i+(n-2)][j] += jac * weight * (dy[i+(n-2)]*coeff[0]*dy[j] + dx[i+(n-2)]*coeff[2]*dx[j]);
        }
    }
}

double *matrixSolve(double **A, double *B, int size)
{ 
    int i, j, k;
    double factor;
    /* Gauss elimination */
    for (k=0; k < size; k++) {
        if ( A[k][k] == 0 ) Error("zero pivot");
        for (i = k+1 ; i < size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++)
            A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    /* Back-substitution */
    for (i = (size)-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
        factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    return(B);
}

