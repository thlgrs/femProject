#include "newProject.h"


void systemAssemble(femProblem *theProblem)
{
    femSystem* theSystem = theProblem->system;
    femGeo* theGeometry = theProblem->geometry;
    femDiscrete *theSpace = theProblem->space;
    femIntegration *theRule = theProblem->rule;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int map[4],mapX[4],mapY[4];
    int nLocalNode = theGeometry->theElements->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A;
    double *B;
    switch(theSystem->solverType){
        case 0:
            A = theSystem->fullSystem->A;
            B  = theSystem->fullSystem->B;
        case 1:
            A = theSystem->bandSystem->A;
            B  = theSystem->bandSystem->B;
        case 3:
            A = theSystem->fullSystem->A;
            B  = theSystem->fullSystem->B;
    }


    for(int iElem=0; iElem<theGeometry->theElements->nElem; iElem++){
        
        for (int j=0; j < nLocalNode; j++) {
            map[j]  = theGeometry->theElements->elem[iElem*nLocalNode+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theGeometry->theElements->nodes->X[map[j]];
            y[j]    = theGeometry->theElements->nodes->Y[map[j]];
        } 
        for (int iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (int i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (int i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi)/jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta)/jac; }
            for (int i = 0; i < theSpace->n; i++) { 
                for(int j = 0; j < theSpace->n; j++) {
                    switch(theProblem->planarStrainStress<2){
                        case TRUE:
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                                    dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                                    dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                                    dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                                    dphidx[i] * c * dphidx[j]) * jac * weight;
                            
                        case FALSE:
                            A[mapX[i]][mapX[j]] += jac * weight * ( dphidx[i] * a * x[i] * dphidx[j] + 
                                                                    dphidy[i] * c * x[i] * dphidy[j] +
                                                                    phi[i] * ( b * dphidx[j] + a * phi[j]/x[i]) +
                                                                    dphidx[i] * b * phi[j]);                                                                                            
                            A[mapX[i]][mapY[j]] += jac * weight * ( dphidx[i] * b * x[i] * dphidy[j] + 
                                                                    dphidy[i] * c * x[i] * dphidx[j] +
                                                                    phi[i]* b * dphidy[j]);                                                                                           
                            A[mapY[i]][mapX[j]] += jac * weight * ( dphidy[i] * b * x[i] * dphidx[j] + 
                                                                    dphidx[i] * c * x[i] * dphidy[j] +
                                                                    dphidy[i] * b * phi[j]);                                                                                            
                            A[mapY[i]][mapY[j]] += jac * weight * ( dphidy[i] * a * x[i] * dphidy[j] + 
                                                                    dphidx[i] * c * x[i] * dphidx[j]);
                    }
                }
            }
            for(int i = 0; i < theSpace->n; i++) {
                switch(theProblem->planarStrainStress<2){
                    case TRUE:
                        B[mapY[i]] -= rho * g * phi[i] * jac * weight;
                    case FALSE:
                        B[mapX[i]] -= rho * g * sin(y[i]) * phi[i] * jac * weight * x[i];
                        B[mapY[i]] -= rho * g * cos(y[i]) * phi[i] * jac * weight * x[i];
                }
            }
        }

    }
}

femSystem* systemConstrain(femSystem *theSystem, femProblem *theProblem)
{
    int* theConstrainedNodes = theProblem->constrainedNodes;
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            femBoundaryCondition* condition = theProblem->conditions[theConstrainedNodes[i]];
            double value = condition->value;
            switch (condition->type<=3)
            {
            case TRUE:
                switch(theSystem->solverType){
                    case 0: dirichlet(theSystem->fullSystem->A, theSystem->fullSystem->B, theSystem->size, i, value);
                    case 1: dirichlet(theSystem->bandSystem->A, theSystem->bandSystem->B, theSystem->size, i, value);
                    case 2: dirichlet(theSystem->fullSystem->A, theSystem->fullSystem->B, theSystem->size, i, value);
                }
                break;
            case FALSE:
                switch(theSystem->solverType){
                    case 0: neumann(theSystem->fullSystem->B, i, ((3-sqrt(3))/3)*value*domainLength(condition->domain)/2);
                    case 1: neumann(theSystem->bandSystem->B, i, ((3-sqrt(3))/3)*value*domainLength(condition->domain)/2);
                    case 2: neumann(theSystem->fullSystem->B, i, ((3-sqrt(3))/3)*value*domainLength(condition->domain)/2);
                }
            } 
        }
    }
    return theSystem;
}

double *systemSolve(femProblem *theProblem)
{
    systemAssemble(theProblem);
    theProblem->system = systemConstrain(theProblem->system, theProblem);
    return femBandSystemEliminate(theProblem->system->bandSystem);
}
