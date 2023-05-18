#include "fem.h"
#include "myfem.h"



int nActive(femMesh *theElements, int **disparus, int **nouveaux){
    femNodes *theNodes = theElements->nodes;
    int nElem = theElements->nElem;
    int *elem = theElements->elem;
    int nLocal = theElements->nLocalNode;
    int idElem1, idElem2, keep, new, nActive = nLocal;
    for(int k = 0; k<nElem; ++k){
        idElem1 = elem[k-1]; idElem2 = elem[k];
        keep = 0; new = 0;
        for(int n=0; n<nLocal; n++){
            double x1 = theNodes->X[idElem1*nLocal+n];
            double y1 = theNodes->Y[idElem1*nLocal+n];
            for(int m=0; m<nLocal; m++){
                double x2 = theNodes->X[idElem2*nLocal+m];
                double y2 = theNodes->Y[idElem2*nLocal+m];

                if(x1 != x2 && y1 != y2){
                    disparus[idElem1][n] = idElem1*nLocal+n;}
                else{
                    disparus[idElem1][n] = -1;
                    keep++;
                }

                for(int i=0; i<nLocal; i++){
                    double x3 = theNodes->X[idElem1*nLocal+i];
                    double y3 = theNodes->Y[idElem1*nLocal+i];
                    if(x2 != x3 || y2 != y3){
                        nouveaux[idElem1][m] = idElem2*nLocal+m;
                        new++;
                    }
                    else{
                        nouveaux[idElem1][m] = -1;
                    }
                }
            }       
        }
        nActive = nActive > keep+new ? nActive : keep+new;
    }
    return(2*nActive);
};

void femFrontalSolverInit(femFrontalSolver *mySolver){

};

void femAssembleLocal(double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal){
    for(int i=0; i<2*nLocal; i++){
        for(int j=0; j<2*nLocal; j++){
            ALoc[i][j] = AGlob[2*map[i/2]][2*map[j/2]];
        }
        BLoc[i] = BGlob[2*map[i/2]];
    }
};

void femInjectGlobal(double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal){
    for(int i=0; i<2*nLocal; i++){
        for(int j=0; j<2*nLocal; j++)
            AGlob[2*map[i/2]][2*map[j/2]+1] = ALoc[i][j];
        BGlob[2*map[i/2]] = BLoc[i];
    }
};

double* femFrontalSolve(femProblem *theProblem){
    femFullSystem *theSystem = theProblem->system->fullSystem;
    femFrontalSolver *theSolver = theProblem->system->frontSolver;
    femMesh *theElements = theProblem->geometry->theElements;
    double** Aglob = theSystem->A;
    double* Bglob = theSystem->B;
    int nElem = theElements->nElem;
    int nLocal = theElements->nLocalNode;
    double** active = theSolver->active;
    double*  b = theSolver->b;

};


femSystem* femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femGeo *theGeometry){
    femSystem* system = malloc(sizeof(femSystem));
    system->solverType = iSolver;
    int nLocal = theGeometry->theElements->nLocalNode;
    system->local = femFullSystemCreate(nLocal);
    
    switch(iSolver) {
        case FEM_FULL :
            system->fullSystem = femFullSystemCreate(size); break;
        case FEM_BAND : 
            system->fullSystem = femFullSystemCreate(size);
            femRenumberNodes(theGeometry, iRenum);
            int band = femMeshComputeBand(theGeometry->theElements);
            system->bandSystem = femBandSystemCreate(size,band); break;
        case FEM_FRONT :
            system->fullSystem = femFullSystemCreate(size);
            //femRenumberElem(theGeometry, iRenum);
            int nLoc = theGeometry->theElements->nLocalNode;
            int** disparus = malloc(theGeometry->theElements->nElem*sizeof(int*));
            int** nouveaux = malloc(theGeometry->theElements->nElem*sizeof(int*));
            for (int i=0; i < theGeometry->theElements->nElem; i++) {
                disparus[i] = calloc(nLoc,sizeof(int));
                nouveaux[i] = calloc(nLoc,sizeof(int));
            }
            int nAct = nActive(theGeometry->theElements,disparus,nouveaux);
            system->frontSolver = femFrontalSolverCreate(size,nAct); 
            system->frontSolver->disparus = disparus;
            system->frontSolver->nouveaux = nouveaux;
            break;
    }
    return system;
}

double *femElasticitySolve(femProblem *theProblem)
{   
    
    
    switch(theProblem->system->solverType) {
        case FEM_FULL : 
            femGlobalSystemAssemble(theProblem);
            return femFullSystemEliminate(theProblem->system->fullSystem);
            break;
        case FEM_BAND :
            femGlobalSystemAssemble(theProblem);
            return femBandSystemEliminate(theProblem->system->bandSystem);
            break;
        case FEM_FRONT :
            femGlobalSystemAssemble(theProblem);  
            femFrontalSolve(theProblem);
            double *soluce;
            return soluce;
            break;
        default :
            Error("Unexpected solver type");
    };
}


void femGlobalSystemAssemble(femProblem *theProblem){
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femMesh        *theElements = theGeometry->theElements;
    femNodes       *theNodes = theElements->nodes;
    femBoundaryCondition *condition;
    femBoundaryType type;

    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4],coeff[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4],ctr[4],globMap[4];

    double **A   = theProblem->system->local->A;
    double *Aloc = theProblem->system->local->A[0];
    double *Bloc = theProblem->system->local->B;
    
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    coeff[0] = a; coeff[1] = b; coeff[2] = c; coeff[3] = rho*g;

    int nLocal = theElements->nLocalNode;
    for (iElem = 0; iElem < theElements->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        for (j=0; j < nLocal; j++) {
            map[j]  = theElements->elem[iElem*nLocal+j]; 
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
            ctr[j]  = theProblem->constrainedNodes[map[j]]; 
            map[j] = nodeIndex(x[j],y[j],theGeometry->theNodes); //num noeuds de l'élément iElem
        }
        for (iInteg = 0; iInteg < theRule->n; iInteg++){
            double weight = theRule->weight[iInteg];     
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);                  //initialise les fonctions de formes dans phi
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);   //initialise les dérivées des fonctions de formes
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < nLocal; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; } 
                      
            if(theProblem->planarStrainStress == AXISYM)
            {
                femLocalAxsym(theSpace, A, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight);
            }
            else if(theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
            {
                femLocalPlan(theSpace, A, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight);
            }
            for (i = 0; i < theSpace->n; i++)
            { 
                if (ctr[i] == 1){
                    condition = theProblem->conditions[map[i]]; 
                    type = condition->type;
                    if(type == DIRICHLET_X || type == DIRICHLET_Y){
                        constrain(type,theProblem->system->local,i,condition->value,1,1);
                    }
                    else if(type == NEUMANN_X || type == NEUMANN_Y){
                        
                        double dx, dy;
                        if(i==0) {
                            dx = x[i]-x[theSpace->n-1];
                            dy = y[i]-y[theSpace->n-1];}
                        else {
                            dx = x[i]-x[i-1];
                            dy = y[i]-y[i-1];}
                        constrain(type,theProblem->system->local,i,condition->value, dx, dy);
                    }
                       
                }
            }
            
        }
        switch(theProblem->system->solverType){
            case FEM_FULL:
                femFullSystemAssemble(theProblem->system->fullSystem,Aloc,Bloc,map,nLocal); break;
            case FEM_BAND:
                femBandSystemAssemble(theProblem->system->bandSystem,Aloc,Bloc,map,nLocal); break;
            case FEM_FRONT:
                femFullSystemAssemble(theProblem->system->fullSystem,Aloc,Bloc,map,nLocal); break;
        }
    }
};


femFrontalSolver* femFrontalSolverCreate(int size, int nActive){
    femFrontalSolver *mySolver = malloc(sizeof(femFrontalSolver));
    mySolver->nActive = nActive;
    mySolver->b = calloc(nActive, sizeof(double));
    mySolver->active = malloc(sizeof(double*)*nActive);
    for(int i=0; i<nActive; i++){
        mySolver->active[i] = calloc(nActive, sizeof(double));
    }
    mySolver->activeMap = calloc(nActive, sizeof(int));
    mySolver->activeState = calloc(nActive, sizeof(int));
    mySolver->pivots = calloc(size, sizeof(int));
    mySolver->stock = malloc(size*sizeof(double*));
    for(int i=0; i<size; i++){
        mySolver->stock[i] = calloc(nActive, sizeof(double));
    }
    return(mySolver);
};

void femFrontalSolverFree(femFrontalSolver *mySolver){
    free(mySolver->activeMap);
    free(mySolver->activeState);
    free(mySolver->pivots);
    free(mySolver->b);
    for(int i=0; i<mySolver->nActive; i++){
        free(mySolver->active[i]);
    }
    free(mySolver->active);
    for(int i=0; i<sizeof(mySolver->stock)/sizeof(double*); i++){
        free(mySolver->stock[i]);
        free(mySolver->disparus[i]);
        free(mySolver->nouveaux[i]);
    }
    free(mySolver->stock);
    free(mySolver->disparus);
    free(mySolver->nouveaux);
    free(mySolver);
};

int nodeIndex(double x, double y, femNodes* theNodes){
    for(int node=0; node < theNodes->nNodes; node++){
        if (x == theNodes->X[node] && y == theNodes->Y[node]) return node;
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

double integrate(double *x, double *y, double *dphidx, double*dphidy, femIntegration *rule, femDiscrete *space){
    double* phi = calloc(space->n,sizeof(double));
    double* dphidxsi = calloc(space->n,sizeof(double));
    double* dphideta = calloc(space->n,sizeof(double));
    double dxdxsi = 0.0;
    double dxdeta = 0.0;
    double dydxsi = 0.0; 
    double dydeta = 0.0;
    for(int i=0; i<rule->n; i++){
        femDiscretePhi2(space,rule->xsi[i],rule->eta[i],phi);
        femDiscreteDphi2(space,rule->xsi[i],rule->eta[i],dphidxsi,dphideta);
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


