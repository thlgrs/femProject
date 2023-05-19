#include "fem.h"
#include "myfem.h"

double *GLOBAL_COORD;

/*PRIMARY FUNCTIONS*/
femProblem* femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType iSolver, femRenumType iRenum)
{
    FILE* file = fopen(filename,"r");
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    theProblem->constrainedNodes = malloc(size*sizeof(int));
    for (int i=0; i < size; i++) 
        theProblem->constrainedNodes[i] = -1;
    
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }


    char theLine[MAXNAME];
    char theDomain[MAXNAME];
    char theArgument[MAXNAME];
    double value;
    double typeCondition;
    
    while (!feof(file)) {
        ErrorScan(fscanf(file,"%19[^\n]s \n",(char *)&theLine));
        if (strncasecmp(theLine,"Type of problem     ",19) == 0) {
            ErrorScan(fscanf(file,":  %[^\n]s \n",(char *)&theArgument));
            if (strncasecmp(theArgument,"Planar stresses",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRESS; 
            if (strncasecmp(theArgument,"Planar strains",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRAIN; 
            if (strncasecmp(theArgument,"Axi-symetric problem",13) == 0)
               theProblem->planarStrainStress = AXISYM; }
        if (strncasecmp(theLine,"Young modulus       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->E)); }
        if (strncasecmp(theLine,"Poisson ratio       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->nu)); }
        if (strncasecmp(theLine,"Mass density        ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->rho)); }
        if (strncasecmp(theLine,"Gravity             ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->g)); }
        if (strncasecmp(theLine,"Boundary condition  ",19) == 0) {
            ErrorScan(fscanf(file,":  %19s = %le : %[^\n]s\n",(char *)&theArgument,&value,(char *)&theDomain));
            if (strncasecmp(theArgument,"Dirichlet-X",19) == 0)
                typeCondition = DIRICHLET_X;
            if (strncasecmp(theArgument,"Dirichlet-Y",19) == 0)
                typeCondition = DIRICHLET_Y;                
            if (strncasecmp(theArgument,"Neumann-X",19) == 0)
                typeCondition = NEUMANN_X;
            if (strncasecmp(theArgument,"Neumann-Y",19) == 0)
                typeCondition = NEUMANN_Y;                
            femElasticityAddBoundaryCondition(theProblem,theDomain,typeCondition,value); }
        ErrorScan(fscanf(file,"\n")); 
    }
 
    int iCase = theProblem->planarStrainStress;
    double E = theProblem->E;
    double nu = theProblem->nu;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }

    femSystemCreate(size, iSolver, iRenum, theProblem);

    fclose(file);
    return theProblem;
}

void femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femProblem *theProblem){
    femGeo *theGeometry = theProblem->geometry;
    femSystem* system = malloc(sizeof(femSystem));
    system->solverType = iSolver;
    int nLocal = theGeometry->theElements->nLocalNode;
    system->local = femFullSystemCreate(nLocal);
    
    switch(iSolver) {
        case FEM_FULL :
            system->fullSystem = femFullSystemCreate(size); break;
        case FEM_BPOST :
            renumberNodes(theProblem, iRenum);
            system->bandSystem = femPostBandSystemCreate(size); break;
        case FEM_BAND : 
            renumberNodes(theProblem, iRenum);
            int band = computeBand(theGeometry->theElements);
            system->bandSystem = femBandSystemCreate(size,band); break;
        case FEM_FRONT :
            system->fullSystem = femFullSystemCreate(size);
            renumberElem(theGeometry, iRenum);
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
    theProblem->system = system;
}

void femElasticitySet(femProblem *theProblem){
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femMesh        *theElements = theGeometry->theElements;
    femNodes       *theNodes = theElements->nodes;
    femBoundaryCondition *condition;
    femBoundaryType      type;

    double x[4],y[4],dx,dy,phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4],coeff[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4],constrained[4],globMap[4];

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
            constrained[j]  = theProblem->constrainedNodes[map[j]]; 
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
                      
            switch(theProblem->planarStrainStress){
                case PLANAR_STRESS:
                    localPlannar(theSpace, A, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight); break;
                case PLANAR_STRAIN:
                    localPlannar(theSpace, A, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight); break;
                case AXISYM:
                    localAxisym(theSpace, A, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight); break;
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
    for (i = 0; i < theSpace->n; i++){ 
        if (constrained[i] == 1){
            condition = theProblem->conditions[map[i]]; 
            switch (condition->type<2){
                case TRUE:
                    constrain(type,theProblem->system->local,i,condition->value,1,1); break;
                case FALSE:
                    if(i==0) {dx = x[i]-x[theSpace->n-1]; dy = y[i]-y[theSpace->n-1];}
                    else {dx = x[i]-x[i-1]; dy = y[i]-y[i-1];}
                    constrain(type,theProblem->system->local,i,condition->value, dx, dy); break;
            }                       
        }
    }
};

double *femElasticitySolve(femProblem *theProblem)
{   

    switch(theProblem->system->solverType) {
        case FEM_FULL : 
            femElasticitySet(theProblem);
            return femFullSystemEliminate(theProblem->system->fullSystem);
            break;
        case FEM_BAND :
            femElasticitySet(theProblem);
            return femBandSystemEliminate(theProblem->system->bandSystem);
            break;
        case FEM_BPOST:
            femPostBandSet(theProblem);
            setBand(theProblem->system->bandSystem);
            return femBandSystemEliminate(theProblem->system->bandSystem);
            break;
        case FEM_FRONT :
            femElasticitySet(theProblem);  
            femFrontalSolve(theProblem);
            double *soluce;
            return soluce;
            break;
        default :
            Error("Unexpected solver type");
    };
}

void femElasticityFree(femProblem *theProblem)
{
    femSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    free(theProblem->conditions);
    free(theProblem->constrainedNodes);
    free(theProblem);
}

void femSystemFree(femSystem* mySystem){
    switch(mySystem->solverType) {
        case FEM_FULL : femFullSystemFree(mySystem->fullSystem); break;
        case FEM_BAND :
            femFullSystemFree(mySystem->local); 
            femBandSystemFree(mySystem->bandSystem); break;
        case FEM_FRONT : 
            femFullSystemFree(mySystem->fullSystem);
            femFrontalSolverFree(mySystem->frontSolver); break;
    }
    free(mySystem);
}



/*SECONDARY FUNCTIONS*/

int computeBand(femMesh *theElements)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theElements->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theElements->nElem; iElem++) {
        for (j=0; j < nLocal; ++j){
            int elemNumber = theElements->elem[iElem*nLocal+j];
            map[j] = 2*elemNumber;
        } 
            map[j] = theElements->elem[iElem*nLocal+j];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}

int compare(const void *i1, const void *i2){
    int *iOne = (int *)i1;
    int *iTwo = (int *)i2;
    double diff = GLOBAL_COORD[*iOne] - GLOBAL_COORD[*iTwo];

    return ( diff < 0) - ( diff > 0) ;
};

void renumberNodes(femProblem *theProblem, femRenumType renumType){
    int i;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes = theGeometry->theNodes;
    double *newX = malloc(sizeof(double)*theNodes->nNodes);
    double *newY = malloc(sizeof(double)*theNodes->nNodes);
    femMesh *theElements = theGeometry->theElements;
    int *newElem = malloc(sizeof(int)*theElements->nElem);
    femMesh *theEdges = theGeometry->theEdges;
    int *newEdges = malloc(sizeof(int)*theEdges->nElem);
    int *constrained = theProblem->constrainedNodes;
    int* newConstrained = malloc(sizeof(int)*theNodes->nNodes);
    
    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM :
            GLOBAL_COORD = theNodes->X; 
            break;
        case FEM_YNUM : 
            GLOBAL_COORD = theNodes->Y; 
            break;
        default : Error("Unexpected renumbering option"); 
    }
    
    int *map = malloc(sizeof(int)*theNodes->nNodes);
    int *inverse = malloc(sizeof(int)*theNodes->nNodes);
    for (i = 0; i < theNodes->nNodes; i++){map[i] = i;}

    qsort(map, theNodes->nNodes, sizeof(int), compare);
    
    for (i=0; i<theNodes->nNodes; i++){
        inverse[map[i]] = i;
        newX[i] = theNodes->X[map[i]];
        newY[i] = theNodes->Y[map[i]];
        newConstrained[i] = constrained[map[i]];
    }
    for(i=0; i<theElements->nElem; i++){
        newElem[i] = inverse[theElements->elem[i]];
    }
    for(i=0; i<theEdges->nElem; i++){
        newEdges[i] = inverse[theEdges->elem[i]];
    }

    for(i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] = newX[i];
        theNodes->Y[i] = newY[i];
        constrained[i] = newConstrained[i];
    }
    for(i=0; i<theElements->nElem; i++){
        theElements->elem[i] = newElem[i];
    }
    for(i=0; i<theEdges->nElem; i++){
        theEdges->elem[i] = newEdges[i];
    }
    free(map);
    free(newX);
    free(newY);
    free(newElem);
    free(newEdges);
    free(inverse);
};

void renumberElem(femGeo *theGeometry, femRenumType renumType){

    int i, j, nLocalNode;
    femMesh *theElements = theGeometry->theElements;
    nLocalNode = theElements->nLocalNode;
    double* slice = malloc(sizeof(double)*nLocalNode);
    int *oldElements = theGeometry->theElements->elem;
    int *newElements = malloc(sizeof(int)*theElements->nElem);
    for (i = 0; i < theElements->nElem; i++){
        newElements[theElements->nElem-i-1] = i;
    }
    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM :
            GLOBAL_COORD = realloc(GLOBAL_COORD, sizeof(double)*theElements->nElem);
            for(i = 0; i < theElements->nElem; i++){
                for(j = 0; j < nLocalNode; j++)
                    slice[j] = theElements->nodes->X[i*nLocalNode+j];
                GLOBAL_COORD[i] = femMin(slice, nLocalNode);
            }             
            break;
        case FEM_YNUM :
            GLOBAL_COORD = realloc(GLOBAL_COORD, sizeof(double)*theElements->nElem);
            for(i = 0; i < theElements->nElem; i++){
                for(j = 0; j < nLocalNode; j++)
                    slice[j] = theElements->nodes->Y[i*nLocalNode+j];
                GLOBAL_COORD[i] = femMin(slice, nLocalNode);
            }
            break;
        default : Error("Unexpected renumbering option"); 
    }
    qsort(newElements, theElements->nElem, sizeof(int), compare);
    for(i = 0; i < theElements->nElem; i++){
        oldElements[i] = newElements[i];
    }
    free(newElements);
    free(slice);
    free(GLOBAL_COORD);

};

int nodeIndex(double x, double y, femNodes* theNodes){
    for(int node=0; node < theNodes->nNodes; node++){
        if (x == theNodes->X[node] && y == theNodes->Y[node]) return node;
    }

}

void localPlannar(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight){
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
    for(int i = 0; i < n; i++){
        B[i+1] -= x[i+1] * coeff[3] * jac * weight;
    }

    
}

void localAxisym(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight){
    int n = space->n;
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            A[i][j] += jac * weight * ( dx[i] * coeff[0] * x[i] * dx[j] + 
                                        dy[i] * coeff[2] * x[i] * dy[j] +
                                        phi[i]* (coeff[1] *dx[j] + coeff[0] * phi[j]/x[i]) +
                                        dx[i] *coeff[1] * phi[j]);
        }
        for(int j=(n-2); j<n; j++){
            A[i][j] += jac * weight * ( dx[i] * coeff[1] * x[i] * dy[j] + 
                                        dy[i] * coeff[2] * x[i] * dx[j] +
                                        phi[i]* coeff[1] * dy[j]);
        }
        for(int j=0; j<2; j++){
            A[i+(n-2)][j] += jac * weight * (dy[i+(n-2)] * coeff[1] * x[i+(n-2)] * dx[j] + 
                                             dx[i+(n-2)] * coeff[2] * x[i+(n-2)] * dy[j] +
                                             dy[i+(n-2)] * coeff[1] * phi[j]);
        }
        for(int j=(n-2); j<n; j++){
            A[i+(n-2)][j] += jac * weight * (dy[i+(n-2)] * coeff[0] * x[i+(n-2)] * dy[j] + 
                                             dx[i+(n-2)] * coeff[2] * x[i+(n-2)] * dx[j]);
        }
    }
    for(int i = 0; i < n; i++){
        B[i+1] -= x[i+1] * phi[i+1] * coeff[3] * jac * weight;
    }
}

void constrain(femBoundaryType type, femFullSystem* system, int myNode, double value, double dx, double dy){
    double** A = system->A;
    double* B = system->B;
    int size = system->size;

    double length = sqrt(dx*dx+dy*dy);
    double cos = dx/length;
    double sin = dy/length;
    switch (type){
        case DIRICHLET_X:
            dirichlet(A,B,size,myNode,value); break;
        case DIRICHLET_Y:  
            dirichlet(A,B,size,myNode,value); break;
        case NEUMANN_X: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value);   //line integral approximation 2point gauss
            break;                        // value = derivative of the flux on boundary
        case NEUMANN_Y: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value); 
            break;
        /*A IMPOSER AVEC TABLEAU METHODE DES DEPLACEMENTS*/
        case NEUMANN_N:
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value*sin);   //en X
            neumann(B,myNode,value*cos); //en Y
            break;
        case NEUMANN_T:
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value*cos);   //en X
            neumann(B,myNode,value*sin); //en Y
            break;
        default: break;
    }
}

void constrainb(femBoundaryType type, femBandSystem* system, int myNode, double value, double dx, double dy){
    double** A = system->A;
    double* B = system->B;
    int size = system->size;

    double length = sqrt(dx*dx+dy*dy);
    double cos = dx/length;
    double sin = dy/length;
    switch (type){
        case DIRICHLET_X:
            dirichlet(A,B,size,myNode,value); break;
        case DIRICHLET_Y:  
            dirichlet(A,B,size,myNode,value); break;
        case NEUMANN_X: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value);   //line integral approximation 2point gauss
            break;                        // value = derivative of the flux on boundary
        case NEUMANN_Y: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value); 
            break;
        /*A IMPOSER AVEC TABLEAU METHODE DES DEPLACEMENTS*/
        case NEUMANN_N:
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value*sin);   //en X
            neumann(B,myNode,value*cos); //en Y
            break;
        case NEUMANN_T:
            value = ((3-sqrt(3))/3)*value*length/2;
            neumann(B,myNode,value*cos);   //en X
            neumann(B,myNode,value*sin); //en Y
            break;
        default: break;
    }
}

void dirichlet(double **A, double *B, int size, int myNode, double myValue){
    int i;
    for(i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0.0; }
    for(i=0; i < size; i++){
        A[myNode][i] = 0.0;
    }
    A[myNode][myNode] = 1.0;
    
    B[myNode] = myValue;
    

}

void neumann(double *B, int myNode, double myValue){
    B[myNode] += myValue;
}

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


/*BAND FUNCTIONS*/

femBandSystem* femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;    
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);

}

void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = 2*map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = 2*map[j];
            if (myCol >= myRow ){ 
                myBandSystem->A[2*map[i]][2*map[j]]     += Aloc[2*i*nLoc+2*j];
                myBandSystem->A[2*map[i]][2*map[j]+1]   += Aloc[2*nLoc*i+2*j+1];
                myBandSystem->A[2*map[i]+1][2*map[j]]   += Aloc[(2*i+1)*nLoc+2*j];
                myBandSystem->A[2*map[i]+1][2*map[j]+1] += Aloc[(2*i+1)*nLoc+2*j+1];
            }
        }
        myBandSystem->B[myRow] += Bloc[2*i];
        myBandSystem->B[myRow+1] += Bloc[2*i+1]; }
    
}

double* femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Check for isolated node */

    for (k = 0; k < size/2; k++) {
        if ((A[2*k][2*k] == 0) && (A[2*k+1][2*k+1] == 0)) {
            printf("Warning : disconnected node %d\n", k);
            A[2*k][2*k] += 1;
            A[2*k+1][2*k+1] += 1; }}
    
    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}

/*POST BAND*/
femBandSystem* femPostBandSystemCreate(int size){
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    int i;  
    myBandSystem->A = malloc(sizeof(double*) * size); 
    myBandSystem->B = calloc(size,sizeof(double)); 
    myBandSystem->size = size;
    myBandSystem->band = 0;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = calloc(size,sizeof(double));
    return(myBandSystem);
};

void femPostBandSet(femProblem* theProblem){
    femSystem      *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->bandSystem->A;
    double *B  = theSystem->bandSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            switch (theProblem->planarStrainStress<2){
                case TRUE:
                    for (i = 0; i < theSpace->n; i++) {
                        B[mapY[i]] -= phi[i] * g * rho * jac * weight;
                        for(j = 0; j < theSpace->n; j++) {
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                                    dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                                    dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                                    dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                                    dphidx[i] * c * dphidx[j]) * jac * weight;
                        }
                    }
                case FALSE:
                    for (i = 0; i < theSpace->n; i++) {
                        B[mapY[i]] -= x[i] * phi[i] * rho * g * jac * weight;
                        for(j = 0; j < theSpace->n; j++) {
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * x[i] * dphidx[j] + 
                                                    dphidy[i] * c * x[i] * dphidy[j] +
                                                    phi[i] * (b *dphidx[j] + a * phi[j]/x[i]) +
                                                    dphidx[i] *b * phi[j]) * jac * weight;                                                                                          
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * x[i] * dphidy[j] + 
                                                    dphidy[i] * c * x[i] * dphidx[j] +
                                                    phi[i] * b * dphidy[j]) * jac * weight;                                                                                           
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * x[i] * dphidx[j] + 
                                                    dphidx[i] * c * x[i] * dphidy[j] +
                                                    dphidy[i] * b * phi[j]) * jac * weight;                                                                                            
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * x[i] * dphidy[j] + 
                                                    dphidx[i] * c * x[i] * dphidx[j]) * jac * weight;
                        }
                    }
            }
        }
    }    
    femBoundaryCondition* condition;
    femBoundaryType type;
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (i = 0; i < theNodes->nNodes; i++){ 
        if (theConstrainedNodes[i] == 1){
            condition = theProblem->conditions[i]; 
            type = condition->type;
            if(type == DIRICHLET_X || type == DIRICHLET_Y){
                constrainb(type,theSystem->bandSystem,i,condition->value,1,1);
            }
            else if(type == NEUMANN_X || type == NEUMANN_Y){
                double dx, dy;
                if(i==0) {
                    dx = x[i]-x[theSpace->n-1];
                    dy = y[i]-y[theSpace->n-1];}
                else {
                    dx = x[i]-x[i-1];
                    dy = y[i]-y[i-1];}
                constrainb(type,theSystem->bandSystem,i,condition->value, dx, dy);
            }    
        }
    }                   
}

void setBand(femBandSystem* myBandSystem){
    int band = myBandSystem->band;
    printf("band = %d\n", myBandSystem->band);
    for(int i = 0; i < myBandSystem->size; i++){
        for(int j = i; j < myBandSystem->size; j++){
            if(myBandSystem->A[i][j] != 0){
                band = ( band > abs(i-j)) ? band : abs(i-j);
            }
        }
    }
    myBandSystem->band = band+1;
    printf("band = %d\n", myBandSystem->band);
};


/*FRONTAL SOLVER not-implemented*/

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

void femFrontalSolverInit(femFrontalSolver *mySolver){};

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


/*DUMP*/

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

void femBoundaryConstrain(femProblem *theProblem, double **A, double *B){
    femBoundaryCondition **conditions = theProblem->conditions;
    femBoundaryCondition *theBoundary;
    femBoundaryType type;
    femDomain *theDomain;
    femMesh *theEdge;
    double h = theProblem->geometry->h;

    int i,node,iElem, nElem;
    double value,length, dx, dy;
    int* map  = calloc(1, 4);
    int* mapX = calloc(1, 4);
    int* mapY = calloc(1, 4);
    double* x = calloc(1, 4);
    double* y = calloc(1, 4);

    
    for(i=0; i<theProblem->nBoundaryConditions; i++){
        theBoundary = conditions[i];
        type = theBoundary->type;
        value = theBoundary->value;
        theDomain = theBoundary->domain;
        theEdge = theDomain->mesh;
        nElem = theDomain->nElem;
        x    = realloc(x, nElem*sizeof(double));
        y    = realloc(y, nElem*sizeof(double));
        map  = realloc(map, nElem*sizeof(int));
        mapX = realloc(mapX, nElem*sizeof(int));
        mapY = realloc(mapY, nElem*sizeof(int));

        map[0] = theEdge->elem[0];
        mapX[0] = 2*map[0];
        mapY[0] = 2*map[0]+1;
        x[0] = theEdge->nodes->X[map[0]];
        y[0] = theEdge->nodes->Y[map[0]];

        
        for(node=0; node<nElem; ++node){
            map[node] = theEdge->elem[node];
            mapX[node] = 2*map[node];
            mapY[node] = 2*map[node]+1;
            x[node] = theEdge->nodes->X[map[node]];
            y[node] = theEdge->nodes->Y[map[node]];
            dx = x[node]-x[node-1]; dy = y[node]-y[node-1];
            length += sqrt(dx*dx+dy*dy);

            switch (type){
                case DIRICHLET_X:
                    dirichlet(A,B,sizeof(A[0])/sizeof(double),mapX[node],value); break;
                case DIRICHLET_Y:  
                    dirichlet(A,B,sizeof(A[0])/sizeof(double),mapY[node],value); break;
                case NEUMANN_X:
                    neumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); //line integral approximation 2point gauss
                    break;                                                   // value = derivative of the flux on boundary
                case NEUMANN_Y:
                    neumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                /*A IMPOSER AVEC TABLEAU METHODE DES DEPLACEMENTS*/
                case NEUMANN_N:
                    neumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); 
                    neumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                case NEUMANN_T:
                    neumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); 
                    neumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                default: break;
            }       
        }

    }
    free(map);
    free(mapX);
    free(mapY);
    free(x);
    free(y);
}
