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



/*--------------------------------------------OK----------------------------------------------------------*/
double *theGlobalCoord;

int femCompare(const void *i1, const void *i2){
    int *iOne = (int *)i1;
    int *iTwo = (int *)i2;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];

    return ( diff < 0) - ( diff > 0) ;
};

void femRenumberNodes(femGeo *theGeometry, femRenumType renumType){
    int i;
    femNodes *oldNodes = theGeometry->theNodes;
    double *newX = malloc(sizeof(double)*oldNodes->nNodes);
    double *newY = malloc(sizeof(double)*oldNodes->nNodes);
    
    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM :
            theGlobalCoord = oldNodes->X; 
            break;
        case FEM_YNUM : 
            theGlobalCoord = oldNodes->Y; 
            break;
        default : Error("Unexpected renumbering option"); 
    }
    
    int *inverse = malloc(sizeof(int)*oldNodes->nNodes);
    for (i = 0; i < oldNodes->nNodes; i++){
        inverse[oldNodes->nNodes-i-1] = i;
    }
    
    qsort(inverse, oldNodes->nNodes, sizeof(int), femCompare);
    
    for (i=0; i<oldNodes->nNodes; i++){
        newX[i] = oldNodes->X[inverse[i]];
        newY[i] = oldNodes->Y[inverse[i]];
    }
    for(i=0; i<oldNodes->nNodes; i++){
        oldNodes->X[i] = newX[i];
        oldNodes->Y[i] = newY[i];
    }

    free(inverse);
    free(newX);
    free(newY);
};

void femRenumberElem(femGeo *theGeometry, femRenumType renumType){

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
            theGlobalCoord = realloc(theGlobalCoord, sizeof(double)*theElements->nElem);
            for(i = 0; i < theElements->nElem; i++){
                for(j = 0; j < nLocalNode; j++)
                    slice[j] = theElements->nodes->X[i*nLocalNode+j];
                theGlobalCoord[i] = femMin(slice, nLocalNode);
            }             
            break;
        case FEM_YNUM :
            theGlobalCoord = realloc(theGlobalCoord, sizeof(double)*theElements->nElem);
            for(i = 0; i < theElements->nElem; i++){
                for(j = 0; j < nLocalNode; j++)
                    slice[j] = theElements->nodes->Y[i*nLocalNode+j];
                theGlobalCoord[i] = femMin(slice, nLocalNode);
            }
            break;
        default : Error("Unexpected renumbering option"); 
    }
    qsort(newElements, theElements->nElem, sizeof(int), femCompare);
    for(i = 0; i < theElements->nElem; i++){
        oldElements[i] = newElements[i];
    }
    free(newElements);
    free(slice);
    free(theGlobalCoord);

};

int femMeshComputeBand(femMesh *theMesh)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->elem[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}

femBandSystem*  femBandSystemCreate(int size, int band)
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
            A[2*k][2*k] = 1;
            A[2*k+1][2*k+1] = 1; }}

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

void femBandSystemAssemble(femBandSystem* myBandSystem, double **Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow) {
                myBandSystem->A[2*myRow][2*myCol] += Aloc[2*i][2*j];
                myBandSystem->A[2*myRow][2*myCol+1] += Aloc[2*i][2*j+1];
                if(2*myCol >= 2*myRow+1) 
                    myBandSystem->A[2*myRow+1][2*myCol] += Aloc[2*i+1][2*j];
                myBandSystem->A[2*myRow+1][2*myCol+1] += Aloc[2*i+1][2*j+1];
            }
        }
        myBandSystem->B[myRow] += Bloc[i]; }
}


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

    theProblem->system = femSystemCreate(size, iSolver, iRenum, theGeometry);

    fclose(file);
    return theProblem;
}

femSystem* femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femGeo *theGeometry){
    femSystem* system = malloc(sizeof(femSystem));
    system->solverType = iSolver;

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
            femRenumberElem(theGeometry, iRenum);
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
    femGlobalSystemAssemble(theProblem, theProblem->system->fullSystem->A, theProblem->system->fullSystem->B);
    switch(theProblem->system->solverType) {
        case FEM_FULL : 
            return femFullSystemEliminate(theProblem->system->fullSystem);
            break;
        case FEM_BAND :
            return femBandSystemEliminate(theProblem->system->bandSystem);
            break;
        case FEM_FRONT :
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
            femFullSystemFree(mySystem->fullSystem); 
            femBandSystemFree(mySystem->bandSystem); break;
        case FEM_FRONT : 
            femFullSystemFree(mySystem->fullSystem);
            femFrontalSolverFree(mySystem->frontSolver); break;
    }
    free(mySystem);
}

void femGlobalSystemAssemble(femProblem *theProblem, double **A, double *B){
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femMesh        *theElements = theGeometry->theElements;
    femNodes       *theNodes = theElements->nodes;

    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4],coeff[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4],ctr[4];

    double **Aloc = theProblem->system->Aloc;
    double *Bloc = theProblem->system->Bloc;
    
    
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    coeff[0] = a; coeff[1] = b; coeff[2] = c; coeff[3] = rho*g;

    int nLocal = theElements->nLocalNode;
    for (iElem = 0; iElem < theElements->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      
            Bloc[i] = 0;
            for (j = 0; j < theSpace->n; j++) 
                Aloc[i] = 0;
        for (j=0; j < nLocal; j++) {
            map[j]  = theElements->elem[iElem*nLocal+j]; //num noeuds de l'élément iElem
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
            ctr[j]  = theProblem->constrainedNodes[map[j]];
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
                femLocalAxsym(theSpace, Aloc, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight);
            }
            else if(theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
            {
                femLocalPlan(theSpace, Aloc, Bloc, x, phi, dphidx, dphidy, coeff, jac, weight);
            }
            for (i = 0; i < theSpace->n; i++)
            {
                if (ctr[i] == 1){
                    femBoundaryCondition *condition = theProblem->conditions[map[i]];
                    femBoundaryType type = condition->type;
                    if(type == DIRICHLET_X || type == DIRICHLET_Y){
                        constrain(type,Aloc,Bloc,i,condition->value,1,1);
                    }
                    else if(type == NEUMANN_X || type == NEUMANN_Y){
                        double dx, dy;
                        if(i==0) {
                            dx = x[i]-x[theSpace->n-1];
                            dy = y[i]-y[theSpace->n-1];}
                        else {
                            dx = x[i]-x[i-1];
                            dy = y[i]-y[i-1];}
                        constrain(type,Aloc,Bloc,i,condition->value, dx, dy);
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

void constrain(femBoundaryType type, double** A, double* B, int myNode, double value, double dx, double dy){
    double length = sqrt(dx*dx+dy*dy);
    double cos = dx/length;
    double sin = dy/length;
    switch (type){
        case DIRICHLET_X:
            femDirichlet(A,B,sizeof(A[0])/sizeof(double),2*myNode,value); break;
        case DIRICHLET_Y:  
            femDirichlet(A,B,sizeof(A[0])/sizeof(double),2*myNode+1,value); break;
        case NEUMANN_X: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            femNeumann(B,2*myNode,value); //line integral approximation 2point gauss
            break;                                                   // value = derivative of the flux on boundary
        case NEUMANN_Y: //value = ((3-sqrt(3))/3)*value*length/2
            value = ((3-sqrt(3))/3)*value*length/2;
            femNeumann(B,2*myNode+1,value); 
            break;
        /*A IMPOSER AVEC TABLEAU METHODE DES DEPLACEMENTS*/
        case NEUMANN_N:
            value = ((3-sqrt(3))/3)*value*length/2;
            femNeumann(B,2*myNode,value*sin); 
            femNeumann(B,2*myNode+1,value*cos); 
            break;
        case NEUMANN_T:
            value = ((3-sqrt(3))/3)*value*length/2;
            femNeumann(B,2*myNode,value*cos); 
            femNeumann(B,2*myNode+1,value*sin); 
            break;
        default: break;
    }
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
                    femDirichlet(A,B,sizeof(A[0])/sizeof(double),mapX[node],value); break;
                case DIRICHLET_Y:  
                    femDirichlet(A,B,sizeof(A[0])/sizeof(double),mapY[node],value); break;
                case NEUMANN_X:
                    femNeumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); //line integral approximation 2point gauss
                    break;                                                   // value = derivative of the flux on boundary
                case NEUMANN_Y:
                    femNeumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                /*A IMPOSER AVEC TABLEAU METHODE DES DEPLACEMENTS*/
                case NEUMANN_N:
                    femNeumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); 
                    femNeumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                case NEUMANN_T:
                    femNeumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); 
                    femNeumann(B,mapY[node],((3-sqrt(3))/3)*value*length/2); 
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

void femDirichlet(double **A, double *B, int size, int myNode, double myValue){
    int i;
    for(i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0.0; }
    for(i=0; i < size; i++)
        A[myNode][i] = 0.0;
    A[myNode][myNode] = 1.0;
    B[myNode] = myValue;

}

void femNeumann(double *B, int myNode, double myValue){
    B[myNode] += myValue;
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

void femLocalPlan(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight){
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

void femLocalAxsym(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight){
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