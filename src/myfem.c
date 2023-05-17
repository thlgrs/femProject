#include "fem.h"
#include "myfem.h"


femFrontalSolver* femFrontalSolverCreate(int size, int nLoc){
    femFrontalSolver *mySolver = malloc(sizeof(femFrontalSolver));
    
    mySolver->size = size;
    mySolver->nLoc = nLoc;
    mySolver->old = calloc(size,sizeof(int*));
    mySolver->A = calloc(size,sizeof(double*));
    mySolver->B = calloc(size,sizeof(double));
    mySolver->ALoc = calloc(2*nLoc, sizeof(double*));    
    mySolver->BLoc = calloc(2*nLoc, sizeof(double));
    
    int i;
    for(i=0; i<size; i++){
        mySolver->old[i] = calloc(nLoc,sizeof(int));
        mySolver->A[i] = calloc(size,sizeof(double));
    }
    for(i=0; i<2*nLoc; i++)
        mySolver->ALoc[i] = calloc(2*nLoc, sizeof(double));
        mySolver->BLoc[i] = 0;
        
    return(mySolver);
};

void femFrontalSolverFree(femFrontalSolver *mySolver){
    int i;
    for(i=0; i<2*mySolver->nLoc; i++)
        free(mySolver->ALoc[i]);
    for(i=0; i<mySolver->size; i++){
        free(mySolver->A[i]);
        free(mySolver->old[i]);
    }
    free(mySolver->ALoc);
    free(mySolver->BLoc);
    free(mySolver->A);
    free(mySolver->B);
    free(mySolver->old);
    free(mySolver);
    
};

void femFrontActivity(femGeo *theGeometry, int **old){
    femMesh *theElements = theGeometry->theElements;
    femNodes *theNodes = theGeometry->theElements->nodes;
    int nElem = theElements->nElem;
    int *elem = theElements->elem;
    int nLocal = theElements->nLocalNode;
    
    int found;
    for(int iElem=0; iElem<nElem; ++iElem){
        int first = elem[iElem-1];
        int second = elem[iElem];
        for (int i=0; i<nLocal; i++){
            found = 0;
            double xfirst = theNodes->X[first*nLocal+i];
            double yfirst = theNodes->Y[first*nLocal+i];
            for(int j=0; j<nLocal; j++){
                double xsecond = theNodes->X[second*nLocal+j];
                double ysecond = theNodes->Y[second*nLocal+j];
                if(xfirst == xsecond && yfirst == ysecond){
                    found++;
                    break;
                }
            }
            old[iElem][i] = found;
        }
        for(int i=0; i<nLocal; i++)
            printf("%d ", old[iElem][i]);
            printf("\n");
    }
};

void femGetMap(int* elem, int iElem, int *map, int nLocal){
    for(int l=0; l < nLocal; l++)
        map[l] = elem[iElem*nLocal+l];
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
    femFrontalSolver *theSolver = theProblem->system->frontSolver;
    femMesh *theElements = theProblem->geometry->theElements;
    int nElem = theElements->nElem;
    int nLocal = theSolver->nLoc;
    int* map = malloc(nLocal*sizeof(int));
    int i, iElem;
    double **AGlob = theSolver->A;
    double *BGlob = theSolver->B;
    double **ALoc = theSolver->ALoc;
    double *BLoc = theSolver->BLoc;
    int **old = theSolver->old;

    femFrontActivity(theProblem->geometry, theSolver->old);
    for(iElem = 0; iElem < nElem; iElem++){
        for(i=0; i<nLocal; i++)
            map[i] = theElements->elem[iElem*nLocal+i];
        femAssembleLocal(ALoc, BLoc, AGlob, BGlob, map, nLocal);
        for(i=0; i<sizeof(theSolver->old[iElem])/sizeof(int); i++){
            femGausFrontal(ALoc, BLoc, old[iElem][i]);
        }
        femInjectGlobal(ALoc, BLoc, AGlob, BGlob, map, nLocal);
    }
    free(map);
    return BGlob;
};

void femGlobalSystemAssemble(femProblem *theProblem, double **A, double *B){
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theElements = theGeometry->theElements;
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;

    int nLocal = theElements->nLocalNode;
    for (iElem = 0; iElem < theElements->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theElements->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];}
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
                for (i = 0; i < theSpace->n; i++) { 
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
                                                dphidx[i] * c * x[i] * dphidx[j]) * jac * weight; }}
                
                for (i = 0; i < theSpace->n; i++){
                    B[mapY[i]] -= x[i] * phi[i] * rho * g * jac * weight;
                }
            }
            else if(theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
            {
                for (i = 0; i < theSpace->n; i++) { 
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
                }
                for (i = 0; i < theSpace->n; i++){
                    B[mapY[i]] -= phi[i] * rho * g * jac * weight; 
                }
        }   
            
    }

};

void femBoundaryConstrain(femProblem *theProblem, double **A, double *B){
    femBoundaryCondition **conditions = theProblem->conditions;
    femIntegration *theRule           = theProblem->rule;
    femDiscrete    *theSpace          = theProblem->space;
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
                    femNeumann(B,mapX[node],((3-sqrt(3))/3)*value*length/2); 
                    break;
                case NEUMANN_Y:
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

void femFulltoBand(double **Aglob, double **A, int size, int band){
    int i,j,k,jend;
    for(i=0; i<size; i++){
        printf("loop1\n");
        for(j=i; j<jend; j++){
            //printf("loop2\n");
            if(fabs(A[i][j])>1e-10){
                for(k=0; k<band; k++){
                    if(j+k<size){
                        printf("loop3\n");
                        A[i][k] = Aglob[i][j];
                    }
                }
            }
        }
    }
}



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
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    myBandSystem->Aglob = malloc(sizeof(double*)*size);
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;    
    myBandSystem->Aglob[0] = elem + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
        myBandSystem->Aglob[i] = myBandSystem->Aglob[i-1] + size;
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
    A    = myBand->Aglob;
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

void femElasticityFree(femProblem *theProblem)
{
    femSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    free(theProblem->conditions);
    free(theProblem->constrainedNodes);
    free(theProblem);
}

double *femElasticitySolve(femProblem *theProblem)
{   
    switch(theProblem->system->solverType) {
        
        case FEM_FULL : 
            femGlobalSystemAssemble(theProblem, theProblem->system->fullSolver->A, theProblem->system->fullSolver->B);
            femBoundaryConstrain(theProblem, theProblem->system->fullSolver->A, theProblem->system->fullSolver->B);
            return femFullSystemEliminate(theProblem->system->fullSolver);
            break;
        case FEM_BAND :
            femGlobalSystemAssemble(theProblem, theProblem->system->bandSolver->Aglob, theProblem->system->bandSolver->B);
            femBoundaryConstrain(theProblem, theProblem->system->bandSolver->Aglob, theProblem->system->bandSolver->B);
            femFulltoBand(theProblem->system->bandSolver->Aglob, theProblem->system->bandSolver->A, theProblem->system->bandSolver->size, theProblem->system->bandSolver->band);
            return femBandSystemEliminate(theProblem->system->bandSolver);
            break;
        case FEM_FRONT :
            femGlobalSystemAssemble(theProblem, theProblem->system->frontSolver->A, theProblem->system->frontSolver->B);
            femBoundaryConstrain(theProblem, theProblem->system->frontSolver->A, theProblem->system->frontSolver->B);
            return femFrontalSolve(theProblem);
            break;
        default :
            Error("Unexpected solver type");
    };
}

femSystem* femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femGeo *theGeometry){
    femSystem* system = malloc(sizeof(femSystem));
    system->solverType = iSolver;
    
    switch(iSolver) {
        case FEM_FULL :
            system->fullSolver = femFullSystemCreate(size); break;
        case FEM_BAND : 
            femRenumberNodes(theGeometry, iRenum);
            int band = femMeshComputeBand(theGeometry->theElements);
            system->bandSolver = femBandSystemCreate(size,band);break;
        case FEM_FRONT :
            femRenumberElem(theGeometry, iRenum);
            int nLocalNode = theGeometry->theElements->nLocalNode;
            system->frontSolver = femFrontalSolverCreate(size,nLocalNode); break;
    }
    return system;
}

void femSystemFree(femSystem* mySystem){
    switch(mySystem->solverType) {
        case FEM_FULL : femFullSystemFree(mySystem->fullSolver); break;
        case FEM_BAND : femBandSystemFree(mySystem->bandSolver); break;
        case FEM_FRONT : femFrontalSolverFree(mySystem->frontSolver); break;
    }
    free(mySystem);
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