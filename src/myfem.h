#include "fem.h"

/*MAIN FUNCTIONS*/
femProblem* femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType iSolver, femRenumType iRenum);
femSystem*  femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femGeo *theGeometry);
void        femElasticitySet(femProblem *theProblem);
double*     femElasticitySolve(femProblem *myProblem);
void        femElasticityFree(femProblem *theProblem);


/*SECUNDARY FUNCTIONS*/
void    femSystemFree(femSystem* mySystem);
int     computeBand(femMesh *theElements);
int     compare(const void *i1, const void *i2);
void    renumberNodes(femGeo *theGeometry, femRenumType renumType);
void    renumberElem(femGeo *theGeometry, femRenumType renumType);
int     nodeIndex(double x, double y, femNodes* theNodes);
void    localPlannar(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);
void    localAxisym(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);
void    constrain(femBoundaryType type, femFullSystem* system, int myNode, double value, double dx, double dy);
void    dirichlet(double **A, double *B, int size, int myNode, double myValue);
void    neumann(double *B, int myNode, double myValue);

/*BAND SYSTEM*/
femBandSystem*      femBandSystemCreate(int size, int band);
void                femBandSystemFree(femBandSystem* myBand);
void                femBandSystemInit(femBandSystem *myBand);
double*             femBandSystemEliminate(femBandSystem *myBand);
void                femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);

/*FRONTAL SOLVER*/
femFrontalSolver*   femFrontalSolverCreate(int size, int nLoc);
void                femFrontalSolverFree(femFrontalSolver* mySolver);
double*             femFrontalSolve(femProblem *thePorblem);
int                 nActive(femMesh *theElements, int **disparus, int **nouveaux);
void                femAssembleLocal(double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);
void                femInjectGlobal (double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);




void                femBoundaryConstrain(femProblem *theProblem, double **A, double *B);