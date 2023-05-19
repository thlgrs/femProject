/*
myfem.c

*/

#include "fem.h"

/*MAIN FUNCTIONS*/
femProblem* femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType iSolver, femRenumType iRenum);
void        femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femProblem *theProblem);
void        femElasticitySet(femProblem *theProblem);
void        femBoundaryConstrain(femProblem *theProblem, double **A, double *B);
double*     femElasticitySolve(femProblem *theProblem);
void        femElasticityFree(femProblem *theProblem);
void        femSystemFree(femSystem* theSystem);

/*SECUNDARY FUNCTIONS*/
int     computeBand(femMesh *theElements);
int     compare(const void *i1, const void *i2);
void    renumberNodes(femProblem *theProblem, femRenumType renumType);
void    renumberElem(femGeo *theGeometry, femRenumType renumType);
int     nodeIndex(double x, double y, femNodes* theNodes);
void    localPlannar(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);
void    localAxisym(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);
void    constrain(femBoundaryType type, femFullSystem* system, int myNode, double value, double dx, double dy);
void    constrainb(femBoundaryType type, femBandSystem* system, int myNode, double value, double dx, double dy);
void    dirichlet(double **A, double *B, int size, int myNode, double myValue);
void    neumann(double *B, int myNode, double myValue);
double  longueur(double* x, double* y, int size);

/*BAND SYSTEM*/
femBandSystem* femBandSystemCreate(int size, int band);
void           femBandSystemFree(femBandSystem* myBand);
void           femBandSystemInit(femBandSystem *myBand);
double*        femBandSystemEliminate(femBandSystem *myBand);
void           femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);

/*POST BAND*/
femBandSystem* femPostBandSystemCreate(int size);
void           femPostBandSet(femProblem* theProblem);
void           setBand(femBandSystem* myBandSystem);


/*FRONTAL SOLVER not-implemented*/
femFrontalSolver* femFrontalSolverCreate(int size, int nLoc);
void              femFrontalSolverFree(femFrontalSolver* mySolver);
double*           femFrontalSolve(femProblem *thePorblem);
int               nActive(femMesh *theElements, int **disparus, int **nouveaux);
void              femAssembleLocal(double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);
void              femInjectGlobal (double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);