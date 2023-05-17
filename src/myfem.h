#include "fem.h"

int                 femCompare(const void *i1, const void *i2);
void                femRenumberNodes(femGeo *theGeometry, femRenumType renumType);
void                femRenumberElem(femGeo *theGeometry, femRenumType renumType);
int                 femMeshComputeBand(femMesh *theMesh);

femBandSystem*      femBandSystemCreate(int size, int band);
void                femBandSystemFree(femBandSystem* myBand);
void                femBandSystemInit(femBandSystem *myBand);
double*             femBandSystemEliminate(femBandSystem *myBand);
void                femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);

femFrontalSolver*   femFrontalSolverCreate(int size, int nLoc);
void                femFrontalSolverFree(femFrontalSolver* mySolver);
double*             femFrontalSolve(femProblem *thePorblem);
void                femFrontActivity(femGeo *theGeometry, int **old);
void                femGetMap(int* elem, int iElem, int *map, int nLocal);
void                femAssembleLocal(double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);
void                femInjectGlobal (double **ALoc, double *BLoc, double **AGlob, double *BGlob, int *map, int nLocal);
double*             femGauss(double **A, double *B, int size);


femSystem*          femSystemCreate(int size, femSolverType iSolver, femRenumType iRenum, femGeo *theGeometry);
void                femSystemFree(femSystem* mySystem);

femProblem*         femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType iSolver, femRenumType iRenum);
void                femElasticityFree(femProblem *theProblem);
void                femGlobalSystemAssemble(femProblem *theProblem, double **A, double *B);
void                femBoundaryConstrain(femProblem *theProblem, double **A, double *B);
void                femDirichlet(double **A, double *B, int size, int myNode, double myValue);
void                femNeumann(double *B, int myNode, double myValue);
double*             femElasticitySolve(femProblem *theProblem);