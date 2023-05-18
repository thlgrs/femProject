#include "fem.h"

int                 femCompare(const void *i1, const void *i2);
void                femRenumberNodes(femGeo *theGeometry, femRenumType renumType);
void                femRenumberElem(femGeo *theGeometry, femRenumType renumType);


femBandSystem*      femBandSystemCreate(int size, int band);
int                 femComputeBand(femMesh *theElements);
void                femBandSystemFree(femBandSystem* myBand);
void                femBandSystemInit(femBandSystem *myBand);
double*             femBandSystemEliminate(femBandSystem *myBand);
void                femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
void                femLocalPlan(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);
void                femLocalAxsym(femDiscrete* space, double **A, double* B, double* x, double* phi, double *dx, double *dy, double* coeff, double jac, double weight);

femProblem*         femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType iSolver, femRenumType iRenum);
void                femElasticityFree(femProblem *theProblem);
void                femGlobalSystemAssemble(femProblem *theProblem);
void                constrain(femBoundaryType type, femFullSystem* system, int myNode, double value, double dx, double dy);
void                femDirichlet(double **A, double *B, int size, int myNode, double myValue);
void                femNeumann(double *B, int myNode, double myValue);


double*             femElasticitySolve(femProblem *myProblem);