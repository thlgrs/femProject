#include "fem.h"
#include "myfem.h"

void       systemAssemble(femProblem *theProblem);
femSystem* systemConstrain(femSystem *theSystem, femProblem *myProblem);
double*    systemSolve(femProblem *theProblem);