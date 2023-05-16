/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.c"
#include "myfem.c"

int main(void)
{   
    femSolverType solverType = FEM_FULL;
    femRenumType renumType = FEM_YNUM;

    femGeo* theGeometry = geoGetGeometry();  
    geoMeshRead("../data/mesh.txt");
    
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt", solverType, renumType);
    femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem); 
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}

 
