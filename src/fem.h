
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,DIRICHLET_N,DIRICHLET_T,NEUMANN_X,NEUMANN_Y,NEUMANN_N,NEUMANN_T} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;
typedef enum {FEM_FULL,FEM_BAND,FEM_FRONT, FEM_BPOST} femSolverType;


typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;  // forces nodales globales
    double **A; // matrice pleine globale
    int size;
} femFullSystem;

typedef struct {
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    double* b; // forces nodales
    double** active; // matrice active
    int nActive; // nombre de noeuds actifs
    int* activeMap; // coordonnées globales des noeuds actifs
    int* activeState; // 0: occupé, 1: libre
    int** disparus; //disparu[iElem][iLoc] = 0 si le noeud iLoc de l'élément iElem est actif, globalCoord si à éliminer
    int** nouveaux; //nouveau[iElem][iLoc] = globalCoord si le noeud iLoc de l'élément iElem est à ajouter
    double** stock; //stock[nElem][nElem] = stock des éléments éliminés
    double* pivots; //pivots[nElem] = stock des pivots
} femFrontalSolver;

typedef struct {
    int size;
    femFullSystem *local;
    femRenumType renum;
    femSolverType solverType;
    femFullSystem *fullSystem;
    femBandSystem *bandSystem;
    femFrontalSolver *frontSolver;
} femSystem;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;


typedef struct {
    double E,nu,rho,g;
    double A,B,C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femSystem *system;
} femProblem;


void                geoInitialize();
femGeo*             geoGetGeometry();
double              geoSize(double x, double y);
double              geoSizeDefault(double x, double y);
void                geoSetSizeCallback(double (*geoSize)(double x, double y));
void                geoMeshGenerate();
void                geoMeshImport();
void                geoMeshPrint();
void                geoMeshWrite(const char *filename);
void                geoMeshRead(const char *filename);
void                geoSetDomainName(int iDomain, char *name);
int                 geoGetDomain(char *name);
void                geoFinalize();
void                geoFree();


void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
double*             femElasticitySolve(femProblem *theProblem);
void                femElasticityWrite(femProblem *theProbconst, const char *filename);


void                femFieldWrite(int size, int shift, double* value, const char *filename);
int                 femFieldRead(int* size, int shift, double* value, const char *filename);

femIntegration*     femIntegrationCreate(int n, femElementType type);
void                femIntegrationFree(femIntegration *theRule);

femDiscrete*        femDiscreteCreate(int n, femElementType type);
void                femDiscreteFree(femDiscrete* mySpace);
void                femDiscretePrint(femDiscrete* mySpace);
void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
void                femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
double              fmin(double a, double b);    
double              fmax(double a, double b);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);


#endif
