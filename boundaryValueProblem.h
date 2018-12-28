#ifndef boundaryValueProblem_h
#define boundaryValueProblem_h


#include <iostream>
using namespace std;

#include "capd/capdlib.h"
#include "utils.h"
#include "localVField.h"
#include "eigenvalues.h"
#include "localManifold.h"

using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

class boundaryValueProblem 
{
private:
  IMap *pf;
  IMap *pf_minus;

  bool SUCCESS;
  IVector verified_XY_pt;
  IVector verified_XY_nbd;
  
  localManifold *pStable;
  localManifold *pUnstable;
  int dimension;
  int order;
  int frozen; 
public:
  boundaryValueProblem(IMap &pf_,IMap &pf_minus_, localManifold &pStable_,localManifold &pUnstable_, int order_, int frozen_){pf = &pf_;pf_minus = &pf_minus_; pStable = &pStable_; pUnstable = &pUnstable_;order = order_; frozen = frozen_;dimension = (*pStable).dim();SUCCESS=0;}


  IVector Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) ; 
  IMatrix DGxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE);  
  
  IMatrix DG_combine(IVector G, IMatrix DGX, IMatrix DGY,IVector X_pt, IVector X_nbd);//TODO REMOVE FROZEN //TODO
  IVector calcG(IVector X, IVector Y, interval T);
  IVector calcDG(IVector X, IVector Y, interval T);
  interval FindTime( IVector X_mid, interval T);
  IVector NormBound( IVector XY,interval T); // TODO pt/nbd Form
  IVector localNormBound( IVector XY,interval T, bool STABLE ); // TODO pt/nbd Form
  
  IVector NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T) ;//TODO REMOVE FROZEN //TODO
  
  bool Verify( IVector XY_pt, IVector XY_nbd  ,interval T){return 0;};// TODO Write this function 
  bool checkProof(void){return SUCCESS;};
  
  vector < IVector > breakUpXY( IVector XY);
  vector < IVector > breakUpXY_gen( IVector XY);
};


#endif
