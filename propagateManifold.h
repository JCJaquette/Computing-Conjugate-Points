#ifndef propagateManifold_h
#define propagateManifold_h


#include <iostream>
using namespace std;

// #ifdef _OPENMP
#include <omp.h>
// #endif
#include "capd/capdlib.h"
#include "utils.h"
#include "localVField.h"
#include "eigenvalues.h"
#include "localManifold.h"
#include "topFrame.h"

using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

class propagateManifold 
{
private:
  IMap *pf;
  
  vector <IMap> list_of_maps; 

  IVector InitialConditions_local;
  IVector XY_pt;
  IVector XY_nbd;
  
  int step_size;
  
//   IMatrix A_lin;
  
  localManifold *pUnstable;
  localManifold *pStable;
  int dimension;
  int order;

public:
  propagateManifold(IMap &pf_, localManifold &pUnstable_, localManifold &pStable_, IVector XY_pt_,IVector XY_nbd_,int order_,int step_size_){pf = &pf_; pUnstable = &pUnstable_;pStable = &pStable_;XY_pt =XY_pt_; XY_nbd =XY_nbd_; order = order_; dimension = (*pUnstable).dim();InitialConditions_local = XY_pt+XY_nbd;step_size=step_size_;}
  
  vector<IMatrix> computeTotalTrajectory(int eigenvector_NUM, interval T, int grid);
  
  
  
  IMatrix construct_A_lin(void);
  vector <IVector> construct_InitCondU(int eigenvector_NUM);
  
  int frameDet(interval T, interval L_plus, int grid,IVector endPoint_LPlus);
  
  bool lastEuFrame(topFrame &A_frame , IVector endPoint_LPlus);
  
  bool checkL_plus( IMatrix U_coord,interval eps_0,IVector eigenvalues );
  
  bool checkL_plus_local( IMatrix Gamma, IMatrix Beta,interval eps_0,interval nu_1 , interval nu_n);
  
  interval eps_beta( IMatrix Gamma, IMatrix Beta);

};


#endif
