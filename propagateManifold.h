#ifndef propagateManifold_h
#define propagateManifold_h


#include <iostream>
using namespace std;


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
//   Data 
  IMap *pf;         // The system of the 1st variation
  
  int manifold_subdivision = 15; // Used in **construct_Manifold_at_LPlus**

  
  IVector XY_pt;    // the CENTER of the unstable (X) and stable(Y) manifold coordinates of the heteroclinic connection  (in local coordinates)
  IVector XY_nbd;   // the neighborhood enclosure of the heteroclinic connection (in local coordinates).
  
  bool MAKE_PLOT = false;
  
  int step_size; // step size for the solver is set to 2^-step_size
  
//   IMatrix A_lin;
  
  localManifold_Eig *pUnstable;
  localManifold_Eig *pStable;
  int dimension;
  int order;    // Order of the taylor solver
  
// // // // // // // // // // // // // 
//   Methods 
  bool lastEuFrame(topFrame &A_frame , IVector endPoint_LPlus);  
  bool checkL_plus( IMatrix U_coord,interval eps_0,IVector eigenvalues , interval E_norm , IMatrix U_coord_pt, IMatrix U_coord_nbd);  
  bool checkL_plus_local( IMatrix Gamma, IMatrix Beta,interval eps_0,interval nu_1 , interval nu_n, interval EE_norm);
  interval compute_epsilon_beta( IMatrix Gamma, IMatrix Beta, interval EE_norm);
  
  vector<IMatrix> computeTotalTrajectory(int eigenvector_NUM, interval T, int grid);  
  IMatrix construct_A_lin(void);
  vector <IVector> construct_InitCondU(int eigenvector_NUM);
  
  localManifold_Eig construct_Manifold_at_LPlus( const IMatrix &last_Frame, IVector endPoint_LPlus);
  
  vector<IMatrix> projectionGammaBeta(  IMatrix& last_Frame  , const IMatrix & EFunction_Error );
  
public:
  propagateManifold(IMap &pf_, localManifold_Eig &pUnstable_, localManifold_Eig &pStable_, IVector XY_pt_,IVector XY_nbd_,int order_,int step_size_){pf = &pf_; pUnstable = &pUnstable_;pStable = &pStable_;XY_pt =XY_pt_; XY_nbd =XY_nbd_; order = order_; dimension = (*pUnstable).dim();step_size=step_size_;}
  
  void plotting(bool MAKE_PLOT_){MAKE_PLOT = MAKE_PLOT_;};
  

  
  int frameDet( interval L_minus, int grid,IVector endPoint_LPlus);
  
  

};


#endif
