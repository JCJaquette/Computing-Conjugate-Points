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
  IFunction *p_energy_proj;

  bool SUCCESS;
  IVector verified_XY_pt;
  IVector verified_XY_nbd;
  
  localManifold *pStable;
  localManifold *pUnstable;
  int dimension;
  int num_middle_points;  // This NEEDS to be \geq 1 
  int order;
public:
  boundaryValueProblem(IMap &pf_,IMap &pf_minus_, IFunction &p_energy_proj_,localManifold &pStable_,localManifold &pUnstable_, int order_, int num_middle_points_){pf = &pf_;pf_minus = &pf_minus_;p_energy_proj=&p_energy_proj_; pStable = &pStable_; pUnstable = &pUnstable_;order = order_; num_middle_points = num_middle_points_;dimension = (*pStable).dim();SUCCESS=0;}

  void Integrate_point( IVector coord_pt, IVector coord_nbd,interval T, bool FORWARD,IVector &vector_out, IMatrix &derivative) ; 
  IVector Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) ; 
  IMatrix DGxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE);  
  
  IMatrix DG_combine(IMatrix DGX, IMatrix DGY,IVector X_pt, IVector X_nbd); 
  IVector calcG(IVector X, IVector Y, interval T);
  IVector calcDG(IVector X, IVector Y, interval T);
  interval FindTime( IVector X_mid, interval T);
  IVector NormBound( IVector XY,interval T); // TODO pt/nbd Form
  IVector localNormBound( IVector XY,interval T, bool STABLE ); // TODO pt/nbd Form
  
  IVector Compute_G(const vector <IVector> &points,interval T);
  IMatrix Compute_DG(const vector <IVector> &points,const vector <IVector> &neighborhoods,interval integration_time);
  IVector Construct_G( vector < IVector > G_forward, vector < IVector > G_backwards, vector <IVector> points);
  IMatrix Construct_DG(const vector <IMatrix> &DG_forward, const vector <IMatrix> &DG_backwards,const  vector <IVector> &points, const vector <IVector> &neighborhoods,const vector < IVector> &time_derivatives);
  IVector Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods);// TODO REMOVE
  IVector Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods,interval integration_time, bool ADD); 
  IVector Construct_Initial_Vector(vector <IVector> points, interval integration_time);// TODO REMOVE
  vector <IVector> Deconstruct_Output_Vector(IVector initial_vector);
  
  
  IVector NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T) ; 
  vector <IVector> NewtonStep(vector <IVector> &XY_pt, vector <IVector> &XY_nbd  ,interval &T, interval time_nbd);   
  bool Verify( IVector XY_pt, IVector XY_nbd  ,interval T){return 0;};// TODO Write this function 
  bool checkProof(void){return SUCCESS;};
  void setMiddlePoints(int shots){num_middle_points = shots;};
  
  vector < IVector > breakUpXY( IVector XY);
  vector < IVector > breakUpXY_gen( IVector XY);
};


#endif
