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
    
protected:
  IMap *pf;
  IMap *pf_minus;
  

  bool SUCCESS;
  
  localManifold *pStable;
  localManifold *pUnstable;
  
  int dimension;
  
  int order;
  
//   IVector verified_XY_pt;
//   IVector verified_XY_nbd;
  
public:
  boundaryValueProblem(IMap &pf_,IMap &pf_minus_,localManifold &pStable_,localManifold &pUnstable_, int order_){pf = &pf_;pf_minus = &pf_minus_; pStable = &pStable_; pUnstable = &pUnstable_;order = order_;dimension = (*pStable).dim();SUCCESS=0;}

  
  IVector Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) ; 
  IMatrix DGxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE);  
  
  
  interval FindTime( IVector X_mid, interval T);
  IVector NormBound( IVector XY,interval T);            // TODO pt/nbd Form
  IVector localNormBound( IVector XY,interval T, bool STABLE ); // TODO pt/nbd Form
  
  IMatrix DG_combine(IMatrix DGX, IMatrix DGY,IVector X_pt, IVector X_nbd); 
  
  
  IVector NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T) ; 
  
  
  bool checkProof(void){return SUCCESS;};
  
  
  vector < IVector > breakUpXY( IVector XY);
  vector < IVector > breakUpXY_gen( IVector XY);
};


// simple derived class
class bvpMultipleShooting : public boundaryValueProblem
{
private:
    int num_middle_points;  // This NEEDS to be \geq 1 
    IFunction *p_energy_proj;
public:
    bvpMultipleShooting(IMap &pf_,IMap &pf_minus_, IFunction &p_energy_proj_,localManifold &pStable_,localManifold &pUnstable_, int order_, int num_middle_points_);
    
    void Integrate_point( IVector coord_pt, IVector coord_nbd,interval T, bool FORWARD,IVector &vector_out, IMatrix &derivative) ;    
    vector <IVector> NewtonStep(vector <IVector> &XY_pt, vector <IVector> &XY_nbd  ,interval &T, interval time_nbd);       
    
    IVector Compute_G(const vector <IVector> &points,interval T);                                                         //     <<>> Multiple Shooting Function <<>>
    IMatrix Compute_DG(const vector <IVector> &points,const vector <IVector> &neighborhoods,interval integration_time);   //     <<>> Multiple Shooting Function <<>>
    IVector Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods,interval integration_time, bool ADD); //     <<>> Multiple Shooting Function <<>> 
    vector <IVector> Deconstruct_Output_Vector(IVector initial_vector);
    
    IVector Construct_G( vector < IVector > G_forward, vector < IVector > G_backwards, vector <IVector> points);
    IMatrix Construct_DG(const vector <IMatrix> &DG_forward, const vector <IMatrix> &DG_backwards,const  vector <IVector> &points, const vector <IVector> &neighborhoods,const vector < IVector> &time_derivatives);    
    
    bool Verify( IVector XY_pt, IVector XY_nbd  ,interval T){return 0;};// TODO Write this function 
    
    void setMiddlePoints(int shots);
};
    



#endif
