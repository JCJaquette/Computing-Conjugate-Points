#ifndef localManifold_h
#define localManifold_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
#include "utils.h"
#include "eigenvalues.h"
#include "localVField.h"
// #include "eigenvalues.h"


using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

class localManifold 
{
private:
    
protected:
  localVField *pF;  
  interval L; // This cooresponds to \vartheta in the paper. It is the angle of our cones. 
  IMatrix DW;
  int dimension;   //  dimension = 2*n
  bool stable; 	
  int subdivisionNUM ; // Used to bounding the derivative of the function. Runtime scales like subdivNUM^(n)  
  IVector U_flat_global; // The size of the manifold  
  interval xi;
  
  friend class boundaryValueProblem;
  
  
public:
  localManifold(localVField &pF_, IVector U_flat_global_,interval L_, bool stable_, int subdivisionNUM_){pF = &pF_; U_flat_global = U_flat_global_; L=L_; stable = stable_;dimension=(*pF).A.numberOfRows();subdivisionNUM = subdivisionNUM_;constructDW( );}

  int dim( ){return dimension;}
  IVector constructU( IVector U_flat);
  IVector getPoint( IVector X_pt);  
  
  IMatrix boundDFU( IVector U); 			
  
  void constructDW( );   
  bool checkRateCondition(IMatrix DFU);
  
  vector < IVector > getPointNbd( IVector XY_pt, IVector XY_nbd);
  IVector projectPoint( IVector XY_pt);

  bool checkIsolatingBlock( IMatrix DFU, const IVector &U);   
  interval getRadius( void ){ return getMax( abs(U_flat_global));};  // Previously used in "bvpMultipleShooting" class
  bool checkConditions( void );  
  
  IVector containmentRatios( IVector point_test);  
  void updateSize( IVector U_flat_new){U_flat_global = U_flat_new;};
  
  void printPF(void){cout << " pF = " << pF << endl;};
  
  IVector getEigenvector( int eigen_NUM){ return getColumn( (*pF).A, dimension, eigen_NUM);}; // Returns non-validated eigenvector.

};

// // // // // // // // // //
//   localManifold_Eig     //
// // // // // // // // // // 
//      This derived class is used in computing and bounding eigenfunctions.  
//      It is used in/by the  ""propagateManifold"" class
class localManifold_Eig : public localManifold
{
private:
    IVector eigenvalues;
    interval K_store;   
    interval eps_unscaled;
    IMatrix Eigenvector_Error; // Error for unit norm eigenvectors;
    IMatrix Eu_m_Error_Final;  // Error for the unstable eigenfunctions at minus infinity 
    IMatrix Eigenfunction_Error_plus_infty;  // Error for the all eigenfunctions at plus infinity 
    
    friend class propagateManifold;
    
public:
    localManifold_Eig(localVField &pF_, IVector U_flat_global_,interval L_, bool stable_, int subdivisionNUM_);
    localManifold_Eig(const localManifold & lM_);
    

    
        
    interval computeK( void );  
    
    void ErrorEigenfunction( void);  
    void computeEigenError_minus_infty( void);
    void computeEigenError_plus_infty( void);

    IVector getEigenError_minus_infty(int columnNumber);

};

#endif
