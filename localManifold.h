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
  localVField *pF;  
  interval L; // For rate conditions 
  IMatrix DW;
  int dimension;
  bool stable; 	
  int subdivisionNUM ;
  IVector U_flat_global; // TODO Replace this in how class is called
  interval xi;
  
  IVector eigenvalues;
  
  interval K_store;
  interval eps_unscaled;
  IMatrix Eigenvector_Error; // Error for unit norm eigenvectors;  
  IMatrix Eu_m_Error_Final;  // Error for the unstable eigenfucntions at minus infinity
  
  IMatrix Eigenfunction_Error_plus_infty;  // Error for the all eigenfucntions at plus infinity
  
  friend class boundaryValueProblem;
  friend class propagateManifold;
  
public:
  localManifold(localVField &pF_, IVector U_flat_global_,interval L_, bool stable_, int subdivisionNUM_){pF = &pF_; U_flat_global = U_flat_global_; L=L_; stable = stable_;dimension=(*pF).A.numberOfRows();subdivisionNUM = subdivisionNUM_;constructDW( );}

  int dim( ){return dimension;}
  IVector constructU( IVector U_flat);
  IVector getPoint( IVector X_pt);  
  
  
  
  IMatrix boundDFU( IVector U); 			
  
  interval boundDFU_proj( IVector U);
  
  void constructDW( ); 
  
  bool checkRateCondition(IVector U_flat );
  bool checkRateCondition(IMatrix DFU);
  
  vector < IVector > getPointNbd( IVector XY_pt, IVector XY_nbd);
  IVector projectPoint( IVector XY_pt);
  
  IVector getEigenvector( int eigen_NUM){ return getColumn( (*pF).A, dimension, eigen_NUM);};

  bool checkIsolatingBlock( IMatrix DFU, const IVector &U); 
  
  interval ErrorEigenfunction( void);  
  void ErrorEigenfunctionTotal_minus_infty( void);
  void ErrorEigenfunctionTotal_plus_infty( void);
  void computeEigenError_minus_infty( void){ErrorEigenfunction( ); ErrorEigenfunctionTotal_minus_infty( );};
  void computeEigenError_plus_infty( void){ErrorEigenfunction( ); ErrorEigenfunctionTotal_plus_infty( );};
  
  IVector getEigenError_minus_infty(int columnNumber);
  
  interval getRadius( void ){ return getMax( abs(U_flat_global));};
  
  bool checkConditions( IVector U_flat );
  
  interval computeK( void );
  
};

#endif
