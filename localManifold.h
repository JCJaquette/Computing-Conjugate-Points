#ifndef localManifold_h
#define localManifold_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
#include "utils.h"
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
  friend class boundaryValueProblem;
  friend class propagateManifold;
  
public:
  localManifold(localVField &pF_, IVector U_flat_global_,interval L_, bool stable_){pF = &pF_; U_flat_global = U_flat_global_; L=L_; stable = stable_;dimension=(*pF).A.numberOfRows();subdivisionNUM = 8;}

  int dim( ){return dimension;}
  IVector constructU( IVector U_flat);
  IVector getPoint( IVector X_pt);  
  
  
  
  IMatrix boundDFU( IVector U); 			
  
  void constructDW( ); 
  
  bool checkRateCondition(IVector U_flat );
  bool checkRateCondition(IMatrix DFU);
  
  vector < IVector > getPointNbd( IVector XY_pt, IVector XY_nbd);
  IVector projectPoint( IVector XY_pt);
  
  IVector getEigenvector( int eigen_NUM){ return getColumn( (*pF).A, dimension, eigen_NUM);};

  bool checkIsolatingBlock( IMatrix DFU, const IVector &U); 
  
  interval ErrorEigenfunction( void);
  
  interval getRadius( void ){ return getMax( abs(U_flat_global));};
  
  bool checkConditions( IVector U_flat );
};

#endif
