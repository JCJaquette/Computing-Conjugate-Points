#ifndef localVField_h
#define localVField_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
// #include <capdAlg/include/capd/vectalg/vectalgLib.h>

using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

class localVField
{
private:
  IMap *f;
  IMatrix A;
  IVector p;
  IMatrix Ainv;
  friend class localManifold;
  friend class localManifold_Eig;
  friend class boundaryValueProblem;
  friend class propagateManifold;
public:
//     We define psi(w) := p0 + A_0 w
//                 f(w) :=    - A_0^-1 * J \grad H ( \psi( w)) 
//           where    J  = [0 -I; I 0]
//     In the program, f(x) = - J \grad H ( x )
  localVField(IMap &f_,IMatrix A_,IVector p_){f=&f_; A=A_; p=p_; Ainv=gaussInverseMatrix(A);}
  localVField(){}
  IVector image(IVector w){return Ainv*(*f)(p+A*w);}
  IMatrix derivative(IVector w){return Ainv*(*f)[p+A*w]*A;}
  
  IVector operator()(IVector w){return image(w);}
  IMatrix operator[](IVector w){return derivative(w);}
};

#endif
