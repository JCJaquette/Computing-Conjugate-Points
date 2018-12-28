#ifndef localVField_h
#define localVField_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"

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
  friend class boundaryValueProblem;
  friend class propagateManifold;
public:
  localVField(IMap &f_,IMatrix A_,IVector p_){f=&f_; A=A_; p=p_; Ainv=gaussInverseMatrix(A);}
  localVField(){}
  IVector image(IVector w){return Ainv*(*f)(p+A*w);}
  IMatrix derivative(IVector w){return Ainv*(*f)[p+A*w]*A;}
  
  IVector operator()(IVector w){return image(w);}
  IMatrix operator[](IVector w){return derivative(w);}
};

#endif