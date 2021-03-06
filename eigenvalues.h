 
#ifndef eigenvalues_h
#define eigenvalues_h

#include <iostream>
using namespace std;
#include "capd/capdlib.h"
// #include <capdAlg/include/capd/vectalg/vectalgLib.h>
#include "utils.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

IVector boundEigenvalues(IMatrix B);
vector < IVector >  boundSingleEigenvector(IMatrix A,const IVector v, interval lambda);
vector < IVector > krawczykEigenvector(IMatrix A,const IVector V, interval lambda , IVector H_vec );


IVector F_eigenvector(IMatrix A, IVector v, interval lambda );
IMatrix DF_eigenvector(IMatrix A, IVector v, interval lambda);
vector < IMatrix > boundEigenvectors(IMatrix A, IMatrix Q, IVector Lambda);


#endif
