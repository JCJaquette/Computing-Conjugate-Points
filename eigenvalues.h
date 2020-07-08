 
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
vector < IVector >  boundSingleEigenvector(IMatrix A, IVector v, interval lambda, interval local_norm_sq);
vector < IVector > krawczykEigenvector(IMatrix A, IVector V, interval lambda , IVector H_vec, interval local_norm_sq);


IVector F_eigenvector(IMatrix A, IVector v, interval lambda, interval local_norm_sq);
IMatrix DF_eigenvector(IMatrix A, IVector v, interval lambda);
vector < IMatrix > boundEigenvectors(IMatrix A, IMatrix Q, IVector Lambda);


#endif
