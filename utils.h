#ifndef utils_h
#define utils_h

#include <stdlib.h>
#include <utility>
#include <iostream>
using namespace std;

#include "capd/capdlib.h"

using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

DMatrix coordinateChange(DMatrix Df);
void 	bubbleSortEigenvectors(DVector &v, DMatrix &A);
void 	bubbleSort(DVector &v);
void 	swapValues(DVector &v,int i, int j);
void 	swapColumns(DMatrix &A,int i,int j);
IVector toInterval(DVector x);
IMatrix toInterval(DMatrix A);
interval part( interval x, int k, int N);
void 	plot(interval x, interval y,ofstream &file);
vector<IVector> getTrajectory(C0Rect2Set &s,interval T,int grid,ITimeMap &timeMap,IOdeSolver &solver);
vector<IMatrix> getTotalTrajectory(C0Rect2Set &s,interval T,int grid,ITimeMap &timeMap,IOdeSolver &solver);
void 	constructMatrix(const interval &time,const IVector &left_end,const IVector &right_end,const IVector &v_value,const IVector &v_deriv, IMatrix &matrix_mod);
IVector getColumn(const IMatrix &A, int dimension,int column);

interval omega( IVector V, IVector W, int n);

IMatrix identityMat( int n);

interval ml(const IMatrix &A);

 IVector getSubdivision(const IVector &U, const vector < int > &index_list , const vector < int > &part_list , const int &subdivisionNUM);
 int int_pow( int base, int exp);
 vector < IMatrix > blockDecompose( const IMatrix M , int dimension);
 interval getMax(const IVector &V);
 int getMaxIndex(const IVector &V);
 void print( vector <IVector> vectors );
 void print( vector <IMatrix> matrix_list);
 
 IMatrix symplecticNormalization(IMatrix A, int dimension );
 interval tensorNorm( IHessian DDDG , int dimension);
 IHessian compressTensor( IHessian DDDG , int dimension);
  
#endif
