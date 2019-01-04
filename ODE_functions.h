#ifndef ODE_functions_h
#define ODE_functions_h

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"

DVector fixedPoint(int i,int dimension) ;
DMatrix coordinateChange(int i,int dimension, vector < double> All_parameters) ;
IVector knownSolution( double a, interval x);
vector < IVector >  multipleShootingGuess( int shots, interval T ,int dimension, vector < double > All_parameters);
IVector initialGuessGlobal(int dimension, vector <double> All_parameters, interval T, bool STABLE);
vector <IFunction> constructEnergy(int dimension,  vector < double > All_parameters);
vector < IMap > constructFunctions( int dimension, vector < double > All_parameters);


#endif
