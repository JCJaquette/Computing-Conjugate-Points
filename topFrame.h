#ifndef topFrame_h
#define topFrame_h


#include <iostream>
using namespace std;

#include "capd/capdlib.h"
#include "utils.h"


using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

class topFrame 
{
private:
  int num_trajectories; // de facto variable used for defining the dimension;    num_trajectories = n (symplectic);     2*num_trajectories = dimension (as used elsewhere). 
  int series_length;
  
  vector <vector < IMatrix> > * p_traject_list;
// 1st level is  # of trajectories 
// 2nd level is  time series  

//   The Matrix contains in its columns:
  //   0 -- Time Interval +++++ ONLY IN ENTRY (0,0)
  //   1 -- Left Endpoint
  //   2 -- Right Endpoint
  //   3 -- Bound on function Values
  //   4 -- Bound on function Derivative  
  
  vector < interval>   time_series;
//   we pull the time component out of traject_list
  
  
  vector < vector< IMatrix> >   frame_series;
// 1st level is  time series
// 2nd level is  type:
  //   0 -- Left Endpoint
  //   1 -- Right Endpoint
  //   2 -- Bound on function Values
  //   3 -- Bound on function Derivative
  
//   The matrix then represented is the top frame
  
  


  vector < vector <interval> >   det_series;
// 1st level is  time series
// 2nd level is  type
  
//   The interval is:
//   0 -- Left Endpoint
//   1 -- Right Endpoint
//   2 -- Bound on function Values
//   3 -- Bound on function Derivative ------ REMOVED!!

  
    
  void constructTimeSeries( void);
  void constructDetSeries( void);
  
  void constructFrameSeries( void);
  
  
  interval calculateDerivative( const IMatrix &A , const IMatrix &A_prime);
  void adjugate( const IMatrix &A , IMatrix &matrix_out);
  interval minorDet( const IMatrix &A ,int i_hat , int j_hat );
  
  IMatrix getFirstFrame( int j_hat);

public:
  topFrame(vector <vector < IMatrix> >  &p_traject_list_){p_traject_list  = &p_traject_list_;num_trajectories =(*p_traject_list).size();series_length =(*p_traject_list)[0].size();}
  void initialize( void);
  
  
  void makePlot( void);
  vector<int> countZeros( void);
  
  IMatrix getLastFrame( void);
  
  bool checkFirstFrame( int j_hat);
  
  
 
};


#endif
