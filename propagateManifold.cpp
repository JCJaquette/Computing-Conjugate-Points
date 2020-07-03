#include "propagateManifold.h"
// #include <../../home/jonathan/capd-5.0.59/capdDynSys4/include/capd/dynset/dynsetLib.h>
// #include <../../home/jonathan/capd-5.0.59/capdAlg/include/capd/vectalg/vectalgLib.h>


vector <IVector> propagateManifold::construct_InitCondU(int eigenvector_NUM)
{
// We get the point on the unstable manifold. Point is in global coord. 
  vector < IVector >  global_pt_nbd = (*pUnstable).getPointNbd(  XY_pt,  XY_nbd);
  IVector p_i = global_pt_nbd[0];
  IVector Uxy = global_pt_nbd[1];
  
  
//   We construct the point for linearized system
  IVector lin_init_pt(dimension*2);
  
  
//   We add the point into the larger coords.
  for(int i = 0 ; i< dimension;i++)
  {
    lin_init_pt[i] = p_i[i];
  }
  
//   We choose the particular eigenvector 
  IVector EigenVector_component = (*pUnstable).getEigenvector( eigenvector_NUM);
  IMatrix A_i = (*(*pUnstable).pF).A;
  
//   cout << "Eigenvector = " << EigenVector_component << endl;
//   cout << "A           = " << A_i[1] << endl;
  
  for (int i = 0 ; i< dimension;i++)
  {
    lin_init_pt[i+dimension]=EigenVector_component[i]; 
  }
  
  
//   We construct the nbd for linearized system  
  IVector lin_init_nbd(dimension*2);
//   We add the nbd into the larger coords.
  for(int i = 0 ; i< dimension;i++)
  {
    lin_init_nbd[i] = Uxy[i];
  }



    IVector Eigenfunction_Error = (*pUnstable).getEigenError_minus_infty(eigenvector_NUM);

// cout << " Eigenfunction_Error = " << Eigenfunction_Error << endl;

  
    for (int i =0;i< dimension;i++)
  {
    lin_init_nbd[i+dimension]  = Eigenfunction_Error[i];
  }

    cout << " Error = " << lin_init_nbd  << endl;
    
  vector < IVector > output;
  output.push_back(lin_init_pt);
  output.push_back(lin_init_nbd);
  return output;
}

IMatrix propagateManifold::construct_A_lin(void)
{
  IMatrix A_u = (*(*pUnstable).pF).A; 
  
  IMatrix A_u_lin(dimension*2,dimension*2);
  for(int i = 0; i<dimension;i++)
  {
    for(int j = 0; j<dimension;j++)
    {
      A_u_lin[i][j] = A_u[i][j];
      A_u_lin[i+dimension][j+dimension] = A_u[i][j];
    }
  }
  return A_u_lin;
}


vector <IMatrix> propagateManifold::computeTotalTrajectory(int eigenvector_NUM, interval T, int grid)
{

  
    int thread_id =omp_get_thread_num();
    if (thread_id ==1 ){
        cout << "  Using multiple processors " << endl;
        abort();
    }
  
  interval timeStep = interval(pow(2,-step_size)); 
    //   We Create our solvers 
  ITaylor lin_solver(list_of_maps[thread_id ],order);
  
  lin_solver.setStep(timeStep);
  ITimeMap lin_Phi(lin_solver);

  
  
  vector <IVector> lin_init 	= construct_InitCondU(eigenvector_NUM); 
  IVector U_lin_init =lin_init [0];
  
  
  
  IVector lin_init_pt  = lin_init [0];
  IVector lin_init_nbd = lin_init [1];
  
//   cout << "Initial Eigenvector pt  = " << lin_init_pt << endl;
//   cout << "Initial Eigenvector nbd = " << lin_init_nbd << endl;
  
  
  
  IMatrix A_lin 	= construct_A_lin();
  
  
  
  C0Rect2Set S_X_init_0(lin_init_pt,A_lin,lin_init_nbd);  

  
    
  // Forward image of Xinit w/ derivatives 
  IMatrix X_lin_deriv(dimension,dimension);

  
  
  vector<IMatrix> local_trajectory = getTotalTrajectory(S_X_init_0, T, grid,lin_Phi,lin_solver);
    
  IVector out_set =S_X_init_0;
//   cout << endl << " Final Set " << out_set << endl;
  
    
  return local_trajectory;
  
};



int propagateManifold::frameDet(interval T,int grid)
{
  

//   We compute the eigenfunction error;
    (*pUnstable).computeEigenError_minus_infty();
  
  vector < vector< IMatrix> > List_of_Trajectories(dimension/2);
  
  
  
  int max_threads = omp_get_max_threads();
  if (max_threads > dimension/2)
      max_threads = dimension/2;
  
  list_of_maps.resize(max_threads );
  for( int i = 0 ; i < max_threads  ; i ++ ) { list_of_maps[i]=(*pf);}
  
// // //   /// I am trying to parrelize this
// // // //   #pragma omp parallel for  
  for (int i = 0 ; i<dimension/2;i++)
  {
    List_of_Trajectories[i] = computeTotalTrajectory(i,  T,  grid);
  }

  
  
  
  cout << "Checking Stability ... " << endl;
  
  topFrame A_frame(List_of_Trajectories);
//   pA_frame = &A_frame;

//   (*pA_frame) = topFrame(List_of_Trajectories); 
 
  A_frame.initialize();
  
//   A_frame.makePlot();
  
  vector<int> conjugate_points = A_frame.countZeros();

  
  
  lastEuFrame( A_frame);
  
  if ( conjugate_points[1]>1 )
  {
//       We could not verify the conjugate conjugate_points
      return -3;
  }
  else 
      return conjugate_points[0];
  
  
 
}

void propagateManifold::lastEuFrame(topFrame &A_frame)
{
  
  IMatrix A_s = (*(*pStable).pF).A;
  
  
  IMatrix last_Frame = A_frame.getLastFrame();
//   cout <<endl<< "Last frame = " << last_Frame << endl; 
  
  for (int i = 0 ; i < dimension /2 +1;i++)
  {
    IVector col = getColumn(last_Frame,dimension,i);
    IVector vec_in_local_coord = gauss(A_s,col);
    vec_in_local_coord = vec_in_local_coord/getMax(abs(vec_in_local_coord));
    
    if (i < dimension /2)
      cout << " w_" << i << "   = ";
    else if (i == dimension /2)
      cout << " phi'  = ";
    cout << vec_in_local_coord << endl;
  }
  
//    We get the last column to output the final point
  IVector col = getColumn(last_Frame,dimension,dimension /2 +1);
  IVector stable_point = (*(*pStable).pF).p ;
  for (int i = 0;i<dimension;i++){ stable_point[i]=stable_point[i]-col[i];}
  
  cout << "Dist from stable equilibrium = " << stable_point << endl;
  
  
}


