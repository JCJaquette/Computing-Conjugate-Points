 #include "boundaryValueProblem.h"
 
 


// BEGIN        SINGLE SHOOTING CLASS FUNCTIONS


IVector boundaryValueProblem::Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) 
{ 
    
//   Input: 
//     XY_pt    - A vector of length n=dim/2, in local coords on the (un)stable manifold 
//     XY_nbd   - The error of the initial point, in local coords, within the (un)stable manifold. 
//     T        - The amount of time to integrate forward. 
//     STABLE   - Specifies whether the input is on stable or unstable manifold. 
//   Output: 
//     A vector, - a validated enclosure of the *IMAGE* of the time-T map of the input point under the flow.
    
    
  //   We Create our validated solvers, using high order  Taylor method
  ITaylor* solver;
  localManifold *pManifold;
  ITaylor solver_minus((*pf_minus),order);
  ITaylor solver_plus((*pf),order);
  
  if (STABLE)
  {
    solver = &solver_minus;
    pManifold = pStable;
  }
  else
  {
    pManifold = pUnstable; 
    solver = &solver_plus;
  }
  
  ITimeMap Phi((*solver));
    
  
  //    We get the local linearization
  IMatrix A_i = (*(*pManifold).pF).A; 
  
//   Constructs an enclosure of XY_pt/XY_nbd on the manifold with error bounds. 
  vector < IVector >  global_pt_nbd = pManifold -> getPointNbd(  XY_pt,  XY_nbd);
  IVector p_i = global_pt_nbd[0];
  IVector Uxy = global_pt_nbd[1];

  
  C0Rect2Set S_XY(p_i,A_i,Uxy);  
  
    
  // Forward/backwards image of Xinit & Yinit  w/ derivatives 
  IVector XY_new(dimension);
  
  XY_new = Phi(T,S_XY); 
  
  return XY_new ;
}

IMatrix boundaryValueProblem::DGxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) 
{
  
//   Input: 
//     XY_pt    - A vector of length n=dim/2, in local coords on the (un)stable manifold 
//     XY_nbd   - The error of the initial point, in local coords, within the (un)stable manifold. 
//     T        - The amount of time to integrate forward. 
//     STABLE   - Specifies whether the input is on stable or unstable manifold. 
//   Output: 
//     A matrix, - a validated enclosure of the *DERIVATIVE* of the time-T map of the input point under the flow.
    
  //   We Create our solvers 
  ITaylor* solver;
  localManifold *pManifold;
  ITaylor solver_minus((*pf_minus),order);
  ITaylor solver_plus((*pf),order);
  
  if (STABLE)
  {
    solver = &solver_minus;
    pManifold = pStable;
  }
  else
  {
    pManifold = pUnstable; 
    solver = &solver_plus;
  }
    ITimeMap Phi((*solver));
    

      //    We get the local linearization
  IMatrix A_i = (*(*pManifold).pF).A;   
  
  //   Constructs an enclosure of XY_pt/XY_nbd on the manifold with error bounds. 
  vector < IVector >  global_pt_nbd = (*pManifold).getPointNbd(  XY_pt,  XY_nbd);
  IVector p_i = global_pt_nbd[0];
  IVector Uxy = global_pt_nbd[1];
  


    
  C1Rect2Set S_XY(p_i,A_i,Uxy);  
    
  // Forward/backwards image of Xinit/Yinit  w/ derivatives 
  IVector XY_new(dimension);
  IMatrix XY_deriv(dimension,dimension);

  XY_new = Phi(T,S_XY,XY_deriv); 
  
//   Apply chain rule
  IMatrix DG_XY = XY_deriv*A_i*( pManifold -> DW );
  
  
  return DG_XY ;
}


interval boundaryValueProblem::FindTime( IVector XY_pt, interval T)
{
    //     <<>> Single Shooting Function <<>>
    
//     For fixed points x and y on the unstable/stable manifolds, 
//     this function finds T such that the distance || Phi(x,T) - Phi(y,T) || is minimized 
  
//     We use 10 steps of a Newton's method. 
    int repetitions = 10;
//     If Newton's method does not decrease the distance, we return the original T.
    interval T_original = T; 
    
//     We get the input points on the manifolds. 
  vector < IVector > XY_vect_pt  = breakUpXY( XY_pt );
  
  IVector X_pt = XY_vect_pt[0];
  IVector Y_pt = XY_vect_pt[1];
  
  
//   IVector X_mid = XY_vect[2];
  
//   Flowed points
  IVector phi_T_x ;
  IVector phi_T_x_deriv ;
//   First Derivative
  IVector phi_T_y ;
  IVector phi_T_y_deriv ;
//   Second Derivative
  IMatrix phi_T_x_DD;
  IMatrix phi_T_y_DD;
  
  interval distance_D0 =0;
  interval distance_D1 =0;
  interval distance_D2 =0;
  
  IVector zero(dimension/2);
  interval original_dist;
  
  for (int i = 0 ; i< repetitions;i++)
  {
//   We try to minimize   || Phi(x,T)-Phi(y,-T) ||^2 /2
//   let Gx = Phi(x,T)   and Gy = Phi(y,-T)
      
//   This is equivalent to finding a zero of 
//      F(T) =  \sum_{1<=i<=2n} ( Gx _i - Gy _i) ( f( Gx )_i -  f( Gy)_i )
      
//   We use a Newton's method, having defined F'(T) as 
//      F'(T) = \sum_{1<=i<=2n} ( f( Gx )_i -  f( Gy)_i )^2 + ( Gx _i - Gy _i) ( [ DF( Gx)f( Gx )]_i -  [ DF( Gy)f( Gy )]_i )

//       We get the flowed points and their 1st and 2nd derivatives. 
    phi_T_x         = Gxy( X_pt, zero,T, 0);
    phi_T_x_deriv   = (*pf)(phi_T_x);
    phi_T_x_DD      = (*pf)[phi_T_x];
    
    phi_T_y         = Gxy( Y_pt,zero, T, 1);
    phi_T_y_deriv   = -(*pf)(phi_T_y); // Negative, because we flow backwards
    phi_T_y_DD      = -(*pf)[phi_T_y]; // Negative, because we flow backwards
    
//     We take the difference between the the Flowed points 
    
    IVector diff_Val    = phi_T_x   -phi_T_y;
    IVector diff_D1     = phi_T_x_deriv - phi_T_y_deriv;
    IVector diff_D2     = phi_T_x_DD*phi_T_x_deriv - phi_T_y_DD*phi_T_y_deriv;  
    
    distance_D0 =0;
    distance_D1 =0;
    distance_D2 =0;
    
    for (int i = 0; i<dimension;i++)
    {
        distance_D0 += sqr(diff_Val[i]);
        distance_D1 += diff_Val[i]*diff_D1[i];
        distance_D2 += sqr(diff_D1[i])+ diff_Val[i]*diff_D2[i];
    }
    
    if (i==0){
        original_dist = distance_D0;}
    
    T = T - (1/distance_D2.mid())*distance_D1.mid();
    
    T=T.mid();

  }
  
  if (original_dist < distance_D0)
    {
      return T_original;
    }
  else
    {
      return T;
    }
  
}


vector < IVector > boundaryValueProblem::breakUpXY( IVector XY)
{
//   This function takes a vector XY --- XY_pt OR XY_nbd, doesn't matter;
//   	and turns it into two vectors X & Y of half the length
  
  //   We construct X and Y 
  IVector X(dimension/2);
  IVector Y(dimension/2);
  
  for(int i =0;i<dimension/2;i++)
  {
    X[i]=XY[i];
    Y[i]=XY[i+dimension/2];
  }
  
  vector < IVector > XY_out;
  
  XY_out.push_back(X);
  XY_out.push_back(Y);  
  
  return XY_out;
}

IVector  boundaryValueProblem::NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T , interval r_u_sqr ) 
{
//     <<>> Single Shooting Function <<>>
    
  //    We assume we fail
  SUCCESS =0; 
  
  
// //   We start the new things 
  vector < IVector > XY_vect_pt  = breakUpXY( XY_pt);
  vector < IVector > XY_vect_nbd = breakUpXY( XY_nbd);
  
  IVector X_pt = XY_vect_pt[0];
  IVector Y_pt = XY_vect_pt[1];
  
  IVector X_nbd = XY_vect_nbd[0];
  IVector Y_nbd = XY_vect_nbd[1];
  
  
  
//    NOTE THIS SHOULD BE A ZERO NBD
  IVector zero_nbd(dimension/2);
  
  IVector G_x = Gxy( X_pt, zero_nbd, T, 0); 
  IVector G_y = Gxy( Y_pt, zero_nbd, T, 1); 
  
//   cout << " Phi_T (X) = " << G_x << endl;
//   cout << " Phi_T (Y) = " << G_y << endl;
  
  
  IVector G = G_x - G_y; 
  
//   cout << " X = " << X << endl;
//   cout << " Y = " << Y << endl;  
//   cout << " Phi_T (X) - Phi_T (Y) = " << G << endl << endl;
  
//       cout << " X_pt = " << X_pt << endl;
//   cout << " Y_pt = " << Y_pt << endl;
//   cout << " X_nbd = " << X_nbd << endl;
//   cout << " Y_nbd = " << Y_nbd << endl;
  
  IMatrix  DG_x =  DGxy( X_pt, X_nbd, T, 0);
  IMatrix  DG_y = -DGxy( Y_pt, Y_nbd, T, 1);
  
  
  IMatrix DG = DG_combine(DG_x,DG_y,X_pt, X_nbd);// TODO Make pt/nbd version // Maybe?
  
   
//   cout <<  " DG = " << DG << endl;
// // // //     We fix the radius of the point on the unstable manifold 

    
    interval radius = r_u_sqr;
    
  interval x_radius_sqr = 0;
  for(int i = 0 ; i<dimension/2;i++){x_radius_sqr += sqr(X_pt[i]);}
  
//   cout << " x_radius^2 = " << x_radius_sqr << endl;
//   cout << "   radius^2 = " << radius << endl;
// // //    We replace the last term of G with point^2 - radius^2 
  G[dimension-1] = x_radius_sqr - radius; //TODO Remove --- why?


  // //   NOTE Make a Krawczyk version of this. 
  IVector XY_out_nbd = gauss(DG,G); 
  
// //   NOTE Make a Krawczyk version of this. 
// //   
  //     Get approximate inverse
    IMatrix  ApproxInverse = midMatrix(gaussInverseMatrix(midMatrix(DG)));
//     cout << "ApproxInverse = " << ApproxInverse  << endl;
// //     Define identity matrix
    IMatrix eye(dimension,dimension); // TODO Replace by identity matrix function 
    for (int i =0;i<dimension+1;i++){ eye[i][i]=1;}
    IVector new_Nbd     = - ApproxInverse*G + ( eye - ApproxInverse*DG )*XY_nbd;
    IVector kraw_image  = XY_pt + new_Nbd ;

  
  
//   IVector XY_out = XY_pt - XY_out_nbd;
  IVector XY_out = kraw_image;
  XY_out_nbd = new_Nbd;
  
  
//     cout << "output nbd " << XY_out_nbd << endl;
//     cout << "output nbd!" << new_Nbd << endl;
  
//   We check to see if we have a proof of existence/uniqueness
// //  We check that the image is in the interior of the domain (Except in the frozen variable)
  bool verify = 1;
  bool verify_local;
  
    for (int i = 0 ; i< dimension;i++)
    {
        verify_local = subsetInterior(-XY_out_nbd[i],XY_nbd[i]);
        if (verify_local ==0)
            verify=0;
    }
//     cout << " verify = " << verify << endl;
//     verify = subsetInterior(-XY_out_nbd,XY_nbd);
//     cout << " verify = " << verify << "  (New Method) " << endl;
  if (verify ==1)
  {
    SUCCESS = 1;
  }
  
//   cout << "Domain = " << XY_nbd << endl;
//   cout << " Image = " << -XY_out_nbd << endl;

    
  return XY_out;
}


IMatrix boundaryValueProblem::DG_combine( IMatrix DGX, IMatrix DGY, IVector X_pt, IVector X_nbd)
{
    

  IMatrix DG(dimension,dimension);
  for (int i = 0; i< dimension;i++)
  {
    for (int j =0;j<dimension;j++)
    {
      if (j<dimension/2) // 	Add unstable
      {
        DG[i][j] = DGX[i][j];
      }
      else// 	Add stable
      {
        DG[i][j] = DGY[i][j-dimension/2];
      }
    }
  }
  
//   Adds derivative from moving the point on the unstable manifold to the last row
    for ( int i = 0 ; i< dimension; i++)
    {
        if (i < dimension/2)
            DG[dimension-1][i]=2 * (X_pt[i]+X_nbd[i]);
        else
            DG[dimension-1][i] = 0;                
    }
  
  return DG;
}



IVector boundaryValueProblem::ComponentBound( IVector XY,interval T)
{
//     <<>> Auxillary Function <<>>
    
// NOTE This would be improved by rewriting it to use the pt/nbd form
    
//     !!! Vet getTrajectory and getTotalTrajectory first !!! 
    
//   We break up XY into X and Y
  
  vector < IVector > XY_vect = breakUpXY( XY);
  
  IVector X = XY_vect[0];
  IVector Y = XY_vect[1];

//   We get the bound from integrating forward, 
//   then we get the bound from integrating backwards
  
  IVector bound_forward 	= localComponentBound( X,T,0);
  IVector bound_backwards 	= localComponentBound( Y,T,1);
  
 
  
  IVector Bounds_out = intervalHull(bound_forward,bound_backwards);

  
  return Bounds_out;

}

IVector boundaryValueProblem::localComponentBound( IVector XY,interval T, bool STABLE )
{
    //     <<>> Auxillary Function <<>>
    
    // NOTE This would be improved by rewriting it to use the pt/nbd form
    
//     !!! Vet getTrajectory and getTotalTrajectory first !!! 
    //   We Create our solvers 
  ITaylor* solver;
  localManifold *pManifold;
  ITaylor solver_minus((*pf_minus),order);
  ITaylor solver_plus((*pf),order);
  
  if (STABLE)
  {
    solver = &solver_minus;
    pManifold = pStable;
  }
  else
  {
    pManifold = pUnstable; 
    solver = &solver_plus;
  }
    ITimeMap Phi((*solver));
    

  IVector Uxy = (*pManifold).constructU(XY);  
  IMatrix A_i = (*(*pManifold).pF).A; 
  IVector p_i = (*(*pManifold).pF).p;    
  
  C0Rect2Set S_XY(p_i,A_i,Uxy);  
    
  // Forward/backwards image of Xinit & Yinit  w/ derivatives 
  IVector XY_new(dimension);
//   IMatrix XY_deriv(dimension,dimension);

//   XY_new = Phi(T,S_XY,XY_deriv); 
  int grid = 64;
  
  vector<IVector> local_trajectory = getTrajectory(S_XY, T, grid,Phi,(*solver));

//    We compute the bound 
  IVector total_bound = local_trajectory[0]; 
  int length =  local_trajectory.size();
  
  for (int i = 1;i< length ; i++)
  {
    total_bound = intervalHull(total_bound,local_trajectory[i]);
  }
  
  return total_bound ;
}


// END          SINGLE SHOOTING CLASS FUNCTIONS


// // // // // // // // // // // // // // // // // // // // // // // // // // // // 

// // // // // // // // // // // // // // // // // // // // // // // // // // // // 

// // // // // // // // // // // // // // // // // // // // // // // // // // // // 

// // // // // // // // // // // // // // // // // // // // // // // // // // // // 


// BEGIN        MULTIPLE SHOOTING CLASS FUNCTIONS


// Constructor
bvpMultipleShooting::bvpMultipleShooting(IMap &pf_,IMap &pf_minus_, IFunction &p_energy_proj_,localManifold &pStable_,localManifold &pUnstable_, int order_, int num_middle_points_)
:boundaryValueProblem(pf_, pf_minus_,pStable_,pUnstable_, order_)
{
    num_middle_points = num_middle_points_;
    p_energy_proj=&p_energy_proj_;
}


 void bvpMultipleShooting::Integrate_point( IVector coord_pt, IVector coord_nbd,interval T, bool FORWARD,IVector &vector_out, IMatrix &derivative_out)  
 {
     //     <<>> Multiple Shooting Function <<>>
     
//      we assume that coord_pt & coord_nbd are IVectors of size (dimension-1)
//      vector_out and derivative_out are the output, should not have dimension set 
//      NOTE  Not sure what to do about derivative!!
     
  //   We Create our solvers 
  ITaylor* solver;
  ITaylor solver_minus((*pf_minus),order);
  ITaylor solver_plus((*pf),order);
  
  if (FORWARD)
    solver = &solver_plus;
  else
    solver = &solver_minus;
  ITimeMap Phi((*solver));
  
//   We get global (size dimension) vectors for the point and the neighborhood 
  IVector global_pt(dimension);
  IVector global_nbd(dimension);
  for(int i =0;i<dimension-1;i++)
  {
      global_pt[i]  = coord_pt[i];
      global_nbd[i] = coord_nbd[i];
  }
  global_pt[dimension-1]=(*p_energy_proj)(coord_pt);
//   We add the last velocity which gets us on the zero-energy surface
//   NOTE  Turn this into a point, and add the thick interval to the nbd

//   We compute the gradient of the projection in the neighborhood of our point. 
  IVector projection_grad = (*p_energy_proj).gradient(coord_pt+coord_nbd); // Maybe optimize this using second derivative
  
//   cout << " grad = " << projection_grad <<endl;
  
//   NOTE  Make this a thin matrix
//   We compute the local frame for this energy section
  IMatrix A_energy(dimension,dimension);
  for(int i = 0;i<dimension-1;i++)
  {
      A_energy[i][i] =1;
      A_energy[dimension-1][i] = projection_grad[i];
  }
  A_energy[dimension-1][dimension-1]=1;

  
  IMatrix A_energy_error = A_energy - midVector(A_energy);
  A_energy = midVector(A_energy);  
  global_nbd = global_nbd +  gauss(A_energy,A_energy_error*global_nbd) ;
  
  
//   We create the set we will integrate  
  C1Rect2Set S_XY(global_pt,A_energy,global_nbd);  
    
  // Forward/backwards image of our point  w/ derivatives 
  IVector XY_new(dimension);

  IMatrix monodromy_matrix(dimension,dimension);
  vector_out = Phi(T,S_XY,monodromy_matrix);   
  
  
//   cout << "monodromy_matrix = " << monodromy_matrix - midVector(monodromy_matrix)<< endl;
//   We pre-multiply the monodromy_matrix by the derivative of the 0-energy section chart
//     cout << "monodromy_matrix = " << monodromy_matrix - midVector(monodromy_matrix)<< endl;
    IMatrix Test1 = monodromy_matrix * A_energy;
    IMatrix Test2 = monodromy_matrix * (A_energy + A_energy_error);
    monodromy_matrix = monodromy_matrix *A_energy + monodromy_matrix* A_energy_error;
    
//     cout << "monodromy_matrix1 = " << Test1 - midVector(Test1)<< endl;
//     cout << "monodromy_matrix2 = " << Test2 - midVector(Test2)<< endl;
//     cout << "monodromy_matrix0 = " << monodromy_matrix - midVector(monodromy_matrix)<< endl;
//   monodromy_matrix = monodromy_matrix *A_energy + monodromy_matrix* A_energy_error;
//   We get rid of the last column & row
//   * * * X
//   * * * X
//   * * * X
//   X X X X
IMatrix output_matrix(dimension-1,dimension-1);
for( int i =0;i<dimension-1;i++)
{
    for( int j =0;j<dimension-1;j++)
    {
        output_matrix[i][j] = monodromy_matrix[i][j];
    }
}
// cout << " output_matrix = " << output_matrix - midVector(output_matrix)<< endl;
derivative_out = output_matrix;

     return;
 }


vector <IVector> bvpMultipleShooting::NewtonStep(vector <IVector> &points, vector <IVector> &neighborhoods ,interval &integration_time, interval time_nbd) 
{ 
    //     <<>> Multiple Shooting Function <<>>
    
//     We Assume that we verify the newton step, and update this if we fail.
    bool verify = 1;
    
  IVector G = Compute_G(points,integration_time); 
//   cout <<  "___ G = " << G  << endl;  
  
  
//   sum_G >> is just used to see if we want to do this non-rigorously
  interval sum_G =0;
  for(unsigned i =0;i<G.dimension();i++) 
      sum_G += abs(G[i]);
  sum_G = sum_G.right();
  cout << " |G| = " << sum_G << endl;
  
  IMatrix DG = Compute_DG(points, neighborhoods , integration_time+time_nbd) ;
  
//   cout <<  "___ DG = " << DG<< endl;  
  cout <<  "___ det " << det(DG) << endl;  
  
  cout << " size( DG)  = " << DG.dimension( ) << endl;
  cout << " size( G)  = " << G.dimension( ) << endl;
//   cout <<  "___ DGw= " << DG - midVector(DG)<< endl;  
  
//   If |G| is large, we do a non-rigorous newton step
  if (sum_G > pow(10,-4))
  {
    G = midVector(G);
    DG = midVector(DG);
    verify =0;
  }
  IVector XY_out_nbd = gauss(DG,G); 
  

//   NOTE Reimpliment
  IVector initial_vector = Construct_Initial_Vector(points ,neighborhoods,integration_time,0);
  IVector initial_nbd    = Construct_Initial_Vector(points ,neighborhoods,(integration_time+time_nbd),1);

  //   We perform the newton step
  IVector out_vector = initial_vector - XY_out_nbd; 
   
  
vector < IVector > output_regions = Deconstruct_Output_Vector(out_vector);



integration_time = out_vector[(1+num_middle_points)*(dimension-1)+1];

//   BEGIN Verify the newton step

  
//   NOTE  Implement Verification
//   We check to see if we have a proof of existence/uniqueness
// //  We check that the image is in the interior of the domain (Except in the frozen variable)
  
  bool verify_local;
  
    for (int i = 0 ; i< dimension;i++)
    {
        verify_local = subsetInterior(out_vector[i],initial_nbd[i]);
        cout <<  verify_local ;
        if (verify_local ==0)
            verify=0;
    }
    
    cout << endl;
    
  if (verify ==1)
  {
    SUCCESS = 1;
  }
  
//   cout << "Domain = " << XY_nbd << endl;
//   cout << " Image = " << -XY_out_nbd << endl;
//     END
     
  return output_regions;
}



IVector bvpMultipleShooting::Compute_G(const vector <IVector> &points,interval integration_time) 
{
    //     <<>> Multiple Shooting Function <<>>
    
  bool STABLE    = 1;
  bool UNSTABLE  = 0;
  bool FORWARD   = 1;
  bool BACKWARDS = 0;

  
//   int points_length = points.size();
//   int num_middle_points = points_length-2; // The number of points between our stable/unstable coordinates.
  
  //   We create a list of the output from integrating forward / backwards
  vector <IVector> G_forward(num_middle_points+1);
  vector <IVector> G_backwards(num_middle_points+1);
    
// //   We get the stable / unstable coordinates 
  IVector X_pt = points[0];
  IVector Y_pt = points.back();
  
  
  //    NOTE This is for integrating the center point
  IVector zero_nbd(dimension/2);
  IVector zero_nbd_middle(dimension-1);
  
  interval time_zero =0;
  
  //   We integrate from the points on the stable/unstable manifold
  G_forward[0]       = Gxy( X_pt, zero_nbd, integration_time , UNSTABLE); 
  G_backwards.back() = Gxy( Y_pt, zero_nbd, time_zero  , STABLE);
  
  
      // We initialize output for integating the middle points. 
    IVector forward_vector;
    IMatrix forward_derivative;
    IVector backwards_vector;
    IMatrix backwards_derivative;
  
  //   We integrate from the middle points.
  for ( int i =0 ; i< num_middle_points;i++)
  {
    Integrate_point(points[i+1], zero_nbd_middle ,integration_time,FORWARD,forward_vector,forward_derivative);
    G_forward[i+1] = forward_vector;
    G_backwards[i] = points[i+1];
  }
  
  
  
  IVector G  = Construct_G(  G_forward, G_backwards, points);
  cout << " G_forward   = " << G_forward << endl;
  cout << " G_backwards = " << G_backwards << endl;
  
//   We compute the local coordinates of the stable point
  
  
  
  
    return G; 
    
}

IMatrix bvpMultipleShooting::Compute_DG(const vector <IVector> &points, const vector <IVector> &neighborhoods , interval integration_time) 
{
//     <<>> Multiple Shooting Function <<>>
//     
     // We initialize output for integating the middle points. 
    IVector forward_vector;
    IMatrix forward_derivative;
    IVector backwards_vector;
    IMatrix backwards_derivative;

  bool STABLE    = 1;
  bool UNSTABLE  = 0;
  bool FORWARD   = 1;
  bool BACKWARDS = 0;
  
   
//   int points_length = points.size();
//   int num_middle_points = points_length-2; // The number of points between our stable/unstable coordinates.
  
// //   We start the new things 
  IVector X_pt = points[0];
  IVector Y_pt = points.back();
  
  IVector X_nbd = neighborhoods[0];
  IVector Y_nbd = neighborhoods.back();
  
//   NOTE We get the derivative for integrating the CUBE
  //   We create a list of the output from integrating forward / backwards
  vector <IMatrix> DG_forward(num_middle_points+1);
  vector <IMatrix> DG_backwards(num_middle_points+1);
  
  vector <IVector> time_derivatives(num_middle_points+1);
  
//   We integrate from the points on the stable/unstable manifold  
// NOTE These matricies still have an extra bottom row.   
  interval time_zero =0;
    DG_forward[0]        =  DGxy( X_pt, X_nbd, integration_time, UNSTABLE);
    DG_backwards.back()  = -DGxy( Y_pt, Y_nbd, time_zero , STABLE);
    time_derivatives[0]                   =         (*pf)(Gxy( X_pt, X_nbd, integration_time , UNSTABLE)); 
    time_derivatives[num_middle_points+0] +=  -(*pf_minus)(Gxy( Y_pt, Y_nbd, time_zero , STABLE)); 
    
    

//     We construct an identity matrix
    
  IMatrix eye(dimension-1,dimension-1);
  for(int i =0;i<dimension-1;i++){      eye[i][i]=1;    }

  
  //   We integrate from the middle points.
  for ( int i =0 ; i< num_middle_points;i++)
  {
      //       We integrate forward, or backwards 
          Integrate_point(points[i+1], neighborhoods[i+1] ,integration_time,FORWARD,forward_vector,forward_derivative);
          DG_forward[i+1] = forward_derivative;
          time_derivatives[i+1] += (*pf)(forward_vector);
          DG_backwards[i] = -eye;

  }

  
  IMatrix DG = Construct_DG( DG_forward, DG_backwards, points,neighborhoods,time_derivatives);//NOTE  Add G_forward
  
  
  return DG;
  
}



IVector bvpMultipleShooting::Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods, interval integration_time, bool ADD)
{
    //     <<>> Multiple Shooting Function <<>>
    
//     We construct one big vector, adding the points to the neighborhoods, and including the integration time. 
    
    int num_regions = points.size();
    
    if (ADD)
    {
    //     We add the points and the neighborhoods together
        for (int i=0;i<num_regions;i++)
            points[i] = points[i] + neighborhoods[i];
    }
    
//     We erase the middle, middle point NOTE this is not very efficient :( 
//     points.erase( points.begin()+ 1+floor(num_middle_points/2));
    
    IVector initial_vector( (num_middle_points+1)*(dimension-1)+2);
    
    initial_vector[(num_middle_points+1)*(dimension-1)+1]=integration_time;

    int counter = 0;
    
    for (int i=0;i<num_middle_points+2;i++)
    {
        int region_length = points[i].dimension();
        for (int j = 0;j<region_length ;j++)
        {
            initial_vector[counter] = points[i][j];
            counter ++;
        }
    }
    
    cout << " initial_vector <" << initial_vector << endl;

    return initial_vector;   
}

vector <IVector>  bvpMultipleShooting::Deconstruct_Output_Vector(IVector vector_in)
{
    //     <<>> Multiple Shooting Function <<>>
    
  vector <IVector> points_out(num_middle_points+2); //THIS IS THE OUTPUT
  
//     We add the unstable coordinates
  IVector X(dimension/2); 
  for (int i =0;i<dimension/2;i++)
    X[i]=vector_in[i];
  points_out[0]=X;
        

        

  for (int i =0;i< num_middle_points;i++)
  {
      IVector middle_point(dimension-1);

      //     We add the middle points coordinates // NOTE  This is messy
//       if ( i !=   floor(num_middle_points/2))
//       {
        int start_index=0;
//         if (i < floor(num_middle_points/2))
            start_index = dimension/2+ i*(dimension-1);
//         else if (i > floor(num_middle_points/2))
//             start_index = dimension/2+ (i-1)*(dimension-1);
        
        for (int j = 0 ; j < dimension-1;j++)
        {
            middle_point[j] = vector_in[start_index+j];
        }
//       }
      
      points_out[i+1]=middle_point;
  }
  
  
  
  //     We add the stable coordinates
  IVector Y(dimension/2);
  int index_start = (num_middle_points+0)*(dimension-1)+(dimension/2);
//   cout << " vect_in_ size =  "<< vector_in.dimension() << endl;
//   cout << " index_start " << index_start << endl;
  for (int i =0; i< dimension/2;i++)
      Y[i] = vector_in[i+index_start];
  
  points_out[num_middle_points+1]=Y;
  
  cout << " points_out " << points_out << endl;
//   cout << " Y = " << endl;
//   cout << Y << endl;
  
//   cout << " Output !!!!" << endl;
//   print(points_out);
  
    return points_out;   
}

IVector bvpMultipleShooting::Construct_G( vector < IVector > G_forward, vector < IVector > G_backwards, vector <IVector> points) // NOTE  replace with "COMPUTE G" 
{
    //     <<>> Multiple Shooting Function <<>>
    
    int N = G_forward.size();
    IVector G(N*(dimension-1)+2);
    
    for (int i =0;i<N;i++)
    {
        for(int j=0;j<dimension-1;j++)
        {
            G[i*(dimension-1)+j]= G_forward[i][j] - G_backwards[i][j];
        }
    }
    
//     cout << "G = " << G << endl;
    
    //     We getRadiustry to fix the radius of the unstable coordinate
    IVector X_pt = points[0];
    IVector Y_pt = points.back();

//     NOTE This depends on all radii of the (un)stable manifolds being the same width.
//          this should get changed to take the radius for greater relibability. 
    interval radius_x = sqr((*pUnstable).getRadius())*dimension/2;
    interval radius_y = sqr((*pStable).getRadius())*dimension/2;
  
    interval x_radius_sqr = 0;
    interval y_radius_sqr = 0;
    for(int i = 0 ; i<dimension/2;i++){x_radius_sqr += sqr(X_pt[i]);}
    for(int i = 0 ; i<dimension/2;i++){y_radius_sqr += sqr(Y_pt[i]);}
// // //    We replace the last term of G with point^2 - radius^2 
  G[N*(dimension-1)] = x_radius_sqr - radius_x; 
  G[N*(dimension-1)+1] = y_radius_sqr - radius_y; 
            
  
    return G;   
}

IMatrix bvpMultipleShooting::Construct_DG( const vector <IMatrix> &DG_forward, const  vector <IMatrix> &DG_backwards, const  vector <IVector> &points, const vector <IVector> &neighborhoods, const vector < IVector> &time_derivatives )//NOTE  Revise
{
    //     <<>> Multiple Shooting Function <<>>
    
    
    int N =  (num_middle_points+1)*(dimension-1)+2;
    
    
    
    IMatrix DG(N,N);
    
//     Add the DG from the unstable manifold
    for ( int i =0;i<dimension/2;i++)  // column
    {
        for (int j=0;j<dimension-1;j++) // row
        {
            DG[j][i] = DG_forward[0][j][i];
        }
    }
    
    //     Add the DG from the stable manifold
    for ( int i =0;i<dimension/2;i++)  // column
    {
        for (int j=0;j<dimension-1;j++) // row
        {
            DG[N-(dimension-1)-1+j-1][N-(dimension/2)+i-1] = DG_backwards.back()[j][i];
        }
    }

    
//  Add the middle points

    for ( int i =0;i <num_middle_points;i++)
    {
        int point_index = i*(dimension-1);
        for (int j =0;j<dimension-1;j++) // row
        {
            for (int k =0;k<dimension-1;k++) // column
            {
                DG[j+point_index][k+point_index+dimension/2]                =DG_backwards[i][j][k];   // on top 
                DG[j+point_index+dimension-1][k+point_index+dimension/2]    =DG_forward[i+1][j][k];   // on bottom 
            }
        }
    }
    
//     We add the derivative from fixing the unstable radius -- to penultimate row
    for (int i =0;i<dimension/2; i++)
    {
        DG[N-2][i]=2*(points[0][i]+neighborhoods[0][i]);
    }

    //     We add the derivative from fixing the stable radius -- to ultimate row
    for (int i =0;i<dimension/2; i++)
    {
        DG[N-1][N-dimension/2-1+i]=2*(points.back()[i]+neighborhoods.back()[i]);
    }
    


    //     We add the derivative with respect to time to the last column
    for (int i = 0 ; i < num_middle_points +1;i++)
    {
//         IVector local_derivative = (*pf)(G_forward[i]);
        for ( int j = 0 ; j < dimension -1; j++)
        {
            DG[i*(dimension-1)+j][N-1] = time_derivatives[i][j];
        }
    }

  return DG;   
}


void bvpMultipleShooting::setMiddlePoints(int shots){num_middle_points = shots;};

// END        MULTIPLE SHOOTING CLASS FUNCTIONS





// // // // // // // // // // // // // // // // //   

//  Below is the code from heteroclinic.cpp used to run the multiple shooting feature. 

// // // //   BEGIN OLD CODE  :::  Testing integration of middle points
// // //   
// // // //     We define the XY_pt that will get used in the single-shooting newton's method. 
// // //   IVector XY_pt(dimension);  
// // //   
// // //   boundaryValueProblem BVP(f,f_minus,localStable,localUnstable,order); // Single Shooting   
// // //   bvpMultipleShooting BVP_ms(f,f_minus,energy_projection,localStable,localUnstable,order ,shots ); // Multiple Shooting 
// // //   
// // //   
// // //   
// // //   
// // //     
// // //     IVector XY_pt_T(dimension);  
// // //     for (int i = 0 ; i < dimension / 2 ; i++)
// // //     {
// // //         XY_pt_T[i] = mid(points[0][i]);
// // //         XY_pt_T[i+dimension/2] = mid( points.back()[i]);
// // //     }
// // //     interval T_new = BVP.FindTime(  XY_pt_T,  T);
// // //     
// // //     
// // //     cout << " T old = " << T << endl; 
// // //     cout << " T new = " << T_new << endl; 
// // //     
// // //     
// // //     T = T_new;
// // //     integration_time = 2*T/(shots+1);
// // //     //   Guess = Guess_pt_nbd( dimension, All_parameters, T, localUnstable,  localStable, shots);
// // //     //   points         = Guess[0];
// // //     //   neighborhoods  = Guess[1];
// // //     
// // //     
// // //     if (USE_MULTIPLE_SHOOTING==1){
// // //     
// // //     
// // //     //    We do the newton method
// // //     vector < IVector > regions;
// // //     interval time_nbd =0;
// // //     for (int i = 0 ; i<multiple_newton_steps ;i++)
// // //     {
// // //         integration_time = integration_time.mid();
// // //         cout << " integration_time " << integration_time <<endl;
// // //         regions = BVP_ms.NewtonStep(points, neighborhoods ,integration_time, time_nbd ) ;
// // //         for (unsigned j=0;j<regions.size();j++)
// // //             points[j] = midVector(regions[j]);
// // //     }
// // //     
// // //     interval multiplier = 0.;
// // //     time_nbd = 0 * (integration_time - integration_time.mid());
// // //     
// // //     cout << " Time nbd " << time_nbd <<endl;
// // //     integration_time = integration_time.mid();
// // //     for (unsigned j=0;j<regions.size();j++){
// // //         neighborhoods[j] = multiplier *(regions[j] -  points[j]);
// // //     }
// // // 
// // //     for (unsigned i = 0 ; i < regions.size();i++)
// // //             {
// // //                 cout << " neighborhoods["<<i<<"] = " << neighborhoods[i] << endl;
// // //             }
// // //     cout << " Verifying .... " << endl;      
// // //     regions = BVP_ms.NewtonStep(points, neighborhoods ,integration_time, time_nbd ) ;
// // //     cout << " Time nbd " << integration_time - integration_time.mid() <<endl;
// // //     
// // //         if (multiple_newton_steps>0)
// // //         {
// // //             cout << endl << "Final Guess " << endl;
// // //             
// // //             for (unsigned i = 0 ; i < regions.size();i++)
// // //             {
// // //     //             cout << " Region["<<i<<"] = " << regions[i] << endl;
// // //             }
// // //             
// // //             cout << endl;
// // //             
// // //                 for (unsigned i = 0 ; i < regions.size();i++)
// // //             {
// // //                 cout << " Region["<<i<<"] width = " << regions[i] - midVector(regions[i]) << endl;
// // //             }
// // //             cout << "done testing " << endl;
// // //         }
// // //     //     return -5;
// // //         
// // //     cout << " Old T = " << T << endl;
// // //     T = (shots+1)*integration_time/2;
// // //     cout << " New T = " << T << endl;
// // //   
// // //   
// // //     
// // //     
// // //     for (int i = 0 ; i < dimension / 2 ; i++)
// // //     {
// // //         XY_pt[i] = mid(regions[0][i]);
// // //         XY_pt[i+dimension/2] = mid( regions.back()[i]);
// // //     }
// // //     
// // //     cout << " XY_pt = " << XY_pt << endl;
// // // 
// // //   
// // //     
// // //   }// Multiple Shooting 
// // // //   END Testing integration of middle points
  
  
