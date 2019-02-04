 #include "boundaryValueProblem.h"


 
 void boundaryValueProblem::Integrate_point( IVector coord_pt, IVector coord_nbd,interval T, bool FORWARD,IVector &vector_out, IMatrix &derivative_out)  
 {
//      we assume that coord_pt & coord_nbd are IVectors of size (dimension-1)
//      vector_out and derivative_out are the output, should not have dimension set 
//      TODO Not sure what to do about derivative!!
     
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
//   We add the last velocity which gets us on the zero-energy surface
  global_pt[dimension-1]=(*p_energy_proj)(coord_pt);
//   TODO Turn this into a point, and add the thick interval to the nbd

//   We compute the gradient of the projection in the neighborhood of our point. 
  IVector projection_grad = (*p_energy_proj).gradient(coord_pt+coord_nbd); // Maybe optimize this using second derivative
  
//   cout << " grad = " << projection_grad <<endl;
  
//   TODO Make this a thin matrix
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
    cout << "monodromy_matrix = " << monodromy_matrix - midVector(monodromy_matrix)<< endl;
    IMatrix Test1 = monodromy_matrix * A_energy;
    IMatrix Test2 = monodromy_matrix * (A_energy + A_energy_error);
    monodromy_matrix = monodromy_matrix *A_energy + monodromy_matrix* A_energy_error;
    
    cout << "monodromy_matrix1 = " << Test1 - midVector(Test1)<< endl;
    cout << "monodromy_matrix2 = " << Test2 - midVector(Test2)<< endl;
    cout << "monodromy_matrix0 = " << monodromy_matrix - midVector(monodromy_matrix)<< endl;
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

IVector boundaryValueProblem::Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) 
{ 
  
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
  
  
  vector < IVector >  global_pt_nbd = (*pManifold).getPointNbd(  XY_pt,  XY_nbd);
  IVector p_i = global_pt_nbd[0];
  IVector Uxy = global_pt_nbd[1];

  
  C0Rect2Set S_XY(p_i,A_i,Uxy);  
  
    
  // Forward/backwards image of Xinit & Yinit  w/ derivatives 
  IVector XY_new(dimension);
  
  XY_new = Phi(T,S_XY); 
  
  return XY_new ;
}


IMatrix boundaryValueProblem::DGxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) // TODO Point separated copy
{
  
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
  
  vector < IVector >  global_pt_nbd = (*pManifold).getPointNbd(  XY_pt,  XY_nbd);
  IVector p_i = global_pt_nbd[0];
  IVector Uxy = global_pt_nbd[1];
  


    
  C1Rect2Set S_XY(p_i,A_i,Uxy);  
    
  // Forward/backwards image of Xinit & Yinit  w/ derivatives 
  IVector XY_new(dimension);
  IMatrix XY_deriv(dimension,dimension);

  XY_new = Phi(T,S_XY,XY_deriv); 
  
  IMatrix DG_XY = XY_deriv*A_i*((*pManifold).DW);
  
  
  return DG_XY ;
}


interval boundaryValueProblem::FindTime( IVector XY_pt, interval T)
{
  int repetitions = 10;
    
  vector < IVector > XY_vect_pt  = breakUpXY_gen( XY_pt );
  
  
  IVector X_pt = XY_vect_pt[0];
  IVector Y_pt = XY_vect_pt[1];
  
  
//   IVector X_mid = XY_vect[2];
  
  IVector phi_T_x ;
  IVector phi_T_x_deriv ;
  IVector phi_T_y ;
  IVector phi_T_y_deriv ;
  
  IVector zero(dimension/2);
  
  for (int i = 0 ; i< repetitions;i++)
  {
//     cout << " T = " << T << endl;
    phi_T_x = Gxy( X_pt, zero,T, 0);
    phi_T_x_deriv = (*pf)(phi_T_x);
    
    phi_T_y = Gxy( Y_pt,zero, T, 1);
    phi_T_y_deriv = (*pf)(phi_T_y);
    
//     cout << "Phi_T(x) = " << phi_T_x << endl;
//     cout << "Phi_T(y) = " << phi_T_y << endl;
    
    interval sum = phi_T_x[0] + phi_T_x[1];
    interval D_sum = phi_T_x_deriv[0] + phi_T_x_deriv[1];
    
    interval sum1=0;// =  phi_T_x[0] + phi_T_x[2];
    interval sum2 =0;//  phi_T_y[0] + phi_T_y[2];
    interval D_sum1 =0;// phi_T_x_deriv[0] + phi_T_x_deriv[2];
    interval D_sum2 =0;// phi_T_y_deriv[0] + phi_T_y_deriv[2];
    
    for (int i = 0; i<dimension;i++)
    {
      sum1 	= sum1 	 + phi_T_x[i];
      sum2 	= sum2 	 + phi_T_y[i];
      D_sum1 	= D_sum1 + phi_T_x_deriv[i];
      D_sum2 	= D_sum2 + phi_T_y_deriv[i];
    }
      
    
    interval diff   = sum1 - sum2;
    interval D_diff = D_sum1 + D_sum2;
    
//     cout << " Diff = " << diff << endl;
//     cout << " D_diff = " << D_diff << endl;
    
//     T = T - (1/D_sum)*(sum-1);
    T = T - (1/D_diff)*(diff);
    T=T.mid();


  }
//       cout << "Phi_T (x)      = " << phi_T_x << endl;
  cout << " T   = " << T << endl;
  return T;
}


vector < IVector > boundaryValueProblem::breakUpXY( IVector XY)
{
//   This function takes a vector XY and turns it into vectors X & Y of half the length
//   It also makes a copy of the midpoint objects
  
  //   We construct X and Y 
  IVector XY_mid(dimension);
  XY_mid = midVector(XY);
  

  IVector X(dimension/2);
  IVector Y(dimension/2);
  IVector X_mid(dimension/2);
  IVector Y_mid(dimension/2);
  
  for(int i =0;i<dimension/2;i++)
  {
    X[i]=XY[i];
    X_mid[i]=XY_mid[i];
    Y[i]=XY[i+dimension/2];
    Y_mid[i]=XY_mid[i+dimension/2];
  }

  vector < IVector > XY_out;
  
  XY_out.push_back(X);
  XY_out.push_back(Y);
  XY_out.push_back(X_mid);
  XY_out.push_back(Y_mid);
  XY_out.push_back(XY_mid);
  
  
  return XY_out;
}


vector < IVector > boundaryValueProblem::breakUpXY_gen( IVector XY)
{
//   This function takes a vector XY_pt OR XY_nbd
//   	and turns it into vectors X & Y of half the length

  
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



IVector  boundaryValueProblem::NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T) 
{

  //    We assume we fail
  SUCCESS =0; 
  
  
// //   We start the new things 
  vector < IVector > XY_vect_pt  = breakUpXY_gen( XY_pt);
  vector < IVector > XY_vect_nbd = breakUpXY_gen( XY_nbd);
  
  IVector X_pt = XY_vect_pt[0];
  IVector Y_pt = XY_vect_pt[1];
  
  IVector X_nbd = XY_vect_nbd[0];
  IVector Y_nbd = XY_vect_nbd[1];
  
  
  
//    TODO THIS SHOULD BE A ZERO NBD
  IVector zero_nbd(dimension/2);
  
  IVector G_x = Gxy( X_pt, zero_nbd, T, 0); 
  IVector G_y = Gxy( Y_pt, zero_nbd, T, 1); 
  
//   cout << " Phi_T (X) = " << G_x << endl;
//   cout << " Phi_T (Y) = " << G_y << endl;
  
  
  IVector G = G_x - G_y; 
  
//   cout << " X = " << X << endl;
//   cout << " Y = " << Y << endl;  
//   cout << " Phi_T (X) - Phi_T (Y) = " << G << endl << endl;
  

  
  IMatrix  DG_x =  DGxy( X_pt, X_nbd, T, 0);
  IMatrix  DG_y = -DGxy( Y_pt, Y_nbd, T, 1);
  
  
  IMatrix DG = DG_combine(DG_x,DG_y,X_pt, X_nbd);// TODO Make pt/nbd version // Maybe?
  
   
//   cout <<  " DG = " << DG << endl;
// // //     We fix the radius of the point on the unstable manifold 
    interval radius = sqr((*pUnstable).getRadius())*dimension/2;
  interval x_radius_sqr = 0;
  for(int i = 0 ; i<dimension/2;i++){x_radius_sqr += sqr(X_pt[i]);}
  
//   cout << " x_radius^2 = " << x_radius_sqr << endl;
//   cout << "   radius^2 = " << radius << endl;
// // //    We replace the last term of G with point^2 - radius^2 
  G[dimension-1] = x_radius_sqr - radius; //TODO Remove --- why?


  IVector XY_out_nbd = gauss(DG,G); 
  
  
  
  IVector XY_out = XY_pt - XY_out_nbd;
  
  
//   cout << "output nbd " << XY_out_nbd << endl;
  
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
  if (verify ==1)
  {
    SUCCESS = 1;
  }
  
//   cout << "Domain = " << XY_nbd << endl;
//   cout << " Image = " << -XY_out_nbd << endl;

    
  return XY_out;
}



vector <IVector> boundaryValueProblem::NewtonStep(vector <IVector> &points, vector <IVector> &neighborhoods ,interval &integration_time, interval time_nbd) 
{ 
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
  
  cout <<  "___ DG = " << DG<< endl;  
  cout <<  "___ det " << det(DG) << endl;  
  
//   cout <<  "___ DGw= " << DG - midVector(DG)<< endl;  
  
//   If |G| is large, we do a non-rigorous newton step
  if (sum_G > pow(10,-4))
  {
    G = midVector(G);
    DG = midVector(DG);
    verify =0;
  }
  IVector XY_out_nbd = gauss(DG,G); 
  
  
//   TODO Reimpliment
  IVector initial_vector = Construct_Initial_Vector(points ,neighborhoods,integration_time,0);
  IVector initial_nbd    = Construct_Initial_Vector(points ,neighborhoods,(integration_time+time_nbd),1);

  //   We perform the newton step
  IVector out_vector = initial_vector - XY_out_nbd; 
   
  
vector < IVector > output_regions = Deconstruct_Output_Vector(out_vector);


integration_time = out_vector[num_middle_points*(dimension-1)+1];

//   BEGIN Verify the newton step

  
//   TODO Implement Verification
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



IVector boundaryValueProblem::Compute_G(const vector <IVector> &points,interval integration_time) 
{
  bool STABLE    = 1;
  bool UNSTABLE  = 0;
  bool FORWARD   = 1;
  bool BACKWARDS = 0;

  
//   int points_length = points.size();
//   int num_middle_points = points_length-2; // The number of points between our stable/unstable coordinates.
  
  //   We create a list of the output from integrating forward / backwards
  vector <IVector> G_forward(num_middle_points);
  vector <IVector> G_backwards(num_middle_points);
    
// //   We get the stable / unstable coordinates 
  IVector X_pt = points[0];
  IVector Y_pt = points.back();
  
  
  //    NOTE This is for integrating the center point
  IVector zero_nbd(dimension/2);
  IVector zero_nbd_middle(dimension-1);
  
  //   We integrate from the points on the stable/unstable manifold
  interval time_zero =0;
  G_forward[0]       = Gxy( X_pt, zero_nbd, integration_time , UNSTABLE); 
  G_backwards.back() = Gxy( Y_pt, zero_nbd, integration_time , STABLE);
  
  
      // We initialize output for integating the middle points. 
    IVector forward_vector;
    IMatrix forward_derivative;
    IVector backwards_vector;
    IMatrix backwards_derivative;
  
  //   We integrate from the middle points.
  for ( int i =0 ; i< num_middle_points;i++)
  {
//       We integrate forward, or backwards 
      if ( i < floor(num_middle_points/2))
      {
          Integrate_point(points[i+1], zero_nbd_middle ,integration_time,FORWARD,forward_vector,forward_derivative);
          G_forward[i+1] = forward_vector;
          G_backwards[i] = points[i+1];
      }
      else if ( i > floor(num_middle_points/2) )
      {
          Integrate_point(points[i+1], zero_nbd_middle ,integration_time,BACKWARDS,backwards_vector,backwards_derivative);
          G_forward[i] = points[i+1];
          G_backwards[i-1] = backwards_vector;
      }
  }
  
  IVector G  = Construct_G(  G_forward, G_backwards, points);
  
  
//   We compute the local coordinates of the stable point
  
//   cout << " G_forward   = " << G_forward << endl;
//   cout << " G_backwards = " << G_backwards << endl;
  
  
    return G; 
    
}

IMatrix boundaryValueProblem::Compute_DG(const vector <IVector> &points, const vector <IVector> &neighborhoods , interval integration_time) 
{
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
  vector <IMatrix> DG_forward(num_middle_points);
  vector <IMatrix> DG_backwards(num_middle_points);
  
  vector <IVector> time_derivatives(num_middle_points);
  
//   We integrate from the points on the stable/unstable manifold  
// NOTE These matricies still have an extra bottom row.   
  interval time_zero =0;
    DG_forward[0]        =  DGxy( X_pt, X_nbd, integration_time, UNSTABLE);
    DG_backwards.back()  = -DGxy( Y_pt, Y_nbd, integration_time , STABLE);
    time_derivatives[0]                   =         (*pf)(Gxy( X_pt, X_nbd, integration_time , UNSTABLE)); 
    time_derivatives[num_middle_points-1] +=  -(*pf_minus)(Gxy( Y_pt, Y_nbd, integration_time , STABLE)); 
    
    

//     We construct an identity matrix
    
  IMatrix eye(dimension-1,dimension-1);
  for(int i =0;i<dimension-1;i++){      eye[i][i]=1;    }

  
  //   We integrate from the middle points.
  for ( int i =0 ; i< num_middle_points;i++)
  {
      //       We integrate forward, or backwards 
      if ( i < floor(num_middle_points/2))
      {
          Integrate_point(points[i+1], neighborhoods[i+1] ,integration_time,FORWARD,forward_vector,forward_derivative);
          DG_forward[i+1] = forward_derivative;
          time_derivatives[i+1] += (*pf)(forward_vector);
          DG_backwards[i] = -eye;
      }
      else if ( i > floor(num_middle_points/2) )
      {
          Integrate_point(points[i+1], neighborhoods[i+1] ,integration_time,BACKWARDS,backwards_vector,backwards_derivative);
          DG_forward[i] = eye;
          DG_backwards[i-1] = -backwards_derivative;
          time_derivatives[i-1] += -(*pf_minus)(backwards_vector);
      }

  }

  
  IMatrix DG = Construct_DG( DG_forward, DG_backwards, points,neighborhoods,time_derivatives);//TODO Add G_forward
  
  
  return DG;
  
}

IVector boundaryValueProblem::Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods) // TODO I think we can delete this function
{
    
    
    int num_regions = points.size();
    int num_middle_points = num_regions-2; // The number of points between our stable/unstable coordinates.
    
//     We add the poitn and the neighborhood together
    for (int i=0;i<num_regions;i++)
        points[i] = points[i] + neighborhoods[i];

    IVector initial_vector( num_middle_points*(dimension-1)+dimension);
    
    int counter = 0;
    
    for (int i=0;i<num_regions;i++)
    {
        int region_length = points[i].dimension();
        for (int j = 0;j<region_length ;j++)
        {
            initial_vector[counter] = points[i][j];
            counter ++;
        }
    }

    return initial_vector;   
}

IVector boundaryValueProblem::Construct_Initial_Vector(vector <IVector> points,const vector <IVector> &neighborhoods, interval integration_time, bool ADD)
{
//     We construct one big vector, adding the points to the neighborhoods, and including the integration time. 
    
    int num_regions = points.size();
    
    if (ADD)
    {
    //     We add the points and the neighborhoods together
        for (int i=0;i<num_regions;i++)
            points[i] = points[i] + neighborhoods[i];
    }
    
//     We erase the middle, middle point TODO this is not very efficient :( 
    points.erase( points.begin()+ 1+floor(num_middle_points/2));
    
    IVector initial_vector( num_middle_points*(dimension-1)+2);
    
    initial_vector[num_middle_points*(dimension-1)+1]=integration_time;

    int counter = 0;
    
    for (int i=0;i<num_middle_points+1;i++)
    {
        int region_length = points[i].dimension();
        for (int j = 0;j<region_length ;j++)
        {
            initial_vector[counter] = points[i][j];
            counter ++;
        }
    }
    

    return initial_vector;   
}
IVector boundaryValueProblem::Construct_Initial_Vector(vector <IVector> points, interval integration_time)
{
//     We construct one big vector,  including the integration time. 
    
//     We erase the middle, middle point TODO this is not very efficient :( 
    points.erase( points.begin()+ 1+floor(num_middle_points/2));
    
    IVector initial_vector( num_middle_points*(dimension-1)+2);
    
    initial_vector[num_middle_points*(dimension-1)+1]=integration_time;

    int counter = 0;
    
    for (int i=0;i<num_middle_points+1;i++)
    {
        int region_length = points[i].dimension();
        for (int j = 0;j<region_length ;j++)
        {
            initial_vector[counter] = points[i][j];
            counter ++;
        }
    }
    

    return initial_vector;   
}
vector <IVector>  boundaryValueProblem::Deconstruct_Output_Vector(IVector vector_in)
{
    
  vector <IVector> points_out(num_middle_points+2); //THIS IS THE OUTPUT
  
//     We add the unstable coordinates
  IVector X(dimension/2); 
  for (int i =0;i<dimension/2;i++)
    X[i]=vector_in[i];
  points_out[0]=X;
        

        

  for (int i =0;i< num_middle_points;i++)
  {
      IVector middle_point(dimension-1);

      //     We add the middle points coordinates // TODO This is messy
      if ( i !=   floor(num_middle_points/2))
      {
        int start_index=0;
        if (i < floor(num_middle_points/2))
            start_index = dimension/2+ i*(dimension-1);
        else if (i > floor(num_middle_points/2))
            start_index = dimension/2+ (i-1)*(dimension-1);
        
        for (int j = 0 ; j < dimension-1;j++)
        {
            middle_point[j] = vector_in[start_index+j];
        }
      }
      
      points_out[i+1]=middle_point;
  }
  
  //     We add the stable coordinates
  IVector Y(dimension/2);
  int index_start = (num_middle_points-1)*(dimension-1)+(dimension/2);
//   cout << " vect_in_ size =  "<< vector_in.dimension() << endl;
//   cout << " index_start " << index_start << endl;
  for (int i =0; i< dimension/2;i++)
      Y[i] = vector_in[i+index_start];
  
  points_out[num_middle_points+1]=Y;
  
  
//   cout << " Y = " << endl;
//   cout << Y << endl;
  
//   cout << " Output !!!!" << endl;
//   print(points_out);
  
    return points_out;   
}

IVector boundaryValueProblem::Construct_G( vector < IVector > G_forward, vector < IVector > G_backwards, vector <IVector> points) // TODO replace with "COMPUTE G" 
{
    int N = G_forward.size();
    IVector G(N*(dimension-1)+2);
    
    for (int i =0;i<N;i++)
    {
        for(int j=0;j<dimension-1;j++)
        {
            G[i*(dimension-1)+j]= G_forward[i][j] - G_backwards[i][j];
        }
    }
    
    //     We try to fix the radius of the unstable coordinate
    IVector X_pt = points[0];
    IVector Y_pt = points.back();

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

IMatrix boundaryValueProblem::Construct_DG( const vector <IMatrix> &DG_forward, const  vector <IMatrix> &DG_backwards, const  vector <IVector> &points, const vector <IVector> &neighborhoods, const vector < IVector> &time_derivatives )//TODO Revise
{
    
    
    int N =  num_middle_points*(dimension-1)+2;
    
    
    
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

    for ( int i =0;i <num_middle_points-1;i++)
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
    for (int i = 0 ; i < num_middle_points ;i++)
    {
//         IVector local_derivative = (*pf)(G_forward[i]);
        for ( int j = 0 ; j < dimension -1; j++)
        {
            DG[i*(dimension-1)+j][N-1] = time_derivatives[i][j];
        }
    }

  return DG;   
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


IVector boundaryValueProblem::NormBound( IVector XY,interval T)
{

//   We break up XY into X and Y
  
  
  vector < IVector > XY_vect = breakUpXY( XY);
  
  IVector X = XY_vect[0];
  IVector Y = XY_vect[1];

//   We get the bound from integrating forward, 
//   then we get the bound from integrating backwards
  
  IVector bound_forward 	= localNormBound( X,T,0);
  IVector bound_backwards 	= localNormBound( Y,T,1);
  
 
  
  IVector Bounds_out = intervalHull(bound_forward,bound_backwards);

  
  return Bounds_out;

}

IVector boundaryValueProblem::localNormBound( IVector XY,interval T, bool STABLE )
{

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
  
//   vector<IVector> getTrajectory(C0Rect2Set &s,interval T,int grid,ITimeMap &timeMap,IOdeSolver &solver)
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

