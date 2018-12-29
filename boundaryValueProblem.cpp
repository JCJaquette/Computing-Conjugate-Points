 #include "boundaryValueProblem.h"


IVector boundaryValueProblem::Gxy( IVector XY_pt, IVector XY_nbd,interval T, bool STABLE) // TODO  point separated copy
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
  int repetitions = 5;
    
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



IVector  boundaryValueProblem::NewtonStep( IVector XY_pt, IVector XY_nbd  ,interval T) // TODO Make pt/nbd version 
{
//   NOTE: WE FIX X[frozen] AS A SINGLE NUMBER 

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
  
  
  IMatrix DG = DG_combine(G,DG_x,DG_y,X_pt, X_nbd);// TODO Make pt/nbd version // Maybe?
  
   
//   cout <<  " DG = " << DG << endl;
// // //     We fix the radius of the point on the unstable manifold 
    interval radius = sqr((*pUnstable).getRadius())*dimension/2;
  interval x_radius_sqr = 0;
  for(int i = 0 ; i<dimension/2;i++){x_radius_sqr += sqr(X_pt[i]+X_nbd[i]);}
  
//   cout << " x_radius^2 = " << x_radius_sqr << endl;
//   cout << "   radius^2 = " << radius << endl;
// // //    We replace the last term of G with point^2 - radius^2 
  G[dimension-1] = x_radius_sqr - radius; //TODO Remove
  
  
//   //TODO Delete down 
//   IVector G_trunc(dimension-1);
//   for (int i =0 ; i<dimension-1;i++)
//   {
//     G_trunc[i]=G[i];
//   }
//   //TODO Delete up

//   IVector invDG_G = gauss(DG,G_trunc); // TODO  invDG_G -->> XY_out_nbd ;  G_trunc -->> G
  IVector XY_out_nbd = gauss(DG,G); // TODO  invDG_G -->> XY_out_nbd ;  G_trunc -->> G
  
  
  
//   //TODO Delete down 
// //   We put our output in the right format
//   IVector XY_out_nbd(dimension);
//   for (int i = 0; i < dimension;i++)
//   {
//     if (i<frozen)
//       XY_out_nbd[i] = invDG_G[i];
//     else if (i>frozen)
//       XY_out_nbd[i] = invDG_G[i-1];
//   }
//   //TODO Delete up
  
  
  
  IVector XY_out = XY_pt - XY_out_nbd;
  
//    We don't change the coordinate we've frozen
//   XY_out[frozen]=XY_pt[frozen]; //NOTE REMOVED FROZEN
  
  
//   cout << "output nbd " << XY_out_nbd << endl;
  
//   We check to see if we have a proof of existence/uniqueness
// //  We check that the image is in the interior of the domain (Except in the frozen variable)
  //TODO REMOVE FROZEN
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
  
  cout << "Domain = " << XY_nbd << endl;
  cout << " Image = " << -XY_out_nbd << endl;
    
     

    
  return XY_out;
}




IMatrix boundaryValueProblem::DG_combine(IVector G, IMatrix DGX, IMatrix DGY, IVector X_pt, IVector X_nbd)
{
  
//    We create the 'frozen' version of DG i=[0,dim-1] j = [1,dim]
  

//   | *  *  *
//   | *  *  *
//   | *  *  *
//   X -- -- --

//   or maybe 
  
//   * | *  *
//   * | *  *
//   * | *  *
//   --X -- --
    
    
    //TODO REMOVE FROZEN
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
  
/*
    
    //TODO REMOVE FROZEN
  IMatrix DG(dimension-1,dimension-1);
  for (int i = 0; i< dimension-1;i++)
  {
    for (int j =0;j<dimension;j++)
    {
      if (j<dimension/2) // 	Add stable
      {
	if (j<frozen)
	  DG[i][j] = DGX[i][j];
	else if ( j> frozen)
	  DG[i][j-1] = DGX[i][j];
      }
      else
      {
	if (j < frozen)
	  DG[i][j] = DGY[i][j-dimension/2];
	else if ( j > frozen)
	  DG[i][j-1] = DGY[i][j-dimension/2];
      }
    }
  }*/
  
 
  
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

