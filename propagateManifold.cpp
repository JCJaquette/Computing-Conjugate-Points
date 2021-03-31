#include "propagateManifold.h"
// #include <../../home/jonathan/capd-5.0.59/capdDynSys4/include/capd/dynset/dynsetLib.h>
// #include <../../home/jonathan/capd-5.0.59/capdAlg/include/capd/vectalg/vectalgLib.h>


vector <IVector> propagateManifold::construct_InitCondU(int eigenvector_NUM)
{
//     Outputs the center approximation for the eigenvector, and the error it is away from the eigenfunction.
//     Error is in local coordinates.
    
    
// We get the point on the unstable manifold and put it into global coord. 
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
  
//   cout << " Uxy = " << Uxy << endl; // The error on the heteroclinic orbit.
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

//     cout << " Eigenfunction_Error = " << Eigenfunction_Error << endl;

  
    for (int i =0;i< dimension;i++)
  {
    lin_init_nbd[i+dimension]  = Eigenfunction_Error[i];
  }

    cout << "E-fun error = " << lin_init_nbd  << endl;
    
  vector < IVector > output;
  output.push_back(lin_init_pt);
  output.push_back(lin_init_nbd);
  return output;
}

IMatrix propagateManifold::construct_A_lin(void)
{
//  Used in specifying the initial condition, and the directions in which its error goes.
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
// INPUT
//  T       --  Total time to integrate
//  grid    --  equal subdivision (1/grid) used to get tighter bounds when computing '3 function values' and '4 function derivative'.
//  eigenvector_NUM --  which eigenfunction to integrate


  interval timeStep = interval(pow(2,-step_size)); 
    //   We Create our solvers 
  ITaylor lin_solver((*pf),order);
  
  lin_solver.setStep(timeStep);
  ITimeMap lin_Phi(lin_solver);

  
//   Constructs two vectors with initial condition for the coupled heteroclinic orbit / eigenfunction system
  vector <IVector> lin_init 	= construct_InitCondU(eigenvector_NUM); 
 
  IVector lin_init_pt  = lin_init [0];  // the center of the vector
  IVector lin_init_nbd = lin_init [1];  // the error of the approximation
  
//   cout << "Initial Eigenvector pt  = " << lin_init_pt << endl;
//   cout << "Initial Eigenvector nbd = " << lin_init_nbd << endl;
  
//   Creates the matrix along which our initial error is set
  IMatrix A_lin 	= construct_A_lin(); 
  
  C0Rect2Set S_X_init_0(lin_init_pt,A_lin,lin_init_nbd);  
  
  vector<IMatrix> local_trajectory = getTotalTrajectory(S_X_init_0, T, grid,lin_Phi,lin_solver);
    
//   IVector out_set =S_X_init_0;
//   cout << endl << " Final Set " << out_set << endl;
    
  return local_trajectory;
  
};



int propagateManifold::frameDet( interval L_minus, int grid,IVector endPoint_LPlus)
{
//  Input
//     L_minus          -   The time amount of time we need to integrate the frame matrix forward.  
//     grid             -   used in 'getTotalTrajectory'
//     endPoint_LPlus   -   The point phi(L_+) - phi(+\infty)  on the heteroclinic orbit. 
//  Output
//     non-negative number -- number of validated conjugate poitns
//     -3 ---   failed to count conjugate poitns
//     -4 ---   failed L_+ condition
//     -5 ---   failed no conjugate points below L_- condition

  

//   We compute the eigenfunction error;
    cout << endl<<"Computing the eigenfunction error at minus infinity ..." << endl;
    (*pUnstable).computeEigenError_minus_infty();
    
    
//  Check Prop 2.5 / (2.12) that there are no conjugate points below -L_-. 
  bool ConjugatePointsBelowLminus = (*pUnstable).checkConjugatePointsBelowLminus();
  if ( ConjugatePointsBelowLminus == 0){
    cout << "We could not verify the no conjugate points below L_- condition " << endl;
    return -5;
  }   

//   We flow each eigen-vector (with error) forward. 
  vector < vector< IMatrix> > List_of_Trajectories(dimension/2);
  for (int i = 0 ; i<dimension/2;i++)
  {
    List_of_Trajectories[i] = computeTotalTrajectory(i,  L_minus,  grid);
  }
  
  
  cout << endl<< "Checking Stability ... " << endl << endl;
  
  topFrame A_frame(List_of_Trajectories);
 
  A_frame.initialize();
  
  if (MAKE_PLOT){
    A_frame.makePlot();
  }

//   We check the L_+ condition
  bool L_PLUS = lastEuFrame( A_frame,endPoint_LPlus);
  
  if ( L_PLUS == 0){
//       We could not verify the L_+ condition
      return -4;
  }      
  
  //   We count all of the conjugate points
  vector<int> conjugate_points = A_frame.countZeros();
  
  if ( conjugate_points[1]>0 )  
  {
//       We could not verify the conjugate conjugate_points
      return -3;
  }
  else 
      return conjugate_points[0];
 
}



bool propagateManifold::lastEuFrame(topFrame &A_frame , IVector endPoint_LPlus)
{
//     INPUT
//         endPoint_LPlus -- the point ( \phi(L_+) - p_s  ) 
//     OUTPUT
//         true/false     -- whether the L_+ condition is satisfied.   
    


    
//      We validate the stable manifold in a larger nbd    
    IMatrix last_Frame = A_frame.getLastFrame();
    //   cout <<endl<< "Last frame = " << last_Frame << endl; 
        
    localManifold_Eig localStableBig = construct_Manifold_at_LPlus(endPoint_LPlus);
   
    bool conditions_S = localStableBig.checkConditions(  ); 
    if( !(conditions_S) )
    {
        cout << "Failed to validate manifold at L_+  " << endl;
        return 0;
    }  
    
        
    cout << endl<<"Computing the eigenfunction error at plus  infinity ..." << endl;
    localStableBig.computeEigenError_plus_infty();
    IMatrix EFunction_Error = localStableBig.Eigenfunction_Error_plus_infty ;    
    
//     cout <<" Eigenfunction_Error_plus_infty  " << EFunction_Error  << endl;
    
    // We compute <<epsilon_0>>     
    // Each component of 'EFunction_Error' is already scaled by the norm of the e-vectors. Thus, we take the union of all these error bounds.
    interval eps_0;
    for (int i = 0 ; i < dimension;i++){
        eps_0 = intervalHull(eps_0,EFunction_Error[0][i]);
    }
//     cout << "eps0 = "<<  eps_0 << endl;
    
    // Validated eigenvalues
    IVector eigenvalues = localStableBig.eigenvalues;
    //     cout << " eigenvalues= " << eigenvalues<< endl;
    
    IMatrix U_coord = projectionGammaBeta(  last_Frame , EFunction_Error );
    
    cout << endl<<"Checking L_+ conditions ..." << endl ;
    bool L_PLUS = checkL_plus(A_frame, U_coord,eps_0,eigenvalues);

    return L_PLUS;
}



localManifold_Eig propagateManifold::construct_Manifold_at_LPlus( IVector endPoint_LPlus)
{
//     Input
//         endPoint_LPlus  -- the point ( \phi(L_+) - p_s  )     , in global coordinates 
//     Output
//          localStableBig -- a stable manifold containing    \phi(L_+)  .
    
    
    IMatrix A_s = (*(*pStable).pF).A;
    
    
//     We validate the stable manifold in a larger nbd   
    cout << "--Dist from stable equilibrium         = " << endPoint_LPlus << endl;
//     We put ourselfs into local coordinates
    endPoint_LPlus = gauss(A_s,endPoint_LPlus);
    cout << "--Dist from stable equilibrium (eigen) = " << endPoint_LPlus << endl;
    
    IVector vec_Lplus_s(dimension/2);
    IVector vec_Lplus_u(dimension/2);
    for (int i = 0 ; i< dimension/2; i++){
        vec_Lplus_u[i] = endPoint_LPlus[i];
        vec_Lplus_s[i] = endPoint_LPlus[i+dimension/2];
    }
    

    cout << " vec_Lplus_u  = " << vec_Lplus_u  << endl;
    cout << " vec_Lplus_s  = " << vec_Lplus_s  << endl;

//     We set the nbd on which we validate the manifold so that it 
//      1) includes the point vec_Lplus_s
//      2) includes the origin, 
//      3) has +/- 1e-10 wiggle room
    IVector U_flat_new(dimension/2);
    for (int i =0;i<dimension/2;i++){ 
        U_flat_new[i]=vec_Lplus_s[i]*interval(-1e-10,1+1e-10) ;
    }
    
    
//     Recall, L here is \vartheta in the paper. 
    interval L_new = euclNorm(vec_Lplus_u)/euclNorm(vec_Lplus_s);
    L_new = L_new.right();
    //   cout << " L_new  = " << L_new  << endl;

    
    //   We create the local manifold object, which encloses our final trajectory;
    bool STABLE = 1;
    localManifold_Eig localStableBig((*(*pStable).pF),U_flat_new, L_new, STABLE,manifold_subdivision);            
    
//     We update/increase L until the manifold is verified.
    int adjust_L = 10;
    bool conditions_S;
    for (int i = 0 ; i< adjust_L;i++){
        conditions_S = localStableBig.checkConditions(  ); 
        if (conditions_S )
            break;
        else
            localStableBig.L = localStableBig.L*2;
        
//         cout << " localStableBig.L = " << localStableBig.L << endl;
    }
    
    return localStableBig;
}



IMatrix propagateManifold::projectionGammaBeta(  IMatrix &last_Frame ,const IMatrix & EFunction_Error )
{
//   INPUT
//      EFunction_Error -- Error of eigenfunctions at L_+, in global coordinates 
//      last_Frame      -- the frame matrix U(x) at L_+, in global coordinates
//   OUTPUT 
//      U_coord         -- An enclosure of the stacked matrices Gamma and Beta.
    
    
    IMatrix A_s = (*(*pStable).pF).A;
    IMatrix eye = identityMat(dimension);
      
    
// Old output for the projection without error.     
// // BEGIN We project the final endpoint into the large stable manifold coordinates 
//     for (int i = 0 ; i < dimension /2 +1;i++)
//     {
//         IVector col = getColumn(last_Frame,dimension,i);
//         IVector vec_in_local_coord = gauss(A_s,col);
//         vec_in_local_coord = vec_in_local_coord/getMax(abs(vec_in_local_coord));
//         
//         if (i < dimension /2)
//             cout << " w_" << i << "   = ";
//         else if (i == dimension /2)
//             cout << " phi'  = ";
//         cout << vec_in_local_coord << endl;
//     }
// //     END
    
    cout << "Final Frame, in Gamma Beta coordinates" << endl;
    
    
    IMatrix U_coord(dimension,dimension/2);

    
    // BEGIN We project the final endpoint into the large stable manifold coordinates 
    for (int j = 0 ; j < dimension /2 +1;j++)
    {
        IVector col = getColumn(last_Frame,dimension,j);
        IVector vec_in_local_coord = gauss(A_s,col);
        
//      NOTE The MAJOR source of error here is inverting "eye+krawczykInverse(A_s)*EFunction_Error " 
//         This wouldn't be such a problem if 'EFunction_Error' wasn't so big. Maybe in the future this could be improved. 
        vec_in_local_coord = gauss(eye+krawczykInverse(A_s)*EFunction_Error ,vec_in_local_coord);  // EigenfunctionError
        vec_in_local_coord = vec_in_local_coord/getMax(abs(vec_in_local_coord));
        
//      store and output the vector.
        if (j < dimension /2){
            cout << "  U_" << j << "   = ";
            for (int i = 0 ; i  < dimension ; i++){
                U_coord[i][j]    = vec_in_local_coord[i];
            }
        }
        else if (j == dimension /2){
            cout << "  phi'  = ";
        }
        cout << vec_in_local_coord << endl;
    }    

//     END 


    return U_coord;
}


bool propagateManifold::checkL_plus(topFrame &A_frame, IMatrix U_coord,  interval eps_0,IVector eigenvalues  ){
//     Different L_+ estimates may be obtained by removing one of the unstable eigenfunctions, and doing that computation.
//     This function creates the Gamma & Beta matrices for each possibility, and then checks whether the L_+ condition is satisfied. 
//     
//   INPUT
//      U_coord     -- The stacked Gamma / Beta matrix, thick interval matrices
//      eps_0       -- eps_0  from the paper 
//      eigenvalues -- eigenvalues 
//   OUTPUT
//         true/false     -- whether the L_+ condition is satisfied. 

//     cout << "U_coord = " << U_coord << endl;
        
    vector < IMatrix > Gamma_List;
    vector < IMatrix > Beta_List;
    
    IMatrix Gamma_local(dimension/2,dimension/2-1);
    IMatrix Beta_local(dimension/2,dimension/2-1);

    for (int k =0;k<dimension/2;k++)                //remove k, the hat-index
    {
        int k_adjust =0;
        for (int j = 0 ; j<dimension/2;j++)         // COLUMNS
        {
            if(j==k){k_adjust =1;continue;}
            for (int i = 0 ; i < dimension/2;i++){  // ROWS
                Gamma_local[i][j-k_adjust]  = U_coord[i][j];
                Beta_local[i ][j-k_adjust]  = U_coord[i+dimension/2][j];
            }
        }
        Gamma_List.push_back(Gamma_local);
        Beta_List.push_back(Beta_local);
    }
    
    interval nu_1 = eigenvalues[0];                 // Largest 
    interval nu_n = eigenvalues[dimension/2-1];     // Smallest
    
//     cout << " nu_1 = " << nu_1 << endl;
//     cout << " nu_n = " << nu_n << endl;
    
    bool L_PLUS = 0;
    for (int k = 0;k<dimension/2;k++){
        cout << endl;
        //  Check to see if we have a full rank frame matrix for E^u_- if we use (*) \varphi' and (*) all the unstable eigenfunctions except U_k . 
        if ( false == A_frame.checkFirstFrame(k) ){
            cout << " Frame matrix defined with U_{ hat{k} } NOT of full rank   (k=" << k <<")" << endl; 
            continue;
        }


        //  We check the L_+ condition, breaking if successful. 
        L_PLUS = checkL_plus_local( Gamma_List[k], Beta_List[k],eps_0, nu_1, nu_n );
        if (L_PLUS ==1)
            break; 
    }
    
    return L_PLUS;
}

bool propagateManifold::checkL_plus_local( IMatrix Gamma, IMatrix Beta,interval eps_0, interval nu_1 , interval nu_n ){
//   INPUT
//      Gamma       -- A Gamma interval matrix, for a particular  hat{i} // hat{k} 
//      Beta        -- A Beta  interval matrix, for a particular  hat{i} // hat{k} 
//      eps_0       -- eps_0  from the paper 
//      nu_1        -- Largest  eigenvalue
//      nu_n        -- Smallest eigenvalue
//   OUTPUT
//      true/false  -- for the choice of Gamma, Beta, whether the L_+ condition is satisfied. cf. proposition 2.9
    
    interval epsilon_beta = compute_epsilon_beta( Gamma, Beta );
    
//     If mu^* could not be bounded above 0, the L_+ condition FAILS
    if (epsilon_beta < 0)
        return 0;
    
    epsilon_beta = epsilon_beta .right();
    eps_0  = eps_0 .right();
    
    cout << " eps_0 = " << eps_0 << endl; 
    cout << " epsilon_beta  = " << epsilon_beta  << endl;
    
//     NOTE The example we've considered is with the diffusion matrix D equal to the identity. 
    interval d_min =1;
    interval d_max =1;
    interval n = dimension/2;   
    
    interval C_P    = ( 2* sqrt( 2*n* nu_1 * d_max)) / ( 1 - eps_0 *  sqrt( 2*n*nu_1*d_max )) ;
    interval C_Q    = sqrt(2*n/(nu_n*d_min))  +  sqrt( 2*n*nu_1 * d_max)  +  2* eps_0*n  ;
    
    C_P = C_P.right();
    C_Q = C_Q.right();
    
    interval C_M1   = eps_0 * C_Q * ( 2 + eps_0 * C_Q); 
    interval C_M2   = eps_0 * C_P * ( 2 + eps_0 * C_P) +  2* epsilon_beta*(1+ eps_0 *C_P)  + sqr(epsilon_beta); 
    
    interval sum_for_neumann_test =  eps_0 * sqrt( 2 * n * nu_1 * d_max);
//     cout << " sum_for_neumann_test = " << sum_for_neumann_test << endl; 
    
    bool NEUMANN_SERIES_TEST;
    if (( C_M1 <1 ) && ( C_M2 < 1) && (sum_for_neumann_test <1))
        NEUMANN_SERIES_TEST =1;
    else{
        cout << " Unable to apply Neumann series " << endl << endl;
        NEUMANN_SERIES_TEST =0;
    }
    
    
    cout << " C_P   = " << C_P << endl; 
    cout << " C_Q   = " << C_Q << endl; 
    cout << " C_M1  = " << C_M1 << endl; 
    cout << " C_M2  = " << C_M2 << endl << endl; 
    
    interval sum_1 = 2*eps_0*(C_Q + C_P);
    interval sum_2 = sqr(eps_0)* ( sqr(C_Q) + sqr(C_P));
    interval sum_3 = sqr(epsilon_beta) + 2*epsilon_beta*(1+eps_0*C_P);
    interval sum_4 = sqr(1+eps_0*C_Q)*C_M1/abs(1-C_M1);
    interval sum_5 = sqr(1+eps_0*C_P+epsilon_beta)*C_M2/abs(1-C_M2);
    
    interval total_sum = sum_1 + sum_2 + sum_3 + sum_4 + sum_5;
    
    cout << " sum_1 = " << sum_1 << endl; 
    cout << " sum_2 = " << sum_2 << endl; 
    cout << " sum_3 = " << sum_3 << endl; 
    cout << " sum_4 = " << sum_4 << endl; 
    cout << " sum_5 = " << sum_5 << endl; 
    
    cout << "total_sum = " << total_sum << endl; 
    
    bool L_PLUS = (NEUMANN_SERIES_TEST && ( total_sum <1));
    
    return L_PLUS ;
}

interval propagateManifold::compute_epsilon_beta( IMatrix Gamma, IMatrix Beta ){
//   Uses methods in Section 3.4 to bound the constant \epsilon_{b}
// 
//   INPUT
//      Gamma       -- A Gamma interval matrix, for a particular  hat{i} // hat{k} 
//      Beta        -- A Beta  interval matrix, for a particular  hat{i} // hat{k} 
//   OUTPUT
//      true/false  -- for the choice of Gamma, Beta, whether the L_+ condition is satisfied, 
    
    capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm; // Matrix norm
    
// STEP 1: We compute a bound on \mu^* 
    
    IMatrix Gamma_center = midMatrix(Gamma);
    IMatrix Gamma_delta  = Gamma - Gamma_center;

//     cout << " Gamma_center = " << Gamma_center << endl;
//     cout << " Gamma_delta = " << Gamma_delta << endl;
    
    IMatrix GcT_by_Gc = transpose(Gamma_center)*Gamma_center    ;
    IMatrix GdT_by_Gc = transpose(Gamma_delta)*Gamma_center     ;
    IMatrix GdT_by_Gd = transpose(Gamma_delta)*Gamma_delta      ;
    
//     cout << " Gc^T * Gc     = " << GcT_by_Gc <<endl;
//     cout << " Gd^T * Gc     = " << GdT_by_Gc <<endl;
//     cout << " Gd^T * Gd     = " << GdT_by_Gd <<endl;
        
//      Recall logarithmic minimum is   ml(A) :=  \min \{ \lambda \in \sigma( (A^T+A)/2 \} 
//      Note also GcT_by_Gc is symmetric     
    interval mu_c = ml( GcT_by_Gc );
//     cout << " mu_center  = " << mu_c << endl;
    
    interval mu = mu_c - interval(2)*euclNorm(GdT_by_Gc) - euclNorm(GdT_by_Gd); 
    cout << " mu^* = " << mu << endl;

    
//  This will cause  the parent function  "checkL_plus_local" to return FAILURE
    if (mu.left() < 0){
        cout << " Failure: Lower bound on mu^* is negative and equals " << mu << endl << endl;
        return interval(-1);
    }
    
// STEP 2: We compute a bound on ||Beta|| and define \epsilon_b 
    cout << " ||B|| = " << euclNorm(Beta) << endl;

    interval epsilon_beta = euclNorm(Beta) / sqrt(mu);
    
    
    return epsilon_beta ;    
}

