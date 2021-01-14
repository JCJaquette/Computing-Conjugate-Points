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
  
  cout << " Uxy = " << Uxy << endl;
  
  cout << "Eigenvector = " << EigenVector_component << endl;
  cout << "A           = " << A_i[1] << endl;
  
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

  
//     int thread_id =omp_get_thread_num();
//     if (thread_id ==1 ){
//         cout << "  Using multiple processors " << endl;
//         abort();
//     }
  
  interval timeStep = interval(pow(2,-step_size)); 
    //   We Create our solvers 
//   ITaylor lin_solver(list_of_maps[thread_id ],order);
  ITaylor lin_solver((*pf),order);
  
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



int propagateManifold::frameDet(interval T, interval L_plus, int grid,IVector endPoint_LPlus)
{
  

//   We compute the eigenfunction error;
    (*pUnstable).computeEigenError_minus_infty();
  
  vector < vector< IMatrix> > List_of_Trajectories(dimension/2);
  

  for (int i = 0 ; i<dimension/2;i++)
  {
    List_of_Trajectories[i] = computeTotalTrajectory(i,  T + L_plus,  grid);
  }
  
  
  cout << "Checking Stability ... " << endl;
  
  topFrame A_frame(List_of_Trajectories);
//   pA_frame = &A_frame;

//   (*pA_frame) = topFrame(List_of_Trajectories); 
 
  A_frame.initialize();
  
//   A_frame.makePlot();
  
  
  bool L_PLUS = lastEuFrame( A_frame,endPoint_LPlus);

  vector<int> conjugate_points = A_frame.countZeros();
  
  if ( L_PLUS == 0){
//       We could not verify the L_+ condition
      return -4;
  }      
  
  
  

  
  if ( conjugate_points[1]>1 )
  {
//       We could not verify the conjugate conjugate_points
      return -3;
  }
  else 
      return conjugate_points[0];
  
  
 
}


bool propagateManifold::lastEuFrame(topFrame &A_frame , IVector endPoint_LPlus)
{
    int STABLE = 1;
    int manifold_subdivision = 15;
  
    IMatrix A_s = (*(*pStable).pF).A;
//     BEGIN We validate the stable manifold in a larger nbd
    
    IMatrix last_Frame = A_frame.getLastFrame();
    //   cout <<endl<< "Last frame = " << last_Frame << endl; 
    
    cout << "--Dist from stable equilibrium         = " << endPoint_LPlus << endl;
    endPoint_LPlus = gauss(A_s,endPoint_LPlus);
    cout << "--Dist from stable equilibrium (eigen) = " << endPoint_LPlus << endl;
    
    IVector vec_Lplus_s(dimension/2);
    IVector vec_Lplus_u(dimension/2);
    for (int i = 0 ; i< dimension/2; i++){
        vec_Lplus_u[i] = endPoint_LPlus[i];
        vec_Lplus_s[i] = endPoint_LPlus[i+dimension/2];
    }
    
//     cout << " endPoint_LPlus  = " << endPoint_LPlus  << endl;
//     cout << " vec_Lplus_u  = " << vec_Lplus_u  << endl;
//     cout << " vec_Lplus_s  = " << vec_Lplus_s  << endl;
    
  interval L_angle_new = euclNorm(vec_Lplus_u)/euclNorm(vec_Lplus_s);
  
//   cout << " L_angle_new  = " << L_angle_new  << endl;
  
  
    
    //    We get the last column to output the final point
//     IVector col = getColumn(last_Frame,dimension,dimension /2 +1);
//     IVector stable_point = (*(*pStable).pF).p ;
//     for (int i = 0;i<dimension;i++){ stable_point[i]= col[i] - stable_point[i];}
    

    
    IVector U_flat_new(dimension/2);
    for (int i =0;i<dimension/2;i++){ U_flat_new[i]=vec_Lplus_s[i]*interval(-1e-10,1+1e-10);}
    
    
    localVField F_s_new     = (*(*pStable).pF);
    interval    L_new       = L_angle_new.right();//5*euclNorm(U_flat_new).right();
    
    
      
    localManifold localStableBig((*(*pStable).pF),U_flat_new, L_new, STABLE,manifold_subdivision);            //   We create the local manifold object, which encloses our final trajectory;
    
    int adjust_L = 10;
    bool conditions_S;
    for (int i = 0 ; i< adjust_L;i++){
        conditions_S = localStableBig.checkConditions( U_flat_new ); 
        if (conditions_S )
            break;
        else
            localStableBig.L = localStableBig.L*2;
        
//         cout << " localStableBig.L = " << localStableBig.L << endl;
    }
    
    if(! (conditions_S))
    {
        cout << "Failed to validate L_+ Conditions " << endl;
        return 0;
    }  
//     TODO There is a problem with the normalization of the eigenvectors .!! 
    
    localStableBig.computeEigenError_plus_infty();
    IMatrix EFunction_Error = localStableBig.Eigenfunction_Error_plus_infty ;
    
    cout <<" Eigenfunction_Error_plus_infty  " << EFunction_Error  << endl;
    
    
    
    
    IMatrix eye = identityMat(dimension);

    
//     cout.precision(16);
//     IMatrix Kronic(dimension,dimension);
//     for (int i = 0 ; i< dimension;i++){
//         IVector V_left = getColumn(A_s,dimension,i);
//         V_left = V_left/euclNorm(V_left);
//         for (int j = 0 ; j<dimension ; j++){
//             IVector V_right = getColumn(A_s,dimension,j);
//             V_right = V_right/euclNorm(V_right);
//             Kronic[i][j] = omega(V_left,V_right,dimension/2);
//         }
//     }
//     cout << " omega product  = " << Kronic << endl;
            
    
    
//     END
    
//     TODO Make sure that this projection onto coordinates is up to code!!
    
// BEGIN We project the final endpoint into the large stable manifold coordinates 
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
    
    cout << " With additional error" << endl;
    
    IMatrix U_coord(dimension,dimension/2);
    
//     IMatrix rightMat = EFunction_Error;
//     for (int i = 0 ; i < dimension;i++){
//         for( int j =0; j<dimension;j++){
//             rightMat[i][j] = abs(rightMat[i][j]).right();
//         }
//     }
//     IMatrix max_inverse_bound = krawczykInverse(eye-rightMat) - eye;
//     max_inverse_bound = interval(-1,1)* max_inverse_bound;
    
//     IMatrix max_inverse_bound = krawczykInverse(eye + EFunction_Error) - eye;
    
//     cout << rightMat << endl;
    
//     IMatrix eye_plus_Eu_inv = boundEyeInverseDefect(EFunction_Error,dimension);

    IMatrix U_coord_pt(dimension,dimension/2);
    IMatrix U_coord_nbd(dimension,dimension/2);
    
    IVector vec_in_local_coord_pt;
    IVector vec_in_local_coord_nbd;


    
    // BEGIN We project the final endpoint into the large stable manifold coordinates 
    for (int j = 0 ; j < dimension /2 +1;j++)
    {
        IVector col = getColumn(last_Frame,dimension,j);
        IVector vec_in_local_coord = gauss(A_s,col);
        
        vec_in_local_coord_pt = vec_in_local_coord ; // NEW 
        
//         vec_in_local_coord = vec_in_local_coord/getMax(abs(vec_in_local_coord));
        
//         cout << " New Bound = " << vec_in_local_coord + max_inverse_bound*vec_in_local_coord  << endl;        
//         cout << " Old Bound = " << gauss(eye+EFunction_Error ,vec_in_local_coord) << endl;  // This is better than krawczykInverse 

                
        vec_in_local_coord = gauss(eye+EFunction_Error ,vec_in_local_coord);  // EigenfunctionError
        vec_in_local_coord = vec_in_local_coord/getMax(abs(vec_in_local_coord));
        
        
        if (j < dimension /2){
            cout << " w_" << j << "   = ";
            for (int i = 0 ; i  < dimension ; i++){
                U_coord[i][j] = vec_in_local_coord[i];
                U_coord_pt[i][j] = vec_in_local_coord[i];
            }
        }
        else if (j == dimension /2){
            cout << " phi'  = ";
        }
        cout << vec_in_local_coord << endl;
    }
    
//     cout << " U = " << U_coord << endl;
//     cout << " U_pt = " << U_coord_pt << endl;
//     cout << " |U_pt| = " << euclNorm(U_coord_pt) << endl;
    
    
    
//     END 

    interval eps_0;
    for (int i = 0 ; i < dimension;i++){
        eps_0 = intervalHull(eps_0,EFunction_Error[0][i]);
    }
    
    interval V_inverse_norm = euclNorm(  krawczykInverse(A_s));
    interval E_norm = V_inverse_norm.right() * eps_0.right() * sqrt( dimension )  ; 
    E_norm = E_norm /(1-E_norm );
    
    
    
    cout << "eps0 = "<<  eps_0 << endl;
    cout << "E_norm = "<<  E_norm << endl;
    
    IVector eigenvalues = localStableBig.eigenvalues;
//     cout << " eigenvalues= " << eigenvalues<< endl;

    bool L_PLUS = checkL_plus(U_coord,eps_0,eigenvalues,E_norm,U_coord_pt,U_coord_nbd);

    return L_PLUS;
}

bool propagateManifold::checkL_plus( IMatrix U_coord,  interval eps_0,IVector eigenvalues , interval E_norm , IMatrix U_coord_pt, IMatrix U_coord_nbd){
    
//     cout << "U_coord = " << U_coord << endl;
    
    vector < IMatrix > Gamma_List;
    vector < IMatrix > Beta_List;
//     vector < IMatrix > VinvU_List;
    vector < IMatrix > VinvU_List_pt;
    vector < IMatrix > VinvU_List_nbd;
    IVector EE_norm_list(dimension/2);
    IMatrix Gamma_local(dimension/2,dimension/2-1);
    IMatrix Beta_local(dimension/2,dimension/2-1);
    IMatrix VinvU_local_pt(dimension,dimension/2-1);
    IMatrix VinvU_local_nbd(dimension,dimension/2-1);
    for (int k =0;k<dimension/2;k++) //remove k
    {
        int k_adjust =0;
        for (int j = 0 ; j<dimension/2;j++) 
        {
            if(j==k){k_adjust =1;continue;}
            for (int i = 0 ; i < dimension/2;i++){
                Gamma_local[i][j-k_adjust] = U_coord[i][j];
                Beta_local[i ][j-k_adjust] = U_coord[i+dimension/2][j];
                
                VinvU_local_pt[i][j-k_adjust] = U_coord_pt[i][j];
                VinvU_local_pt[i+dimension/2][j-k_adjust] = U_coord_pt[i+dimension/2][j];
                
                VinvU_local_nbd[i][j-k_adjust] = U_coord_nbd[i][j];
                VinvU_local_nbd[i+dimension/2][j-k_adjust] = U_coord_nbd[i+dimension/2][j];
            }
            
        }
        
        
        Gamma_List.push_back(Gamma_local);
        Beta_List.push_back(Beta_local);
        
        VinvU_List_pt.push_back(VinvU_local_pt);
        VinvU_List_nbd.push_back(VinvU_local_nbd);
        
//         cout << " Vinv U = " << VinvU_local_pt << endl;
        
        cout << " ||E|| .||V^-1 U|| = " << E_norm * euclNorm(VinvU_local_pt ) << endl ;
        
        EE_norm_list[k] = E_norm * euclNorm(VinvU_local_pt ) ;
        
        
    }
    
    
    
    interval nu_1 = eigenvalues[0];                 // Largest 
    interval nu_n = eigenvalues[dimension/2-1];     // Smallest
    
//     cout << " nu_1 = " << nu_1 << endl;
//     cout << " nu_n = " << nu_n << endl;
    
    bool L_PLUS = 0;
    for (int k = 0;k<dimension/2;k++){
        L_PLUS = checkL_plus_local( Gamma_List[k], Beta_List[k],eps_0, nu_1, nu_n , EE_norm_list[k]);
        if (L_PLUS ==1)
            break;
    }
    
    return L_PLUS;
}

bool propagateManifold::checkL_plus_local( IMatrix Gamma, IMatrix Beta,interval eps_0, interval nu_1 , interval nu_n , interval EE_norm){
    
    interval epsilon_beta = compute_epsilon_beta( Gamma, Beta, EE_norm);
    
    
// // //     A_s = (*(*pStable).pF).A;
// // //     krawczykInverse( A_s )*
    
    
    if (epsilon_beta < 0)
        return 0;
    
    epsilon_beta = epsilon_beta .right();
    eps_0  = eps_0 .right();
    
    cout << " eps_0 = " << eps_0 << endl; 
    cout << "epsilon_beta  = " << epsilon_beta  << endl;
//     epsilon_beta =0;
    
    interval d_min =1;
    interval d_max =1;
    int n = dimension/2;
    
    
    IMatrix A_s = (*(*pStable).pF).A;
    
    IMatrix pi_1_Vs(dimension/2,dimension/2);
    IMatrix pi_1_Vu(dimension/2,dimension/2);
    IMatrix pi_2_Vs(dimension/2,dimension/2);
    IMatrix pi_2_Vu(dimension/2,dimension/2);
    
    for (int i = 0 ; i < dimension/2 ; i++){
        for (int j = 0 ; j < dimension/2 ; j ++){
            pi_1_Vu[i][j] = A_s[i][j];
            pi_1_Vs[i][j] = A_s[i][j+dimension/2];
            pi_2_Vu[i][j] = A_s[i+dimension/2][j];
            pi_2_Vs[i][j] = A_s[i+dimension/2][j+dimension/2];
        }
    }

    
//     cout << " A_s  = " << A_s << endl;
//     
//     cout << " pi_1_Vs  = " << pi_1_Vs << endl;
//     cout << " pi_1_Vu  = " << pi_1_Vu << endl;
//     cout << " pi_2_Vs  = " << pi_2_Vs << endl;
//     cout << " pi_2_Vu  = " << pi_2_Vu << endl;
//     
//     cout << " |pi_1_Vs|  = " << euclNorm(pi_1_Vs) << endl;
//     cout << " |pi_1_Vu|  = " << euclNorm(pi_1_Vu) << endl;
//     cout << " |pi_2_Vs|  = " << euclNorm(pi_2_Vs) << endl;
//     cout << " |pi_2_Vu|  = " << euclNorm(pi_2_Vu) << endl;
    
    interval C_Q_new= sqrt(n)*( euclNorm(pi_1_Vs) + euclNorm(pi_1_Vu) + euclNorm(pi_2_Vs) + euclNorm(pi_2_Vu) + 2* eps_0*sqrt(n));
    
    interval C_P    = ( 2*sqrt(n) * sqrt( 2* nu_1 * d_max)) / ( 1 - eps_0 * sqrt(n) * sqrt( 2*nu_1*d_max )) ;
    interval C_Q    = sqrt(n)*( sqrt(2/(nu_n*d_min))+ sqrt( 2* nu_1 * d_max) + 2* eps_0*sqrt(n));
    
    C_P = C_P.right();
    C_Q = C_Q.right();
    
    
    
    
    interval C_M1   = eps_0 * C_Q * ( 2 + eps_0 * C_Q); 
    interval C_M2   = eps_0 * C_P * ( 2 + eps_0 * C_P) +  2* epsilon_beta*(1+ eps_0 *C_P)  + sqr(epsilon_beta); 
    
    interval sum_for_neumann_test =  eps_0 * n * sqrt( 2 * nu_1 * d_max);
//     cout << " sum_for_neumann_test = " << sum_for_neumann_test << endl; 
    
    bool NEUMANN_SERIES_TEST;
    if (( C_M1 <1 ) && ( C_M2 < 1) && (sum_for_neumann_test <1))
        NEUMANN_SERIES_TEST =1;
    else
        NEUMANN_SERIES_TEST =0;
    
    
    cout << " C_P  = " << C_P << endl; 
    cout << " C_Q  = " << C_Q << endl; 
    cout << " C_Q!!= " << C_Q_new << endl; 
    cout << " C_M1 = " << C_M1 << endl; 
    cout << " C_M2 = " << C_M2 << endl << endl; 
    
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
    
    cout << " total_sum = " << total_sum << endl; 
    
    bool L_PLUS = (NEUMANN_SERIES_TEST && ( total_sum <1));
    
    return L_PLUS ;
}

interval propagateManifold::compute_epsilon_beta( IMatrix Gamma, IMatrix Beta, interval EE_norm){
    
    
    IMatrix Gamma_center = midMatrix(Gamma);
    IMatrix Gamma_delta  = Gamma - Gamma_center;
    
//     cout << " Gamma_center = " << Gamma_center << endl;
//     cout << " Gamma_delta = " << Gamma_delta << endl;
    
//     cout << " Gc^T * Gc     = " << transpose(Gamma_center)*Gamma_center <<endl;
//     cout << " 2* Gd^T * Gc  = " << interval(2)*((transpose(Gamma_delta)*Gamma_center)) <<endl;
//     cout << " Gd^T * Gd     = " << transpose(Gamma_delta)*Gamma_delta  <<endl;
    
    
    IMatrix GcT_by_Gc = transpose(Gamma_center)*Gamma_center    ;
    IMatrix GdT_by_Gc = transpose(Gamma_delta)*Gamma_center     ;
    IMatrix GdT_by_Gd = transpose(Gamma_delta)*Gamma_delta      ;
    
    
    interval mu_c = ml( transpose(Gamma_center)*Gamma_center );
    cout << "mu_c = " << mu_c << endl;
    
    interval mu_Rayleigh = mu_c - interval(2)*euclNorm(GdT_by_Gc) - euclNorm(GdT_by_Gd); 
    
    
    
//     cout << "mu_Rayleigh       = " << mu_Rayleigh << endl;
    
    cout << "mu_Rayleigh ( -E) = " << mu_Rayleigh - interval(2)*EE_norm*(euclNorm(Gamma_center)+euclNorm(Gamma_delta)) - sqr(EE_norm)<< endl;
    
    interval mu_old = ml( transpose(Gamma)*Gamma);
    

    
//     cout << " Gamma = " << Gamma << endl;
//     cout << " Gamma^t*Gamma = " << transpose(Gamma)*Gamma << endl;
    
//     cout << "mu_old = " << mu_old << endl;
    
    interval mu = mu_Rayleigh;
    
    if (mu.left() < 0)
        return interval(-1);
    
//     cout << " ||B|| (old) = " << euclNorm(Beta) << endl;
    cout << " ||B|| ( +E) = " << euclNorm(Beta) +EE_norm<< endl;
    interval epsilon_beta = euclNorm(Beta) / sqrt(mu);
    
    
    
//     epsilon_beta = interval(0.0000001);
    
    
    return epsilon_beta ;
    
}



