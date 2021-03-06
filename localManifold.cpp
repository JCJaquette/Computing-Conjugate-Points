#include "localManifold.h"

// BEGIN Base class: "localManifold" methods
IVector localManifold::constructU( IVector U_flat)
{
//     Input U_flat is in local coords, of length  n=dim/2
//     Returns an interval vector in local coords 
//      which encloses the cooresponding point/nbd on the stable/unstable manifold.
//      (with error bounds)
  IVector U(dimension);
  
//   Compute the Euclidean norm of U_flat
  interval norm = 0;
  for (int i=0;i<dimension/2;i++)
  {
    norm += U_flat[i]^2;
  }
  norm = sqrt(norm);
//   Defines the error as    +/- L* || U_flat|| 
  interval error = L*interval(-1,1)*norm; 
  
  for (int i = 0;i<dimension;i++)
  {  
    if (stable)
    {
      if( i < dimension/2)
        U[i]=error; 
      else
        U[i]=U_flat[i-dimension/2]; 
    }
    else // Unstable case
    {
      if( i < dimension/2)
        U[i]=U_flat[i]; 
      else
        U[i]=error;     
    }
  }
  return U;
//   cout << "U = " << U << endl;
}

IVector localManifold::getPoint( IVector X_pt)  
{
//   This funtion takes input X_pt (length dimension /2), 
//   It returns pt_out (length dimension)  putting X_pt into the (LOCAL) stable/unstable coordinates
  
  IVector pt_out(dimension);
 
  
  for (int i = 0;i<dimension;i++)
  {  
    if (stable)
    {
        if( i < dimension/2)
            pt_out[i]=0; 
        else
            pt_out[i]=X_pt[i-dimension/2]; 
    }
    else // Unstable case
    {
        if( i < dimension/2)
            pt_out[i]=X_pt[i]; 
        else
            pt_out[i]=0;     
    }
  }
  
//    We return the output in LOCAL coordinates 
  return pt_out;
}

vector < IVector > localManifold::getPointNbd( IVector XY_pt, IVector XY_nbd)
{
//     Receives as input a point and neighborhood of lengths n=dim/2 , in local coords
//     Outputs the point in global coordinates, and the neighborhood in local.

//      (with error bounds)
  vector < IVector >  output;
  
  //    We get the local linearization
  IMatrix A_i = (*pF).A;   
  //    We get the local point offset
  IVector local_pt = getPoint(XY_pt);
  //    We define p_i to be the equilibrium + offset
  IVector p_i = (*pF).p + A_i*local_pt;   
  
  // 	We need to translate this neighborhood to get the right errorbounds. 
    IVector Uxy =constructU(XY_pt+XY_nbd);  
  //     We translate back to enclose zero
    Uxy 	= Uxy - local_pt;
  
//     cout << " A_i   " << A_i << endl;
//     cout << " Point " << p_i << endl;

//    We return this output in GLOBAL coordinates 
    output.push_back(p_i);
//    We return this output in LOCAL  coordinates     
    output.push_back(Uxy);
    
  return output;
}


IVector localManifold::projectPoint( IVector XY_pt)
{
//     Takes an input point in ambient coordinates. 
//     If at the stable (unstable) manifold, we project the point >> ( pt - equilibrium ) << onto the stable (unstable) eigenspace.
//     Returns a vector of length  n = dimension/2
  IVector local_coord(dimension/2);
  
  if (stable)
    XY_pt = XY_pt- (*pF).p;
  
  IVector eigen_vector;
  for (int i = 0 ; i< dimension/2;i++)
  {
//     Get Eigenvector 
    if (stable)
      eigen_vector = getEigenvector( i + dimension/2);
    else
      eigen_vector = getEigenvector( i );
    
//    Take scalar product of XY_pt with  eigenvector i 
    interval a_i = eigen_vector * XY_pt;
//    Store in local_coord
    local_coord[i]=a_i;
//    Subtract from XY_pt
    XY_pt = XY_pt - a_i * eigen_vector;
  }
  return local_coord;
}

void localManifold::constructDW( )
{
//   For the chart of the unstable manifold,   W(x)= (x , \omega^u (x)
//     this function computes the derivative DW. 
//     (And does the analogous thing for the stable manifold.) 

    
  int k= dimension;
  DW = IMatrix(k,k/2);
  
  IMatrix identity = IMatrix(k/2,k/2);
  IMatrix errorMat = IMatrix(k/2,k/2);
  
  interval error = L*interval(-1,1);
  for (int i = 0;i<k/2;i++)
  {
    identity[i][i] = 1;
  }

  
    for (int i = 0;i<k;i++)
    {
        for (int j=0;j<k/2;j++)
            {
            if (i< k/2)
            {
                if (stable) //     We put the error on the top for stable	    
                    DW[i][j] = error;
                else//     We put the identity on the top for unstable
                    DW[i][j] = identity[i][j];
            }
            else
            {
                if (stable)//     We put the identity on the bottom for stable
                    DW[i][j] = identity[i-k/2][j];
                else //     We put the error on the bottom for unstable
                    DW[i][j] = error;
            }
        }
    }
//   cout << "DW" << DW << endl;
}
  



bool localManifold::checkIsolatingBlock( IMatrix DFU, const IVector &U)
{
// We use the same technique to check the isolating block condition as in Capinski, Wasieczko-Zajac 2017 (see also their code & comments)

//    This check is the same for stable / unstable manifold
  capd::vectalg::EuclNorm<IVector,IMatrix> Norm;
  
  interval r1 = abs(U[0]).right();              // r_u
  interval r2 = abs(U[dimension/2]).right();    // r_s
  
  
  vector < IMatrix > ABCD = blockDecompose( DFU,dimension );
  
  IMatrix Fxx = ABCD[0]; // A
  IMatrix Fxy = ABCD[1]; // B
  IMatrix Fyx = ABCD[2]; // C
  IMatrix Fyy = ABCD[3]; // D
  
  IMatrix C11(dimension/2,dimension/2);  // mid( diag( Fxx ) )    -- matrix
  IMatrix C22(dimension/2,dimension/2);  // mid( diag( Fyy ) )    -- matrix
  
  IVector lambda(dimension/2);           // mid( diag( Fxx ) )    -- vector
  IVector beta(dimension/2);             // mid( diag( Fyy ) )    -- vector
  
  for(int i =0;i<dimension/2;i++)
  {
    C11[i][i] = Fxx[i][i].mid();
    lambda[i] = Fxx[i][i].mid();
    C22[i][i] = Fyy[i][i].mid();
    beta[i]   = Fyy[i][i].mid();
  }
  
//   Define 'defect' terms
  IMatrix R11 = Fxx - C11;
  IMatrix R22 = Fyy - C22;
  IMatrix R12 = Fxy;
  IMatrix R21 = Fyx;
  
  interval lambda_min = - getMax(-lambda);  // approx slow unstable eigen
  interval beta_max = getMax(beta);         // approx slow stable   eigen
  
//   The local vector field has already moved the equilibrium to zero. 
  IVector F_of_zero = (*pF)(IVector(dimension)) ;
  
  
  interval unstable_bound = lambda_min - Norm(R11) * r1 - Norm(F_of_zero)- Norm(R12)*r2;
  interval stable_bound   = beta_max +   Norm(R21) * r1 + Norm(F_of_zero)+ Norm(R22)*r2;
  
  
//   cout << "Unstable Bound = " << unstable_bound << endl;
//   cout << "  Stable Bound = " << stable_bound << endl;
  
  bool output;
  if (( unstable_bound > 0) && ( stable_bound < 0))
    output = 1;
  else
    output =0;
  
  return output;
}


bool localManifold::checkRateCondition(IMatrix DFU) 
{
    
  vector < IMatrix > ABCD = blockDecompose( DFU,dimension );
//   DFU = [ A B ]
//         [ C D ]
  
//   We block decompose our function
  IMatrix Fxx ;// A
  IMatrix Fxy ;// B
  IMatrix Fyx ;// C
  IMatrix Fyy ;// D
//   This proof is constructed to validate unstable manifolds. 
//   To validate stable manifolds, we must change F \mapsto -F, and put the (formerly) stable Fyy component into the Fxx block
  if (stable)
  {
    Fxx = -ABCD[3]; // D
    Fxy = -ABCD[2]; // C
    Fyx = -ABCD[1]; // B
    Fyy = -ABCD[0]; // A
  }
  else
  {
    Fxx = ABCD[0]; // A
    Fxy = ABCD[1]; // B
    Fyx = ABCD[2]; // C
    Fyy = ABCD[3]; // D
  }
  
//   cout << " A = " << Fxx << endl;
//   cout << " B = " << Fxy << endl;
//   cout << " C = " << Fyx << endl;
//   cout << " D = " << Fyy << endl;
  
  
  capd::vectalg::EuclLNorm<IVector,IMatrix> l;
  capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm;
  
//   NOTE   The calculation here for mu and xi could be improved (maybe by 30%?) 
//          if we were to subdivide Fxx, Fxy,Fyx,Fyy before taking the norm.
  
  interval mu = l(Fyy)+euclNorm(Fyx)/L;
  
  // xi - is a class variable. 
  xi = ml(Fxx) - L*euclNorm(Fxy); 
  
  
  cout << "  mu  = " << mu << endl ;
  cout << "  xi  = " << xi << endl ;
  
  bool conditions ;
  
  bool condition_1 = (mu <0);
  bool condition_2 = (0< xi);
   
  if ( condition_1 && condition_2 )
    conditions = 1;
  else 
    conditions = 0;
  
  
 return conditions; 
}


 bool  localManifold::checkConditions( void  )
{

  IVector U = constructU( U_flat_global );
  
  
  cout << " U = " << U << endl;

  //   M  is DF[U]
// // //   IMatrix M_old = (*pF)[U];
  IMatrix M = boundDFU(U);
  
  bool Isolating_block = checkIsolatingBlock( M, U);
  bool Rate_Conditions = checkRateCondition(M );
  
  
  cout << " Rate  Condition " << Rate_Conditions << endl;  
  cout << " Block Condition " << Isolating_block << endl;  
   
  return (Isolating_block && Rate_Conditions);
   
}
  
  


IVector localManifold::containmentRatios( IVector point_test){
    IVector ratios(dimension/2);
// Returns, a vector which, for each coordinate, 
//     is the ratio of the input point_test, and the neighborhood **U_flat** of the local manifold
    for (int i = 0 ; i<dimension/2;i++){
        if ( point_test[i] > 0 )
            ratios[i] = point_test[i]/(U_flat_global[i].right());
        else
            ratios[i] = point_test[i]/(U_flat_global[i].left());
    }
    return ratios;
}


IMatrix localManifold::boundDFU( IVector U)
{
    //   Divides U into N^n pieces and bounds DF on it.
    //      N = subdivisionNUM
    //      n = dimension/2
    //   If part of the stable manifold, only divides the stable part / parity for unstable manifold.
    
  int subdivision_dim = dimension/2;
  
//   We create the index list ( either all the stable coords or all the unstable coords) 
  vector < int > index_list(dimension/2);

//     We subdivide the stable part if we are looking at the stable manifold
  if (stable)
  {
    for (int i=0;i<dimension/2;i++)
    {
      index_list[i] = i +dimension/2;
    }
  }
//     We subdivide the unstable part if we are looking at the unstable manifold
  else
  {
    for (int i=0;i<dimension/2;i++)
    {
      index_list[i] = i;
    }
  }

//    We get a bound on the first part
  vector < int > part_list(dimension/2);
  
  IVector U_part = getSubdivision(U, index_list , part_list , subdivisionNUM);
  IMatrix A = (*pF)[U_part];
  
  
//    We make a single for loop to go through all the subdivision dimensions
  int sum;
  int N = subdivisionNUM;
  int n = subdivision_dim;
  int N_to_i;
  for (int j = 0;j< (int_pow(N,n));j++)
  {
//     We index the part we want
    sum = 0;
    N_to_i = 1;
    for (int i = 0 ; i<n;i++)
    {
      part_list[i] = ( (j-sum) / N_to_i ) % subdivisionNUM;
      if (i < n-1)
      {
        sum += part_list[i]*N_to_i;
        N_to_i = N_to_i * N;
      }
    }
//     We get the part 
    U_part = getSubdivision(U, index_list , part_list , subdivisionNUM);
    A = intervalHull(A,(*pF)[U_part]);

  }
  
  return A;
}

// END

// // // // // // // // // // // // // 
// BEGIN Derived class: " localManifold_Eig "  methods
// // // // // // // // // // // // // 

localManifold_Eig::localManifold_Eig(localVField &pF_, IVector U_flat_global_,interval L_, bool stable_, int subdivisionNUM_)
    : localManifold( pF_, U_flat_global_,L_, stable_, subdivisionNUM_){}

localManifold_Eig::localManifold_Eig(const localManifold & lM_)
    : localManifold( lM_ ){}
        



void localManifold_Eig::ErrorEigenfunction( void)
{  
//     This function uses the the estimate in Section 3 of the paper to bound the constant 
//     
//              \lambda_{T_-} := (1/xi) * K_- * C_G * r_u * ||A_0|| * sqrt{  1+ vartheta^2 }
//     
//     Recall "vartheta" in the paper is "L" in the code. 
//     Then, it defines the class variable.
//     
//              eps_unscaled  := lambda/(1-lambda)
//     
//     which, as in Proposition 2.3, is used to bound the error of the eigenfunction.  
  

    
// BEGIN Compute C_G 

  IVector U = constructU(U_flat_global);

  cout << "U = " << U << endl;
    
//   Note G(v) is a function from R^n to R^n, and does not depend on velocity, hence the need for pi_1
  IMatrix pi_1(dimension,dimension);
  for (int i =0;i<dimension/2;i++){ pi_1[i][i]=1;}
  
  IMatrix pi1_A = pi_1 * (*pF).A;
    
//   We compute the third derivative
  (*(*pF).f).setDegree(2);
  
  IMatrix Df(dimension,dimension);
  IHessian Hf(dimension,dimension);
  
  // simultaneous computation of value, derivative and normalized hessian
  // NOTE Hf contains second order Taylor coefficients of f at x, i.e. normalized derivatives.
  IVector y = (*(*pF).f)( (*pF).p +  pi1_A*U   ,Df,Hf); 
      
  (*(*pF).f).setDegree(1);
  
  // Most of the components in the tensor are zero; this is somewhat related to precomposing with pi_1.
  IHessian DDDG_small = compressTensor(  Hf , dimension);     
  
  interval tensor_norm = tensorNorm( DDDG_small , dimension/2);  
  interval C_G = right( tensor_norm);  

//   END
  
  interval K = computeK();   
  interval eta = xi; // This needs xi to already have been computed. 
  
  capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm;   // Use matrix norm. 
  interval norm_A0 = euclNorm((*pF).A);
  interval r_u = euclNorm(abs(U_flat_global)) ; 
  
  
//   cout << " C_G      = " << C_G << endl; 
//   cout << " eta      = " << eta << endl;  
//   cout << " norm_A0  = " << norm_A0 << endl;  
//   cout << " sqrt(1+L^2)  = " << sqrt(1 + sqr(L) ) << endl;  
//   cout << "        L^2   = " << sqr(L) << endl;  
  
  interval lambda = K* C_G * r_u *norm_A0 *sqrt(1 + sqr(L) )  /eta;
  
  cout << " lambda = " << lambda<< endl;
  lambda = lambda.right(); // This reduces wrapping effect. 
  
  interval error = lambda/(1-lambda);
  eps_unscaled = error;
  
  cout << " error = " << error<< endl;
  
}


interval localManifold_Eig::computeK( void )
{    
//     Computes the constant K_-
//     Stores validated bounds on the eigenvalues   -- K_store
//     Stores error bounds on the eigenvectors.     -- Eigenvector_Error
    
    IMatrix A_u = (*pF).A;                                  // Eigenvectors
    IMatrix A_infty = (*(*pF).f)[(*pF).p];                  // Asymptotic Matrix
    
    IMatrix Lambda = gaussInverseMatrix(A_u)*A_infty*A_u;   // Eigenvalues     
    
//     Store class variable "eigenvalues", a validated enclosure
    eigenvalues = boundEigenvalues( Lambda);
//     cout << " eigenvalues =" << eigenvalues << endl;
    
//     Develop a bound on the eigenvectors, with out put as center + error
    IVector Lambda_vec(dimension) ;
    for (int i =0;i<dimension;i++){ Lambda_vec[i]=Lambda[i][i];}

    vector < IMatrix >  output= boundEigenvectors( A_infty, A_u , Lambda_vec); 
    
    IMatrix Q_center = output[0];  
    IMatrix Q_error  = output[1];
    
    IMatrix Q = Q_center +Q_error;    
    

    capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm;     // Use matrix norm.
    interval K =euclNorm(Q)* euclNorm( krawczykInverse(Q));     
    
//     cout << " K      = " << K << endl;
    
//     Store eigenvector error
    Eigenvector_Error = Q_error;
//     Store K constant.
    K_store = K;
    
    return K;
}

bool localManifold_Eig::checkConjugatePointsBelowLminus( void){
// checks whether there are NO Conjugate Points Below L_- ; see Prop 2.5 / eq (2.12). That is, whether 
//     \lambda_T- / (1-\lambda_T-)    <    1 /  \sqrt{ n * ( 1 + ||\mathcal{M}_-|| ) }  
    
//     This function has to be called AFTER computeEigenError_minus_infty.
//     If it is called beforehand, the test automatically fails, as  *eps_unscaled* will be undefined/=0
    
//     OUTPUT
//          true    ---  Everything is GOOD! There are NO Conjugate Points Below L_-
//          false   ---  Everything is BAD!! There MAY BE Conjugate Points Below L_-
    
//     We compute the 2n by 2n matrix containing $\mathcal{M}_-$ in the lower left block.
    IMatrix M_large = (*(*pF).f)[ pF -> p];
    
//     We get the lower left block.
    vector < IMatrix > M_vec = blockDecompose( M_large, dimension);
    IMatrix M_minus = M_vec[2];
    
//     We bound its matrix norm 
    capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm;
    interval M_bound = euclNorm(M_minus); 
//     cout << " M_bound = " << M_bound << endl;
    
    interval n = interval(dimension/2);    
    interval RHS_2p12 = 1/ sqrt( n*(1+ M_bound) );
    
//     Recall, as computed in *computeEigenError_minus_infty* that 
//          eps_unscaled =  \lambda_T- / (1-\lambda_T-)
    bool output = ( eps_unscaled  < RHS_2p12 );
    
    return output;
}


void localManifold_Eig::computeEigenError_minus_infty( void){
//     Computes the Eigenfunction error, and puts it all in the stable modes. 
//     Output/Effects:
//      -- Creates class matrix "Eu_m_Error_Final" , the distance W is from our approximate eigenvector basis
//     NOTE It would have been more efficient to combine the function "computeEigenError_minus_infty" and "computeEigenError_plus_infty". Alas.
    
    
//     We compute the error in our eigenfunction approximation. This computes ''eps_unscaled''
    ErrorEigenfunction( );
    
    IMatrix Error_mat(dimension,dimension/2);
    
    IMatrix A_u = (*pF).A;      // Approximate Eigenvectors
    
//  We compute the error on the frame matrix W. 
//  Suppose that \tilde{V}  is the true asymptotic eigenvector, and \bar{V} is our numerically approximate one, and \eps = eps_unscaled. Then
//          || W - \tilde{V} || <= \eps * || \tilde{V} ||      &  ( \tilde{V}  - \bar{V} )  \in  Eigenvector_Error 
//  Hence
//       W - \bar{V}  \in   Eigenvector_Error + [-1,1] * ( || \bar{V} || + ||  Eigenvector_Error|| )
//  This collected error term, we define as Error_mat, and it is in general ambient coordinates.  Later we put this into eigencoordinates.
    IVector V_col(dimension);
    IVector Error_col(dimension);
    for (int j =0;j<dimension/2;j++){
        V_col = getColumn(A_u,dimension,j);
        Error_col = getColumn(Eigenvector_Error,dimension,j);
        interval asymptotic_evect_norm = euclNorm(V_col) + euclNorm(Error_col ) ;
        interval eps_local = eps_unscaled *interval(-1,1);
        
        for (int i=0;i<dimension;i++){
            Error_mat[i][j] = Eigenvector_Error[i][j] + eps_local ;
        }        
    }
//     cout << " Error_mat " << Error_mat<< endl;
    
//  We put our error into eigen-coordinates   We define the unstable(top) and stable(bottom) part of the error in our unstable eigenvectors. 
//     The error here goes from e-4 / e-5 in each component equally to e-3 / e-4 in the eigen-directions.
    IMatrix E_u(dimension/2,dimension/2);
    IMatrix E_s(dimension/2,dimension/2);
    IVector V_sol(dimension);

    for (int j =0;j<dimension/2;j++){
        V_col = getColumn(Error_mat,dimension,j);
        
        V_sol = gauss(A_u,V_col);
        
        for (int i=0;i<dimension/2;i++){
            E_u[i][j] = V_sol[i];
            E_s[i][j] = V_sol[i+dimension/2];
        }
    }
    IMatrix eye = identityMat(dimension/2);
    
//  We do a change of coordinates to put all the error in the stable coordinates.    
    Eu_m_Error_Final = E_s*krawczykInverse(eye+E_u);
    
//         cout << "E_u = " << E_u << endl;
//         cout << "E_s = " << E_s << endl;
//         cout << "Eu_m_Error_Final= " << Eu_m_Error_Final<< endl;
    
}

void localManifold_Eig::computeEigenError_plus_infty( void){
//     Computes the Eigenfunction error at plus infty 
//     
//     Output/Effects:
//      -- Creates class matrix "Eigenfunction_Error_plus_infty", the distance W is from our approximate eigenvector basis in global coordinates 
//     NOTE It would have been more efficient to combine the function "computeEigenError_minus_infty" and "computeEigenError_plus_infty". Alas.
    
//     We compute the error in our eigenfunction approximation.
    ErrorEigenfunction( );
    
//  We compute the error on the frame matrix W. 
//  Suppose that \tilde{V}  is the true asymptotic eigenvector, and \bar{V} is our numerically approximate one, and \eps = eps_unscaled. Then
//          || W - \tilde{V} || <= \eps * || \tilde{V} ||      &  ( \tilde{V}  - \bar{V} )  \in  Eigenvector_Error 
//  Hence
//       W - \bar{V}  \in   Eigenvector_Error + [-1,1] * ( || \bar{V} || + ||  Eigenvector_Error|| )
//  NOTE This collected error term, we define as Error_mat, and it is in general ambient coordinates. 
    IMatrix Error_mat(dimension,dimension);     // All eigenfunction error
    
    IMatrix A_s = (*pF).A;                                  // Eigenvectors
    
    IVector V_col(dimension);
    IVector Error_col(dimension);
    for (int j =0;j<dimension;j++){
        V_col = getColumn(A_s,dimension,j);
        Error_col = getColumn(Eigenvector_Error,dimension,j);
        interval asymptotic_evect_norm = euclNorm(V_col) + euclNorm(Error_col ) ;
        interval eps_local = eps_unscaled *interval(-1,1);
        
        for (int i=0;i<dimension;i++){
            Error_mat[i][j] = Eigenvector_Error[i][j] + eps_local ;
        }        
    }
    
    Eigenfunction_Error_plus_infty = Error_mat;
    
}


IVector localManifold_Eig::getEigenError_minus_infty(int columnNumber){
    IVector v_out(dimension);
//     Outputs the error of the unstable eigenfunctions, all put into the stable eigendirections
    for (int i = 0 ; i<dimension/2;i++){
        v_out[i+dimension/2] = Eu_m_Error_Final[i][columnNumber];
    }
    
    return v_out;
}


// END
