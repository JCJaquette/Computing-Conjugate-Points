 
#include "localManifold.h"

IVector localManifold::constructU( IVector U_flat)
{
  IVector U(dimension);
  
  interval norm = 0 ;
  for (int i=0;i<dimension/2;i++)
  {
    norm += U_flat[i]^2;
  }
  norm = sqrt(norm);
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
//   It returns pt_out (length dimension)  putting X_pt into the stable/unstable coordinates
  
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
//   
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
  

IMatrix localManifold::boundDFU( IVector U)
{
  int subdivision_dim = dimension/2;
  
//   We Create the index list
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


/////////////// -- From Capinski, Wasieczko-Zajac 2017
// This function computes the bound on the logarithmic min-norm of a matrix
interval ml(const IMatrix &A)
{
	capd::vectalg::EuclLNorm<IVector,IMatrix> l;
	return -l(-A);
}

bool localManifold::checkIsolatingBlock( IMatrix DFU, const IVector &U)
{
// We use the same technique to check the isolating block condition as in Capinski, Wasieczko-Zajac 2017 (see also their code & comments)

//    This check is the same for stable / unstable manifold
  capd::vectalg::EuclNorm<IVector,IMatrix> Norm;
  
  interval r1 = abs(U[0]).right();
  interval r2 = abs(U[dimension/2]).right();
  
  
  vector < IMatrix > ABCD = blockDecompose( DFU,dimension );
  
  IMatrix Fxx = ABCD[0]; // A
  IMatrix Fxy = ABCD[1]; // B
  IMatrix Fyx = ABCD[2]; // C
  IMatrix Fyy = ABCD[3]; // D
  
  IMatrix C11(dimension/2,dimension/2);
  IMatrix C22(dimension/2,dimension/2);
  
  IVector lambda(dimension/2);
  IVector beta(dimension/2);
  
  for(int i =0;i<dimension/2;i++)
  {
    C11[i][i] = Fxx[i][i].mid();
    lambda[i] = Fxx[i][i].mid();
    C22[i][i] = Fyy[i][i].mid();
    beta[i]   = Fyy[i][i].mid();
  }
  
  IMatrix R11 = Fxx - C11;
  IMatrix R22 = Fyy - C22;
  IMatrix R12 = Fxy;
  IMatrix R21 = Fyx;
  
  interval lambda_min = - getMax(-lambda);
  interval beta_max = getMax(beta);
  
  
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

bool localManifold::checkRateCondition(IVector U_flat ) //TODO Separate this from the isolating block computuation
{
//    We calculate mu_1, mu_2, xi

  IVector U = constructU( U_flat);
  
  cout << " U = " << U << endl;

  //   M  is DF[U]
  IMatrix M = boundDFU(U);
  
  return checkRateCondition(M);
}


bool localManifold::checkRateCondition(IMatrix DFU) //TODO Separate this from the isolating block computuation
{
    
  vector < IMatrix > ABCD = blockDecompose( DFU,dimension );
  
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
  
  
  interval mu = l(Fyy)+euclNorm(Fyx)/L;
//   interval mu_2 = l(Fyy)+euclNorm(Fxy)/L;
  
  xi = ml(Fxx) - L*euclNorm(Fxy); //This is a class variable. 
  
//   cout << " l(Fyy) = " << l(Fyy)<< endl ;
//   cout << " euclNorm(Fyx) = " << euclNorm(Fyx)<< endl ;
  
  cout << "  mu  = " << mu << endl ;
  cout << "  xi  = " << xi << endl ;
  
  
  bool conditions ;
  
  bool condition_1 = (mu <0);
  bool condition_2 = (0< xi);
//   bool condition_3 = (mu_2< xi);
  
  if ( condition_1 && condition_2 )
    conditions = 1;
  else 
    conditions = 0;
  
  
 return conditions; 
}


 bool  localManifold::checkConditions( IVector U_flat )
{

  IVector U = constructU( U_flat);
  
  
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
  
  


interval localManifold::ErrorEigenfunction( void)
{
//  TODO This is a rough calculation and needs to be properly implemented
   
  
  
  interval K =1; // TODO
  
  capd::vectalg::EuclNorm <IVector,IMatrix> euclNorm;
  
//   We estimate A(U) -A(p) --- this is NOT correct
  
  IVector U_flat(dimension/2);
  interval scale = 0.00001;
  for (int i =0;i<dimension/2;i++){ U_flat[i]=scale*interval(-1,1);}
  IVector U = constructU(U_flat);
  
  
  
  
  interval c = sqrt(1 + sqr(L) );
  IMatrix guess = (*pF)[U] - (*pF)[(*pF).p];
  
  
  interval C = c*sqr(euclNorm(guess)); // We should really be taking a second derivative ..... and calculating this in a different class. 
//   cout << " C = " << C<< endl;
  interval eta = xi; // This needs xi to already have been computed. 
  interval lambda = K*C/eta;
  
  lambda = lambda.right(); // This reduces wrapping effect. 
  
  interval error = lambda/(1-lambda);
  error = error/10000;
  return error;
}





