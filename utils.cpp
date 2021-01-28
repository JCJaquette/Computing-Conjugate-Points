#include "utils.h"

interval part( interval x, int k, int N)
{
//   For an interval subdivided into N parts, we return, the k-th part
//     0  \leq  k  \leq  N-1
  
  interval y = x.left() + k* ( x.right() - x.left())/N +(x-x.left())/N;
  return y;
}

IMatrix toInterval(DMatrix A)
{
        int n=A.numberOfRows();
        int k=A.numberOfColumns();
        IMatrix B(n,k);
        for(int i=0;i<n;i++)
        {
                for(int j=0;j<k;j++) B[i][j]=A[i][j];
        }
        return B;
}


IVector toInterval(DVector x)
{
        int n=x.dimension();
        IVector y(n);
        for(int i=0;i<n;i++)
        {
             y[i]=x[i];
        }
        return y;
}


void swapColumns(DMatrix &A,int i,int j)
{
    double a;
    int k=A.numberOfRows();
    for(int n=0;n<k;n++)
    {
	a=A[n][i];
	A[n][i]=A[n][j];
	A[n][j]=a;
    }
}


void swapValues(DVector &v,int i, int j)
{
    double a=v[i]; 
    v[i]=v[j];
    v[j]=a;
}

void bubbleSort(DVector &v)
{
    int n=v.dimension();
    for(int i=n-1;i>=0;i--)
    {
      for(int j=0;j<i;j++) if(v[j]<v[j+1]) swapValues(v,j,j+1);
    }
}


void bubbleSortEigenvectors(DVector &v, DMatrix &A)
{
    int n=v.dimension();
    for(int i=n-1;i>=0;i--)
    {
      for(int j=0;j<i;j++) if(v[j]<v[j+1]) 
      {
	swapValues(v,j,j+1);
	swapColumns(A,j,j+1);
      }
    }
}

DMatrix coordinateChange(DMatrix Df)
{
//     Returns eigenvectors sorted  by eigenvalue (from largest (+) to smallest (-) )
  int n=Df.numberOfColumns();
  DVector rE(n), iE(n);         	// real and imaginary parts of eigenvalues
  DMatrix rVec(n,n), iVec(n,n); 	// real and imaginary parts of eigenvectors
  
//   The program is built for Matrices with real eigenvalues, and so we only return the real eigenvectors. 
//   If the input matrix has complex eigenvalues, things will fail to be validated later in the program.

  computeEigenvaluesAndEigenvectors(Df,rE,iE,rVec,iVec);
  bubbleSortEigenvectors(rE,rVec);

  return rVec;
}

void plot(interval x, interval y,ofstream &file)
{
//     Given two intervals x & y, writes to a file 
//     the average left end points, and 
//     the average right end points.
    file << (x.rightBound()+x.leftBound())/2. <<" ";
    file << (y.rightBound()+y.leftBound())/2. <<" ";
    file << (x.rightBound()-x.leftBound())/2. <<" ";
    file << (y.rightBound()-y.leftBound())/2. <<endl;
}



void constructMatrix(const interval &time,const IVector &left_end,const IVector &right_end,const IVector &v_value,const IVector &v_deriv, IMatrix &matrix_mod)
{
//   The Matrix matrix_out  will contain in its columns:
//   0 -- Time Interval +++++ ONLY IN ENTRY (0,0)
//   1 -- Left Endpoint
//   2 -- Right Endpoint
//   3 -- Bound on function Values
//   4 -- Bound on function Derivative
  matrix_mod[0][0] = time ;
  
  int dim = left_end.dimension();
  
  for (int i = 0; i < dim ; i++)
  {
    matrix_mod[i][1] = left_end[i];
    matrix_mod[i][2] = right_end[i];
    matrix_mod[i][3] = v_value[i];
    matrix_mod[i][4] = v_deriv[i];
  }

  
}

IVector getColumn(const IMatrix &A, int num_rows, int column)
{
 IVector output(num_rows);
 
 for(int i =0;i<num_rows;i++)
 { 
   output[i] = A[i][column];
 }
 
 return output;
}


vector<IMatrix> getTotalTrajectory(C0Rect2Set &s,interval T,int grid,ITimeMap &timeMap,IOdeSolver &solver)
{
  //   The Matrix contains in its columns:
  //   0 -- Time Interval +++++ ONLY IN ENTRY (0,0)
  //   1 -- Left Endpoint
  //   2 -- Right Endpoint
  //   3 -- Bound on function Values
  //   4 -- Bound on function Derivative  
  
  
//   NOTE This does not seem to work if solver has been used before
//   (in the sense that it does not make intermediate steps, not that it gives a wrong anwer)
    timeMap.stopAfterStep(true);
    interval prevTime(0.);
    
    
    vector<IMatrix> Matrix_List(0);
    int dim = s.dimension();
    
    do 
    {
      timeMap(T,s);
      interval stepMade = solver.getStep();
      const IOdeSolver::SolutionCurve& curve = solver.getCurve();
      interval domain = interval(0,1)*stepMade;


      
//       We get bounds on the value and derivative      
//       IVector value_bound(dim);// = vector_value; // TODO Write a function to do this
//       IVector deriv_bound(dim);// = vector_value; //
      
      for(int i=0;i<grid;++i)
      {
	
	
	
        interval subsetOfDomain = interval(i,i+1)*stepMade/grid;

        intersection(domain,subsetOfDomain,subsetOfDomain);
	
	
	
	//       We get the left and right end points
	IVector vector_left = curve(subsetOfDomain.left());
	IVector vector_right = curve(subsetOfDomain.right());
	
	
	
        IVector vector_value = curve(subsetOfDomain);
	IVector vector_deriv = curve.timeDerivative(subsetOfDomain); 
	
	

	interval localTime = prevTime + subsetOfDomain;
	
	
	  //       We output our results
	IMatrix local_matrix(dim,5);
  
	
	constructMatrix(localTime,vector_left,vector_right,vector_value,vector_deriv, local_matrix);
	
	Matrix_List.push_back(local_matrix);
	
	
      }      
      
      
//       cout << "Matrix added to list = " << local_matrix << endl;
      
// 	We update the time
      prevTime = timeMap.getCurrentTime();

    }while(!timeMap.completed());
    
    

    return Matrix_List;
}

vector<IVector> getTrajectory(C0Rect2Set &s,interval T,int grid,ITimeMap &timeMap,IOdeSolver &solver)
{
//     Only used in          boundaryValueProblem::localNormBound
//     
//     Not used for verified computing
    
    timeMap.stopAfterStep(true);
    interval prevTime(0.);
    vector<IVector> V(0);
    
    int dim = s.dimension();
    
    do 
    {
      timeMap(T,s);
      interval stepMade = solver.getStep();
      const IOdeSolver::SolutionCurve& curve = solver.getCurve();
      interval domain = interval(0,1)*stepMade;
      for(int i=0;i<grid;++i)
      {
        interval subsetOfDomain = interval(i,i+1)*stepMade/grid;

        intersection(domain,subsetOfDomain,subsetOfDomain);
        IVector v = curve(subsetOfDomain);
        cout << " v = " << v << endl;
	
        interval localTime = prevTime + subsetOfDomain;


        IVector v_aug(dim+1);
        
        for(int j = 0 ; j<dim+1;j++)	{	  v_aug[j]=v[j];	}
        
        v_aug[dim]=localTime;

        V.push_back(v_aug); // add v to the stack
        
      }
      prevTime = timeMap.getCurrentTime();

    }while(!timeMap.completed());
    return V;
}


IVector getSubdivision(const IVector &U, const vector < int > &index_list , const vector < int > &part_list , const int &subdivisionNUM)
{  
  IVector U_out = U;
  
  int length = index_list.size(); 
  
  for (int i =0;i<length;i++)
  {
    U_out[index_list[i]] = part(U[index_list[i]],part_list[i],subdivisionNUM);
  }   
    
  return U_out;
}


int int_pow(int base, int exp)
{
    // Computes the integer valued (base)^(exp)
 if (exp <0)
 {
   cout << " int_pow :: error :: negative exponent";
   abort();   
 }
 else
 {  
  int out = 1;
  for (int i = 0; i < exp;i++)
  {
    out = out * base;
  }
  return out;
 }
}

vector < IMatrix > blockDecompose( const IMatrix M , int dimension)
{
//    We decompose a square, even-dimensional matrix M = [A,B;C,D]
  vector < IMatrix > ABCD_out ; 

  
  IMatrix A(dimension/2,dimension/2) ;
  IMatrix B(dimension/2,dimension/2) ;
  IMatrix C(dimension/2,dimension/2) ;
  IMatrix D(dimension/2,dimension/2) ;
  
  for (int i = 0 ; i < dimension/2 ; i++)
  {
    for (int j = 0;j < dimension/2;j++)
    {
      A[i][j] = M[i][j];
      B[i][j] = M[i][j+dimension/2];
      C[i][j] = M[i+dimension/2][j];
      D[i][j] = M[i+dimension/2][j+dimension/2];
    }
  }
  
  
  ABCD_out.push_back(A);
  ABCD_out.push_back(B);
  ABCD_out.push_back(C);
  ABCD_out.push_back(D);
  
// //   cout << " A = " << A << endl;
// //   cout << " B = " << B << endl;
// //   cout << " C = " << C << endl;
// //   cout << " D = " << D << endl;
  
  
  return ABCD_out;
}

interval getMax(const IVector &V)
{
  int dim = V.dimension();
  interval max = V[0].right();
  
  for (int i =1;i<dim;i++)
  {
    if (V[i].right() > max )
      max = V[i].right();
  }
  return max;
}



int getMaxIndex(const IVector &V)
{
  int index =0;
  int dim = V.dimension();
  interval max = V[0].right();
  
  for (int i =1;i<dim;i++)
  {
    if (V[i].right() > max )
    {
      max = V[i].right();
      index =i;
    }
  }
  return index;
}


void print( vector <IVector> vector_list )
{
    int length = vector_list.size();
    for ( int i =0;i <length;i++)
    {
        cout << " i = " << i << "  --  ";
        cout << vector_list[i] << endl;
    }   
}

void print( vector <IMatrix> matrix_list)
{
    int length = matrix_list.size();
    for ( int i =0;i <length;i++)
    {
        cout << " i = " << i << "  --  ";
        cout << matrix_list[i] << endl;
    }
}

interval omega( IVector V, IVector W, int n){
    IVector JW(2*n);
    for (int i =0;i<n;i++){
        JW[i]=-W[i+n];
    }
    for (int i =0;i<n;i++){
        JW[i+n]=W[i];
    }
    
    interval prod = V*JW;
    return prod;
}

/////////////// -- From Capinski, Wasieczko-Zajac 2017
// This function computes the bound on the logarithmic min-norm of a matrix
interval ml(const IMatrix &A)
{
	capd::vectalg::EuclLNorm<IVector,IMatrix> l;
	return -l(-A);
}



IMatrix symplecticNormalization(IMatrix A_s, int dimension){
  interval normalization_factor;
  
  IVector local_vec_half_p(dimension/2);
  IVector local_vec_half_q(dimension/2);  
  
  IVector local_vec_full_p(dimension);
  IVector local_vec_full_q(dimension);
  
  for (int j = 0 ; j < dimension/2 ; j ++ ) {
      local_vec_full_p = getColumn(A_s,dimension,j);
      local_vec_full_q = getColumn(A_s,dimension,dimension-j-1);
      for ( int i =0; i < dimension/2 ; i++){
          local_vec_half_p[i] = local_vec_full_p[i];
          local_vec_half_q[i] = local_vec_full_q[i];
      }
//       Normalize the top 'z' part
      local_vec_full_p = local_vec_full_p/euclNorm(local_vec_half_p);
      local_vec_full_q = local_vec_full_q/euclNorm(local_vec_half_q);
      normalization_factor = omega(local_vec_full_p,local_vec_full_q,dimension/2) ;
//       Flip sign if necessary
      if (normalization_factor<0){
          normalization_factor = - normalization_factor;
          local_vec_full_q = -local_vec_full_q;
      }
      
      local_vec_full_p = local_vec_full_p/sqrt(normalization_factor);
      local_vec_full_q = local_vec_full_q/sqrt(normalization_factor);
      for (int i=0;i<dimension;i++){
          A_s[i][j]             = local_vec_full_p[i];
          A_s[i][dimension-j-1] = local_vec_full_q[i];
      }
//             cout << " local_vec_full_p + local_vec_full_q = " << local_vec_full_p + local_vec_full_q << endl;
//             cout << " local_vec_full_p - local_vec_full_q = " << local_vec_full_p - local_vec_full_q << endl;
            
//       cout << " || local_vec_full_p || = " << euclNorm(local_vec_full_p) << endl;
//       cout << " || local_vec_full_q || = " << euclNorm(local_vec_full_q) << endl;
  }
  return A_s;
}

// // // // IMatrix boundEyeInverseDefect(IMatrix E,int n){
// // // //     
// // // //     IMatrix rightMat = E;
// // // //     IMatrix eye = identityMat(n);
// // // //     
// // // //     for (int i = 0 ; i < n;i++){
// // // //         for( int j =0; j<n;j++){
// // // //             rightMat[i][j] = abs(rightMat[i][j]).right();
// // // //         }
// // // //     }
// // // //     cout << " rightMat " <<  rightMat <<endl;
// // // //     cout << " max_inverse_bound 0 " <<  krawczykInverse(eye+E) - eye   <<endl;
// // // //     IMatrix max_inverse_bound = krawczykInverse(eye+rightMat) - eye;
// // // //     cout << " max_inverse_bound 1 " <<  max_inverse_bound <<endl;
// // // //     max_inverse_bound = krawczykInverse(eye-rightMat) - eye;
// // // //     cout << " max_inverse_bound 2 " <<  max_inverse_bound <<endl;
// // // //     
// // // //     max_inverse_bound = interval(-1,1)* max_inverse_bound;   
// // // //     
// // // //     return max_inverse_bound;    
// // // // }

IMatrix identityMat( int n){
    IMatrix eye(n,n);
    for (int i = 0 ; i  < n;i++){
        eye[i][i] = 1;
    }
    return eye;
}

interval tensorNorm( IHessian DDDG , int n){
    
//     Matrix list 
    vector < IMatrix > Matrix_list;
    IMatrix A_local(n,n);
    
    for (int k = 0 ; k  < n;k++){// Outer loop
        for (int i = 0 ; i  < n;i++){ 
            for (int j = 0 ; j  < n;j++){
                A_local[i][j] = DDDG(i,j,k);
            }
        }
        Matrix_list.push_back(A_local);
    }
    
    IVector sliced_norms(n);
    IMatrix A_combo(n,n);
    
    for (int i = 0 ; i  < n;i++){ 
        sliced_norms[i] = euclNorm(Matrix_list[i]);
        A_combo += abs(Matrix_list[i]);
        
    }
    //     A_combo = A_combo.right();
//     cout << " A_combo  = " << A_combo  << endl;
    
    interval out = euclNorm(sliced_norms);
    
    cout << " Slice estimate  = " << out << endl;
//     cout << " Combo estimate  = " << euclNorm(A_combo) << endl;
    return out;
}

IHessian compressTensor( IHessian DDDG , int dimension){    
    IHessian DDDG_out(dimension/2,dimension/2);
    
    for (int i = 0 ; i  < dimension/2;i++){
        for (int j = 0 ; j  < dimension/2;j++){
            for (int k = 0 ; k  < dimension/2;k++){
                DDDG_out(i,j,k) = DDDG(i+dimension/2,j,k);
            }
        }
    }

    return DDDG_out;
}
