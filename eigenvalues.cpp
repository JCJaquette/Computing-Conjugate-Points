#include "eigenvalues.h" 
// #include <capdAlg/include/capd/vectalg/vectalgLib.h>

IVector boundEigenvalues(IMatrix B)
{
//     Returns the diagonal of the matrix with bounds from Gershgorin theorem 
    int k=B.numberOfRows();
    IVector eigenvalues(k);
    for(int i = 0;i<k;i++)
    {
      interval sum =0;
      for (int j = 0;j<k;j++)
      {
        if( i!=j) sum += abs(B[i][j]); 
      }
      eigenvalues[i] = B[i][i] + interval(-1,1)*sum;
    }
    return eigenvalues;
}

vector < IVector > boundSingleEigenvector(IMatrix A,const IVector V, interval lambda )
{
//     Returns an approximate eigenvector and validated enclosure 
    
//     Does not update either the approximate eigenvector or eigenvalue. 
    
//     If verification fails, the entire program is aborted  
    
    int dimension = A.numberOfRows();
    
    int Newton_Steps = 5;
    
//     Initial neighborhood enclosing eigenvector
    interval scale = 3e-14;
    
    IVector H_vec(dimension+1);
    for (int i =0;i<dimension+1;i++){ H_vec[i]=scale*interval(-1,1);}
    
    vector < IVector > kraw_out ;
    IVector kraw_image;
    IVector new_Nbd;
    
    bool verify = 1;
    bool verify_local;
    
    for (int it = 0;it<Newton_Steps;it++){
        kraw_out = krawczykEigenvector(A, V, lambda , H_vec );
        
        kraw_image  = kraw_out[0];
        new_Nbd     = kraw_out[1];
        
//         cout << "kraw_image " << kraw_image << endl;
//         cout << "H_vec      " << H_vec << endl;
//         cout << "new_Nbd    " << new_Nbd << endl;
        
        verify = 1;
        for (int i = 0 ; i< dimension;i++)
        {
            verify_local = subsetInterior(new_Nbd[i],H_vec[i]);
            if (verify_local ==0)
                verify=0;
        }
//         cout << " Validated Eigenvector = " << verify << endl;
        
        if (it == Newton_Steps -1){break;}
        
        for (int i = 0 ; i< dimension+1;i++){
            H_vec[i] = new_Nbd[i]*(1+1e-10); // This is meant to ensure a future inclusion ... maybe unnecessary
        }        
    }
    if (verify == 0 ){
        cout << "Could not validate eigenvectors" << endl;
        throw -344;
    }
    
//     Get the vector part.

    IVector V_err(dimension);
    for (int i=0; i<dimension;i++){
        V_err[i] = H_vec[i];
    }
    
    vector < IVector > output;    
    output.push_back(V);
    output.push_back(V_err);    
    
    return output;
}

vector < IVector > krawczykEigenvector(IMatrix A,const IVector V, interval lambda , IVector H_vec )
{    
//     K = z - C F(z) + ( Id - C DF([z]) )( [z] - z )
    
//     OUTPUTS: [ K(image), K(image)- mid ] 
    
//     e.g. outputs: [ z , - C F(z) + ( Id - C DF([z]) )( [z] - z ) ]
    
//     -- does not really update the middle approximation point. 
    
    int dimension = A.numberOfRows();
    
    
//     Define the enlarged neighborhoods
    interval lambda_nbd = lambda + H_vec[dimension];    
    IVector  V_nbd(dimension);
    for (int i =0;i<dimension;i++){ V_nbd[i] = V[i] + H_vec[i];}
    
//     NOTE I think F_out & DF_out could be improved by passing V & lambda as thin intervals 
//          ... but I don't want to fix what ain't broke, 
//              also, this is not the major source of error when computing eigenfunctions.
    IVector F_out  = F_eigenvector(A,  V,  lambda);
    IMatrix DF_out = DF_eigenvector(A, V,  lambda);
    IMatrix DF_outH = DF_eigenvector(A, V_nbd,  lambda_nbd );
    
//     Define augmented vector
    IVector moveit(dimension+1);
    for (int i =0;i<dimension;i++){ moveit[i]=V[i];}
    moveit[dimension] = lambda;
    
//     Get approximate inverse
    IMatrix  ApproxInverse = midMatrix(gaussInverseMatrix(midMatrix(DF_out))); 
//     Define identity matrix
    IMatrix eye(dimension+1,dimension+1);
    for (int i =0;i<dimension+1;i++){ eye[i][i]=1;}

    
    IVector new_Nbd     = - ApproxInverse*F_out + ( eye - ApproxInverse*DF_outH )*H_vec;
    IVector kraw_image  = moveit + new_Nbd ;

    vector < IVector >  output;
    output.push_back(kraw_image);
    output.push_back(new_Nbd);
    
    return output;
}

IVector F_eigenvector(IMatrix A, IVector v, interval lambda )
{
//     Computes in first n-components   A*v-lambda*v 
//     Computes in n+1 component        || pi_1 (v) ||^2 - 1/( 2 * abs(lambda_i) )    cf. symplectic normalization in (2.18), (2.19). 
    
//     NOTE Old :: Computes in n+1 component        |v|^2-local_norm_sq  NOTE This is the old version . no longer used. 
    
    int dimension = A.numberOfRows();
    
    IVector F_out(dimension +1);
    
    IVector base = (A*v-lambda*v);
    for (int i =0;i<dimension;i++){ F_out[i]=base[i];}
    
    interval normalization = v*v;
    
    IVector pi_1_V(dimension/2); 
    for (int i =0;i<dimension/2;i++){
        pi_1_V[i]=v[i];
    }
    
    interval symplectic_normalizer = pi_1_V*pi_1_V - 1/(2*abs(lambda)) ; 
    
//     cout << " ||pi_1_V||^2  -  1/(2*lambda_i)" << symplectic_normalizer<< endl; 
    
    F_out[dimension] = symplectic_normalizer;

    
    return F_out;
    
}

IMatrix DF_eigenvector(IMatrix A, IVector v, interval lambda)
{
//     Computes derivative of "F_eigenvector"
    int dimension = A.numberOfRows();
    
    IMatrix DF_out(dimension +1,dimension +1);
        
    for(int i = 0;i<dimension;i++)
    {
      for (int j = 0;j<dimension;j++)
      {
        if( i==j) 
            DF_out[i][j]=A[i][j]-lambda;
        else
            DF_out[i][j]=A[i][j];
      }
    }
//     for (int i =0;i<dimension;i++){ DF_out[i][dimension]=-v[i];}     // OLD METHOD
//     for (int j =0;j<dimension;j++){ DF_out[dimension][j]=2*v[j];}    // OLD METHOD
    
//     Derivative in lambda component
    for (int i =0;i<dimension;i++){ DF_out[i][dimension]=-v[i];}
//     Derivative for last row
    for (int j =0;j<dimension/2;j++){ DF_out[dimension][j]=2*v[j];}
    
//     d_lambda F_{n+1} = sgn(lambda) /( 2 lambda^2 )
    if (lambda>0)
        DF_out[dimension][dimension] = 1/( 2 *lambda^2 );
    else
        DF_out[dimension][dimension] = -1/( 2 *lambda^2 );

    
//     cout << " DF_out = " << DF_out << endl;
    return DF_out;
    
    
}

vector < IMatrix > boundEigenvectors(IMatrix A, IMatrix Q, IVector Lambda)
{
//     Returns <0> approximate e-vectors, and <1> rigorous enclosure
//         We do not normalize the vectors passed through here;
//     Returns enclosure of each e-vector having norm equal to that of the approximate e-vectors.     
    
//     Currently, the output <0> will be the input Q unchanged. (but somehow, in the act of returning the matrix, the intervals might get inflated?)
//     Assumes Q is a matrix with columns of e-vectors of A. 

//     Uses the normalization condition described in the paper on equations (2.18) and (2.19).
    
    int dimension = A.numberOfRows();   
    
    IVector  V_col(dimension); // Local e-vector
     
    IVector  V_cen(dimension);
    IVector  V_err(dimension);
    
    IMatrix Q_center(dimension,dimension);
    IMatrix Q_error(dimension,dimension);
    
    
    for (int j = 0 ; j< dimension; j++)
    {
        //  Pull the j^th-column 
        interval lambda = Lambda[j];
        V_col = getColumn(Q, dimension, j);
        
        
        //  Feed to local function
        vector < IVector > output  =  boundSingleEigenvector(A, V_col, lambda);
        V_cen = output[0];
        V_err = output[1];
        
        // Define output;
        for(int i =0;i<dimension;i++)
        { 
            Q_center[i][j] = V_cen[i];
            Q_error[i][j] = V_err[i];
        }
        
    }
    
    vector < IMatrix >  output;
    output.push_back(Q_center);
    output.push_back(Q_error);    
    
    return output;
}

