#include "eigenvalues.h" 


IVector boundEigenvalues(IMatrix B)
{
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

IVector boundSingleEigenvector(IMatrix A, IVector V, interval lambda)
{
    int dimension = A.numberOfRows();
    
    int Newton_Steps = 5;
    
    interval scale = 3e-14;
    
    IVector H_vec(dimension+1);
    for (int i =0;i<dimension+1;i++){ H_vec[i]=scale*interval(-1,1);}
    
    vector < IVector > kraw_out ;
    IVector kraw_image;
    IVector new_Nbd;
    
    bool verify = 1;
    bool verify_local;
    
    for (int it = 0;it<Newton_Steps;it++){
        kraw_out = krawczykEigenvector(A, V, lambda , H_vec);
        
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
            H_vec[i] = new_Nbd[i]*(1+1e-10);
        }        
    }
    if (verify == 0 ){
        cout << "Could not validate eigenvectors" << endl;
        abort();
    }
    
//     Get the vector part.
    IVector V_out(dimension);
    for (int i=0; i<dimension;i++){
        V_out[i] = kraw_image[i];
    }
    
    return V_out;
}

vector < IVector > krawczykEigenvector(IMatrix A, IVector V, interval lambda , IVector H_vec)
{
    // OUTPUTS: [ K(image), K(image)- mid ]
    
    int dimension = A.numberOfRows();
    
    
//     Define the enlarged neighborhoods
    interval lambda_nbd = lambda + H_vec[dimension];    
    IVector  V_nbd(dimension);
    for (int i =0;i<dimension;i++){ V_nbd[i] = V[i] + H_vec[i];}
    
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

IVector F_eigenvector(IMatrix A, IVector v, interval lambda)
{
    int dimension = A.numberOfRows();
    
    IVector F_out(dimension +1);
    
    IVector base = (A*v-lambda*v);
    for (int i =0;i<dimension;i++){ F_out[i]=base[i];}
    
    interval normalization = (euclNorm(v)^2 );
    normalization = normalization-1;

    F_out[dimension] = normalization;

    
    return F_out;
    
}

IMatrix DF_eigenvector(IMatrix A, IVector v, interval lambda)
{
    int dimension = A.numberOfRows();
    
    IMatrix DF_out(dimension +1,dimension +1);
    
    IVector base = (A*v-lambda*v);

    
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
    
    for (int i =0;i<dimension;i++){ DF_out[i][dimension]=2*v[i];}
    for (int j =0;j<dimension;j++){ DF_out[dimension][j]=-lambda;}
    
    
    return DF_out;
    
    
}

IMatrix boundEigenvectors(IMatrix A, IMatrix Q, IVector Lambda)
{
        
    int dimension = A.numberOfRows();   
    
    IVector  V_col(dimension);
    
    for (int j = 0 ; j< dimension; j++)
    {
        //  Pull the j^th-column 
        interval lambda = Lambda[j];
        V_col = getColumn(Q, dimension, j);
        V_col = V_col /euclNorm(V_col );
        
        
        //  Feed to local function
        V_col =  boundSingleEigenvector(A, V_col, lambda);
        
        // Define output;
        for(int i =0;i<dimension;i++)
        { 
            Q[i][j] = V_col[i];
        }
        
    }
    
    
    
    
    return Q;
}

