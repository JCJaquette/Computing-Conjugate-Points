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