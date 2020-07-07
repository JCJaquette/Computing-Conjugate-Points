#include "topFrame.h"


void topFrame::initialize( void)
{
  
  
  
  
    int s_length = (*p_traject_list)[1].size();
    if (s_length != series_length)
    {
      cout << "FATAL ERROR: Trajectories do not have the same length!! Inspect step size! " << endl;
      return;
    }
  
  constructTimeSeries();
  constructFrameSeries();
  constructDetSeries(); //TODO
  improveDetBound();   //TODO
}


void topFrame::constructTimeSeries( void)
{

//   cout << " dimension    " << dimension << endl;
  for (int j =0;j<series_length;j++)
  {
    time_series.push_back((*p_traject_list)[0][j][0][0]);
  }

  for (int j = 0 ; j < series_length;j++)
  {
//     cout<< "Time["<<j<<"] = " << time_series[j] << endl;
  }
  
}

void topFrame::constructFrameSeries( void)
{
//   frame_series
// 1st level is  time series
// 2nd level is  type:
  //   0 -- Left Endpoint
  //   1 -- Right Endpoint
  //   2 -- Bound on function Values
  //   3 -- Bound on function Derivative
  
//   The matrix then represented is the top frame
  
  
  
  
  for (int i_time =0;i_time<series_length;i_time++)
  {
//     for each time unit
    vector< IMatrix> local_frame_info(0);
    for(int j_type=0;j_type<4;j_type++)
    {
//       j is type 
      
      IMatrix local_matrix(num_trajectories,num_trajectories);
      
      
      for (int k_traj =0; k_traj <num_trajectories;k_traj++)
      {
	for (int m_dim=0; m_dim <num_trajectories;m_dim++)
	{
// 	  For the dimension index - m_dim - we need to skip over the heteroclinic orbit, 
// 	  and only include the top frame. 
	  local_matrix[m_dim][k_traj] = (*p_traject_list)[k_traj][i_time][m_dim + 2*num_trajectories][j_type+1];
	}
      }
      
      
      local_frame_info.push_back(local_matrix);
    }
    frame_series.push_back(local_frame_info);    
  }
  
  
}
  

void topFrame::constructDetSeries( void )
{
  
  
//   vector < vector <interval> >   det_series;
// // 1st level is  time series
// // 2nd level is  type
//   
// //   The Matrix contains in its columns:
// //   0 -- Left Endpoint
// //   1 -- Right Endpoint
// //   2 -- Bound on function Values
// //   3 -- Bound on function Derivative   ----- REMOVED
  
//   We calculate the det-value
  for (int i_time =0;i_time<series_length;i_time++)
  {
//     for each time unit
    vector<  interval > local_det_info(0);
    for(int j_type=0;j_type<4;j_type++)
    {
      interval det_quantity;
      if (j_type < 3)
      {
// 	just calculate determinent
	det_quantity = det(frame_series[i_time][j_type]);
      }
//       else
//       {
// // 	Calculate the derivative
// 	det_quantity = calculateDerivative(frame_series[i_time][2], frame_series[i_time][3]);
//       }
      local_det_info.push_back(det_quantity);
    }
    det_series.push_back(local_det_info);
  }
  
  
  
}

interval topFrame::calculateDerivative( const IMatrix &A , const IMatrix &A_prime)
{
//  We use Jacobi's formula:
//   d/dt det A(t) = tr ( adj (A(t) ) * A'(t) )
  
  IMatrix adj_A(num_trajectories,num_trajectories);
  
  adjugate(A,adj_A);
  
  interval output = trace(adj_A*A_prime);
  return output;
  
}

interval topFrame::minorDet( const IMatrix &A ,int i_hat , int j_hat )
{
  interval minor_out;
  
  IMatrix Minor_ij(num_trajectories-1,num_trajectories-1);
  
  for (int i = 0 ; i<num_trajectories-1;i++)
  {
    int i_pull;
    if (i< i_hat)
      i_pull = i;
    else
      i_pull = i +1;
    
    for (int j = 0 ; j<num_trajectories-1;j++)
    {
      int j_pull;
      if (j<j_hat)
	j_pull=j;
      else
	j_pull = j+1;
      
      Minor_ij[i][j] = A[i_pull][j_pull];
    }
  }
  
  minor_out = det(Minor_ij);

  return minor_out;
}


void topFrame::adjugate( const IMatrix &A , IMatrix &matrix_out)
{
  
  for (int i =0;i<num_trajectories;i++)
  {
    for (int j =0;j<num_trajectories;j++)
    {
      int sign ;
      if ((i+j)%2 == 0)
	sign=1;
      else
	sign = -1;
      
      matrix_out[i][j] = sign*minorDet(A,i,j);
    }
  }
  matrix_out = transpose(matrix_out);
  
}

void topFrame::improveDetBound( void)
{
//   If the determinent is bounded away from zero, then we replace the det value by the hull of the endpoints
}


void topFrame::makePlot( void)
{
//   We Plot the derivative of the top frame
  ofstream file;
  file.open("plot_det.txt");
  file.precision(16);

  
  

  for(int i_time=0;i_time<series_length;i_time++)
  {
      plot(time_series[i_time],det_series[i_time][2],file);
      
  }
  
  file.close();
  
}

vector<int> topFrame::countZeros( void)
{
    //   output [ 0 ] =  # of conjugate pts 
    //   output [ 1 ] =  failure count
    vector<int> output;
    int failure_count  =0;
    int zero_count  =0;
    //   We go through each thing in the constructDetSeries
    
    interval derivative;
    interval MVT_deriv;
    interval MVT_bound;
  
    for (int i_time = 0;i_time < series_length;i_time ++)
    {
        //    if left and right endpoints are of the same sign
        if (det_series[i_time][0] * det_series[i_time][1]>0 )
        {
            //   	check that the entire interval is bounded away from zero
            if ((det_series[i_time][2]>0)||(det_series[i_time][2]<0))
                continue;
            else
            {
                //       Alternatively, check if the derivative is bounded away from zero
                derivative = calculateDerivative(frame_series[i_time][2], frame_series[i_time][3]);
                if ((derivative>0)||(derivative<0))
                    continue;
                else
                {
                    // Use the mean-value-theorem
                    MVT_deriv = (time_series[i_time]-time_series[i_time].left() )*derivative;
                    MVT_bound = intervalHull(det_series[i_time][0],det_series[i_time][1]) + MVT_deriv/2;
                    if ((MVT_bound>0)||(MVT_bound<0))
                        continue;
                    else{
                        failure_count ++;
                        cout << "Time is " << time_series[i_time] << endl;
                        cout << "   LEFT        = " << det_series[i_time][0] << "     RIGHT = " << det_series[i_time][1] << endl;
                        cout << "   MVT Bound = " << (time_series[i_time]-time_series[i_time].left())*derivative  << endl;        
                    //         cout << "   Value Bound = " << det_series[i_time][2]  << endl;
                    //         cout << "   Deriv Bound = " << derivative  << endl;
                    }
                }
            }
        }
        //   else if endpoints are of different sign
        else if (det_series[i_time][0] * det_series[i_time][1] < 0 )
        {
        //   	if derivative is bounded away from zero
            derivative = calculateDerivative(frame_series[i_time][2], frame_series[i_time][3]);
            if ((derivative>0)||(derivative<0))
            zero_count++;
            else
            failure_count ++;
        }
        else
        {
            failure_count ++;
            cout << "Cannot determine sign of vertex  --- decrease grid or stepsize" << endl ;
        }
    }
            
    if (failure_count > 0)
    cout << " Failure Count = " <<     failure_count << endl;

    cout << endl << " Unstable Eigenvalues    = " <<     zero_count << endl;
    output.push_back(zero_count);
    output.push_back(failure_count);
    return output;
}

IMatrix topFrame::getLastFrame( void)
{ 
// For the last time recorded, we return the right-endpt num_trajectories + the derivative of \varphi
  
  IMatrix LastFrame(num_trajectories*2,num_trajectories+2);
  
//   Traj / column
  for (int i = 0 ; i < num_trajectories ; i++)
  {
//     cout << "traj[i] = " << (*p_traject_list)[i].back()  << endl;
//     Row / Entry 
    for (int j = 0 ; j < num_trajectories*2 ; j++)
    {
      LastFrame[j][i] = (*p_traject_list)[i].back()[j+num_trajectories*2][2];
    }
  }
  
//   We add \varphi' to the 2nd to last column
  for (int j =0;j< num_trajectories*2;j++)
  {
//       This is the derivative over the entire interval, and could be improved to just take the right one. 
    LastFrame[j][num_trajectories]=(*p_traject_list).back().back()[j][4];
  }
  
  //   We add \varphi to the last column
  for (int j =0;j< num_trajectories*2;j++)
  {
//       This is the derivative over the entire interval, and could be improved to just take the right one. 
    LastFrame[j][num_trajectories+1]=(*p_traject_list).back().back()[j][2];
  }
  
  return LastFrame;
  
}

