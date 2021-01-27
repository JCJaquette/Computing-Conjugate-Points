// INTEGRATION OF ODEs
// This is an example of computation of time shift maps for ODEs in CAPD library
// Below we have two applications: 
// 1. non-rigorous computations in "NonRigorous"
// 2. interval-rigorous computations in "IntervalRigorous"

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"
#include "eigenvalues.h"
#include "ODE_functions.h"
#include "localVField.h"
#include "localManifold.h"
#include "boundaryValueProblem.h"
#include "propagateManifold.h"
#include "topFrame.h"
#include <ctime>
#include <cmath>

int test(int dimension,vector < double > All_parameters)
{
//  //  This program has 5 major parts: 
//  //     (i)      Compute a standing wave \varphi // // // //
//  //     (ii)     Determine L_-                   // // // //
//  //     (iii)    Calculate a frame matrix        // // // //
//  //     (iv)     Prove no conj pts past L_+      // // // //
//  //     (v)      Count conjugate points          // // // //
  clock_t begin = clock();
  
  
  
  
  int order = 20;
  
  int manifold_subdivision = 15;
  int shots = 9;
  int multiple_newton_steps = 1;
  int single_newton_steps = 20;
 
  int grid = 14; // count zeros
  int stepsize = 7;
  interval L_plus = 13.55;
  interval T ;
  if (dimension ==4 ){
    T = 16.2; // n=2
    L_plus = 13;
    stepsize =5;
  }
  else{ 
    T = 20.; // n=3
    L_plus = 13.55;
    L_plus = 14;
    stepsize =6;
  }
  
  bool USE_MULTIPLE_SHOOTING = 0;
  
  bool CHECK_MANIFOLD 		    = 1;
  bool CHECK_CONNECTING_ORBIT 	= 1;
  
//  //     (i)      Compute a standing wave \varphi // // // //


  vector <IMap> functions = constructFunctions(  dimension,All_parameters);
  vector <IFunction> energy_vec   = constructEnergy(dimension,  All_parameters);
  IMap f             = functions[0]; // For unstable manifold
  IMap f_minus       = functions[1]; // For stable manifold
  IMap f_linearize   = functions[2]; // For non-autonomous system
  

  
  
  
  IFunction energy   = energy_vec[0];  // Not sure if this gets used
  IFunction energy_projection   = energy_vec[1];

 
  //   Cone angle
  //   We define the box used for the boundary value problem
  interval initial_box;
  interval L; 
  interval scale;
  
  
  if (dimension ==4)
  {
      initial_box =interval(-.000001,.000001);
      L = interval(.00007); // n=2
      scale = 0.00001;
  }
  else if (dimension ==6)
  {
      L = interval(.000015); // n=3
      scale = 0.000001;
      initial_box =interval(-.000001,.000001);
  }
  

    //   BEGIN we construct the manifolds  
  IVector U_flat(dimension/2);
  for (int i =0;i<dimension/2;i++){ U_flat[i]=scale*interval(-1,1);}
  
//   Create local Unstable Manifold
  int UNSTABLE = 0;

  IVector p_u=toInterval(fixedPoint(UNSTABLE,dimension));                           //   We create the point
  IMatrix A_u = toInterval(coordinateChange(UNSTABLE,dimension,All_parameters));    //   We create the approximate linearization
  A_u = midMatrix(  symplecticNormalization(A_u,dimension)  );                      //   Impose symplecticNormalization  
  localVField F_u(f,A_u,p_u);                                                       //   We create the local vector field object  
  localManifold localUnstable(F_u,U_flat, L, UNSTABLE,manifold_subdivision);        //   We create the local manifold object
  
  //   Create local Stable Manifold
  int STABLE = 1;
  
  IVector p_s=toInterval(fixedPoint(STABLE,dimension));                             //   We create the point
  IMatrix A_s = toInterval(coordinateChange(STABLE,dimension,All_parameters));      //   We Create the approximate linearization
  A_s = midMatrix(  symplecticNormalization(A_s,dimension)  );                      //   Impose symplecticNormalization  
  localVField F_s(f,A_s,p_s);                                                       //   We create the local vector field object  
  localManifold localStable(F_s,U_flat, L, STABLE,manifold_subdivision);            //   We create the local manifold object
  
  bool conditions_U = localUnstable.checkConditions( U_flat );
  cout << endl;
  bool conditions_S = localStable.checkConditions( U_flat ); 
  
  if(! (conditions_S&&conditions_U))
  {
    cout << "Failed to validate Conditions " << endl;
    if (CHECK_MANIFOLD )
      return -1;
  }  
//   END we construct the manifolds
 
  
//   TODO Clean up how we get the initial points 
//  BEGIN  We Get our initial condition 
    vector <vector < IVector > > Guess = Guess_pt_nbd( dimension, All_parameters, T, localUnstable,  localStable, shots);
    vector <IVector> points         = Guess[0];
    vector <IVector> neighborhoods  = Guess[1];
    interval integration_time = 2*T/(shots+1);
//   END we construct an initial guess   
    
    
    
    
// //     // // // // //    NOTE  Output where the heteroclinic orbit intersects the stable/unstable manifolds. 
// //     cout.precision(10);
// //     IVector Zero_n(dimension/2);
// //     vector < IVector > Unstable_Global_pt_Coord =   localUnstable.getPointNbd( Guess[0][0], Zero_n);
// //     vector < IVector > Stable_Global_pt_Coord   =   localStable.getPointNbd( Guess[0][shots+1], Zero_n);
// //     cout << endl;
// //     cout << " Unstable Point = " << Unstable_Global_pt_Coord[0] << endl;
// //     cout << " Stable Point = " << Stable_Global_pt_Coord[0] << endl;
// //     cout << endl;
// //     cout.precision(6);
// //     // // // // //    NOTE 
    
    
  
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
//   BEGIN Testing integration of middle points
  
//     We define the XY_pt that will get used in the single-shooting newton's method. 
  IVector XY_pt(dimension);  
  
  boundaryValueProblem BVP(f,f_minus,localStable,localUnstable,order); // Single Shooting   
  bvpMultipleShooting BVP_ms(f,f_minus,energy_projection,localStable,localUnstable,order ,shots ); // Multiple Shooting 
  
  
  
  
    
    IVector XY_pt_T(dimension);  
    for (int i = 0 ; i < dimension / 2 ; i++)
    {
        XY_pt_T[i] = mid(points[0][i]);
        XY_pt_T[i+dimension/2] = mid( points.back()[i]);
    }
    interval T_new = BVP.FindTime(  XY_pt_T,  T);
    
    
    cout << " T old = " << T << endl; 
    cout << " T new = " << T_new << endl; 
    
    
    T = T_new;
    integration_time = 2*T/(shots+1);
    //   Guess = Guess_pt_nbd( dimension, All_parameters, T, localUnstable,  localStable, shots);
    //   points         = Guess[0];
    //   neighborhoods  = Guess[1];
    
    
    if (USE_MULTIPLE_SHOOTING==1){
    
    
    //    We do the newton method
    vector < IVector > regions;
    interval time_nbd =0;
    for (int i = 0 ; i<multiple_newton_steps ;i++)
    {
        integration_time = integration_time.mid();
        cout << " integration_time " << integration_time <<endl;
        regions = BVP_ms.NewtonStep(points, neighborhoods ,integration_time, time_nbd ) ;
        for (unsigned j=0;j<regions.size();j++)
            points[j] = midVector(regions[j]);
    }
    
    interval multiplier = 0.;
    time_nbd = 0 * (integration_time - integration_time.mid());
    
    cout << " Time nbd " << time_nbd <<endl;
    integration_time = integration_time.mid();
    for (unsigned j=0;j<regions.size();j++){
        neighborhoods[j] = multiplier *(regions[j] -  points[j]);
    }

    for (unsigned i = 0 ; i < regions.size();i++)
            {
                cout << " neighborhoods["<<i<<"] = " << neighborhoods[i] << endl;
            }
    cout << " Verifying .... " << endl;      
    regions = BVP_ms.NewtonStep(points, neighborhoods ,integration_time, time_nbd ) ;
    cout << " Time nbd " << integration_time - integration_time.mid() <<endl;
    
        if (multiple_newton_steps>0)
        {
            cout << endl << "Final Guess " << endl;
            
            for (unsigned i = 0 ; i < regions.size();i++)
            {
    //             cout << " Region["<<i<<"] = " << regions[i] << endl;
            }
            
            cout << endl;
            
                for (unsigned i = 0 ; i < regions.size();i++)
            {
                cout << " Region["<<i<<"] width = " << regions[i] - midVector(regions[i]) << endl;
            }
            cout << "done testing " << endl;
        }
    //     return -5;
        
    cout << " Old T = " << T << endl;
    T = (shots+1)*integration_time/2;
    cout << " New T = " << T << endl;
  
  
    
    
    for (int i = 0 ; i < dimension / 2 ; i++)
    {
        XY_pt[i] = mid(regions[0][i]);
        XY_pt[i+dimension/2] = mid( regions.back()[i]);
    }
    
    cout << " XY_pt = " << XY_pt << endl;

  //   END Testing integration of middle points
  
  }// Multiple Shooting 
  else{
      XY_pt = XY_pt_T;
  }
  
  
  
  
  
//   BEGIN  do the single-shooting newton method // TODO Remove this
  
  
  IVector XY_nbd_ZERO(dimension);
  IVector Newton_out ;
  for (int i = 0 ; i<single_newton_steps ;i++)
  {
    
    Newton_out =BVP.NewtonStep(XY_pt,XY_nbd_ZERO,T);
//     cout << "Answ XY_pt 	= " << XY_pt << endl;
    XY_pt = midVector(Newton_out);
//     cout << "XY_pt " << XY_pt << endl;

  }
  
  
//    We inflate our neighborhood and try to verify
    IVector XY_nbd(dimension);
    
  for (int i = 0 ; i< dimension;i++)
  {
      XY_nbd[i] = scale * initial_box; 
  }
  
  
    for (int i = 0 ; i<single_newton_steps ;i++){
        Newton_out =BVP.NewtonStep(XY_pt,XY_nbd,T);
    //     cout << "Answ XY_pt 	= " << XY_pt << endl;
        XY_pt = midVector(Newton_out);
        XY_nbd = (Newton_out - XY_pt)*interval(1.5);
//         cout << "XY_pt " << XY_pt << endl;
//         cout << "XY_nbd " << XY_nbd << endl;
    }

  Newton_out =BVP.NewtonStep(XY_pt,XY_nbd,T);
  cout << " Answ T       = " << T<< endl;
  cout << " Answ XY pt   = " << XY_pt<< endl;
  cout << " Answ XY nbd  = " << XY_nbd<< endl;
  
  
// // // // // // //    NOTE  Output where the heteroclinic orbit intersects the stable/unstable manifolds. 
// //   cout.precision(10);
// // IVector X_pt_final(dimension/2); 
// // IVector Y_pt_final(dimension/2);
// //   for( int i = 0 ; i < dimension/2;i++){
// //       X_pt_final[i] = XY_pt[i];
// //       Y_pt_final[i] = XY_pt[i+dimension/2];
// //   }
// // Unstable_Global_pt_Coord =   localUnstable.getPointNbd( X_pt_final, Zero_n);
// // Stable_Global_pt_Coord   =   localStable.getPointNbd( Y_pt_final, Zero_n);
// // cout << endl;
// // cout << " Unstable Point = " << Unstable_Global_pt_Coord[0] << endl;
// // cout << " Stable Point = " << Stable_Global_pt_Coord[0] << endl;
// // cout << endl;
// // cout.precision(6);
// // // // // // //    NOTE 

//   TODO We need to make sure that our initial conditions are inside the verified manifold.
  if (BVP.checkProof())
      cout << "We are verified!!!!! :) " << endl <<endl;
  else
  {
    cout << "Cannot verify connecting orbit" << endl;
    if (CHECK_CONNECTING_ORBIT )
      return -2;
  }
    

    
//      END single shooting 



cout << " Done with BVP at line ~350 " << endl;
return -20;





//  //     (ii)     Determine L_-                   // // // //
//  //     (iii)    Calculate a frame matrix        // // // //


//  BEGIN globalize manifold
//    -- We calculate the norm bounds 
// IVector XY = XY_pt + XY_nbd; // TODO Make pt + nbd version 
//   IVector norm_bounds = BVP.NormBound(XY,T);
// //   cout << "Norm Bounds = " << norm_bounds << endl;
  
  
  cout << endl << "Globalizing Manifold ... " << endl;
  
  propagateManifold E_u(f_linearize, localUnstable,localStable, XY_pt,XY_nbd,order,stepsize);
      
//   Get the endpoint
  
  
  IVector Y_pt(dimension/2); 
  IVector Y_nbd(dimension/2);
  
  for (int i = 0 ; i< dimension/2 ; i++){
      Y_pt[i]   = XY_pt[i+dimension/2];
      Y_nbd[i]  = XY_nbd[i+dimension/2]; 
  }
  
  
  IVector endPoint_LPlus = BVP.Gxy( Y_pt, Y_nbd, T-L_plus,  STABLE) ;
  
  
  cout << " endPoint_LPlus = " << endPoint_LPlus -p_s << endl;
  cout << " endpoint eig = " << gauss(A_s,endPoint_LPlus -p_s) << endl;
  
  
  
  
  
  
  
  
//  //     (iv)     Prove no conj pts past L_+      // // // //
//  //     (v)      Count conjugate points          // // // //
  
  int unstable_e_values = E_u.frameDet(T,L_plus,grid,endPoint_LPlus-p_s);
    
  
  
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  
  cout <<endl <<  "Runtime = " << elapsed_secs  << endl;
  
  cout << " Result " << unstable_e_values  << endl;
  
  return unstable_e_values;
  
//   END
  
}

// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   


////////////////////////////////////////////////////////////////

vector < double >  RetrieveParameters( int dimension)
{
  double param_in;
  vector < double > vector_out  ;
  cout << "\nEnter Parameters: ";
  for (int i = 1 ; i< dimension/2+1 ; i++)
  {
    cout << "\n b_"<<i<<" = " ;
    cin >> param_in;
    vector_out.push_back(param_in);
  }
  for (int i = 1 ; i< dimension/2 ; i++)
  {
    cout << "\n c_"<<i<<i+1<<" = " ;
    cin >> param_in;
    vector_out.push_back(param_in);
  }

 return vector_out  ;
 
}

vector < vector < double > > Construct_ParameterList( int subdivisions_circle,int subdivisions_radius, double radius_min, double radius_max)
{
    vector < double > Input;
    Input.push_back(1); // b1
    Input.push_back(.98);// b2 
    Input.push_back(.99);// b3
    
    
    vector < vector < double > > Parameter_list;
    double radius_step = (radius_max-radius_min)/(subdivisions_radius-1);
    
    /*vector < vector < double > > Parameter_list(subdivisions_circle*subdivisions_radius)*/;
    for (int i = 0 ; i< subdivisions_circle;i++)
    {
        double theta  = (i+.5)*2*M_PI/subdivisions_circle;
        
        for (int j = 0 ; j< subdivisions_radius; j ++ ) 
        {
            
            double radius = radius_min + radius_step*j;
            
            vector < double > Local_Input = Input;
            
            double c_12 = radius*cos(theta);
            double c_23 = radius*sin(theta);
            
            Local_Input.push_back(c_12);
            Local_Input.push_back(c_23);            
            
            Parameter_list.push_back(Local_Input);
            theta = theta + M_PI/(subdivisions_circle*subdivisions_radius);
            
        }
    }
    
    return Parameter_list;
}

void SampleParameters( void)
{
    
    
  clock_t begin = clock();
    
    int subdivisions_circle = 24;
    int subdivisions_radius = 4;
    
    double radius_min = .005;
    double radius_max = .03;
    
    vector < vector < double > > Parameter_list = Construct_ParameterList(subdivisions_circle,subdivisions_radius, radius_min, radius_max);
    
    int no_params = Parameter_list.size();
    vector < int > Unstable_EigenValues(no_params);
    
    int dimension =6;
    for (int i = 0 ; i< no_params; i++)
    {
        cout <<endl<<endl<< "Trial #"<<i+1<< " of " << no_params<<endl;
        cout << "Parameters = " << Parameter_list[i][3] << " " ;
        cout << Parameter_list[i][4] << " " << endl << endl;
        
        try
        {
            Unstable_EigenValues[i] =  test(dimension,Parameter_list[i]);
            
        }
        catch(exception& e)
        {
    		cout << "\n\nException caught: "<< e.what() << endl;
            Unstable_EigenValues[i] = -5;
        }
//         We consolidate our error messages;
        if (Unstable_EigenValues[i] <0)
            Unstable_EigenValues[i] = -1;
    }
    
    
    ofstream file;
    file.open("plot_parameters.txt");
    file.precision(8);
  
  
    
    for (int i = 0 ; i < no_params; i++)
    {
        file << Parameter_list[i][3] << " " ;
        file << Parameter_list[i][4] << " " ;
        file << Unstable_EigenValues[i] << endl ;
    }
    file.close();
    
    
      
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cout <<endl <<  "Total time for " << no_params << " runs = " << elapsed_secs  << endl;
}

int main(int argc, char* argv[])
{
//     SampleParameters();
//     return 0;
  	cout.precision(6);
	try
	{
// 	  plotDemo();
	  bool Get_Param = 0;
	  int dimension;
	  vector < double > Input;
	  
	  if (!Get_Param) 
	  {
          
	    dimension=4; // TESTING DIMENSION
	    
        if (dimension ==4)
        {
            Input.push_back(1); // a
            Input.push_back(.97);// b 
            Input.push_back(.05);// c
            test(dimension,Input);
        }
        else if (dimension ==6)
        {            
            Input.push_back(1); // b1
            Input.push_back(.97);// b2              0.98 previous
            Input.push_back(.95);// b3              0.95 previous
                        
            Input.push_back(-.05);// c12                         (.050,.015)                     L_+ = 13.55
            Input.push_back(-.015);// c23         0 - unstable   (-.015,.050) & (-.015,-.050)    L_+ = 14
            
            //  Finding computational parameters which get (+.015,-.050) to work is difficult / not yet successful. 
            //   -- For this, it is probably best to perturb up from the n=2 case.            
            

            test(dimension,Input);
        }
	  }
	  else
	  {
	    cout << "no. of equations:  = " ;
	    cin >> dimension;
	    dimension=dimension*2;
	    Input = RetrieveParameters( dimension);
    
	    test(dimension,Input);
	  }
				
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 
