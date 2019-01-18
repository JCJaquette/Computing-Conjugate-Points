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
  clock_t begin = clock();
  
  int order = 20;
  
  int manifold_subdivision = 8;
  int shots = 5;
  int multiple_newton_steps = 20;
  int single_newton_steps = 20;
 
  int grid = 10; // count zeros
  int stepsize = 4;
  interval L_plus = 8;
  
  bool CHECK_MANIFOLD 		= 0;
  bool CHECK_CONNECTING_ORBIT 	= 0;

    vector <IMap> functions = constructFunctions(  dimension,All_parameters);
  vector <IFunction> energy_vec   = constructEnergy(dimension,  All_parameters);
  IMap f             = functions[0];
  IMap f_minus       = functions[1];
  IMap f_linearize   = functions[2];
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
//   We create the point
  IVector p_u=toInterval(fixedPoint(UNSTABLE,dimension));
//    We Create the approximate linearization
  IMatrix A_u = toInterval(coordinateChange(UNSTABLE,dimension,All_parameters));
//     We create the local vector field object  
  localVField F_u(f,A_u,p_u);  
//   We create the local manifold object
  localManifold localUnstable(F_u,U_flat, L, UNSTABLE,manifold_subdivision);
  
//   cout << "Dfg = " << f[p_u] << endl;
//   cout << "Df  = " << F_u[p_u] << endl;

  
  bool conditions_U = localUnstable.checkConditions( U_flat );
  
  
  
  cout << endl;
  

  //   Create local Stable Manifold
  int STABLE = 1;
  //   We create the point
  IVector p_s=toInterval(fixedPoint(STABLE,dimension));
  
  
//    We Create the approximate linearization
  IMatrix A_s = toInterval(coordinateChange(STABLE,dimension,All_parameters));
//     We create the local vector field object  
  localVField F_s(f,A_s,p_s);  
//   We create the local manifold object
  localManifold localStable(F_s,U_flat, L, STABLE,manifold_subdivision);  
  
  
  bool conditions_S = localStable.checkConditions( U_flat ); 

  
  if(! (conditions_S&&conditions_U))
  {
    cout << "Failed to validate Conditions " << endl;
    if (CHECK_MANIFOLD )
      return -1;
  }
  
  localStable.constructDW();
  localUnstable.constructDW();
  
//   clock_t end_1 = clock();
//   double elapsed_secs_1 = double(end_1 - begin) / CLOCKS_PER_SEC;
//   
//   cout << "Runtime = " << elapsed_secs_1  << endl;
  
  
  
  
//   We Get our initial condition
//   TODO This thing below needs to be made into a function; consider moving to BVP classs??
  
  IVector XY_pt(dimension);  
  interval T ;
  if (dimension ==4 )
    T = 16.3; // n=2
  else 
    T = 20; // n=3
    
    
  IVector guess_U = initialGuessGlobal(dimension, All_parameters, T, 1);
  IVector guess_S = initialGuessGlobal(dimension, All_parameters, T, 0);
  
  IVector local_guess_U = localUnstable.projectPoint(  guess_U);
  IVector local_guess_S = localStable.projectPoint(  guess_S); 
  
  
    interval radius_x = sqr(localUnstable.getRadius())*dimension/2;
    interval radius_y = sqr(localStable.getRadius())*dimension/2;
  interval U_radius_sqr = 0;  // TODO Turn this into a function 
  interval S_radius_sqr = 0;  // TODO Turn this into a function 
  for(int i = 0 ; i<dimension/2;i++){U_radius_sqr += sqr(local_guess_U[i]);}
  for(int i = 0 ; i<dimension/2;i++){S_radius_sqr += sqr(local_guess_S[i]);}
//   cout << "U_radius" << sqrt(U_radius_sqr) << endl;
//   cout << "S_radius" << sqrt(S_radius_sqr) << endl;
  local_guess_U = local_guess_U *sqrt(radius_x/U_radius_sqr);
  local_guess_S = local_guess_S *sqrt(radius_y/S_radius_sqr);
  
    
  for (int i = 0 ; i< dimension; i++)
  {
    if (i<dimension/2)
      XY_pt[i]=local_guess_U[i];
    else
      XY_pt[i]=local_guess_S[i-dimension/2];
  }
  cout << " Old XY_pt    = " << XY_pt << endl;  
  
  
cout << endl;
  
  boundaryValueProblem BVP(f,f_minus,energy_projection,localStable,localUnstable,order ,shots );//TODO REMOVE FROZEN
 
  
//   T = BVP.FindTime(XY_pt,T);
// T = T *.95;
//   END we construct the manifolds
  
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
//   BEGIN Testing integration of middle points
  
    //   We get the initial mid-points for multiple shooting
    vector < IVector > multiple_guess ;
    if (shots >0)
        multiple_guess = multipleShootingGuess(shots,T,dimension,All_parameters);
    
    // We output the energies of our guesses
    for (int i = 0 ; i < shots ; i++)
    {
//         cout << " Energy of guess " << i << " = " << energy(multiple_guess[i]) << endl;
    }
    
    // We go to the zero level set via the gradient 
    for (int i = 0 ; i < shots ; i++)
    {
        for (int j = 0;j<15;j++)
        {
            interval local_energy = energy(multiple_guess[i]);
            IVector local_grad = energy.gradient(multiple_guess[i]);
            multiple_guess[i] -= midVector((local_energy /(local_grad*local_grad))*local_grad);
            multiple_guess[i]=midVector(multiple_guess[i]);
        }
//         cout << " diff " << i << " = " << (local_energy /(local_grad*local_grad))*local_grad << endl;
//         cout << " new Energy " << i << " = " << energy(multiple_guess[i]) << endl;
    }
    
    //   We throw away the last coordinate of each of our guesses  
    for (int i = 0 ; i < shots ; i++)
    {
        IVector new_vector(dimension-1);
        for(int j = 0 ; j < dimension-1;j++)
        new_vector[j]= multiple_guess[i][j];
        multiple_guess[i]=new_vector;
    }
    


//   We make our input points and nbd 
  vector <IVector> points(shots+2);
  vector <IVector> neighborhoods(shots+2);
  
  for (int i =1;i<shots+1;i++)
  {
      points[i]         =multiple_guess[i-1];
//       neighborhoods[i]  = coord_nbd;
      neighborhoods[i]  = IVector(dimension-1);
  }
  points[0]=local_guess_U;
  points.back()=local_guess_S;
  
  neighborhoods[0]    =0*U_flat;
  neighborhoods.back()  =0*U_flat;
  
//   cout << "points = " << points << endl;
  
//   for (unsigned i = 0 ; i < points.size();i++)
//   {
//     cout << "Region["<<i<<"] = " << points[i] << endl;
//   }
  
  interval integration_time = 2*T/(shots+1);
  
//    We do the newton method
  vector < IVector > regions;
  for (int i = 0 ; i<multiple_newton_steps ;i++)
  {
//       cout << " ### = " << i << endl; 
      regions = BVP.NewtonStep(points, neighborhoods ,integration_time) ;
      integration_time = integration_time.mid();
      for (unsigned j=0;j<regions.size();j++)
          points[j] = midVector(regions[j]);
      
//       cout << "New Guess " << endl;
//       for (unsigned i = 0 ; i < regions.size();i++)
//     {
//         cout << " Region["<<i<<"] = " << regions[i] << endl;
//     }
  }
  

  
    if (multiple_newton_steps>0)
    {
        cout << endl << "Final Guess " << endl;
        
        for (unsigned i = 0 ; i < regions.size();i++)
        {
            cout << " Region["<<i<<"] = " << regions[i] << endl;
        }
        
        cout << endl;
        
            for (unsigned i = 0 ; i < regions.size();i++)
        {
            cout << " Region["<<i<<"] width = " << regions[i] - midVector(regions[i]) << endl;
        }
        cout << "done testing " << endl;
    }
    
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
  
  
  
//   return;
  
  BVP.setMiddlePoints(0);
//   BEGIN  do the single-shooting newton method & globalize manifold

  
  IVector XY_nbd_ZERO(dimension);
  IVector Newton_out ;
  for (int i = 0 ; i<single_newton_steps ;i++)
  {
    
    Newton_out =BVP.NewtonStep(XY_pt,XY_nbd_ZERO,T);
//     cout << "Answ XY_pt 	= " << XY_pt << endl;
    XY_pt = midVector(Newton_out);
//     cout << "XY_pt " << XY_pt << endl;
   

  }
  
  
  

//   cout << "Answ XY_pt 	= " << XY_pt << endl;
  
//    We inflate our neighborhood and try to verify
    IVector XY_nbd(dimension);
    
  for (int i = 0 ; i< dimension;i++)
  {
      
//     if ( i != frozen)
      XY_nbd[i] = scale * initial_box; //TODO REMOVE FROZEN
  }
 

//     cout << "XY_nbd = " << XY_nbd << endl;
//   cout << "Y_nbd = " << X_nbd << endl;

  Newton_out =BVP.NewtonStep(XY_pt,XY_nbd,T);
  cout << " Answ XY pt   = " << XY_pt<< endl;
  cout << " Answ XY nbd  = " << XY_nbd<< endl;

//   TODO We need to make sure that our initial conditions are inside the verified manifold.
  if (BVP.checkProof())
      cout << "We are verified!!!!! :) " << endl <<endl;
  else
  {
    cout << "Cannot verify connecting orbit" << endl;
    if (CHECK_CONNECTING_ORBIT )
      return -2;
  }
    
    
  

//    -- We calculate the norm bounds 
// IVector XY = XY_pt + XY_nbd; // TODO Make pt + nbd version 
//   IVector norm_bounds = BVP.NormBound(XY,T);
//   cout << "Norm Bounds = " << norm_bounds << endl;
  
  
//   return;
  
  
  cout << endl << "Globalizing Manifold ... " << endl;
  
  propagateManifold E_u(f_linearize, localUnstable,localStable, XY_pt,XY_nbd,order,stepsize);
      
  
  int unstable_e_values = E_u.frameDet(T+L_plus,grid);
    
  
  
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




int plotDemo()
{
    ofstream file;
    file.open("plot.txt");
    file.precision(16);
    
    cout.precision(12);
    
    // This is the vector field for the Planar Restricted Circular 3 Body Problem
    IMap vectorField("par:mu,mj;var:x,y,dx,dy;fun:dx,dy,x-mj*(x+mu)*sqrt((x+mu)^2+y^2)^(-3)-mu*(x-mj)*sqrt((x-mj)^2+y^2)^(-3)+2*dy,y*(1-mj*sqrt((x+mu)^2+y^2)^(-3)-mu*sqrt((x-mj)^2+y^2)^(-3))-2*dx;");
    // mass ratio
    interval mu=interval(123.0)/interval(10000.0);
    interval mj=1.0-mu;
    // set parameter values
    vectorField.setParameter("mu",mu);
    vectorField.setParameter("mj",mj);
    // The solver uses high order enclosure method to verify the existence of the solution. 
    // The order will be set to 20.
    IOdeSolver solver(vectorField,20);
    ITimeMap timeMap(solver);
    
    
    
    
    
    // This is our initial condition
    IVector x(4);
    x[0]=0.9928634178;
    x[1]=0.0;
    x[2]=0.0;
    x[3]=2.129213043;
    // define a doubleton representation of the interval vector x
    C0Rect2Set s(x);
    // Here we start to integrate. The time of integration is set to T=5.
    interval T=5;
    
//     cout << "Set 1 = " << s << endl;
    IVector initialInput = s;
    cout << "Set 1  = " << initialInput  << endl;
    IVector output = timeMap(T,s);
    
    IVector setOutput = s;
//     cout << "Output = " << output << endl;
    cout << "Set 2  = " << s.show() << endl;
    
//     vector<IVector> V=getTrajectory(s,T,32,timeMap,solver);
    
    vector<IMatrix> Mat_out =getTotalTrajectory(s,T,32,timeMap,solver);
    
//     int k=V.size();
//     for(int i=0;i<k;i++)
//     {
// 	plot(V[i][4],V[i][0],file);
//     }
//     
    file.close();
    return 0;
}  // END


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
          
	    dimension=6;
        if (dimension ==4)
        {
            Input.push_back(1); // a
            Input.push_back(.95);// b 
            Input.push_back(.5);// c
            test(dimension,Input);
        }
        else if (dimension ==6)
        {            
            Input.push_back(1); // b1
            Input.push_back(.95);// b2 
            Input.push_back(.9);// b3
            
//             Input.push_back(.005);// c12
//             Input.push_back(.025);// c23         2 - unstable

//             Input.push_back(-.015);// c12
//             Input.push_back(.02);// c23          1 - unstable
            
            Input.push_back(-.17);// c12
            Input.push_back(.15);// c23         0 - unstable
            

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
