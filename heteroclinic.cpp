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
    
  bool CHECK_MANIFOLD 		    = 1;
  bool CHECK_CONNECTING_ORBIT 	= 1;
  
  clock_t begin = clock();
  
  int order = 20;               //  High order Taylor method of validated integration
  
  int manifold_subdivision = 15;
  int single_newton_steps = 20;
 
  int grid = 14;                //  Grid for counting conjugate points
  int stepsize = 5;             //  Fixed stepsize for counting conjugate points
  interval L_plus = 13.55;      //  Distance to integrate forward when finding unstable eigenspace.
  
  interval scale;                           //  initial size of the (un)stable manifold;
  interval initial_box(-.000001,.000001);   //  validation nbd for BVP problem; to be multiplied by 'scale'
  interval L;                               //  Cone angle  --  \vartheta in the paper -- 
  
  if (dimension ==4 ){
    scale = 0.000016;     // n=2
    L = interval(.00008); // cone angle
    L_plus = 13;
    stepsize =5;
  }
  else{ 
    scale = 0.000001;     // n=3      
    L = interval(.000015); // cone angle
//     L_plus = 13.55;
    L_plus = 14;
    stepsize =6;
  }
    
  
//  //     (i)      Compute a standing wave \varphi // // // //


  vector <IMap> functions = constructFunctions(  dimension,All_parameters);
  vector <IFunction> energy_vec   = constructEnergy(dimension,  All_parameters);
  IMap f             = functions[0]; // For unstable manifold
  IMap f_minus       = functions[1]; // For stable manifold
  IMap f_linearize   = functions[2]; // For non-autonomous system
  
//  We compute an approximate value of T, for determining the initial guess. 
//  This is chosen so that the uncoupled standing front is a distance 'scale' from the equilibria (in global coordinates) at time +/- T . 
  interval T = approxT(dimension, All_parameters,scale);


    //   BEGIN We construct the manifolds  
    //              -- The size of these manifolds may be adjusted after we find the approximate heteroclinic orbit, and verified again. 
  
  IVector U_flat(dimension/2);  // The size of the initial manifolds.
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
  
  bool conditions_U = localUnstable.checkConditions();
  cout << endl;
  bool conditions_S = localStable.checkConditions(); 

  
  if(! (conditions_S&&conditions_U))
  {
    cout << "Failed to validate stable/unstable manifold conditions " << endl;
    if (CHECK_MANIFOLD )
      return -1;
  }  
//   END we construct the manifolds
 
  
//  We produce a guess in the local eigen-coordinates of the boundary points of our heteroclinic orbit. 
//  This is chosen and rescaled so that the guess will just barely fit inside manifolds <>localUnstable<> and <>localStable<>. 
vector < IVector > Guess = getLocalGuess( dimension, All_parameters, T, localUnstable,  localStable);
    
    
// BEGIN    TESTING This block may be used to find better guesses for obstinant parameters. 

//     // // // // //      Output where the heteroclinic orbit intersects the stable/unstable manifolds. 
//     cout.precision(16);
//     IVector Zero_n(dimension/2);
//     vector < IVector > Unstable_Global_pt_Coord =   localUnstable.getPointNbd( Guess[0], Zero_n);
//     vector < IVector > Stable_Global_pt_Coord   =   localStable.getPointNbd( Guess[1], Zero_n);
//     cout << endl;
//     cout << " Unstable Point = " << Unstable_Global_pt_Coord[0] << endl;
//     cout << " Stable Point = " << Stable_Global_pt_Coord[0] << endl;
//     cout << endl;
// //     cout.precision(6);
// END     // // // // //    
    
  
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
//   BEGIN Testing integration of guess 
  
//     We define the XY_pt that will get used in the single-shooting newton's method. 
    IVector XY_pt(dimension);  
    for (int i = 0 ; i < dimension / 2 ; i++)
    {
        XY_pt[i] = mid(Guess[0][i]);
        XY_pt[i+dimension/2] = mid( Guess.back()[i]);
    }
    cout << " XY pt   = " << XY_pt<< endl;
    
//     When we solve the BVP, we fix the radius away from the equilibrium, of the point on the unstable manifold. 
//     This value we choose for the radius r_u is the initial radius of our guess on the unstable manifold. 
//     To be more precise, we fix the square of the radius. 
//     In any case, this value is fixed and not readjusted. If necessary, the manifold validation nbd is increased to contain the boundary values of the heteroclinic.
    interval r_u_sqr;
    for (int i = 0 ; i < dimension / 2 ; i++){
        r_u_sqr += sqr(XY_pt[i]);
    }
    r_u_sqr=r_u_sqr.left();
    
//      We create the object of the BVP class
    boundaryValueProblem BVP(f,f_minus,localStable,localUnstable,order); 

//  We adjust the integration time to minimize the distance between our initial points. 
//  This time is then fixed for the rest of the program. 
    interval T_new = BVP.FindTime(  XY_pt,  T);
    
    
    cout << " T old = " << T << endl; 
    cout << " T new = " << T_new << endl;     
    T = T_new;
  
  
//   BEGIN  Perform a single-shooting Newton-like method 
  
  
//     Here, we perform a non-rigorous Newton method; we use intervals of zero-width.
  IVector XY_nbd_ZERO(dimension);
  IVector Newton_out ;
  for (int i = 0 ; i<single_newton_steps ;i++)
  {
    Newton_out =BVP.NewtonStep(XY_pt,XY_nbd_ZERO,T,r_u_sqr);
//     cout << "Answ XY_pt 	= " << XY_pt << endl;
    XY_pt = midVector(Newton_out);
  }
  
    
// BEGIN  We update the validity size of the neighborhoods of the manifolds
//      This is done such that the final 'approximate' solution is just inside the neighborhood 
//      'just inside' here being 0.1% larger than the approximate solution. 

  IVector X_pt_final(dimension/2); 
  IVector Y_pt_final(dimension/2);
  for( int i = 0 ; i < dimension/2;i++){
      X_pt_final[i] = XY_pt[i];
      Y_pt_final[i] = XY_pt[i+dimension/2];
  }
  IVector x_Ratios = localUnstable.containmentRatios(X_pt_final);
  IVector y_Ratios = localStable.containmentRatios(Y_pt_final);
  
  IVector U_flat_U_new(dimension/2);
  IVector U_flat_S_new(dimension/2);
  for (int i =0;i<dimension/2;i++){
      U_flat_U_new[i]=U_flat[i]*x_Ratios[i]*1.001;
      U_flat_S_new[i]=U_flat[i]*y_Ratios[i]*1.001;
  }
  localUnstable.updateSize(U_flat_U_new);
  localStable.updateSize(U_flat_S_new);
  
//   After the manifolds have been updated, we re-check the manifold conditions, and then validate the BVP. 
  conditions_U = localUnstable.checkConditions(  );
  cout << endl;
  conditions_S = localStable.checkConditions(  ); 
  
  if ( ! ( conditions_U && conditions_S) ){
      cout << " Cannot verify stable/unstable manifolds after resizing for approximate solution." << endl;
      return -1;
  }
//   END
  
  cout << endl<< "Validating BVP ... " << endl << endl;
  
//   BEGIN We validate the solution to our BVP using the Krawczyk method. 

//    First, we obtain a neighborhood about our approximate solution, which should map into itself under the Krawczyk operator. 
  IVector XY_nbd(dimension);
  for (int i = 0 ; i< dimension;i++)
  {
      XY_nbd[i] = scale * initial_box; 
  }
  
//    We map this neighborhood several times under the Krawczyk operator. This should get a near optimal nbd to be verified. 
    for (int i = 0 ; i<single_newton_steps ;i++){
        Newton_out =BVP.NewtonStep(XY_pt,XY_nbd,T,r_u_sqr);
        XY_pt = midVector(Newton_out);
        XY_nbd = (Newton_out - XY_pt);
//         cout << "XY_pt " << XY_pt << endl;
//         cout << "XY_nbd " << XY_nbd << endl;
    }    
//     We slightly inflate the neighborhood so that it will ensure proper containment. 
    XY_nbd = XY_nbd *interval(1.5);

    //     We perform a final Newton/Krawczyk step and check for proper containment. 
  Newton_out =BVP.NewtonStep(XY_pt,XY_nbd,T,r_u_sqr);
  cout << " Answ XY pt   = " << XY_pt<< endl;
  cout << " Answ XY nbd  = " << XY_nbd<< endl;
  
//   BEGIN We check to see that the boundry value of our heteroclinic orbits lie inside the validated neighborhood of the manifold. 
  IVector X_pt_nbd(dimension/2);
  IVector Y_pt_nbd(dimension/2);
  
  for (int i = 0 ; i< dimension/2 ; i++){
      X_pt_nbd[i] = XY_pt[i]+XY_nbd[i];
      Y_pt_nbd[i] = XY_pt[i+dimension/2]+XY_nbd[i+dimension/2];
  }
  
  bool X_pt_INSIDE = subsetInterior(X_pt_nbd,U_flat_U_new);
  bool Y_pt_INSIDE = subsetInterior(Y_pt_nbd,U_flat_S_new);
  bool INSIDE_BOX = X_pt_INSIDE && Y_pt_INSIDE;
  
  if ( !(INSIDE_BOX)  ){
    cout << " Endpoints of heteroclinic lie outside the (un)stable manifolds. " << endl ; 
    return -1;
  }
//   END

  if (BVP.checkProof())
    cout << endl <<  "We have verified the connecting orbit!!  " << endl <<endl;
  else
  {
    cout << "Cannot verify connecting orbit" << endl;
    if (CHECK_CONNECTING_ORBIT )
      return -2;
  }

//      END single shooting 



// BEGIN This block may be used to find better guesses for obstinant parameters. 
//     // // // // //    NOTE  Output where the heteroclinic orbit intersects the stable/unstable manifolds. 
//     cout.precision(10);

//     Unstable_Global_pt_Coord =   localUnstable.getPointNbd( X_pt_final, Zero_n);
//     Stable_Global_pt_Coord   =   localStable.getPointNbd( Y_pt_final, Zero_n);
//     cout << endl;
//     cout << " Unstable Point = " << Unstable_Global_pt_Coord[0] << endl;
//     cout << " Stable Point = " << Stable_Global_pt_Coord[0] << endl;
//     cout << endl;
//     cout.precision(6);
//  END   // // // // //    


cout << " Done with BVP at line ~330 " << endl;
return -20;





//  //     (ii)     Determine L_-                   // // // //
//  //     (iii)    Calculate a frame matrix        // // // //


//  BEGIN globalize manifold
//    -- We calculate the norm bounds 
//   interval T_nn = 10;
//   IVector XY = XY_pt + XY_nbd; // TODO Make pt + nbd version 
//   IVector component_bounds = BVP.ComponentBound(XY,T_nn);
//   cout << "Component Bounds = " << component_bounds << endl;
  
//   return -4;
  
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
  
    if (!(INSIDE_BOX ))
      cout << endl << " Endpoints of heteroclinic lie outside the (un)stable manifolds :( " << endl ; 
  
    
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
    
    int dimension =4;
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
  	cout.precision(16);
	try
	{
// 	  plotDemo();
	  bool Get_Param = 0;
	  int dimension;
	  vector < double > Input;
	  
	  if (!Get_Param) 
	  {
          
	    dimension=6; // TESTING DIMENSION
	    
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
            Input.push_back(.98);// b2              0.98 previous
            Input.push_back(.97);// b3              0.95 previous
                        
            Input.push_back(.05);// c12                         (.050,.015)                     L_+ = 13.55
            Input.push_back(.05);// c23         0 - unstable   (-.015,.050) & (-.015,-.050)    L_+ = 14
            
            //  Finding computational parameters which get (+.015,-.050) to work is difficult / not yet successful. 
            //   -- For this, it is probably best to perturb up from the n=2 case.   
            
            //  -- NOTE UPDATE! This works for b = ( 1, .98, .97 ) and c = ( -0.05, 0.015 ) .   ~~~~ Not sure about this anymore. 
            
            // 2/2/21 ~ It also looks like a damped Newton's method will help it converge for more difficult parameters. 
            

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
