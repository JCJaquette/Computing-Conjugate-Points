// Validated Computation of Conjugate Points

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
#include <string>

int computeFrontAndConjugatePoints(int dimension,vector < double > All_parameters, vector < interval > All_parameters_interval) 
{
// This function attempts to construct a computer assisted proof for 
//     -- The existence of a standing front in the PDE, the equation for which is specified in the ODE_functions.cpp file, and
//     -- A rigorous count of the number of conjugate points, and hence the spectral stability of the front. 
//     
//     Note, we pass both a double version and an interval version of the parameters because there is no conversion from interval to double, and the function for computing eigenvalues and eigenvectors is just for double type matrices. 
//     
//     
// INPUT
//      dimension               -- either 4 or 6 -- the dimension of the spatial dynamics, a Hamiltonian ODE.
//      All_parameters          -- type *double* - the center of the interval valued parameters. 
//      All_parameters_interval -- interval enclosure of the parameters 
//     
// OUTPUT
//      [non-negative int ] --   The computer assisted proof IS     successful. The output corresponds to a rigorous count of the number of conjugate points.  
//      [    negative int ] --   The computer assisted proof IS NOT successful. The output corresponds to a particular way the proof failed. The different error codes are not standardized 

//  As described in Section 3 of the paper, this program can be roughly divided into 5 parts: 
//     
//     (i)      Compute a standing wave \varphi 
//     (ii)     Determine a valid L_- 
//     (iii)    Calculate a frame matrix 
//     (iv)     Prove there are no conj pts past L_+ 
//     (v)      Count conjugate points 
//     
//  About half of the program implements part (i), and the other half implements (ii)-(v). 
    
  clock_t begin = clock();

// // // // // // // // // // // // // // // // //   
//   We define our computational parameters     //
// // // // // // // // // // // // // // // // //
  
  bool CHECK_MANIFOLD           = 1;        //  One may choose to override the manifold validation for de-bugging purposes. 
  bool CHECK_CONNECTING_ORBIT   = 1;        //  One may choose to override the connecting orbit validation for de-bugging purposes. 
  
  //  Option to plot the determinent of the top frame. Produces file to create "plot_det.txt". Note: this plots a rough bound of the determinent over each time interval, and does not incorporate extra information from the derivative. 
  bool MAKE_PLOT                = 0;
    
  int order = 20;                           //  High order Taylor method of validated integration
  int manifold_subdivision = 15;            //  Uniform subdivision number for bounding functions when validating the (un)stable manifolds
  int newton_steps = 20;                    //  Max iterates when applying newton's method
  int grid = 14;                            //  When counting conjugate points, each time step the integrator is subdivided for tighter bounds
  int stepsize = 5;                         //  Fixed step size for counting conjugate points, equal to  2^(- stepsize) 
  
  //  When validating the heteroclinic orbit, we integrate on the time interval (-T,T) points on the (un)stable manifolds, meeting in the middle. Selection of T is automatized. 
  //  When computing conjugate points, we integrate on (-L_-,L_+). Nominally we say L_+ is 0.
  //  We choose L_- = L_minus_percent * 2T .
  interval L_minus_percent = 0.84;   
  
  interval scale;                           //  initial size of the (un)stable manifold; later modified based on location of approximate heteroclinic
  interval initial_box(-.000001,.000001);   //  validation nbd for BVP problem; to be multiplied by 'scale'; will get shrunk by Newton's method
  interval L;                               //  Cone angle  --  \vartheta in the paper -- 
  double inflation_ratio = 1.001;           //  After we find an approximate heteroclinic orbit, we redefine our validity region of the manifolds so it just includes the orbit, and then slightly inflate this neighborhood. 
  
//   We adjust some parameters based on the dimension of the system. 
  if (dimension ==4 ){
    scale = 0.000016;     // n=2
    L = interval(.00008); // cone angle           
  }
  else{ 
    scale = 0.000001;     // n=3      
    L = interval(.000015); // cone angle 
  }

  
// // // // // // // // // // // // // // // // //   
//   We begin the computation                   //
// // // // // // // // // // // // // // // // //
  

  //  We create the function objects defining the spatial dynamics, and the linearized systems
  vector <IMap> functions = constructFunctions(  dimension,All_parameters_interval);         // NOTE for a different class of functions, this function will need to be adjusted 
  IMap f             = functions[0]; // For unstable manifold
  IMap f_minus       = functions[1]; // For stable manifold
  IMap f_linearize   = functions[2]; // For non-autonomous system
  
//  We compute an approximate value of T, for determining the initial guess. 
//  This is chosen so that the uncoupled standing front is a distance 'scale' from the equilibria (in global coordinates) at time +/- T . 
  interval T = approxT(dimension, All_parameters,scale);                            // NOTE for a different class of functions, this function may need to be adjusted    
    
    
    //   BEGIN We construct the manifolds  
    //              -- The size of these manifolds will be adjusted after we find the approximate heteroclinic orbit, and verified again. 
  
  IVector U_flat(dimension/2);  
  for (int i =0;i<dimension/2;i++){ U_flat[i]=scale*interval(-1,1);} // The size of the initial manifolds.
  
//   Create local Unstable Manifold
  int UNSTABLE = 0;

  IVector p_u=toInterval(fixedPoint(UNSTABLE,dimension));                           //   We create the equilibrium point
  IMatrix A_u = toInterval(coordinateChange(UNSTABLE,dimension,All_parameters));    //   We create the approximate linearization
  A_u = midMatrix(  symplecticNormalization(A_u,dimension)  );                      //   Impose symplecticNormalization  
  localVField F_u(f,A_u,p_u);                                                       //   We create the local vector field object  
  localManifold localUnstable(F_u,U_flat, L, UNSTABLE,manifold_subdivision);        //   We create the local manifold object
  
  //   Create local Stable Manifold
  int STABLE = 1;
  
  IVector p_s=toInterval(fixedPoint(STABLE,dimension));                             //   We create the equilibrium point
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
vector < IVector > Guess = getLocalGuess( dimension, All_parameters, T, localUnstable,  localStable);       // NOTE for a different class of functions, this function may need to be adjusted 
  
// // // // // // // // // // // // // // // // //   
// // // // // // // // // // // // // // // // //   
//   BEGIN Testing integration of guess 
  
//     We define the XY_pt that will get used in the single-shooting Newton's method. 
    IVector XY_pt(dimension);  
    for (int i = 0 ; i < dimension / 2 ; i++){
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
    cout << " T old = " << T << endl; 
    T = BVP.FindTime(  XY_pt,  T);
    cout << " T new = " << T << endl;     
    
//   END
  
    
//   BEGIN  Perform a single-shooting Newton-like method 
//     Here, we perform a non-rigorous Newton method; we use intervals of zero-width.
  IVector XY_nbd_ZERO(dimension);
  IVector Newton_out ;
  for (int i = 0 ; i<newton_steps ;i++)
  {
    Newton_out =BVP.NewtonStep(XY_pt,XY_nbd_ZERO,T,r_u_sqr);
//     cout << "Answ XY_pt 	= " << XY_pt << endl;
    XY_pt = midVector(Newton_out);
  }
//   END
  

// BEGIN  We update the validity size of the neighborhoods of the manifolds     
//      This is done such that the final 'approximate' solution is just inside the neighborhood 
//      'just inside' here being 0.1% larger than the approximate solution. (see  inflation_ratio defined earlier)
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
      U_flat_U_new[i]=U_flat[i]*x_Ratios[i]*inflation_ratio; 
      U_flat_S_new[i]=U_flat[i]*y_Ratios[i]*inflation_ratio;
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
    for (int i = 0 ; i<newton_steps ;i++){
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
  
  
// We check to see that the boundry value of our heteroclinic orbits lie inside the validated neighborhood of the manifold. 
  IVector X_pt_nbd(dimension/2);
  IVector Y_pt_nbd(dimension/2);
  
  for (int i = 0 ; i< dimension/2 ; i++){
      X_pt_nbd[i] = XY_pt[i]+XY_nbd[i];
      Y_pt_nbd[i] = XY_pt[i+dimension/2]+XY_nbd[i+dimension/2];
  }
  
  bool X_pt_INSIDE = subsetInterior(X_pt_nbd,U_flat_U_new);
  bool Y_pt_INSIDE = subsetInterior(Y_pt_nbd,U_flat_S_new);
  bool INSIDE_BOX = X_pt_INSIDE && Y_pt_INSIDE;
  
  if (CHECK_CONNECTING_ORBIT ){
    if ( !(INSIDE_BOX)  ){
        cout << " Endpoints of heteroclinic lie outside the (un)stable manifolds. " << endl ; 
        return -1;
    }
    if (BVP.checkProof())
        cout << endl <<  "We have verified the connecting orbit!!  " << endl <<endl;
    else
    {
        cout << "Cannot verify connecting orbit" << endl;
        return -2;
    }
  }
  else{
      cout << "Not checking if connecting orbit is verified" << endl;
  }

//      END single shooting    



//  //     (ii)     Determine L_-                   // // // //
//  //     (iii)    Calculate a frame matrix        // // // //


//  BEGIN globalize manifold
  
  cout << endl << "Globalizing Manifold ... " << endl;
  
//   We upcast our manifold objects to allow for eigenfunction calculations
  localManifold_Eig localUnstable_eig(localUnstable);
  localManifold_Eig localStable_eig(  localStable  );
  
  propagateManifold E_u(f_linearize, localUnstable_eig,localStable_eig, XY_pt,XY_nbd,order,stepsize);
  E_u.plotting(MAKE_PLOT);
      
//   Get the endpoint  
  IVector Y_pt(dimension/2); 
  IVector Y_nbd(dimension/2);
  
  for (int i = 0 ; i< dimension/2 ; i++){
      Y_pt[i]   = XY_pt[i+dimension/2];
      Y_nbd[i]  = XY_nbd[i+dimension/2]; 
  }
  
  
// NOTE ON L_- & L_+
// We impose a coordinate system on our heteroclinic orbit, as when solving the BVP, 
//      we formally define the left  point on our unstable manifold to be phi(-L_-) and
//      we formally define the right point on our   stable manifold to be phi(2T-L_-) 
// However we do not need to go all the way to the right to check the 'no more conj pts' condition. 
// For that, we only go from phi(-L_-) on the left to phi(L_+) on the right, in particular having defined L_+ := 0. 
// We obtain the point phi(L_+) by integrating backwards from the right endpoint phi(2T-L_-) 

// We define L_minus 
interval L_minus = sup(2*T*L_minus_percent);
// The amount of time we need to integrate backwards from q_1 to get to L_+ = 0
interval effective_L_plus = 2*T*(1-L_minus_percent);
//     The point: phi(0) = Phi_{ effective_L_plus }(q_1) 
IVector endPoint_LPlus = BVP.Gxy( Y_pt, Y_nbd, effective_L_plus,  STABLE) ;

  
//  //     (iv)     Prove no conj pts past L_+      // // // //
//  //     (v)      Count conjugate points          // // // //
  
  int unstable_e_values = E_u.frameDet(L_minus,grid,endPoint_LPlus-p_s);
//  unstable_e_values -----
//      non-negative int -- validated # of conjugate points 
//          negative int -- the proof failed 
  
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


////////////////////////////////////////////////////////////////
    

vector < string >  RetrieveParameters( int dimension){
//  Provides a terminal interface for the user to select the parameters for the PDE.
//     The parameters are input as strings, and then later cast as doubles / intervals.
  string param_in;
  vector < string > vector_out  ;
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

vector < vector < double > > Construct_ParameterList( int subdivisions_circle,int subdivisions_radius, double radius_min, double radius_max){
//  For the n=3 case, and for fixed values of the 'b' vector, 
//  This function samples different values of 'c = [ c_{12} , c_{23} ]', and returns a list of parameters (as double).
//     
//  We are sampling 2-dimensional vectors from the plane, using polar coordinates 
//      --  We sample **subdivisions_radius** different radii in [ radius_min, radius_max ] 
//      --  We sample  **subdivisions_circle**  distinct angles
    
//     If desired, these values of b_i can be changed 
    vector < double > Input;
    Input.push_back(1);     // b1
    Input.push_back(.975);  // b2 
    Input.push_back(.95);   // b3
    
    
    vector < vector < double > > Parameter_list;
    double radius_step = (radius_max-radius_min)/(subdivisions_radius-1);
    
    for (int i = 0 ; i< subdivisions_circle;i++)
    {
        // We offset the angle, so we don't end up with c_{12} or c_{23} equalling zero. 
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
//      For the n=3 case, this function runs **computeFrontAndConjugatePoints** for a variety of the parameter *c*. 
//      The list of parameters are defined in **Construct_ParameterList**
//      The results of this computation are output to the file "plot_parameters.txt"
    
  clock_t begin = clock();
    
    int subdivisions_circle = 24;
    int subdivisions_radius = 4;
    
    double radius_min = .02;
    double radius_max = .06;
    
    vector < vector < double > > Parameter_list = Construct_ParameterList(subdivisions_circle,subdivisions_radius, radius_min, radius_max);
    
    int no_params = Parameter_list.size();
    vector < int > Unstable_EigenValues(no_params);
    
    int dimension =6;
    for (int i = 0 ; i< no_params; i++)
    {
        cout <<endl<<endl<< "Trial #"<<i+1<< " of " << no_params<<endl;
        cout << "Parameters = " << Parameter_list[i][3] << " " ;
        cout << Parameter_list[i][4] << " " << endl << endl;
        
        //  If, for whatever reason, the function *computeFrontAndConjugatePoints* throws an error, we use the try/catch block so that we can just continue with our parameter list.
        try
        {
            Unstable_EigenValues[i] =  computeFrontAndConjugatePoints(dimension,Parameter_list[i],vector_double2interval(Parameter_list[i]) );
        }
        catch(exception& e){
            cout << "\n\nException caught: "<< e.what() << endl;
            Unstable_EigenValues[i] = -10;
        }
        catch(int error_int)
        {
            Unstable_EigenValues[i] = error_int;
        }
//         We consolidate our error messages, for graphing purposes
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

int main(int argc, char* argv[]){
//  The entire program is called from this function. This function selects the mathematical parameters for the PDE of interest, and then calls **computeFrontAndConjugatePoints** to compute the standing front and compute the conjugates points. 
//  The results of this program are printed onto the terminal screen. 

    cout.precision(16);
    
//  NOTE If you wish to sample many parameters for the n=3 case, uncomment the two lines below
//     SampleParameters(); 
//     return 0;
    
    
//  NOTE If you wish to enter the parameter values in the terminal, set *Get_Param = 1*. Otherwise, set *Get_Param = 1* and enter the parameters below. 
    bool Get_Param = 1;
    
    
    int dimension;
    vector < string > Input_str;        //  The PDE parameters, as a string 
    vector < double > Input;            //  The PDE parameters, as a double
    vector < interval > Input_interval; //  The PDE parameters, as a interval 
    
    if (!Get_Param){  
        dimension = 6; // NOTE This can be either 4 or 6. 
        
        if (dimension ==4)
        {
            Input_str.push_back("1"); // a
            Input_str.push_back(".97");// b 
            Input_str.push_back(".05");// c
        }
        else if (dimension ==6)
        {            
            Input_str.push_back("1");  // b1              1.00
            Input_str.push_back(".98");// b2              0.98
            Input_str.push_back(".96");// b3              0.96
                        
            Input_str.push_back(".04");// c12            +/- .04 
            Input_str.push_back(".02");// c23            +/- .02  
        }
    }
    else{
        //  The user is prompted (in the terminal) to enter the parameters of the PDE.
        cout << "Enter the number of coupled equations, either 2 or 3:" << endl << "      ";
        cin >> dimension;
        if  ( ( dimension != 2) && ( dimension != 3) ) {
            cout << "Invalid number of equations." << endl;
            return 0;
        }
        dimension=dimension*2;
        Input_str = RetrieveParameters( dimension);
    }
    
    //  We convert the string of parameters to a validated interval enclosure, and the double midpoint. 
    Input           = vector_string2double(Input_str);
    Input_interval  = vector_string2interval(Input_str);
    
    
    
    //  We attempt to construct the Computer Assisted Proof for the existence of the standing front and the # of conjugate points. 
    try{
        computeFrontAndConjugatePoints(dimension,Input,Input_interval);
    }
    catch(exception& e){
        cout << "\n\nException caught: "<< e.what() << endl;
    }
    
    
    return 0; // This return line doesn't do anything important
} 
