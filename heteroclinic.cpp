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
#include "localVField.h"
#include "localManifold.h"
#include "boundaryValueProblem.h"
#include "propagateManifold.h"
#include "topFrame.h"
#include <ctime>


DVector fixedPoint(int i,int dimension) 
{
 DVector p(dimension);
 if(i==0) return p;
 for (int j = 0;j<dimension/2;j++) 
 {
  p[j]=1; 
 }
 return p;
}

DMatrix coordinateChange(int i,int dimension, vector < double> All_parameters) 
{
  
  DVector p=fixedPoint(i,dimension);
  DMap f;
  if ( dimension == 4)
  {
    f        = "par:a,b,c;var:u1,u2,v1,v2;fun: v1, v2,  -a*u1*(u1-1/2)*(1-u1)-c*u2*(1-u2)*(u1-1/2),       -b*u2*(u2-1/2)*(1-u2)-c*u1*(1-u1)*(u2-1/2);";  
    
    f.setParameter("a",All_parameters[0]); // a=1
    f.setParameter("b",All_parameters[1]); // b=1
    f.setParameter("c",All_parameters[2]); // c=.5
    
  }
  else if (dimension == 6 )
  {
      f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));"; 
      
    f.setParameter("b1", All_parameters[0]); // 
    f.setParameter("b2", All_parameters[1]); // 
    f.setParameter("b3", All_parameters[2]); // 
    f.setParameter("c12",All_parameters[3]); // 
    f.setParameter("c23",All_parameters[4]); //       
  }

  
  DMatrix Df=f[p];
  return coordinateChange(Df);
}


IVector knownSolution( double a, interval x)
{
  IVector solution(2);
  
  double rt_a_2 = sqrt(a/2);
  solution[0] = 1/(1+exp(-rt_a_2 *x)); // u
  solution[1] = rt_a_2/( 2 + exp(rt_a_2 *x) + exp(-rt_a_2 *x)); // v
  
  return solution; 
}

vector < IVector >  multipleShootingGuess( int shots, interval T ,int dimension, vector < double > All_parameters)
{
 vector < IVector > output;
//  Output is of length SHOTS+1
 
 vector < interval > times(shots+1);
 interval time_increment = (2*T /(shots+1));
 times[0]=-T;
 for(int i =1;i<shots+1;i++) 
  times[i]=times[i-1]+time_increment;
 
 for(int i = 0 ; i< shots+1 ; i++)
 {
     IVector local_vector(dimension); 
     for (int j = 0;j<dimension/2;j++)
     {
        IVector local_solution = knownSolution( All_parameters[j], times[i]);
        local_vector[j] = local_solution[0];
        local_vector[j+dimension/2] = local_solution[1];
     }
     output.push_back(local_vector);
 }
 
 cout << "Times " << times <<endl;
 
 return output;
}


IVector initialGuessGlobal(int dimension, vector <double> All_parameters, interval T, bool STABLE)
{
  IVector point(dimension);
  
  interval x;
  if (STABLE)
    x = -T;
  else
    x = T;
   
  for (int i =0;i<dimension/2;i++)
  {
    IVector local_solution = knownSolution(All_parameters[i], x);
    point[i]=local_solution[0];
    point[i+dimension/2]=local_solution[1];
  }
  
  return point;
} 

vector <IFunction> constructEnergy(int dimension,  vector < double > All_parameters)
{
    vector <IFunction> output;
    IFunction energy;
    IFunction energy_projection;
    if(dimension == 4)
    {
      
//       U(u1)    = sqr(u1)*sqr(1-u1)/4
//       U(u2)    = sqr(u2)*sqr(1-u2)/4
//       U(u1,u2) = u1*(1-u1)*u2*(1-u2)/2
        energy             = "par:a,b,c;var:u1,u2,v1,v2;fun:         (sqr(v1)+sqr(v2))/2-a*sqr(u1)*sqr(1-u1)/4-b*sqr(u2)*sqr(1-u2)/4 -c*u1*(1-u1)*u2*(1-u2)/2;";
        energy_projection  = "par:a,b,c;var:u1,u2,v1   ;fun:sqrt(-2*((sqr(v1)        )/2-a*sqr(u1)*sqr(1-u1)/4-b*sqr(u2)*sqr(1-u2)/4 -c*u1*(1-u1)*u2*(1-u2)/2));";
        /*
//         Make this parameter set into a function
        energy.setParameter("a",All_parameters[0]); // a=1
        energy.setParameter("b",All_parameters[1]); // b=1
        energy.setParameter("c",All_parameters[2]); // c= +/- 0.1*/
        
    }
    else if( dimension == 6) 
    {
//       U(u1)    = sqr(u1)*sqr(1-u1)/4
//       U(u2)    = sqr(u2)*sqr(1-u2)/4
//       U(u1,u2) = u1*(1-u1)*u2*(1-u2)/2
        
        energy              = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun:         (sqr(v1)+sqr(v2)+sqr(v3))/2-b1*sqr(u1)*sqr(1-u1)/4-b2*sqr(u2)*sqr(1-u2)/4-b3*sqr(u3)*sqr(1-u3)/4 -c12*u1*(1-u1)*u2*(1-u2)/2 -c23*u2*(1-u2)*u3*(1-u3)/2;";
        energy_projection   = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2   ;fun:sqrt(-2*((sqr(v1)+sqr(v2)        )/2-b1*sqr(u1)*sqr(1-u1)/4-b2*sqr(u2)*sqr(1-u2)/4-b3*sqr(u3)*sqr(1-u3)/4 -c12*u1*(1-u1)*u2*(1-u2)/2 -c23*u2*(1-u2)*u3*(1-u3)/2));";
/*
//         Make this parameter set into a function        
        energy.setParameter("b1", All_parameters[0]); // 
        energy.setParameter("b2", All_parameters[1]); // 
        energy.setParameter("b3", All_parameters[2]); // 
        energy.setParameter("c12",All_parameters[3]); // 
        energy.setParameter("c23",All_parameters[4]); // */
    }
    else
    {
        cout << " Unexpected dimension " << endl;
        abort();
    }
    output.push_back(energy);
    output.push_back(energy_projection);
    
    
    for (int i = 0 ; i < 2;i++)
  {
      if (dimension == 4)
      {
        output[i].setParameter("a",All_parameters[0]); // a=1
        output[i].setParameter("b",All_parameters[1]); // b=1
        output[i].setParameter("c",All_parameters[2]); // c= +/- 0.1
      }
      else if( dimension == 6)
      {
        output[i].setParameter("b1", All_parameters[0]); // 
        output[i].setParameter("b2", All_parameters[1]); // 
        output[i].setParameter("b3", All_parameters[2]); // 
        output[i].setParameter("c12",All_parameters[3]); // 
        output[i].setParameter("c23",All_parameters[4]); // 
      }
  }
    return output;
}


vector < IMap > constructFunctions( int dimension, vector < double > All_parameters)
{
 
  IMap f;
  IMap f_minus;
  IMap f_linearize;
  
  

  if(dimension == 4)
  {
      


    f       = "par:a,b,c;var:u1,u2,v1,v2;fun: v1, v2,  -a*u1*(u1-1/2)*(1-u1)-c*u2*(1-u2)*(u1-1/2),       -b*u2*(u2-1/2)*(1-u2)-c*u1*(1-u1)*(u2-1/2);";  
    f_minus = "par:a,b,c;var:u1,u2,v1,v2;fun:-v1,-v2,-(-a*u1*(u1-1/2)*(1-u1)-c*u2*(1-u2)*(u1-1/2)),    -(-b*u2*(u2-1/2)*(1-u2)-c*u1*(1-u1)*(u2-1/2));";  
    

      
// //  We Construct the linearized problem
  
  // F(u1) = (3*u1*(1-u1)-1/2)
  // F(u2) = (3*u2*(1-u2)-1/2)
  // G(u1) = u1*(1-u1)
  // G(u2) = u2*(1-u2)
  // H(u1,u2) = (-2*(u1-1/2)*(u2-1/2))
//   
    
    f_linearize       = "par:a,b,c;var:u1,u2,v1,v2,U1,U2,V1,V2;fun: v1, v2, -a*u1*(u1-1/2)*(1-u1)-c*u2*(1-u2)*(u1-1/2),    -b*u2*(u2-1/2)*(1-u2)-c*u1*(1-u1)*(u2-1/2),   V1   ,V2   ,(-a*(3*u1*(1-u1)-1/2)-c*u2*(1-u2))*U1 -c*(-2*(u1-1/2)*(u2-1/2))*U2,  -c*(-2*(u1-1/2)*(u2-1/2))*U1 + (-b*(3*u2*(1-u2)-1/2)-c*u1*(1-u1))*U2;";
    

  }
  else if (dimension == 6)
  {
//       f(u1) = u1*(u1-.5)*(1-u1)
//       f(u2) = u2*(u2-.5)*(1-u2)
//       f(u3) = u3*(u3-.5)*(1-u3)      
      
//     g(u1,u2) = u2*(1-u2)*(1-2*u1)
//     g(u2,u1) = u1*(1-u1)*(1-2*u2)
//     g(u2,u3) = u3*(1-u3)*(1-2*u2)
//     g(u3,u2) = u2*(1-u2)*(1-2*u3)      
      
      
 
//       f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*f(u1) - c12*g(u1,u2)), (-b2*f(u2) -c12*g(u2,u1)-c23*g(u2,u3)),(-b3*f(u3)-c23*g(u3,u2));"; 
    f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));"; 
    f_minus        = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: -v1, -v2, -v3, -(-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), -(-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),-(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));"; 
    
    
// //  We Construct the linearized problem
    
  // F(u1) = (3*u1*(1-u1)-1/2)
  // F(u2) = (3*u2*(1-u2)-1/2)
  // F(u3) = (3*u3*(1-u3)-1/2)
  // G(u1) = u1*(1-u1)
  // G(u2) = u2*(1-u2)
  // G(u3) = u3*(1-u3)
    f_linearize       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3,U1,U2,U3,V1,V2,V3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3)),V1,V2,V3,(-b1*(3*u1*(1-u1)-1/2)-c12*u2*(1-u2))*U1-c12*(-2*(u1-1/2)*(u2-1/2))*U2,-c12*(-2*(u1-1/2)*(u2-1/2))*U1+(-b2*(3*u2*(1-u2)-1/2)-c12*u1*(1-u1)-c23*u3*(1-u3))*U2-c23*(-2*(u2-1/2)*(u3-1/2))*U3,-c23*(-2*(u2-1/2)*(u3-1/2))*U2+(-b3*(3*u3*(1-u3)-1/2)-c23*u2*(1-u2))*U3;"; 
  // H(u1,u2) = (-2*(u1-1/2)*(u2-1/2))
  // H(u2,u3) = (-2*(u2-1/2)*(u3-1/2))
    
//     f_linearize       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3,U1,U2,U3,V1,V2,V3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3)),V1,V2,V3,(-b1*F(u1)-c12*G(u2))*U1-c12*H(u1,u2)*U2,-c12*H(u1,u2)*U1+(-b2*F(u2)-c12*G(u1)-c23*G(u3))*U2-c23*H(u2,u3)*U3,-c23*H(u2,u3)*U2+(-b3*F(u3)-c23*G(u2))*U3;"; 
    
    
  }
  else
  {
    cout << " Unexpected dimension " << endl;
    abort();
  }
  vector < IMap > output;
  output.push_back(f);
  output.push_back(f_minus);
  output.push_back(f_linearize);
  
  for (int i = 0 ; i < 3;i++)
  {
      if (dimension == 4)
      {
        output[i].setParameter("a",All_parameters[0]); // a=1
        output[i].setParameter("b",All_parameters[1]); // b=1
        output[i].setParameter("c",All_parameters[2]); // c= +/- 0.1
      }
      else if( dimension == 6)
      {
        output[i].setParameter("b1", All_parameters[0]); // 
        output[i].setParameter("b2", All_parameters[1]); // 
        output[i].setParameter("b3", All_parameters[2]); // 
        output[i].setParameter("c12",All_parameters[3]); // 
        output[i].setParameter("c23",All_parameters[4]); // 
      }
  }
  
  return output;
}

void test(int dimension,vector < double > All_parameters)
{
  
  
  clock_t begin = clock();
  


  
  vector <IMap> functions = constructFunctions(  dimension,All_parameters);
  vector <IFunction> energy_vec   = constructEnergy(dimension,  All_parameters);
  IMap f             = functions[0];
  IMap f_minus       = functions[1];
  IMap f_linearize   = functions[2];
  IFunction energy   = energy_vec[0];  // Not sure if this gets used
  IFunction energy_projection   = energy_vec[1];
  
  
  
  int order = 20;
  int grid = 30; 
  int stepsize = 7;
  interval L_plus = 10;
  
 
  bool CHECK_MANIFOLD 		= 1;
  bool CHECK_CONNECTING_ORBIT 	= 1;

 
    //   Cone angle
  //   We define the box used for the boundary value problem
  interval initial_box;
  interval L; 
  interval scale;
  int shots = 3;
  
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
//       T = 20; // n=2
      
      
      
      initial_box =interval(-.000001,.000001);
  }
      
  
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
  localManifold localUnstable(F_u,U_flat, L, UNSTABLE);
  
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
  localManifold localStable(F_s,U_flat, L, STABLE);  
  
  
  bool conditions_S = localStable.checkConditions( U_flat ); 

  
  if(! (conditions_S&&conditions_U))
  {
    cout << "Failed to validate Conditions " << endl;
    if (CHECK_MANIFOLD )
      return;
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
  int frozen;
  T = 16.3; // n=2
//   T = 20; // n=3
  IVector guess_U = initialGuessGlobal(dimension, All_parameters, T, 1);
  IVector guess_S = initialGuessGlobal(dimension, All_parameters, T, 0);
  
  IVector local_guess_U = localUnstable.projectPoint(  guess_U);
  IVector local_guess_S = localStable.projectPoint(  guess_S); 
  
  
    interval radius = sqr(localUnstable.getRadius())*dimension/2;
  interval U_radius_sqr = 0;  // TODO Turn this into a function 
  for(int i = 0 ; i<dimension/2;i++){U_radius_sqr += sqr(local_guess_U[i]);}
  cout << "U_radius_sqr" << U_radius_sqr << endl;
//   local_guess_U = scale*local_guess_U/getMax(abs(local_guess_U));
  local_guess_U = local_guess_U *sqrt(radius/U_radius_sqr);
  local_guess_S = scale*local_guess_S/getMax(abs(local_guess_S));
  frozen = getMaxIndex( abs(local_guess_U)); // TODO Remove this, and remove frozen from BVP class
    
  for (int i = 0 ; i< dimension; i++)
  {
    if (i<dimension/2)
      XY_pt[i]=local_guess_U[i];
    else
      XY_pt[i]=local_guess_S[i-dimension/2];
  }
  cout << " Old XY_pt    = " << XY_pt << endl;  
  
  
cout << endl;
  
  boundaryValueProblem BVP(f,f_minus,energy_projection,localStable,localUnstable,order ,frozen );//TODO REMOVE FROZEN
 
  
  T = BVP.FindTime(XY_pt,T);

  
//   BEGIN Testing integration of middle points
  
    //   We get the initial mid-points for multiple shooting
    vector < IVector > multiple_guess ;
    if (shots >0)
        multiple_guess = multipleShootingGuess(shots,T,dimension,All_parameters);
    
//     // We output the energies of our guesses
//     for (int i = 0 ; i < shots+1 ; i++)
//     {
//         cout << " Energy of guess " << i << " = " << energy(multiple_guess[i]) << endl;
//     }
    
    // We go to the zero level set via the gradient 
    for (int i = 0 ; i < shots+1 ; i++)
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
    for (int i = 0 ; i < shots +1; i++)
    {
        IVector new_vector(dimension-1);
        for(int j = 0 ; j < dimension-1;j++)
        new_vector[j]= multiple_guess[i][j];
        multiple_guess[i]=new_vector;
    }
    
    
    
    
  
// //   We get a sample neighborhood to integrate along
//   IVector coord_nbd(dimension-1);
//   for(int i = 0;i<dimension-1;i++)
//       coord_nbd[i]=scale*interval(-1,1)/100;


//   We make our input points and nbd 
  vector <IVector> points(shots+3);
  vector <IVector> neighborhoods(shots+3);
  
  for (int i =1;i<shots+2;i++)
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
  
  for (unsigned i = 0 ; i < points.size();i++)
  {
    cout << "Region["<<i<<"] = " << points[i] << endl;
  }
  
//    We do the newton method
  vector < IVector > regions;
  for (int i = 0 ; i<2;i++)
  {
      cout << " ### = " << i << endl; 
      cout << points[0] << endl;
      regions = BVP.NewtonStep(points, neighborhoods ,T) ;
      for (unsigned j=0;j<regions.size();j++)
          points[j] = midVector(regions[j]);
      
      cout << "New Guess " << endl;
      for (unsigned i = 0 ; i < regions.size();i++)
    {
        cout << " Region["<<i<<"] = " << regions[i] << endl;
    }
  }


//   
//   return;
  //   END Testing integration of middle points
  
  
  cout << "done testing " << endl;
  return;
  
//   BEGIN 
   for (unsigned i = 0 ; i < regions.size();i++)
    {
        cout << "Region["<<i<<"] = " << midVector(regions[i]) << endl;
    }
  
  IVector XY_nbd_ZERO(dimension);
  IVector Newton_out ;
  for (int i = 0 ; i<50;i++)
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
      return;
  }
    
    


// 
//    -- We calculate the norm bounds 
//   IVector norm_bounds = BVP.NormBound(XY,T);
//   cout << "Norm Bounds = " << norm_bounds << endl;
  
  
//   return;
  
  
  cout << endl << "Globalizing Manifold ... " << endl;
  
  propagateManifold E_u(f_linearize, localUnstable,localStable, XY_pt,XY_nbd,order,stepsize);
      
  E_u.frameDet(T+L_plus,grid);
    
  
  
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  
  cout <<endl <<  "Runtime = " << elapsed_secs  << endl;
//   END
}


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


int main(int argc, char* argv[])
{
  	cout.precision(6);
	try
	{
// 	  plotDemo();
	  bool Get_Param = 0;
	  int dimension;
	  vector < double > Input;
	  
	  if (!Get_Param) 
	  {
          
	    dimension=4;
        if (dimension ==4)
        {
            Input.push_back(1); // a
            Input.push_back(.95);// b 
            Input.push_back(.05);// c
            test(dimension,Input);
        }
        else if (dimension ==6)
        {            
            Input.push_back(1); // b1
            Input.push_back(.98);// b2 
            Input.push_back(.95);// b3
            Input.push_back(.005);// c12
            Input.push_back(.025);// c23
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
