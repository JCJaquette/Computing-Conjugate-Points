#include "ODE_functions.h"
// #include <../../home/jonathan/capd-5.0.59/capdAlg/include/capd/vectalg/vectalgLib.h>
// #include <../../home/jonathan/capd-5.0.59/capdDynSys4/include/capd/map/mapLib.h>

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
//       f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));";  // NOTE THIS WAS WRONG
      
      f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3,  (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(u1-.5)),  (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(u2-.5)-c23*u3*(1-u3)*(u2-.5)), (-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(u3-.5));"; 

      
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
//  Output is of length SHOTS
 
 vector < interval > times(shots);
 interval time_increment = (2*T /(shots+1));
 times[0]=-T+time_increment;
 for(int i =1;i<shots;i++) 
  times[i]=times[i-1]+time_increment;
 
 for(int i = 0 ; i< shots ; i++)
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
 
//  cout << "Times " << times <<endl;
 
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
      
//     g(u1,u2) = u2*(1-u2)*(1-2*u1) // WRONG
//     g(u2,u1) = u1*(1-u1)*(1-2*u2) // WRONG
//     g(u2,u3) = u3*(1-u3)*(1-2*u2) // WRONG
//     g(u3,u2) = u2*(1-u2)*(1-2*u3) // WRONG
      
//     g(u1,u2) = u2*(1-u2)*(u1-.5)
//     g(u2,u1) = u1*(1-u1)*(u2-.5)
//     g(u2,u3) = u3*(1-u3)*(u2-.5)
//     g(u3,u2) = u2*(1-u2)*(u3-.5)      
 
      //       f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*f(u1) - c12*g(u1,u2)), (-b2*f(u2) -c12*g(u2,u1)-c23*g(u2,u3)),(-b3*f(u3)-c23*g(u3,u2));"; 
      
      
      
      f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3,  (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(u1-.5)),  (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(u2-.5)-c23*u3*(1-u3)*(u2-.5)), (-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(u3-.5));"; 
      f_minus = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun:-v1,-v2,-v3, -(-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(u1-.5)), -(-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(u2-.5)-c23*u3*(1-u3)*(u2-.5)),-(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(u3-.5));"; 
      
      
      
      
 
      
//     f       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));";  // NOTE THIS WAS WRONG
//     f_minus        = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3;fun: -v1, -v2, -v3, -(-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(1-2*u1)), -(-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(1-2*u2)-c23*u3*(1-u3)*(1-2*u2)),-(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(1-2*u3));"; // NOTE THIS WAS WRONG
    
    
// //  We Construct the linearized problem
    
  // F(u1) = (3*u1*(1-u1)-1/2)
  // F(u2) = (3*u2*(1-u2)-1/2)
  // F(u3) = (3*u3*(1-u3)-1/2)
  // G(u1) = u1*(1-u1)
  // G(u2) = u2*(1-u2)
  // G(u3) = u3*(1-u3)
    f_linearize       = "par:b1,b2,b3,c12,c23;var:u1,u2,u3,v1,v2,v3,U1,U2,U3,V1,V2,V3;fun: v1, v2, v3, (-b1*u1*(u1-.5)*(1-u1) - c12*u2*(1-u2)*(u1-.5)), (-b2*u2*(u2-.5)*(1-u2) -c12*u1*(1-u1)*(u2-.5)-c23*u3*(1-u3)*(u2-.5)),(-b3*u3*(u3-.5)*(1-u3)-c23*u2*(1-u2)*(u3-.5)),V1,V2,V3,(-b1*(3*u1*(1-u1)-1/2)-c12*u2*(1-u2))*U1-c12*(-2*(u1-1/2)*(u2-1/2))*U2,-c12*(-2*(u1-1/2)*(u2-1/2))*U1+(-b2*(3*u2*(1-u2)-1/2)-c12*u1*(1-u1)-c23*u3*(1-u3))*U2-c23*(-2*(u2-1/2)*(u3-1/2))*U3,-c23*(-2*(u2-1/2)*(u3-1/2))*U2+(-b3*(3*u3*(1-u3)-1/2)-c23*u2*(1-u2))*U3;"; 
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


vector < IVector > getLocalGuess( int dimension, vector < double> All_parameters, interval T ,  localManifold &localUnstable,  localManifold &localStable)
{

  IVector guess_U = initialGuessGlobal(dimension, All_parameters, T, 1);
  IVector guess_S = initialGuessGlobal(dimension, All_parameters, T, 0);
  
  IVector local_guess_U = localUnstable.projectPoint(  guess_U);
  IVector local_guess_S = localStable.projectPoint(  guess_S); 
  
  
  interval radius_x = sqr(localUnstable.getRadius())*dimension/2;
  interval radius_y = sqr(localStable.getRadius())*dimension/2;
  
  interval U_radius_sqr = 0;   
  interval S_radius_sqr = 0;  
  for(int i = 0 ; i<dimension/2;i++){U_radius_sqr += sqr(local_guess_U[i]);}
  for(int i = 0 ; i<dimension/2;i++){S_radius_sqr += sqr(local_guess_S[i]);}
  local_guess_U = local_guess_U *sqrt(radius_x/U_radius_sqr);
  local_guess_S = local_guess_S *sqrt(radius_y/S_radius_sqr);
    
    vector < IVector > output;
    output.push_back(local_guess_U);
    output.push_back(local_guess_S);
    return output;
    
    
}

vector <vector < IVector > >Guess_pt_nbd( int dimension, vector < double> All_parameters, interval T, localManifold &localUnstable,  localManifold &localStable, int shots)
{
    
    //   We get the initial mid-points for multiple shooting
    vector < IVector > multiple_guess ;
    if (shots >0)
        multiple_guess = multipleShootingGuess(shots,T,dimension,All_parameters);
    
    
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
        neighborhoods[i]  = IVector(dimension-1);
    }
  
    vector < IVector > local_guesses = getLocalGuess( dimension, All_parameters, T ,  localUnstable,  localStable);
    points[0]     = local_guesses[0];  // local_guess_U 
    points.back() = local_guesses[1];  // local_guess_S 
    
    neighborhoods[0]      = IVector(dimension/2);
    neighborhoods.back()  = IVector(dimension/2);
    
    vector <vector < IVector > > output;
    
    output.push_back(points);
    output.push_back(neighborhoods);
    
    return output;
    
}
