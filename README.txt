*********************************************************
README file for Computing-Conjugate-Points Version 1.0

Date: April ## 2021
Code Author: Jonathan Jaquette
Contact: jaquette@math.bu.edu

Summary: C++ code which calculates a standing front in a gradient reaction-diffusion PDE (a connecting orbit in a Hamiltonian ODE) and then computes its stability by counting the number of conjugate points. This code accompanies the paper "Validated spectral stability via conjugate points" by Margaret Beck and Jonathan Jaquette, submitted 2021. This code may be reused according to the MIT copyright License

Dependencies: This code requires the CAPD library to perform interval arithmetic, validated integration, and other functions. 
For more information, see the following website    <<<<<    http://capd.ii.uj.edu.pl/    >>>>>

Getting Started: One should compile the code "heteroclinic.cpp". To that end, included is the ***makefile*** used while developing this code library, which the end user will need to modify accordingly. For example, my CAPD library was located in the folder  *** /home/jonathan/capd-5.0.59/bin/  ***. If your CAPD library is located in a different folder, or you use a different version of the library, you should change this. Furthermore during development, this code was compiled using g++ (Ubuntu 9.3.0-17ubuntu1~20.04). Compatibility with other compilers was not tested. 


*************************
**  heteroclinic.cpp   **
*************************

main
    -- The entire program is called from this function, located at the end of the 'heteroclinic.cpp' file. This function selects the mathematical parameters for the PDE of interest, and then calls **computeFrontAndConjugatePoints** to compute the standing front and compute the conjugates points. 
    The PDE of interest, as described in the accompanying paper, is either 2 or 3 coupled bi-stable equations. 
    The selection of parameters may be: hard coded in the 'main' function; input by the user using the 'RetrieveParameters' function; or a computation may be done for a large number of parameters at once using the 'SampleParameters' function. 
    
computeFrontAndConjugatePoints
    -- The main workhorse of the program. This function attempts to construct a computer assisted proof for 
        ++ The existence of a standing front in the PDE, specified in the ODE_functions.cpp file, and
        ++ A rigorous count of the number of conjugate points, and hence the spectral stability of the front. 
    
RetrieveParameters
    -- Provides a terminal interface for the user to select the parameters for the PDE.
    
SampleParameters
    -- A function which 
    
Construct_ParameterList
    -- Constructs the parameters used in the function 'SampleParameters'


********************************* 
Additional function files 
  
ODE_functions
    -- This file is used to define the specific PDE of interest, the spatial dynamics ODE, and other associated files, such as those used in providing an initial guess for the heteroclinic orbit.  The PDE of interest, as described in the accompanying paper, is either 2 or 3 coupled bi-stable equations.
  
eigenvalues
    -- This file contains functions for validating eigenvalues and eigenvectors (always assumed to be real).

utils
    -- This file contains an assortment of functions. 
  
********************************* 
Class Files


localVField
    -- This class provides a change of coordinates about an equilibrium. NOTE: Throughout, it is assumed that the equilibria are known exactly. If this is not the case, the code could be modified to accommodate this, although this would require going through all the code and adding the extra estimates where necessary. 

localManifold
    -- Contains classes **localManifold** and **localManifold_Eig**. The former encapulates bounds for the stable/unstable manifolds. The latter is a derived class is used in computing and bounding eigenfunctions.  

boundaryValueProblem 
    -- This class solves a boundary value problem between an unstable manifold and a stable manifold. If successful, this provides a computer assisted proof of a heteroclinic orbit in the spatial dynamics, corresponding to a standing front in the PDE.

propagateManifold
    -- After the existence of a standing front has been proven, this class computes frame matrix for the unstable eigenspace of solutions, and checks various conditions needed for computer assisted proof.

topFrame
    -- This class is used to manipulate the frame matrix for the unstable eigenspace of solutions produced by propagateManifold, and provide a rigorous count of the conjugate points. 
 
********************************* 
Miscellaneous Files

plot_det.txt
    --  The zeros of this plot correspond to conjugate points. See variable MAKE_PLOT in **heteroclinic.cpp** function **computeFrontAndConjugatePoints**.
    
plot_parameters.txt
    --  Plot of the # of conjugate points for fronts with a variety of parameters, as produced by the **heteroclinic.cpp** function **SampleParameters**

GNU Plot command
    -- Instructions for graphing the file **plot_det.txt** 
  
Instruction_Bifurcation_Diagram.txt
    -- Instructions for graphing the file **plot_parameters.txt** 


