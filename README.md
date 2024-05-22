# Scattering_Poles_Reconstruction
Numerical reconstruction of scattering poles for the Helmholtz equation in 2 dimensions for Dirichlet/Neumann/Impedance boundary conditions

This code corresponds to the algorithm described in the article "An algorithm for computing scattering poles based on dual characterization to interior eigenvalues", Proceedings A, by Cakoni - Haddar - Zilberberg. 
It consists of 4 Matlab files that rely on the use of Gypsilab : https://github.com/matthieuaussal/gypsilab
as well as auxiliary codes used to find the zeros of the derivative of the Hankel function. 

## --------------------------Main---------------------------------------------   
### main_sc_poles.m :   
	-main programm, produces a array/plot of 1 over the norm of g obtained by  
 	inverting the interior scattering operator. It shows peaks at the values of the scattering poles.  
	-set up the ranges of real and imaginary parts of k (values of k to be tested)  
	-set up of the geometry (shape of the obstacle)  
	-set up of boundary condition (Dirichlet or Impedance)  
	-return graphic (Re k, Im k, R(k))  
	-calls "MatrixF/MatrixF_impedance", "LSM_F"  
  
### MatrixF.m :  
	-return the matrix of the Fourier transform of the interior scattering operator  
 	for a Dirichlet boundary condition by solving the interior problem.  
        -We use Single Layer Potential and finite elements to find the solution to each interior problem   
        -(uses from Gypsilab : femGreenKernel, integral, mesh tools)  
  
### MatrixF_impedance.m :  
	-Definition of the impedance function  
	-return the matrix of the Fourier transform of interior scattering operator   
 	for an impedance boundary condition.  
        -We use Single Layer Potential and finite elements to find the solution to each interior problem   
        -(uses from Gypsilab : femGreenKernel, integral, mesh tools)  
  
### LSM_F.m :  
	-Definition of delta_z (parameter for points z taken at distance delta_z around the domain D)  
	-Inversion of the interior scattering operator matrix using SVD and   
 	compute the norm of the inverse of the interior scattering operator
	-(applied to the incident Fourier coefficients of the 2D fundamental solution of Hemholtz equation)

## -----------------------------Shape_Variation--------------------------  
### shape_variation_kite :   
	code to vary shape parameter and compute pole   
	(must provide a starting point corresponding to a poles)  
   
### shape_variation4 :   
	data corresponding to the run of shape_variation_kite (SavingPoles = poles for a=linspace(0.7,1.2,15))  
  
### kite_shapes :   
	draw the shapes of the kite for various parameter  
	compute the area  
   
### figure_shape_var :  
	make figure from data  
   
### ---------------------------Impedance_Variation-------------------------  
### impedance_variation_kite :  
	code varying impedance and computing the corresponding poles  
	(need a starting point at a pole)  
### impedance_variation_savings :  
	data corresponding to a run of impedance_variation_kite  
### figure_impedance_var :  
	make figure from data  
  
## -----------------------------Cauchy_method_Hankel_functions----------------------------------------    
### Cauchy.m :   
	-compute the Cauchy integral in a given box for a given order  
	-set up : Dirichlet or Impedance (+ impedance function)  
### Cauchy_Newton_Zeros_Hankel.m    
	-compute the poles in case of a circle for Dirichlet/Impedance  
	-calls Cauchy  
## -----------------------------Results-----------------------------------   
-poles_kite_impedance_1e-1k : impedance k/10  
-poles_circle_impedance_1e-1k : impedance k/10, radius 1.3  

