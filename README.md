# Scattering_Poles_Reconstruction
Numerical reconstruction of scattering poles for the Helmholtz equation in 2 dimensions for Dirichlet/Neumann/Impedance boundary conditions

This code corresponds to the algorithm described in the article "An algorithm for computing scattering poles based on dual characterization to interior eigenvalues", Proceedings A, by Cakoni - Haddar - Zilberberg. 
The main algorithm consists of 4 Matlab files that rely on the use of Gypsilab : https://github.com/matthieuaussal/gypsilab  
We then provide codes showing the evolution of one pole when the shape of the domain is changed (in case of a kite) and when the impedence parameter is changed. 
We also give auxiliary codes used to find the zeros of the Hankel function of the first kind, which correspond to the poles of a Dirichlet problem for a disk.   

## --------------------------Main---------------------------------------------   
### main_sc_poles.m :   
	-main programm, produces a array/plot of 1 over the norm of g obtained by inverting the  
 	interior scattering operator. It shows peaks at the values of the scattering poles.  
	-set up the ranges of real and imaginary parts of k (values of k to be tested)  
	-set up of the geometry (shape of the obstacle)  
	-set up of boundary condition (Dirichlet or Impedance)  
	-return graphic (Re k, Im k, R(k))  
	-calls "MatrixF/MatrixF_impedance", "LSM_F"  
  
### MatrixF.m :  
	-return the matrix of the Fourier transform of the interior scattering operator  
 	for a Dirichlet boundary condition by solving the interior problem   
        -We use Single Layer Potential and finite elements to find the solution to each   
	interior problem   
        -(uses from Gypsilab : femGreenKernel, integral, mesh tools)  
  
### MatrixF_impedance.m :  
	-Definition of the impedance function  
	-return the matrix of the Fourier transform of interior scattering operator   
 	for an impedance boundary condition   
        -We use Single Layer Potential and finite elements to find the solution to each   
	interior problem   
        -(uses from Gypsilab : femGreenKernel, integral, mesh tools)  
  
### LSM_F.m :  
	-Definition of delta_z (parameter for points z taken at distance delta_z   
 	around the domain D)  
	-Inversion of the interior scattering operator matrix using SVD and   
 	compute the norm of the inverse of the interior scattering operator
	(applied to the incident Fourier coefficients of the 2D fundamental solution of   
 	Hemholtz equation phi(.,z) for some points z that are delta_z far away from D)

## -----------------------------Shape_Variation--------------------------  
### shape_variation_kite.m :   
	compute the value of a pole for different shapes of the kite
 	the shape of the kite can be modified using its parametrization  
	(must provide a starting point corresponding to a pole)  
   
### shape_variation4 :   
	an example of data corresponding to the run of shape_variation_kite.m
  	(SavingPoles = poles for a=linspace(0.7,1.2,15))  
  
### kite_shapes :   
	draw the shapes of the kite for various parameter  
	and compute its area  
   
### figure_shape_var.m :  
	makes a figure from the data obtained by running shape_variation_kite.m  
   
### ---------------------------Impedance_Variation-------------------------  
### impedance_variation_kite.m :  
	computes the value of a pole for different impedence parameters  
	(must provide a starting point corresponding to a pole)  
### impedance_variation_savings :  
	an example of data corresponding to the run of impedance_variation_kite  
### figure_impedance_var :  
	make figure from data obtained by impedance_variation_kite.m 
  
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

