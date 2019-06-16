# UNSFlow2D
Unstrucured 2d compressible Cell Centered FVM flow solver

Current feature :
	Solver: Euler
        Flux calculation         : roe,vanleer,ausm+ -up,rusanov
        Flux limiter             : Venkatakrishnan    
        Gradient Reconstruction  : Green-Gauss,Least Square,Green-Gauss Diamond Path
        Time Mode                : LUSGS, RK3
        Convergance Accelaration : central implicit residual smoothening
	
Coming soon    :
	Solver               : Laminar, RANS
        Parallelisation      : MPI


Code still in testing and development stage....
