Read me - AUTO code ECCM:
The model code to be used in AUTO consists of two files:
- c.ECCM_AUTO
- ECCM_AUTO.f90

c.ECCM_AUTO
This file contains parameters necessary for running AUTO.
A full description of the parameters is given in the AUTO-07p manual (Doedel et al., 2007; 2021).

ECCM_AUTO.f90
This file contains the model code.

How AUTO is used:
1. Run AUTO in time integration mode to get run the model to steady state.
	--> a = load('ECCM_AUTO')
	--> r1 = run(a)
	
2. Start the first continuation from the last point x of the time integration.
	--> r2 = run(r1(x),ICP=1,IPS=1,DS=0.1,DSMAX=0.1,NPR=10,NMX=1000,SP={'BP0'},ILP=1)
	
	The parameters change certain parameters in c.ECCM_AUTO to be able to run it correctly in continuation mode.
	If necessary, also EPSS, EPSL and EPSU (accuracy parameters) can be added in this sequence. 
	Again, for a full description of the parameters, see the AUTO-07p manual (Doedel et al., 2007; 2021).

3. If the continuation stops at the saddle node on the off or weak AMOC branch, start a second continuation.
This continuation starts at the last point y (i.e. the saddle node) of the first continuation. 
	--> r3 = run(r2(y),UZSTOP={1:1})
	
Important notes:
- PAR(1) represents Ea and is used as control parameter in this study.
- PAR(2) represents Es and is kept constant in this study.
- When different values for lambdaA (PAR(3)) are used, values for air temperatures Tta and Tna need to be changed manually.
- PAR(4) is used to run the model without sea ice (PAR(4)==0) or with sea ice (PAR(4)==1).
- PAR(5), the minimum reduction factor for qn, has to be changed manually to reproduce figure 4d.
	
References:
Doedel, E. J., R. C. Paffenroth, A. C. Champneys, T. F. Fairgrieve, Y. A. Kuznetsov, B. E. Oldeman,
B. Sandstede, and X. J. Wang, 2007: AUTO-07p: Continuation and Bifurcation Software for
Ordinary Differential Equations.

Doedel, E. J., R. C. Paffenroth, A. C. Champneys, T. F. Fairgrieve, Y. A. Kuznetsov, B. E. Oldeman,
B. Sandstede, and X. J. Wang, 2021: auto-07p. URL https://github.com/auto-07p/auto-07p.