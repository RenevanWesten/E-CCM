

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION T0, Tt, Tts, Tn, Ts, Td, S0, St, Sts, Sn, Ss, Sd, rho0, rhot, rhots, rhon, rhos, rhod, alpha, beta
      DOUBLE PRECISION V0, Vt, Vts, Vn, Vs, Vd, VAt, SAt, SAts, SAn, SAs, dt, dts, dn, ds, D, c, qna, qnb, qnmin
      DOUBLE PRECISION d0, LxA, LxS , Ly, A, Ar, t_dim, S_dim, D_dim, lambdaA, Tta, Ttsa, Tna, Tsa
      DOUBLE PRECISION qe, qs, qu, qn, qEk, QN1, QS1, QS2, rn, rs, eta, kappa, Agm, tau, Ea, Es, fs, fn
   

        ! Start with defining the state variables
        St = U(1)         		! Salinity thermocline box
        Sts = U(2)         		! Salinity ts box
        Sn = U(3)         		! Salinity northern box
        Ss = U(4)         		! Salinity southern box
        D = U(5)         		! Thermocline depth
        
        Tt = U(6)				! Temperature thermocline box
        Tts = U(7)				! Temperature ts box
        Tn = U(8)				! Temperature northern box
        Ts = U(9)				! Temperature southern box
        Td = U(10)				! Temperature deep box
        
        t_dim = 1*365*86400	! Dimensionalize time
        D_dim = 1000			! Dimensionalize thermocline depth
        S_dim = 35				! Dimensionalize salinity

        ! Parameters physical model (not from Cimatoribus: dn, ds, Tt, Ts)
        V0 = 3d17       		! Total volume of the basin [m3]
        Vn = 3d15       		! Volume northern box [m3]
        Vs = 9d15       		! Volume southern box [m3]
        A = 1d14        		! Horizontal area of Atlantic pycnocline [m2]
        LxA = 1d7       		! Zonal exten of the Atlantic at its southern end [m]
        LxS = 3d7       		! Zonal extent of the Southern Ocean [m]
        Ly = 1d6        		! Meridional extent of frontal region Southern Ocean [m]
        Ar = 4d14       		! Area of thermocline + box ts [m2]
        Vts = LxA*Ly*D*D_dim/2  ! Volume ts box [m3]
        Vt = A*D*D_dim          ! Volume thermocline box

        Vd = V0-Vn-Vs-Vts-Vt    ! Volume deep box [m3]
        
        SAn = 1d13      		! Surface area northern box [m2]
        SAs = 3d13      		! Surfacer area southern box [m2]
        SAt = 1d14      		! Surdace area thermocline box [m2]
        SAts = 1d13   			! Surface area ts box [m2]
   
        tau = 0.1       		! Average zonal wind stress amplitude [Nm-2]
        Agm = 1700      		! Eddy diffusivity [m2s-1]
        kappa = 1d-5    		! Vertical diffusivity [m2s-1]
        eta = 3d4       		! Hydraulic constant [ms-1]
        fs = -1d-4      		! Coriolis paramter [s-1]
        Es = par(2)*1d6     	! Symmetric freshwater flux [m3s-1]
        rn = 5d6        		! Transport by northern subtropical gyre [m3s-1]
        rs = 1d7        		! Transport by southern subtropical gyre
        Ea = par(1)     		! Asymmetric freshater flux [m3s-1]

        rho0 = 1027.5   		! Reference density [kgm-3]
        S0 = 35         		! Reference salinity [psu]
        T0 = 5          		! Reference temperature [degree C]
        alpha = 2d-4    		! Thermal expansion coefficient [K-1]
        beta = 8d-4     		! Haline contraction coefficient [psu-1]
        
        lambdaA = PAR(3)*(1d-6)	! Heat exchange rate (value stated below) [ms-1]

		! Air temperatures: change for different values of lambdaA
        Tta = 16.67740676229534813 		! Air temperature above thermocline box
        Ttsa = Tta 						! Air temperature above ts box
        Tna = -2.096541343133797142		! Air temperature above northern box
        Tsa = Tna						! Air temperature above southern box

        qnmin = PAR(5)			! Minimum reduction factor for qn [-]; value stated below
		
		! Determine sea-ice fraction
		IF (Tn.lt.0.0) THEN
            fn = 100.
        ELSE IF (Tn.gt.5.0) THEN
            fn = 0.
        ELSE
             fn = 100.*(1.-0.2*Tn)
        ENDIF 
        
        ! Calculate volume fluxes
        qEk = tau*Lxs/(rho0*ABS(fs))
        qe = Agm*D*D_dim*LxA/Ly
        qs = qEk - qe
        qu = kappa*A/(D*D_dim)
        
        ! Calculate AMOC strength (PAR(4) == 1: including sea ice)
        qnb = (1-PAR(4))+PAR(4)*((1.-qnmin)/2*TANH(0.1*(50-fn))+(1.+qnmin)/2)
        qna = eta*(D*D_dim*D*D_dim)*(alpha*(Tts-Tn)+beta*(Sn*S_dim-Sts*S_dim))
		qn = qna*qnb
		
        ! Conservation laws for salinity to determine salinity deep box
        Sd= ((S0*V0-(St*Vt+Sts*Vts+Sn*Vn+Ss*Vs)*S_dim)/Vd)/S_dim 

        ! Heaviside functions for qn and qs
        IF (qn.lt.0) THEN
            QN1 = 0
        ELSE
            QN1 = qn
        ENDIF

        IF (qs.gt.0) THEN
            QS2 = 0
            QS1 = qs
        ELSE
            QS2 = qs
            QS1 = 0
        ENDIF
		
		c = (Tt*(qu+qEk-qe-QN1)*t_dim/(Ar*D_dim))/D  ! (Tt x dVt/dt) to be included in dTt/dt (ODE 6)
				
        ! State the ODEs
		F(1) = ((QS1*Sts+QS2*St)+qu*Sd-QN1*St+rs*(Sts-St)+rn*(Sn-St)+2*Es)*t_dim/Vt-(St*(qu+qEk-qe-QN1)*t_dim/(Ar*D_dim))/D      	! ODE Salinity thermocline box
        F(2) =  (qEk*Ss-qe*Sts-(QS1*Sts+QS2*St)+rs*(St-Sts))*t_dim/Vts-(Sts*(qu+qEk-qe-QN1)*t_dim/(Ar*D_dim))/D                    	! ODE Salinity ts box
        F(3) = ((QN1+rn)*(St-Sn)-(Es+Ea*1e6))*t_dim/Vn            																	! ODE Salinity northern box
        F(4) = ((QS1*Sd+QS2*Ss)+qe*Sts-qEk*Ss-(Es-Ea*1e6))*t_dim/Vs  																! ODE Salinity southern box
        F(5) = (qu+qEk-qe-QN1)*t_dim/(Ar*D_dim)                  																	! ODE Thermocline depth
        
        F(6) = ((QS1*Tts+QS2*Tt)+qu*Td-QN1*Tt+rs*(Tts-Tt)+rn*(Tn-Tt)-lambdaA*SAt*(Tt-Tta))*t_dim/Vt-c 								! ODE Temperature thermocline box
        F(7) = (qEk*Ts-qe*Tts-(QS1*Tts+QS2*Tt)+rs*(Tt-Tts)-lambdaA*SAts*(Tts-Ttsa))*t_dim/Vts-(Tts*(qu+qs-QN1)*t_dim/(Ar*D_dim))/D	! ODE Temperature ts box
        F(8) = ((QN1+rn)*(Tt-Tn)-lambdaA*SAn*(Tn-Tna))*t_dim/Vn																		! ODE Temperature northern box
        F(9) = ((QS1*Td+QS2*Ts)+qe*Tts-qEk*Ts-lambdaA*SAs*(Ts-Tsa))*t_dim/Vs 														! ODE Temperature southern box
        F(10) = (QN1*Tn-qu*Td-((QS1*Td+QS2*Ts)))*t_dim/Vd																						! ODE Temperature deep box
        
      END SUBROUTINE FUNC
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER I 
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      

! Parameters
        PAR(1) = 0.					! Ea [Sv]
        PAR(2) = 0.17				! Es [Sv]
        PAR(3) = 3.5				! lambdaA x 10^6 [ms-1] 
        PAR(4) = 1.					! 0: sea ice excluded; 1: sea ice included [-]
        PAR(5) = 0.02				! Minimum reduction factor qn [-]
 
! Used initial conditions (specifics not relevant) 
        U(1) = 1.0058376418610575   ! Salinity thermocline box
        U(2) = 0.994840665484994    ! Salinity ts box
        U(3) = 0.9984467757381668   ! Salinity northern box
        U(4) = 0.9910742010761921   ! Salinity southern box
        U(5) = 0.7382511281434748   ! Thermocline depth box
        
        U(6) = 10					! Temperature thermocline box
        U(7) = 10					! Temperature ts box
        U(8) = 5					! Temperature northern box
        U(9) = 5					! Temperature southern box
        U(10) = 5					! Temperature deep box

      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND
      END SUBROUTINE BCND

    SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !---------- ----

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
    DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
    

    END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
