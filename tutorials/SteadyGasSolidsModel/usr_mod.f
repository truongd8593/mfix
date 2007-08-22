

      MODULE usr


        Use param
        Use param1


!
!       Include user-defined variables in this module.  To access the variables
!       from a subroutine add the statement "Use usr".  If allocatable arrays
!       are defined in this module   them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!                        a dummy variable to keep the compiler happy.                     
        DOUBLE PRECISION usr_dummy, DUMMY_DP
! for time-averaging flow variables
      DOUBLE PRECISION TS_counter
! for averaging cross correlations
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  U_s_mean,V_s_mean,U_g_mean,V_g_mean,Ep_g_mean, Th_mean
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Th2 , eps2 , eps_Th , eps_Th2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps2 , eps2_Th2 , Th3 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  eps_Th3 , eps2_Th3 , eps3 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps3 , Th2_eps3 , Th3_eps3 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4 , Th1_eps4 , Th2_eps4 , Th3_eps4 , Th4_eps4 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th4 , eps2_Th4 , eps3_Th4 , Th4 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps5 , Th1_eps5 , Th2_eps5 , Th3_eps5 , Th4_eps5 , Th5_eps5 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th5 , eps2_Th5 , eps3_Th5 , eps4_Th5 , Th5 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps6 ,Th1_eps6 ,Th2_eps6 ,Th3_eps6 ,Th4_eps6 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps6 ,Th6_eps6 ,eps1_Th6 ,eps2_Th6 ,eps3_Th6 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4_Th6 ,eps5_Th6 ,Th6 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps7 ,Th1_eps7 ,Th2_eps7 ,Th3_eps7 ,Th4_eps7 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps7 ,Th6_eps7 ,Th7_eps7 ,eps1_Th7 ,eps2_Th7 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps3_Th7 ,eps4_Th7 ,eps5_Th7 ,eps6_Th7 ,Th7 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::    eps8 ,Th1_eps8 ,Th2_eps8 ,Th3_eps8 ,Th4_eps8 , Th5_eps8 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Th6_eps8 ,Th7_eps8 ,Th8_eps8 ,eps1_Th8 ,eps2_Th8 ,eps3_Th8 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  	eps4_Th8 ,eps5_Th8 ,eps6_Th8 ,eps7_Th8 ,Th8 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps9 ,Th1_eps9 ,Th2_eps9 ,Th3_eps9 ,Th4_eps9 ,Th5_eps9 ,Th6_eps9 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th7_eps9 ,Th8_eps9 ,Th9_eps9 ,eps1_Th9 ,eps2_Th9 ,eps3_Th9 ,eps4_Th9 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps5_Th9 ,eps6_Th9 , eps7_Th9 ,eps8_Th9 ,Th9 
! for averaging cross correlations
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  gsSqr, Thg, epsg, Th2g , eps2g , eps_Thg , eps_Th2g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps2g , eps2_Th2g , Th3g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  eps_Th3g , eps2_Th3g , eps3g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps3g , Th2_eps3g , Th3_eps3g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4g , Th1_eps4g , Th2_eps4g , Th3_eps4g , Th4_eps4g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th4g , eps2_Th4g , eps3_Th4g , Th4g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps5g , Th1_eps5g , Th2_eps5g , Th3_eps5g , Th4_eps5g , Th5_eps5g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th5g , eps2_Th5g , eps3_Th5g , eps4_Th5g , Th5g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps6g ,Th1_eps6g ,Th2_eps6g ,Th3_eps6g ,Th4_eps6g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps6g ,Th6_eps6g ,eps1_Th6g ,eps2_Th6g ,eps3_Th6g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4_Th6g ,eps5_Th6g ,Th6g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps7g ,Th1_eps7g ,Th2_eps7g ,Th3_eps7g ,Th4_eps7g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps7g ,Th6_eps7g ,Th7_eps7g ,eps1_Th7g ,eps2_Th7g 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps3_Th7g ,eps4_Th7g ,eps5_Th7g ,eps6_Th7g ,Th7g 
! for averaging cross correlations
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Thg2, epsg2, Th2g2 , eps2g2 , eps_Thg2 , eps_Th2g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps2g2 , eps2_Th2g2 , Th3g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  eps_Th3g2 , eps2_Th3g2 , eps3g2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th_eps3g2 , Th2_eps3g2 , Th3_eps3g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4g2 , Th1_eps4g2 , Th2_eps4g2 , Th3_eps4g2 , Th4_eps4g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th4g2 , eps2_Th4g2 , eps3_Th4g2 , Th4g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps5g2 , Th1_eps5g2 , Th2_eps5g2 , Th3_eps5g2 , Th4_eps5g2 , Th5_eps5g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps1_Th5g2 , eps2_Th5g2 , eps3_Th5g2 , eps4_Th5g2 , Th5g2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps6g2 ,Th1_eps6g2 ,Th2_eps6g2 ,Th3_eps6g2 ,Th4_eps6g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps6g2 ,Th6_eps6g2 ,eps1_Th6g2 ,eps2_Th6g2 ,eps3_Th6g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps4_Th6g2 ,eps5_Th6g2 ,Th6g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps7g2 ,Th1_eps7g2 ,Th2_eps7g2 ,Th3_eps7g2 ,Th4_eps7g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   Th5_eps7g2 ,Th6_eps7g2 ,Th7_eps7g2 ,eps1_Th7g2 ,eps2_Th7g2 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   eps3_Th7g2 ,eps4_Th7g2 ,eps5_Th7g2 ,eps6_Th7g2 ,Th7g2 
! Re shear stresses      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   epsUxUy, UxUy 
! Re normal stresses      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   epsUxUx, UxUx
! G-S drag correlation     
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::   epg2 , Vr2 , epgVr , &
	         Vr3 , epg3 , epg2Vr , epgVr2 , &
		Vr4 , epg4 , epg3Vr , epgVr3 , epg2Vr2 , &
		Vr5 , epg5 , epg4Vr , epgVr4 , epg3Vr2 , epg2Vr3 , &
		Vr6 , epg6 , epg5Vr , epgVr5 , epg4Vr2 , epg2Vr4 ,epg3Vr3 , &
		Vr7 , epg7 , epg6Vr , epgVr6 , epg5Vr2 , epg2Vr5 , &
	        epg4Vr3 , epg3Vr4 
      

      END MODULE usr                                                                             
