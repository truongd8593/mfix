!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PHYSICAL_PROP(DENSITY, PSIZE, SP_HEAT,IER)              C
!  Purpose: Calculate physical properties that vary with time          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI                                                   C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE PHYSICAL_PROP(DENSITY, PSIZE, SP_HEAT, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      USE constant
      USE scalars
      USE compar 
      USE funits 
      USE usr   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      gas constant
      DOUBLE PRECISION, PARAMETER :: RGAS = 1.987207D0  !cal/mol.K
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK, M, N
!
!                      Average molecular weight
      DOUBLE PRECISION MW
!
!                      Flags to tell whether to calculate or not
      LOGICAL          DENSITY(0:DIMENSION_M), PSIZE(0:DIMENSION_M),&
                       SP_HEAT(0:DIMENSION_M)
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG, calc_CpoR, calc_ICpoR
!-----------------------------------------------
      INCLUDE 'species_indices.inc'
      INCLUDE 'usrnlst.inc'
!
      DOUBLE PRECISION TGX, TSX,  DIFF, EP_g2
      DOUBLE PRECISION Sc1o3, UGC, VGC, WGC, USCM, VSCM, WSCM, VREL, Re
      INTEGER IMJK, IJMK, IJKM,I
    
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!  Fluid phase
!

!$omp  parallel do private(IJK, MW, N)  
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
	   TGX = T_g(IJK)
!
!
! 1.1      Density
!
            IF (DENSITY(0)) THEN 
               IF (MW_AVG == UNDEFINED) THEN 
!              Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
!              -----------------------------------------------------
!               MW = CALC_MW(X_g, DIMENSION_3, IJK, NMAX(0), MW_g)
                  if(.not.database_read) call read_database(IER)
                  MW = ZERO
                  N = 1
                  IF (NMAX(0) > 0) THEN
                     MW = MW + SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
                     N = NMAX(0) + 1
                  ENDIF
                  MW = ONE/MAX(MW,OMW_MAX)
!              -----------------------------------------------------
                  MW_MIX_G(IJK) = MW
                  RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ELSE
                  RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF
!	       
	       IF(RO_G(IJK) < ZERO) THEN
		 WRITE(*,1000) I_OF(IJK), J_OF(IJK), K_OF(IJK), &
		               RO_G(IJK), P_G(IJK), T_G(IJK)
		 CALL MFIX_EXIT(myPE)
	       ENDIF
!   
            ENDIF
!
!      Gas specific heat
!       1 Cal = 4.183925 J
!
!           Constant pressure specific heat of air in cal/g.K
            IF (SP_HEAT(0) .AND. C_PG0==UNDEFINED) then
              if(.not.database_read) call read_database(IER)
	      
!	      IF(C(23) == ONE) THEN
	        C_PG(IJK) = ZERO
                DO N = 1, NMAX(0)
	          C_PG(IJK) = C_PG(IJK) + X_g(IJK,N) * calc_CpoR(T_G(IJK), Thigh_g(N), Tlow_g(N), &
		               Tcom_g(N), Ahigh_g(1,N), Alow_g(1,N)) * RGAS / MW_g(N) 
		ENDDO
!	      ELSE
!                C_pg(IJK) =  X_g(IJK,CH4)*CPCH4(TGX)  &
!                       + X_g(IJK,CO2)*CPCO2(TGX)  &
!                       + X_g(IJK,O2)*CPO2(TGX) + X_g(IJK,H2O)*CPH2O(TGX) &
!                       + X_g(IJK,N2)*CPN2(TGX) 
!	      ENDIF
	      
              !to SI, S. Dartevelle
              IF (UNITS == 'SI') C_PG(IJK) = 4183.925D0*C_PG(IJK)    !in J/kg K
	    ENDIF
!
!
         ENDIF 
      END DO 
!
! add by rong
! diameter
        
      IF (CALL_DQMOM) THEN
         DO M=1,MMAX     
!$omp  parallel do private(IJK)
            DO IJK = IJKSTART3, IJKEND3 
               IF (.NOT.WALL_AT(IJK)) THEN 
                  N=phase4scalar(M)
                  IF(EP_s(IJK,N)>small_number) THEN
                     IF(PSIZE(M)) THEN
                        D_p(IJK,M)= Scalar(IJK,N)
                     ELSE
                        D_P(IJK,M)=D_p0(M)
                     ENDIF
                  ELSE
                     D_P(IJK,M)=D_P0(M)
                  ENDIF
               ENDIF
            END DO
         END DO
       ELSE
          DO M=1,MMAX
!$omp  parallel do private(IJK)
             DO IJK = IJKSTART3, IJKEND3
                IF (.NOT.WALL_AT(IJK)) THEN 
                   D_P(IJK,M)=D_P0(M)
                ENDIF
             END DO
          END DO  
       END IF
! 

! add by rong   
!
      DO M = 1, MMAX 
!
!$omp  parallel do private(IJK)  
         DO IJK = IJKSTART3, IJKEND3 
            IF (.NOT.WALL_AT(IJK)) THEN 
!
!             Specific heat of solids (Coal = 0.3 cal/g.K)
!             Perry & Chilton(1973) -- Table 3-201 on page 3-136
!             Specific heat of solids (Ash =  0.310713 cal/g.K)
!             Dobran et al., 1991
              IF (SP_HEAT(M) .AND. C_PS0==UNDEFINED)then
                if(.not.database_read) call read_database(IER)
	        C_PS(IJK,M) = 0.310713
                !to SI, S. Dartevelle
                IF (UNITS == 'SI') C_PS(IJK,M) = 4183.925D0*C_PS(IJK,M)    !in J/kg K
	      ENDIF
!
            END IF
         END DO 
      END DO  

     
      RETURN 
 1000 FORMAT(1X,'Message from: PHYSICAL_PROP',/& 
            'WARNING: Gas density negative in this cell: ', /&
	    'I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
	    'Values of variables: ','RO_g = ', G12.5, 2X, 'P_g = ', G12.5, 2X, 'T_g = ', G12.5, /&
	    'Suggestion: Lower UR_FAC(1) in mfix.dat')
      END SUBROUTINE PHYSICAL_PROP 



!  New functions based on mfix/model/thermochemical/burcat.thr
      DOUBLE PRECISION FUNCTION CPCH4(T)
        IMPLICIT NONE
	DOUBLE PRECISION T, calc_CpoR
	DOUBLE PRECISION, DIMENSION(7) :: Ahigh, Alow
	Ahigh = (/&
   1.65326226000000D0,        1.00263099000000D-002,  -3.31661238000000D-006, &
   5.36483138000000D-010,  -3.14696758000000D-014,  -10009.5936000000D0, &
   9.90506283000000D0/)
	Alow = (/&
   5.14911468000000D0,       -1.36622009000000D-002,   4.91453921000000D-005, &
  -4.84246767000000D-008,   1.66603441000000D-011,  -10246.5983000000D0, &
  -4.63848842000000D0/)
	
        CPCH4 = calc_CpoR(T, 6000D0, 200D0, 1000D0, Ahigh, Alow)*1.987207D0/16.D0
      END FUNCTION CPCH4
      
      DOUBLE PRECISION FUNCTION CPCO2(T)
        IMPLICIT NONE
	DOUBLE PRECISION T, calc_CpoR
	DOUBLE PRECISION, DIMENSION(7) :: Ahigh, Alow
	Ahigh = (/&
   4.63651110000000D0,        2.74145690000000D-003,  -9.95897590000000D-007, &
   1.60386660000000D-010,  -9.16198570000000D-015,  -49024.9040000000D0, &
  -1.93489550000000D0/)
	Alow = (/&
   2.35681300000000D0,        8.98412990000000D-003,  -7.12206320000000D-006, &
   2.45730080000000D-009,  -1.42885480000000D-013,  -48371.9710000000D0, &
   9.90090350000000D0/)
 	
        CPCO2 = calc_CpoR(T, 6000D0, 200D0, 1000D0, Ahigh, Alow)*1.987207D0/44.D0
      END FUNCTION CPCO2
      
      DOUBLE PRECISION FUNCTION CPO2(T)
        IMPLICIT NONE
	DOUBLE PRECISION T, calc_CpoR
	DOUBLE PRECISION, DIMENSION(7) :: Ahigh, Alow
	Ahigh = (/&
   3.66096083000000D0,        6.56365523000000D-004,  -1.41149485000000D-007, &
   2.05797658000000D-011,  -1.29913248000000D-015,  -1215.97725000000D0, &
   3.41536184000000D0/)
	Alow = (/&
   3.78245636000000D0,       -2.99673415000000D-003,   9.84730200000000D-006, &
  -9.68129508000000D-009,   3.24372836000000D-012,  -1063.94356000000D0, &
   3.65767573000000D0/)
 	
        CPO2 = calc_CpoR(T, 6000D0, 200D0, 1000D0, Ahigh, Alow)*1.987207D0/32.D0
      END FUNCTION CPO2
      
      DOUBLE PRECISION FUNCTION CPH2O(T)
        IMPLICIT NONE
	DOUBLE PRECISION T, calc_CpoR
	DOUBLE PRECISION, DIMENSION(7) :: Ahigh, Alow
	Ahigh = (/&
   2.67703890000000D0,        2.97318160000000D-003,  -7.73768890000000D-007, &
   9.44335140000000D-011,  -4.26899910000000D-015,  -29885.8940000000D0, &
   6.88255000000000D0/)
	Alow = (/&
   4.19863520000000D0,       -2.03640170000000D-003,   6.52034160000000D-006, &
  -5.48792690000000D-009,   1.77196800000000D-012,  -30293.7260000000D0, &
  -0.849009010000000D0/)
	
        CPH2O = calc_CpoR(T, 6000D0, 200D0, 1000D0, Ahigh, Alow)*1.987207D0/18.D0
      END FUNCTION CPH2O
      
      DOUBLE PRECISION FUNCTION CPN2(T)
        IMPLICIT NONE
	DOUBLE PRECISION T, calc_CpoR
	DOUBLE PRECISION, DIMENSION(7) :: Ahigh, Alow
	Ahigh = (/&
   2.95257637000000D0,        1.39690040000000D-003,  -4.92631603000000D-007, &
   7.86010195000000D-011,  -4.60755204000000D-015,  -923.948688000000D0, &
   5.87188762000000D0/)
	Alow = (/&
   3.53100528000000D0,       -1.23660988000000D-004,  -5.02999433000000D-007, &
   2.43530612000000D-009,  -1.40881235000000D-012,  -1046.97628000000D0, &
   2.96747038000000D0/)
 	
        CPN2 = calc_CpoR(T, 6000D0, 200D0, 1000D0, Ahigh, Alow)*1.987207D0/28.D0
      END FUNCTION CPN2
      
