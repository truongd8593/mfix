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
      USE compar 
      USE funits 
      USE scalars   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'cp_fun2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Fluid phase
!

!$omp  parallel do private(IJK, MW, N)  
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
!
! 1.1      Density
!
            IF (DENSITY(0)) THEN 
               IF (MW_AVG == UNDEFINED) THEN 
!              Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
!              -----------------------------------------------------
!               MW = CALC_MW(X_g, DIMENSION_3, IJK, NMAX(0), MW_g)
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
            ENDIF
!
!
!       Gas specific heat

!              1 Cal = 4.183925 J
!
!           Constant pressure specific heat of air in cal/g.K
            IF (SP_HEAT(0) .AND. C_PG0==UNDEFINED) then
               C_PG(IJK) = 0.767*CPN2(T_G(IJK)) + 0.233*CPO2(T_G(IJK))
               !to SI, S. Dartevelle
               IF (UNITS == 'SI') C_PG(IJK) = 4183.925*C_PG(IJK)    !in J/kg K
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
               IF (SP_HEAT(M) .AND. C_PS0==UNDEFINED) then
	         C_PS(IJK,M) = 0.310713
                 !to SI, S. Dartevelle
                 IF (UNITS == 'SI') C_PS(IJK,M) = 4183.925*C_PS(IJK,M)    !in J/kg K
	       ENDIF
!
            ENDIF 
         END DO 
      END DO 

     
      RETURN  
      END SUBROUTINE PHYSICAL_PROP 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
