!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_FLUIDBED_P                                         C
!  Purpose: Set the pressure field inside the bed assuming a fluidized C
!           bed with gravity acting the -ve y-direction                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!  Reviewer:M. Syamlal, S. Venkatesan, P. Nicoletti,  Date: 29-JAN-92  C
!           W. Rogers                                                  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for including cylindrical geometry           C
!  Author: M. Syamlal                                 Date:  6-MAR-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 2                                                  C
!  Purpose: Set pressure drop for cyclic boundary condition w/ pr. dropC
!  Author: M. Syamlal                                 Date: 29-APR-94  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:BC_DEFINED, BC_TYPE, P_OUTFLOW, BC_P_g, KMAX1, C
!                       IMAX1, DX, DY, DZ, G, EP_g, MW_AVG, T_g, MMAX, C
!                       ROP_s, IMIN1, KMIN1, JMIN1, X                  C
!                                                                      C
!  Variables modified: P_g, IJK, I, J, K                               C
!                                                                      C
!  Local variables: PJ, BED_WEIGHT, AREA, dAREA                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_FLUIDBED_P 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop
      USE geometry
      USE bc
      USE ic
      USE fldvar
      USE constant
      USE indices
      USE funits 
      USE scales 
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      indices 
      INTEGER          I, J, K, IJK, M 
! 
!                      Local loop counter 
      INTEGER          L 
! 
!                      Gas pressure at the axial location j 
      DOUBLE PRECISION PJ 
! 
!                      Bed weight per unit area 
      DOUBLE PRECISION BED_WEIGHT 
! 
!                      Total area of a x-z plane 
      DOUBLE PRECISION AREA 
! 
!                      x-z plane area of one cell 
      DOUBLE PRECISION dAREA 
! 
!                      Average pressure drop per unit length 
      DOUBLE PRECISION DPoDX, DPoDY, DPoDZ 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG 
!-----------------------------------------------
      INCLUDE 'sc_p_g1.inc'
      INCLUDE 'b_force1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'b_force2.inc'
      INCLUDE 'sc_p_g2.inc'
!
!  IF all initial pressures are specified, then return.
!
      DO L = 1, DIMENSION_IC 
         IF (IC_DEFINED(L)) THEN 
            IF (IC_P_G(L) == UNDEFINED) GO TO 60 
            PJ = IC_P_G(L) 
         ENDIF 
      END DO 
      IF (DO_I .AND. DELP_X/=UNDEFINED) THEN 
         DPODX = DELP_X/XLENGTH 
         PJ = PJ - DPODX*HALF*(DX(IMAX1)+DX(IMAX2)) 
         DO I = IMAX1, IMIN1, -1 
            PJ = PJ + DPODX*HALF*(DX(I)+DX(I+1)) 
            DO K = KMIN1, KMAX1 
               DO J = JMIN1, JMAX1 
                  IJK = FUNIJK(I,J,K) 
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE(PJ) 
               END DO 
            END DO 
         END DO 
      ENDIF 
!
      IF (DO_J .AND. DELP_Y/=UNDEFINED) THEN 
         DPODY = DELP_Y/YLENGTH 
         PJ = PJ - DPODY*HALF*(DY(JMAX1)+DY(JMAX2)) 
         DO J = JMAX1, JMIN1, -1 
            PJ = PJ + DPODY*HALF*(DY(J)+DY(J+1)) 
            DO K = KMIN1, KMAX1 
               DO I = IMIN1, IMAX1 
                  IJK = FUNIJK(I,J,K) 
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE(PJ) 
               END DO 
            END DO 
         END DO 
      ENDIF 
!
      IF (DO_K .AND. DELP_Z/=UNDEFINED) THEN 
         DPODZ = DELP_Z/ZLENGTH 
         PJ = PJ - DPODZ*HALF*(DZ(KMAX1)+DZ(KMAX2)) 
         DO K = KMAX1, KMIN1, -1 
            PJ = PJ + DPODZ*HALF*(DZ(K)+DZ(K+1)) 
            DO J = JMIN1, JMAX1 
               DO I = IMIN1, IMAX1 
                  IJK = FUNIJK(I,J,K) 
                  IF (FLUID_AT(IJK)) P_G(IJK) = SCALE(PJ) 
               END DO 
            END DO 
         END DO 
      ENDIF 
!4/29/94
      RETURN  
!
!  Search for an outflow boundary condition where pressure is specified
!
   60 CONTINUE 
      PJ = UNDEFINED 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L) .AND. BC_TYPE(L)=='P_OUTFLOW') PJ = BC_P_G(L) 
      END DO 
      IF (PJ == UNDEFINED) THEN 
         IF (RO_G0 /= UNDEFINED) THEN 
!
!         If incompressible flow set P_g to zero.
!
            DO IJK = 1, IJKMAX2 
               IF (FLUID_AT(IJK)) P_G(IJK) = ZERO 
            END DO 
            RETURN  
         ELSE 
!
!         Error condition -- no pressure outflow boundary condition is specified
!
            CALL START_LOG 
            WRITE (UNIT_LOG, 1000) 
            CALL MFIX_EXIT 
         ENDIF 
      ENDIF 
!
!  Set an approximate pressure field assuming that the pressure drop
!  balances the weight of the bed, if the initial pressure-field is not
!  specified
!
      DO J = JMAX2, JMIN1, -1 
!
!  Find the average weight per unit area over an x-z slice
!
         BED_WEIGHT = 0.0 
         AREA = 0.0 
         DO K = KMIN1, KMAX1 
            DO I = IMIN1, IMAX1 
               IJK = FUNIJK(I,J,K) 
               IF (FLUID_AT(IJK)) THEN 
                  IF (COORDINATES == 'CARTESIAN') THEN 
                     DAREA = DX(I)*DZ(K) 
                  ELSE IF (CYLINDRICAL) THEN 
                     DAREA = DX(I)*X(I)*DZ(K) 
                  ENDIF 
                  AREA = AREA + DAREA 
                  IF (RO_G0 == UNDEFINED) THEN 
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*EOSG(&
                        MW_MIX_G(IJK),PJ,T_G(IJK))*DAREA 
                  ELSE 
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_G(IJK)*EP_G(IJK)*RO_G0&
                        *DAREA 
                  ENDIF 
                  DO M = 1, MMAX 
                     BED_WEIGHT = BED_WEIGHT - DY(J)*BFY_S(IJK,1)*ROP_S(IJK,M)*&
                        DAREA 
                  END DO 
               ENDIF 
            END DO 
         END DO 
         IF (AREA /= 0.0) BED_WEIGHT = BED_WEIGHT/AREA 
         PJ = PJ + BED_WEIGHT 
         DO K = KMIN1, KMAX1 
            DO I = IMIN1, IMAX1 
               IJK = FUNIJK(I,J,K) 
               IF(FLUID_AT(IJK).AND.P_G(IJK)==UNDEFINED)P_G(IJK)=SCALE(PJ) 
            END DO 
         END DO 
      END DO 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: SET_FLUIDBED_P'/' Message: Outflow ',&
         'pressure boundary condition (P_OUTFLOW) not found.',/&
         'All the initial pressures (IC_P_g) or at least one P_OUTFLOW',/&
         'condition need to be specified',/1X,70('*')/) 
      END SUBROUTINE SET_FLUIDBED_P 
