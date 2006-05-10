!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!     loop counters
      INTEGER :: LC, LCM, M, N, K, KKK
!     
!     Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E
!-----------------------------------------------
!

                   DO LC = 1, PARTICLES 
                      DES_RADIUS(LC) = ZERO
                      PMASS(LC) = ZERO
                      PVOL(LC) = ZERO
                      OMOI(LC) = ZERO
                      RO_Sol(LC) = ZERO 
                         DES_POS_OLD(LC,:) = ZERO
                         DES_POS_NEW(LC,:) = ZERO
                         DES_VEL_OLD(LC,:) = ZERO
                         DES_VEL_NEW(LC,:) = ZERO
                         FC(LC,:) = ZERO
                         FN(LC,:) = ZERO
                         FT(LC,:) = ZERO
                         TOW(LC,:) = ZERO
                         OMEGA_OLD(LC,:) = ZERO
                         OMEGA_NEW(LC,:) = ZERO
                   END DO
                   
                         GRAV(:) = ZERO
                         FNS1(:) = ZERO
                         FTS1(:) = ZERO

                   DO LC = 1, MAXQUADS
                      DO K = 1, NMQD 
                         LQUAD(LC,K) = UNDEFINED_I
                      END DO
                      DO K = 1, NWALLS
                         CQUAD(LC,K) = UNDEFINED
                      END DO
                   END DO

                   DO LC= 1, PARTICLES
                      DO KKK = 2, MAXNEIGHBORS
                         NEIGHBOURS(K,KKK) = -1
                         PN(K,KKK) = -1
                         PV(K,KKK) = 1
                         PFN(LC,KKK,:) = ZERO
                         PFT(LC,KKK,:) = ZERO
                      END DO
                      NEIGHBOURS(LC,1) = 0
                      PN(LC,1) = 0
                      PV(LC,1) = 1
                      PIJK(LC,:) = ZERO
                   END DO

                   DO K = 1, NWALLS
                      DES_WALL_POS(K,:) = UNDEFINED
                      DES_WALL_VEL(K,:) = UNDEFINED
                      WALL_NORMAL(K,:) = UNDEFINED
                   END DO

                   DO LC = 1, DIMENSION_3
                      PINC(LC) = ZERO
                      DO K = 1, MMAX
                         DES_U_s(LC,K) = ZERO
                         DES_W_s(LC,K) = ZERO
                         DES_W_s(LC,K) = ZERO
                         SOLID_DRAG(LC,K,:) = ZERO
                      END DO
                   END DO

		   DO LC = 1, DIMENSION_I
                      XE(LC) = ZERO
                   END DO
		   DO LC = 1, DIMENSION_J
                      YN(LC) = ZERO
                   END DO
		   DO LC = 1, DIMENSION_K
                      ZT(LC) = ZERO
                   END DO
                  
                  RETURN
                  END SUBROUTINE DES_INIT_ARRAYS
