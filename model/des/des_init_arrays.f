!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!>     Purpose: DES - initialize the des-namelist                          
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
      INTEGER :: M, N, K
!     
!     Coefficient of restitution (old symbol)
      DOUBLE PRECISION :: E
!-----------------------------------------------
!

                      DES_RADIUS(:) = ZERO
                      PMASS(:) = ZERO
                      PVOL(:) = ZERO
                      OMOI(:) = ZERO
                      RO_Sol(:) = ZERO 
                      DES_POS_OLD(:,:) = ZERO
                      DES_POS_NEW(:,:) = ZERO
                      DES_VEL_OLD(:,:) = ZERO
                      DES_VEL_NEW(:,:) = ZERO
                      FC(:,:) = ZERO
                      FN(:,:) = ZERO
                      FT(:,:) = ZERO
                      TOW(:,:) = ZERO
                      OMEGA_OLD(:,:) = ZERO
                      OMEGA_NEW(:,:) = ZERO
                      PPOS(:,:) = ZERO
                   
                      GRAV(:) = ZERO
                      FNS1(:) = ZERO
                      FTS1(:) = ZERO

                      LQUAD(:,:) = UNDEFINED_I
                      CQUAD(:,:) = UNDEFINED
		      PQUAD(:) = 0
                      NEIGHBOURS(:,:) = -1
                      PN(:,:) = -1
                      PV(:,:) = 1
                      PFN(:,:,:) = ZERO
                      PFT(:,:,:) = ZERO

                      NEIGHBOURS(:,1) = 0
                      PN(:,1) = 0
                      PV(:,1) = 1
                      PIJK(:,:) = ZERO
                      PN_DIST(:,:) = ZERO
                      PN_RLM(:,:) = ZERO

                      DES_WALL_POS(:,:) = UNDEFINED
                      DES_WALL_VEL(:,:) = UNDEFINED

                      PINC(:) = ZERO
                      DES_U_s(:,:) = ZERO
                      DES_W_s(:,:) = ZERO
                      DES_W_s(:,:) = ZERO
                      SOLID_DRAG(:,:,:) = ZERO

                      XE(:) = ZERO
                      YN(:) = ZERO
                      ZT(:) = ZERO
 
                  RETURN
                  END SUBROUTINE DES_INIT_ARRAYS
