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
      USE des_bc
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I

!-----------------------------------------------

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

      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
        DES_NEIGHBOR_SEARCH .EQ. 3) THEN
          LQUAD(:,:) = UNDEFINED_I
          CQUAD(:,:) = UNDEFINED
          PQUAD(:) = 0
      ENDIF

      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0
      PN(:,:) = -1
      PN(:,1) = 0
      PV(:,:) = 1
      PFN(:,:,:) = ZERO
      PFT(:,:,:) = ZERO

      PINC(:) = ZERO
      PIJK(:,:) = ZERO

      DES_WALL_POS(:,:) = UNDEFINED
      DES_WALL_VEL(:,:) = UNDEFINED

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      SOLID_DRAG(:,:,:) = ZERO

      XE(:) = ZERO
      YN(:) = ZERO
      ZT(:) = ZERO

      PEA(:,:) = .FALSE.
      DO I=1, PARTICLES
         PEA(I,1)=.TRUE.
      ENDDO

      RETURN
      END SUBROUTINE DES_INIT_ARRAYS 


