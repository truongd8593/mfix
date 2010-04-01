!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
    
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      
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

      DES_VEL_OOLD(:,:) = ZERO
      DES_ACC_OLD(:,:) = ZERO

      OMEGA_OLD(:,:) = ZERO
      OMEGA_NEW(:,:) = ZERO
      ROT_ACC_OLD(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      SOLID_DRAG(:,:,:) = ZERO

      FC(:,:) = ZERO
      FN(:,:) = ZERO
      FT(:,:) = ZERO
      TOW(:,:) = ZERO

      PPOS(:,:) = ZERO
      GRAV(:) = ZERO
      DES_WALL_POS(:,:) = UNDEFINED
      DES_WALL_VEL(:,:) = UNDEFINED

      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0
      PN(:,:) = -1
      PN(:,1) = 0
      PV(:,:) = 1
      PFT(:,:,:) = ZERO

      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
        DES_NEIGHBOR_SEARCH .EQ. 3) THEN
          LQUAD(:,:) = UNDEFINED_I
          CQUAD(:,:) = UNDEFINED
          PQUAD(:) = 0
      ENDIF      

      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
         DESGRIDSEARCH_PIJK(:,:) = ZERO
      ENDIF


      PINC(:) = ZERO
      PIJK(:,:) = ZERO

      XE(:) = ZERO
      YN(:) = ZERO
      ZT(:) = ZERO

! J.Musser: DEM inlet/outlet
      PEA(:,:) = .FALSE.
! If RESTART_1, PEA will be read in from the restart file
      IF(RUN_TYPE == 'NEW') THEN
         DO I=1, PARTICLES
            PEA(I,1)=.TRUE.
         ENDDO
      ENDIF
      DES_BC_U_s(:) = ZERO
      DES_BC_V_s(:) = ZERO
      DES_BC_W_s(:) = ZERO

! T.Li : Hertzian collision model
      g_mod(:) = zero
      hert_kn(:,:) = zero
      hert_kwn(:) = zero
      hert_kt(:,:) = zero
      hert_kwt(:) = zero
      
      RETURN
      END SUBROUTINE DES_INIT_ARRAYS 


