!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_MPPIC                                          !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!                                                                      !
!  Purpose: check the initial conditions input section for MPPIC model !
!     - ensure the first IC is defined over the entire domain with     ! 
!        ep_g = 1 when more than one IC has solids                     !
!     - ensure the ICs are non-overlapping                             !
!     - calculate the number of particles needed to initalize the      !
!        MPPIC model                                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_MPPIC
      

! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only : gener_part_config
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Number of particles seeded per phase 
      USE discretelement, only: PART_MPHASE
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles to read from input file.
      USE discretelement, only: PARTICLES

! Implied total number of physical particles 
      USE mfix_pic, only: rnp_pic 
! Total number of computational particles/parcels 
      USE mfix_pic, only: cnp_pic 

      USE mfix_pic, only: cnp_array
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! IC Region solid volume fraction.
      USE ic, only: IC_EP_S
      
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
      USE ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE

! MPPIC specific IC region specification. 
      USE ic, only: IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT
! Maximum number of IC regions and solids phases
      USE param, only: DIMENSION_IC, DIM_M 

! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: xlength, ylength, zlength, dx, dy, dz 

! Constant: 3.14159...
      USE constant, only: PI

      USE mpi_utility
      
! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

      DOUBLE PRECISION :: EP_SM
! Temp variable for storing number of particles.
      DOUBLE PRECISION :: PARTS_TEMP 
! Temp logical variables for checking constant npc and statwt specification
      LOGICAL :: CONST_NPC, CONST_STATWT

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M)
! Volume of the cell 
      DOUBLE PRECISION :: VOLIJK
      INTEGER :: I, J, K, IJK, ICV, M, COUNT_IC
!temporary variable for Dimension_3 
      INTEGER :: DIM3_TEMP
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE '../function.inc'

      IF (.NOT.GENER_PART_CONFIG) RETURN 

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_MPPIC")

! First check if either a constant npc or constant statwt 
! is specified for each IC 
      DO ICV = 1, DIMENSION_IC 
         
         IF(.not.ic_defined(icv)) cycle         

         IF (IC_EP_G(ICV).lt.ONE) THEN 
            DO M = 1, DES_MMAX           
               CONST_NPC    = (IC_PIC_CONST_NPC   (ICV, M) .ne. 0)
               CONST_STATWT = (IC_PIC_CONST_STATWT(ICV, M) .ne. ZERO  )
               IF(CONST_NPC.and.CONST_STATWT.and.ic_ep_s(icv,m).gt.zero) then
                  WRITE(ERR_MSG, 1100) ICV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
               
               IF(.not.CONST_NPC.and.(.not.CONST_STATWT).and. & 
               ic_ep_s(icv,m).gt.zero) then
                  WRITE(ERR_MSG, 1101) ICV, M 
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

 1100 FORMAT('Error 1100: In MPPIC model for IC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'Non zero Values specified for both ', &
      'IC_PIC_CONST_NPC and IC_PIC_CONST_STATWT.', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')


 1101 FORMAT('Error 1101: In MPPIC model for IC # ',i5, & 
      ' and solid phase # ', i5, /, & 
      'A non-zero value not specified for ', &
      'IC_PIC_CONST_NPC or IC_PIC_CONST_STATWT. ', /, &
      'Choose between constant number of parcels per cell or ', &
      'constant statistical weight', /, &
      'See MFIX readme',/'Please correct the data file.')

      
   
      IF(DIMN.EQ.2) THEN 
! require that DZ(1)/ZLENGTH be specified for 2-dimensional case.  
! unclear why this cannot be one - other than the user may be unaware 
! that a depth has been set (a value of one implies default setting) 
         IF (DZ(1) == ONE) THEN
            WRITE(*,'(5X,A,A,/5X,A,A)') &
            'For DIMN = 2, specify a value for DZ(1) or ',&
            'ZLENGTH in mfix.dat which is not',&
            'equal to one. If you want it to be one then ',&
            'set it close to one but not exactly one'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (DZ(1) .NE. ZLENGTH) THEN 
            WRITE(ERR_MSG,'(5X,A,/5x,A,/5X,2(A20,2X,G17.8))')        &
            'For DIMN = 2, DZ(1) and ZLENGTH are used ',         &
            'interchangeably', ' Specify same values for ',      &
            'DZ(1) and ZLENGTH','DZ(1) = ', DZ(1), 'ZLENGTH = ', &
            ZLENGTH
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF


      CALL FINL_ERR_MSG
      
      END SUBROUTINE CHECK_IC_MPPIC
      
