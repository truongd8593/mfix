!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_MODEL_PREREQS                              !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS


! Global Variables:
!---------------------------------------------------------------------//
! Number of ranks.
      use run, only: SOLIDS_MODEL

! Flag: Use DES E-L model
      use discretelement, only: DISCRETE_ELEMENT
! Flag: Use MPPIC E-L model
      use mfix_pic, only: MPPIC
! Flag: Use TFM and DEM solids models.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Flag: gas/solids E-L simulation, otherwise granular flow.
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag: Fluid affects particles, but particles do not impact fluid.
      use discretelement, only: DES_ONEWAY_COUPLED
! Flag: Solving DEM species equations.
      use des_rxns, only: ANY_DES_SPECIES_EQ
! Flag: gas/solids E-L convective heat transfer.
      use des_thermo, only: DES_CONV_EQ
! Flag: Interpolate between gas/solids
      use discretelement, only: DES_INTERP_ON
! Number of discrtete solids phases.
      use discretelement, only: DES_MMAX
! Number of phases specified by the user.
      use physprop, only: MMAX
! Print E-L data.
      use discretelement, only: PRINT_DES_DATA

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids phases.
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M ! Phase index
! Flags for various solids phase models.
      INTEGER :: TFM_SOLIDS = 0
      INTEGER :: DEM_SOLIDS = 0
      INTEGER :: PIC_SOLIDS = 0

! Flag for two-fluid model solids.
      LOGICAL :: CONTINUUM_SOLIDS


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MODEL_PREREQS")


! Check MMAX
      IF (MMAX<0 .OR. MMAX>DIM_M) THEN 
         WRITE(ERR_MSG, 1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF 

 1000 FORMAT('Error 1000: MMAX out of range. Min: 0, Max: ',I2)

! Loop over the phases to see what was specified.
      DO M=1, MMAX
         SOLIDS_MODEL(M) = trim(adjustl(SOLIDS_MODEL(M)))
         SELECT CASE(SOLIDS_MODEL(M))
         CASE ('TFM')
            TFM_SOLIDS = TFM_SOLIDS + 1
         CASE ('DEM')
            TFM_SOLIDS = DEM_SOLIDS + 1
         CASE ('MPPIC')
            PIC_SOLIDS = PIC_SOLIDS + 1
         CASE DEFAULT
            WRITE(ERR_MSG,1001) M, SOLIDS_MODEL(M)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         END SELECT
      ENDDO

 1001 FORMAT('Error 1001: Unknown solids model: ',/3x,'SOLIDS_MODEL(', &
         I2,') = ', A)


! MPPIC and continuum solids don't mix.
      IF((PIC_SOLIDS > 0) .AND. (TFM_SOLIDS > 0))THEN
         WRITE(ERR_MSG, 1002)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1002 FORMAT('Error 1002: MPPIC solids and continuum solids cannot ',&
         'be combined.',/'Please correct the mfix.dat file.')


! MPPIC and DEM solids don't mix.
      IF((PIC_SOLIDS > 0) .AND. (DEM_SOLIDS > 0))THEN
         WRITE(ERR_MSG, 1003)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 1003 FORMAT('Error 1003: MPPIC solids and DES solids cannot be ',     &
         'combined.',/'Please correct the mfix.dat file.')


! Set the number of discrete phases.
      DES_MMAX = DEM_SOLIDS + PIC_SOLIDS


   !****************************************************************!
   !                       -->   WARNING   <---                     !
   !                                                                !
   ! The following checks set the run time flags based on input     !
   ! from the mfix.dat file.  The 'self setting' needs removed once !
   ! these variables are no longer namelist variables.              !
   !                                                                !
   !****************************************************************!

! Set the two-fluid model flag.
      CONTINUUM_SOLIDS = (TFM_SOLIDS > 0)
! Set the DEM runtime flag.
      DISCRETE_ELEMENT = (DEM_SOLIDS > 0) .OR. DISCRETE_ELEMENT
! Set the MMPIC runtime flag.
      MPPIC = (PIC_SOLIDS > 0) .OR. MPPIC
! Set the Hybird flag.
      DES_CONTINUUM_HYBRID = ((DEM_SOLIDS > 0) .AND. (TFM_SOLIDS > 0)) &
         .OR. DES_CONTINUUM_HYBRID


! This is a legacy check that needs removed once this variable
! is removed from being a namelist keyword.
      IF(.NOT.DISCRETE_ELEMENT) DES_CONTINUUM_HYBRID = .FALSE. !<------ TO BE REMOVED

! This is a legacy check that needs removed once this variable
! is removed from being a namelist keyword.
      IF(.NOT.DISCRETE_ELEMENT) MPPIC = .FALSE. !<--------------------- TO BE REMOVED

! This is a legacy check that needs removed once this variable
! is removed from being a namelist keyword.
      IF(MPPIC) DISCRETE_ELEMENT = .TRUE. !<--------------------------- TO BE REMOVED



! Overwrite user settings if no Lagrangian solids
      IF(.NOT.DISCRETE_ELEMENT) THEN
         DES_CONTINUUM_COUPLED = .FALSE.
         DES_INTERP_ON = .FALSE.
         PRINT_DES_DATA = .FALSE.
         DES_ONEWAY_COUPLED = .false. 
         DES_CONV_EQ = .FALSE.
         ANY_DES_SPECIES_EQ = .FALSE.
      ENDIF


! Most likely, the TFM/DEM hybird model is going to break a lot in the 
! process to generalize the input. I'm going to disable it now so it
! is not available.
      IF(DES_CONTINUUM_HYBRID)THEN
         WRITE(ERR_MSG, 9000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 9000 FORMAT('Error 9000: TFM/DEM Hybrid model has been disabled.',/   &
         'This will be restored shortly. Sorry :(')




      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_MODEL_PREREQS
