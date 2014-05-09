!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM

      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE mfix_pic
      use mpi_utility

      use bc

      use error_manager

      IMPLICIT NONE

!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER, PARAMETER :: DIM_BCTYPE = 4

!-----------------------------------------------
! Local variables
!-----------------------------------------------

! Loop counters
      INTEGER BCV, BCV_I
! Solids phase index
      INTEGER M
! tmp variable to calculate solids volume fraction at inlet   
      DOUBLE PRECISION EPs_tmp


!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 

!-----------------------------------------------           


      CALL INIT_ERR_MSG("SET_BC_DEM")




! Set the flag that one or more DEM MI/MO exists.
      DEM_MIO = (DEM_BCMI /= 0 .OR. DEM_BCMO /= 0)

! The variable PARTICLES should already be set by this point if using
! gener_part_config option      
      IF(PARTICLES == UNDEFINED_I .AND. MAX_PIS /= UNDEFINED_I)THEN
         PARTICLES = 0

      ELSEIF(PARTICLES == UNDEFINED_I .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1200 FORMAT('Error 1200: Either PARTICLES or MAX_PIS must specified.',&
         'Please correct the mfix.dat file.')

      ELSEIF(PARTICLES == 0 .AND. MAX_PIS == UNDEFINED_I) THEN
         WRITE(ERR_MSG, 1201)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1201 FORMAT('Error 1201: MAX_PIS must be specified in the mfix.dat ', &
         'file if',/' there are no initial particles (PARTICLES = 0).')

      ENDIF 

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
! Inlet/outlet for MPPIC are based off the regular mfix declarations, 
! and so DEM_BCMI could still be zero.
      IF(PARTICLES == 0 .AND. DEM_BCMI == 0) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG
      ENDIF

 1202 FORMAT('WARNING 1202: The system is initiated with no particles',&
         ' and no',/'solids inlet was detected.')


! Check MAX_PIS requirements
      IF(DEM_BCMI == 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         MAX_PIS = PARTICLES
      ELSEIF(DEM_BCMI /= 0 .AND. MAX_PIS == UNDEFINED_I)THEN
         WRITE(ERR_MSG, 1203)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1203 FORMAT('Error 1203: The maximum number of particles (MAX_PIS) ', &
         'must be',/'spedified if a DEM inlet is specifed. Please ',   &
         'correct the mfix.dat ',/'file.')




      IF(DEM_BCMI > 0) CALL SET_BC_DEM_MI
!      IF(DEM_BCMO > 0) CALL SET_BC_DEM_MO

      write(*,"(   2x,'DEM outlets:',I4)") DEM_BCMO


 !     stop

      CALL FINL_ERR_MSG


      RETURN
      END SUBROUTINE SET_BC_DEM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MI                                           !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE bc
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run

      use error_manager

      IMPLICIT NONE

      INTEGER :: BCV


!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV_I      ! BC loop counter
      INTEGER M, MM           ! Mass phase loop counter
      INTEGER HOLD, I         ! Dummy values
      INTEGER RANGE_TOP, RANGE_BOT ! Dummy values
      INTEGER PHASE_CNT        ! Number of solid phases at bc
      INTEGER PHASE_LIST(DES_MMAX) ! List of phases used in current bc

! the number of particles injected in a solids time step
      DOUBLE PRECISION NPMpSEC(DES_MMAX) ! For solid phase m
      DOUBLE PRECISION NPpSEC
      DOUBLE PRECISION NPpDT        ! Total for BC
      DOUBLE PRECISION SCALED_VAL
      DOUBLE PRECISION MAX_DIA ! Max diameter of incoming particles at bc

      DOUBLE PRECISION :: EPs_ERR
      DOUBLE PRECISION :: VOL_FLOW


      LOGICAL, parameter :: setDBG = .FALSE. ! .FALSE.
      LOGICAL :: dFlag


      LOGICAL :: FATAL

! Temp inlet velocity for solids phase M
      DOUBLE PRECISION VEL_TMP(DIM_M)
      DOUBLE PRECISION EPs_TMP(DIM_M)

! Minimum/maximum solids velocity at inlet.  Also used in the iterative
! steps as the starting and ending velocities 
      DOUBLE PRECISION  MIN_VEL, MAX_VEL

      DOUBLE PRECISION  MINIPV, MAXIPV

      INTEGER :: OCCUPANTS
! jump_here

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL, EXTERNAL :: COMPARE 
!-----------------------------------------------


      CALL INIT_ERR_MSG("SET_BC_DEM_MI")


      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,2x,'DEM inlet count: ',I4)") DEM_BCMI


! Loop over BCs that flagged for DEM mass inflow.
      DO BCV_I = 1, DEM_BCMI

! Get the user defined BC ID.
         BCV = DEM_BC_MI_MAP(BCV_I)

         if(dFlag) write(*,"(2/,'Setting DEM_MI:',I3)") BCV_I

! The number of mass phases at this inlet.  While a system may be
! polydisperse, the inlet could consist of a single mass phase
         PHASE_CNT = 0
! The mass phase indices of incoming particles at this inlet
         PHASE_LIST(:) = -1
! The max diameter of incoming particles at this inlet
         MAX_DIA = ZERO

! Determine if the inlet is mono or polydisperse               
         DO M=1, DES_MMAX
            IF(SOLIDS_MODEL(M) /= 'DEM') CYCLE
            IF(BC_ROP_s(BCV,M) == UNDEFINED) CYCLE
            IF(COMPARE(BC_ROP_s(BCV,M),ZERO)) CYCLE
            PHASE_CNT = PHASE_CNT + 1
            PHASE_LIST(PHASE_CNT) = M
            MAX_DIA = MAX(MAX_DIA,DES_D_P0(M))
         ENDDO

! This subroutine determines the pattern that the particles will need to
! enter the system, if any. This routine only needs to be called if a
! run is new.  If a run is a RESTART_1, all of the setup information
! provided by this subroutine is will be obtained from the *_DES.RES file.
! This is done due to this routine's strong dependence on the 
! RANDOM_NUMBER() subroutine.
         IF(RUN_TYPE == 'NEW') THEN

            SELECT CASE (BC_PLANE(BCV))
            CASE('N','S'); CALL LAYOUT_DEM_MI_NS(BCV, BCV_I, MAX_DIA)
!            CASE('E','W'); CALL LAYOUT_DEM_MI_EW(BCV, BCV_I, MAX_DIA)
!            CASE('T','B'); CALL LAYOUT_DEM_MI_TB(BCV, BCV_I, MAX_DIA)
            END SELECT
         ENDIF


! Initialize
         MAX_VEL = ZERO
         NPMpSEC(:) = ZERO

! Calculate the individual velocities for each solid phase
         DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)

! Pull off the BC velocity normal to the flow plane.
            SELECT CASE(BC_PLANE(BCV))
            CASE('N','S'); VEL_TMP(M) = abs(BC_V_s(BCV,M))
            CASE('E','W'); VEL_TMP(M) = abs(BC_U_s(BCV,M))
            CASE('T','B'); VEL_TMP(M) = abs(BC_W_s(BCV,M))
            END SELECT

! Check for min/max inlet velocity
            MAX_VEL = MAX(ABS(VEL_TMP(M)), MAX_VEL)
! Calculate volumetric flow rate to convert to particle count. BC_AREA
! was already corrected for cut cells and velocity was recalculated
! to ensure user-specified mass or volumetric flow rates.
            VOL_FLOW = VEL_TMP(M) * BC_AREA(BCV) * BC_EP_S(BCV,M)
! Calculate the number of particles of mass phase M are injected per 
! second for each solid phase present at the boundary
            NPMpSEC(M) = VOL_FLOW / (PI/6.d0*DES_D_P0(M)**3)
! Write some debugging information if needed.
            if(dFlag) write(*,1100) M, VEL_TMP(M), NPMpSEC(M)
         ENDDO

! Total number of particles at BCV injected per second
         NPpSEC = sum(NPMpSEC)

 1100 FORMAT(/2x,'Conversion Info: Phase ',I2,/4x,'Velocity: ',g11.5,/ &
         4X,'NPMpSEC = ',F11.1)

         if(dFlag) write(*,"(/2x,'Max Velocity:',3x,g11.5)") MAX_VEL
         if(dFlag) write(*,"( 2x,'NPpSEC:',3x,F11.1)") NPpSEC

! The number of total particles per solid time step DTSOLID
         NPpDT = NPpSEC * DTSOLID

! Inject one particle every PI_FACTOR solids time steps.
         IF(NPpDT .LT. 1)THEN
            PI_COUNT(BCV_I) = 1
            PI_FACTOR(BCV_I) = FLOOR(real(1.d0/NPpDT))
! Inject PI_COUNT particles every soilds time step.
         ELSE
            PI_COUNT(BCV_I) = CEILING(real(NPpDT))
            PI_FACTOR(BCV_I) = 1
         ENDIF

         OCCUPANTS = DEM_MI(BCV_I)%OCCUPANTS

! Calculate the minimum inlet velocity. The cutoff is associated with
! square packing of disks on a plane.
         MINIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * dble(     &
            FLOOR( real(OCCUPANTS)/real(PI_COUNT(BCV_I)))))
! Calculate the velocity needed to ensure that half the inlet is free.
! Inlets with velocities greater than this value can be randomly seeded,
! otherwise, particles are seeded in according to the grid.
         MAXIPV = MAX_DIA/( DTSOLID*dble(PI_FACTOR(BCV_I)) * dble(     &
            FLOOR(CEILING(real(OCCUPANTS)/2.0)/real(PI_COUNT(BCV_I))))) 

         if(dFlag) write(*,"(/2x,'MaxIPV:',3x,g11.5)") MAXIPV
         if(dFlag) write(*,"( 2x,'MinIPV:',3x,g11.5)") MINIPV


         IF(MAX_VEL < MINIPV) THEN
            WRITE(ERR_MSG,1110) BCV, MAX_VEL, MINIPV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1110 FORMAT('Error 1110: Solids velocity for BC ',I3,' is too low ', &
         'to satisfy DEM',/'inlet packing restrictions. Potential ',  &
         'solutions:',//,' > If the BC velocities (BC_U_s, BC_V_s, ', &
         'BC_W_s) are defined, specify',/3x,'a larger value for the ',&
         'velocity normal to the flow plane.',//' > If MASSFLOW or ', &
         'VOLFLOW are defined, decrease the solids volume',/3x,       &
         'fraction to increase solids velocity.',//2x,'Max user-',    &
         'specified BC velocity:   ',g11.5,/2x,'Minimum required ',   &
         'solids Velocity: ',g11.5)


! Set all BC solids velocities to the largest velocity and recalculate
! BC_EP_s to determine the magnitude of the change.
         EPs_ERR = ZERO
         EPs_TMP = ZERO
         DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)
            EPs_TMP(M) = BC_EP_s(BCV,M) * (VEL_TMP(M) / MAX_VEL)
            EPs_ERR = EPs_ERR + (BC_EP_s(BCV,M) - EPs_TMP(M))

! Over-write the current BC value.
            SELECT CASE(BC_PLANE(BCV))
            CASE('N'); BC_V_s(BCV,M) =  VEL_TMP(M)
            CASE('S'); BC_V_s(BCV,M) = -VEL_TMP(M)
            CASE('E'); BC_U_s(BCV,M) =  VEL_TMP(M)
            CASE('W'); BC_U_s(BCV,M) = -VEL_TMP(M)
            CASE('T'); BC_W_s(BCV,M) =  VEL_TMP(M)
            CASE('B'); BC_W_s(BCV,M) = -VEL_TMP(M)
            END SELECT

         ENDDO

! If the net change in solids volume fraction is greatere than 0.01,
! flag this as an error and exit. >> Let the user fix the input.
         IF(.NOT.COMPARE(EPs_ERR,SMALL_NUMBER)) THEN
            IF(EPs_ERR > 0.01) THEN
               WRITE(ERR_MSG,1200) BCV, MAX_VEL, EPs_ERR
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               FATAL = .TRUE.

! Report the amount of changes imposed on the BC in setting a 
! uniform inlet velocity.
            ELSE
               WRITE(ERR_MSG,1205) BCV, MAX_VEL, EPs_ERR
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               FATAL = .FALSE.
            ENDIF

            WRITE(ERR_MSG, 1210)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            DO MM = 1, PHASE_CNT
            M = PHASE_LIST(MM)
               WRITE(ERR_MSG,1211) M, BC_EP_s(BCV,M), EPs_TMP(M), &
                  (BC_EP_s(BCV,M)-EPs_TMP(M))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDDO
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=FATAL)
         ENDIF


 1200 FORMAT('Error 1200: Unable to impose a uniform solids velocity ',&
         'on BC',I3,'.',/'Setting all solids to the highest velocity ',&
         'results in a large error',/'in the solids volume fraction. ',&
         'Please correct the mfix.dat file.',2/5X,'Max Inlet Velocity',&
         ':',1X,ES11.4,/5x,'Total BC_EP_s Error:',ES11.4)


 1205 FORMAT('Warning 1201: Uniform solids velocity imposed on BC',I3, &
         '.',2/,5X,'Uniform Inlet Velocity:',1X,ES11.4,/5X,'Total ',   &
         'BC_EP_s Error:',4X,ES11.4,/' ')

 1210 FORMAT(/,5X,'|',11('-'),'|',3(14('-'),'|'),/5X,'|',3X,'Phase',3X,&
         '|',4X,'BC_EP_s',3X,'|',2X,'Calculated',2X,'|',3X,'ABS ',   &
         'Error',2X,'|',/5X,'|',11('-'),'|',3(14('-'),'|'))

 1211 FORMAT(5X,'|',4X,I2,5X,'| ',1X,ES11.4,1X,'|',1X,ES11.4,2X,'|',   &
         1X,ES11.4,1X,' |',/5X,'|',11('-'),'|',3(14('-'),'|'))


! For polydisperse inlets, construct the DES_POLY_LAYOUT array
         IF(DEM_BC_POLY(BCV_I)) THEN
            RANGE_BOT = 1
            DO MM=1,PHASE_CNT - 1
               M = PHASE_LIST(MM)
               SCALED_VAL = dble(NUMFRAC_LIMIT)*(NPMpSEC(M)/NPpSEC)
               RANGE_TOP = FLOOR(SCALED_VAL) + (RANGE_BOT-1)
               DEM_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:RANGE_TOP) = M 
               RANGE_BOT = RANGE_TOP+1
            ENDDO

            M = PHASE_LIST(PHASE_CNT)
            DEM_BC_POLY_LAYOUT(BCV_I,RANGE_BOT:NUMFRAC_LIMIT) = M
! For monodisperse inlets, store the single mass phase used
         ELSE
            DEM_BC_POLY_LAYOUT(BCV_I,:) = PHASE_LIST(1)
         ENDIF


! Calculate des mass inlet time; time between injection.  If the run
! type is RESTART_1, DES_MI_TIME will be picked up from the restart file
! with an updated value.
         IF(RUN_TYPE == 'NEW') &
            DEM_MI_TIME(BCV_I) = TIME + dble(PI_FACTOR(BCV_I)) * DTSOLID


         WRITE(*,1000) BCV, NPpDT, PI_FACTOR(BCV_I),&
            PI_COUNT(BCV_I), DEM_MI_TIME(BCV_I)

! Flag that the inlet is polydisperse if true
         DEM_MI(BCV_I)%POLYDISPERSE = (PHASE_CNT > 1)


      ENDDO


 1000 FORMAT(2/,2X,'For mass inlet BC: ', I3,/,&
         4X,'No. particles injected per solids time step = ', ES15.8,/,&
         4X,'PI_FACTOR = ', I10,' PI_COUNT = ', I5,/,&
         4X,'start DES_MI_TIME = ', ES15.8)

      RETURN
      END SUBROUTINE SET_BC_DEM_MI




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LAYOUT_DEM_MI_NS                                        !
!                                                                      !
!  Purpose:  This routine determines the layout of the mass inlet as   !
!  either ordered or random based upon the inlet conditions.  This     !
!  routine also verifies that the specified inlet conditions for mass  !
!  or volumetric flow rates along with inlet size (length or area) and !
!  particle inlet velocity will work.  If not an error is flagged and  !
!  the program is exited.                                              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LAYOUT_DEM_MI_NS(BCV, BCV_I, MAX_DIA)

      use bc, only: BC_PLANE, BC_Y_s, BC_AREA
      use bc, only: BC_X_w, BC_X_e
      use bc, only: BC_Z_b, BC_Z_t

      use des_bc, only: DEM_MI

      use compar
      use geometry
      use indices

      use funits, only: DMP_LOG

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I
! Max diameter of incoming particles at bc
      DOUBLE PRECISION, INTENT(IN) :: MAX_DIA
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters.
      INTEGER LL, LC
! Indices for mapping to fluid grid.
      INTEGER IJK, I, J, K
! Local MESH indices
      INTEGER H, W
! Temp variable for double precision
      DOUBLE PRECISION :: TMP_DP
! Temporary variable for integers
      INTEGER :: TMP_INT

      INTEGER, allocatable :: MESH_H(:)
      INTEGER, allocatable :: MESH_W(:)

      DOUBLE PRECISION, allocatable :: MESH_P(:)
      DOUBLE PRECISION, allocatable :: MESH_Q(:)

      INTEGER, allocatable :: RAND_MAP(:)
      INTEGER, allocatable :: FULL_MAP(:,:)

! max number of partitions along length of inlet
      INTEGER :: WMAX, HMAX

! the length of each side of the inlet boundary
      DOUBLE PRECISION :: PLEN, QLEN

! Number of occupied mesh cells
      INTEGER :: OCCUPANTS

      DOUBLE PRECISION :: SHIFT, WINDOW, OFFSET

      LOGICAL, EXTERNAL :: COMPARE 

      LOGICAL, parameter :: setDBG = .FALSE. ! .FALSE.
      LOGICAL :: dFlag


      include 'function.inc'

!-----------------------------------------------

! Initialize the error manager.
      CALL INIT_ERR_MSG('LAYOUT_DEM_MI_NS')

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,'Building DEM_MI: ',I3)") BCV_I

! Store the index that maps back to the user input.
      DEM_MI(BCV_I)%BCV = BCV
      DEM_MI(BCV_I)%PLANE = BC_PLANE(BCV)

      OCCUPANTS = 0

! Calculate number of partitions in first direction.
      PLEN = BC_X_e(BCV) - BC_X_w(BCV)
      WMAX = FLOOR(real(PLEN/MAX_DIA))
      allocate( MESH_W(WMAX) )
      allocate( MESH_P(0:WMAX) )

! Calculate number of partitions in the second direction.
      QLEN = merge(BC_Z_t(BCV) - BC_Z_b(BCV), MAX_DIA, DO_K)
      HMAX = FLOOR(real(QLEN/MAX_DIA))
      allocate( MESH_H(HMAX) )
      allocate( MESH_Q(0:HMAX) )

! Allocate the full map.
      allocate( FULL_MAP(WMAX, HMAX))

! Set the value of the boundary condtion offset value used in the
! placement of new particles.
      CALL CALC_CELL_INTERSECT(ZERO, BC_Y_s(BCV), DY, JMAX, J)
      SHIFT = merge(-ONE, ONE, BC_PLANE(BCV) == 'N')
      DEM_MI(BCV_I)%OFFSET = BC_Y_s(BCV) + MAX_DIA*SHIFT
      DEM_MI(BCV_I)%L = J + int(SHIFT)
      if(dFlag) write(*,"(2x,'Offset: '3x,I4,3x,g11.5)") &
         DEM_MI(BCV_I)%L, DEM_MI(BCV_I)%OFFSET


! Dimension of grid cell for seeding particles; this may be larger than
! than the particle diameter but not smaller: 
      DEM_MI(BCV_I)%WINDOW = MIN(PLEN/WMAX, QLEN/HMAX)
      WINDOW = DEM_MI(BCV_I)%WINDOW
      if(dFlag) write(*,"(2x,'Windows size: ',g11.5)") WINDOW

! Setup the first direction.
      SHIFT = HALF*(PLEN - WMAX*WINDOW)
      MESH_P(0) = BC_X_w(BCV) + SHIFT
      if(dFlag) write(*,8005) 'P', SHIFT, 'P', MESH_P(0)
      DO LC=1,WMAX
         MESH_P(LC) = MESH_P(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_P(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(XMIN, SHIFT, DX, IMAX, MESH_W(LC))
         IF(dFlag)WRITE(*,8006) LC, 'W', MESH_W(LC), 'P', MESH_P(LC)
      ENDDO

! Setup the second direction.
      SHIFT = HALF*(QLEN - HMAX*WINDOW)
      MESH_Q(0) = BC_Z_b(BCV) + SHIFT
      if(dFlag) write(*,8005) 'Q',SHIFT, 'Q',MESH_Q(0)
      DO LC=1,HMAX
         MESH_Q(LC) = MESH_Q(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_Q(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DZ, KMAX, MESH_H(LC))
         IF(dFlag)WRITE(*,8006) LC, 'H', MESH_H(LC), 'Q', MESH_Q(LC)
      ENDDO


! Get the Jth index of the fluid cell
      CALL CALC_CELL_INTERSECT(ZERO, BC_Y_s(BCV), DY, JMAX, J)

! If the computationsl cell adjacent to the DEM_MI mesh cell is a 
! fluid cell and has not been cut, store the ID of the cell owner.
      DO H=1,HMAX
      DO W=1,WMAX
         I = MESH_W(W)
         K = MESH_H(H)
         FULL_MAP(W,H) = 0
         IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE
         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(.NOT.COMPARE(AXZ(IJK),DX(I)*DZ(K))) CYCLE
         OCCUPANTS = OCCUPANTS + 1
         FULL_MAP(W,H) = myPE+1
      ENDDO
      ENDDO

! Sync the full map across all ranks.
      CALL GLOBAL_ALL_SUM(OCCUPANTS)
      CALL GLOBAL_ALL_SUM(FULL_MAP)

! Throw an error and exit if there are no occupants.
      IF(OCCUPANTS == 0) THEN
         WRITE(ERR_MSG, 1100) BCV_I
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: No un-cut fluid cells adjacent to DEM_MI ',  &
         'staging area.',/'Unable to setup the discrete solids mass ', &
         'inlet for BC:',I3)

! Store the number of occupants.
      DEM_MI(BCV_I)%OCCUPANTS = OCCUPANTS

! Display the fill map if debugging
      IF(dFlag) THEN
         WRITE(*,"(2/,2x,'Displaying Fill Map:')")
         DO H=HMAX,1,-1
            WRITE(*,"(2x,'H =',I3)",advance='no')H
            DO W=1,WMAX
               IF(FULL_MAP(W,H) == 0) then
                  WRITE(*,"(' *')",advance='no')
               ELSE
                  WRITE(*,"(' .')",advance='no')
               ENDIF
            ENDDO
            WRITE(*,*)' '
         ENDDO
      ENDIF

! Construct an array of integers randomly ordered.
      if(dFLAG) write(*,"(2/,2x,'Building RAND_MAP:')")
      allocate( RAND_MAP(OCCUPANTS) )
      RAND_MAP = 0

! Only Rank 0 will randomize the layout and broadcast it to the other
! ranks. This will ensure that all ranks see the same layout.
      IF(myPE == 0) THEN
         LL = 1
         DO WHILE (RAND_MAP(OCCUPANTS) .EQ. 0)
            CALL RANDOM_NUMBER(TMP_DP)
            TMP_INT = CEILING(real(TMP_DP*dble(OCCUPANTS)))
            DO LC = 1, LL
              IF(TMP_INT .EQ. RAND_MAP(LC) )EXIT
              IF(LC .EQ. LL)THEN
                 if(dFlag) WRITE(*,"(4x,'LC:',I3,' : ',I3)") LC, TMP_INT
                 RAND_MAP(LC) = TMP_INT
                 LL = LL + 1
              ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL GLOBAL_ALL_SUM(RAND_MAP)

      allocate( DEM_MI(BCV_I)%W(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%P(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%H(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%Q(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%OWNER(OCCUPANTS) )

      if(dFlag) write(*,8010)
! Store the MI layout in the randomized order.
      LC = 0
      DO H=1,HMAX
      DO W=1,WMAX
         IF(FULL_MAP(W,H) == 0) CYCLE
         LC = LC + 1
         LL = RAND_MAP(LC)
         DEM_MI(BCV_I)%OWNER(LL) = FULL_MAP(W,H) - 1

         DEM_MI(BCV_I)%W(LL) = MESH_W(W)
         DEM_MI(BCV_I)%H(LL) = MESH_H(H)

         DEM_MI(BCV_I)%P(LL) = MESH_P(W)
         DEM_MI(BCV_I)%Q(LL) = MESH_Q(H)

         if(dFlag) write(*,8011) DEM_MI(BCV_I)%OWNER(LL), &
            DEM_MI(BCV_I)%W(LL), DEM_MI(BCV_I)%H(LL), DEM_MI(BCV_I)%L, &
            DEM_MI(BCV_I)%P(LL), DEM_MI(BCV_I)%Q(LL), DEM_MI(BCV_I)%OFFSET

      ENDDO
      ENDDO


 8010 FORMAT(2/,2x,'Storing DEM_MI data:',/4X,'OWNER',5X,'W',5X,'H',   &
         5X,'L',7X,'P',12X,'Q',12X,'R')
 8011 FORMAT(4x,I5,3(2X,I4),3(2x,g11.5))


      if(dFlag) write(*,"(2/,2x,'Inlet area sizes:')")
      if(dFlag) write(*,9000) 'mfix.dat: ', PLEN * QLEN
      if(dFlag) write(*,9000) 'BC_AREA:  ', BC_AREA(BCV)
      if(dFlag) write(*,9000) 'DEM_MI:   ', OCCUPANTS * (WINDOW**2)
 9000 FORMAT(2x,A,g11.5)

! House keeping.
      IF( allocated(MESH_H)) deallocate(MESH_H)
      IF( allocated(MESH_W)) deallocate(MESH_W)
      IF( allocated(MESH_P)) deallocate(MESH_P)
      IF( allocated(MESH_Q)) deallocate(MESH_Q)

      IF( allocated(RAND_MAP)) deallocate(RAND_MAP)
      IF( allocated(FULL_MAP)) deallocate(FULL_MAP)

      CALL FINL_ERR_MSG

      RETURN

 8005 FORMAT(2/,2x,'Building MESH_',A1,':',/4x'Shift:',f8.4,/4x,       &
         'MESH_',A1,'(0) = ',f8.4,/)

 8006 FORMAT(4x,'LC = ',I4,3x,A1,' =',I3,3x,A1,' =',f8.4)

      END SUBROUTINE LAYOUT_DEM_MI_NS





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM_MO                                           !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 23-Nov-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM_MO

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE des_bc
      USE bc
      USE discretelement 
      USE funits  
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV, BCV_I      ! BC loop counter
      INTEGER M, MM           ! Mass phase loop counter
      INTEGER HOLD, I         ! Dummy values
      INTEGER RANGE_TOP, RANGE_BOT ! Dummy values
      INTEGER PHASE_CNT        ! Number of solid phases at bc
      INTEGER PHASE_LIST(DES_MMAX) ! List of phases used in current bc

! the number of particles injected in a solids time step
      DOUBLE PRECISION NPMpSEC(DES_MMAX) ! For solid phase m
      DOUBLE PRECISION NPpSEC
      DOUBLE PRECISION NPpDT        ! Total for BC
      DOUBLE PRECISION SCALED_VAL
      DOUBLE PRECISION MAX_DIA ! Max diameter of incoming particles at bc

!-----------------------------------------------
!   External functions
!-----------------------------------------------
      LOGICAL, EXTERNAL :: COMPARE 
!-----------------------------------------------


! Check each discrete mass outlet for necessary data
      IF(DES_BCMO/=0)THEN
         BCV_I = 1
         DO BCV = 1, DIMENSION_BC 
           IF((DES_BC_DEFINED(BCV)) )THEN! .AND. &
!              (DES_BC_TYPE(BCV) == 'MASS_OUTFLOW'))THEN

               DES_BC_MO_ID(BCV_I) = BCV

               CALL DES_MO_CLASSIFY(BCV_I, BCV)


               BCV_I = BCV_I + 1 
            ENDIF
         ENDDO
      ENDIF


 1000 FORMAT(/5X,'For mass inlet BC: ', I3,/,&
         7X,'No. particles injected per solids time step = ', ES15.8,/,&
         7X,'PI_FACTOR = ', I10,' PI_COUNT = ', I5,/,&
         7X,'start DES_MI_TIME = ', ES15.8)

      RETURN
      END SUBROUTINE SET_BC_DEM_MO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MO_CLASSIFY                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date:  5-Oct-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MO_CLASSIFY(BCV_I, BCV)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE des_bc
      USE discretelement 
      USE funits  
      USE geometry
      USE param1
      USE physprop
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV_I, BCV
!-----------------------------------------------


      IF(DIMN == 2)THEN   ! 2D domain

! Check verticle mass outlet: 2D
! ----------------------------------------
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV))THEN
            DES_MO_X = .TRUE.
            IF (DES_BC_X_w(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'XW'
            IF (DES_BC_X_w(BCV) == XLENGTH) DES_MO_CLASS(BCV_I) = 'XE'
         ENDIF

! Check horizontal mass outlet: 2D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV))THEN
            DES_MO_Y = .TRUE.
            IF (DES_BC_Y_s(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'YS'
            IF (DES_BC_Y_s(BCV) == YLENGTH) DES_MO_CLASS(BCV_I) = 'YN'
         ENDIF

     ELSE   !  3D domain

! Check mass outlet on YZ face: 3D 
! ----------------------------------------
! see comments following the yz face section of DEM_MI_CLASSIFY
! for elucidation on this, and subsequent series of if statements
         IF(DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH &
         .OR. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

            DES_MO_X = .TRUE.
            IF (DES_BC_X_w(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'XW'
            IF (DES_BC_X_w(BCV) == XLENGTH) DES_MO_CLASS(BCV_I) = 'XE'

         ENDIF
! End check mass outlet on YZ face: 3D 


! Check mass inlet on XZ face: 3D 
! ----------------------------------------
         IF(DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Z_b(BCV) /= DES_BC_Z_t(BCV) &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Z_b(BCV) /= ZERO .AND. &
            DES_BC_Z_b(BCV) /= ZLENGTH)THEN

            DES_MO_Y = .TRUE.
            IF (DES_BC_Y_s(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'YS'
            IF (DES_BC_Y_s(BCV) == YLENGTH) DES_MO_CLASS(BCV_I) = 'YN'

         ENDIF
! End check mass outlet on XZ face: 3D 

! Check mass inlet on XY face: 3D 
! ----------------------------------------
         IF(DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) /= DES_BC_X_e(BCV) .AND. &
            DES_BC_Y_s(BCV) /= DES_BC_Y_n(BCV) &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_X_w(BCV) == DES_BC_X_e(BCV) .AND. &
            DES_BC_X_w(BCV) /= ZERO .AND. &
            DES_BC_X_w(BCV) /= XLENGTH &
         .OR. &
            DES_BC_Z_b(BCV) == DES_BC_Z_t(BCV) .AND. &
            DES_BC_Y_s(BCV) == DES_BC_Y_n(BCV) .AND. &
            DES_BC_Y_s(BCV) /= ZERO .AND. &
            DES_BC_Y_s(BCV) /= YLENGTH)THEN

            DES_MO_Z = .TRUE.
            IF (DES_BC_Z_b(BCV) == ZERO) DES_MO_CLASS(BCV_I) = 'ZB'
            IF (DES_BC_Z_b(BCV) == ZLENGTH) DES_MO_CLASS(BCV_I) = 'ZT'

         ENDIF
! End check mass outlet on XY face: 3D 

      ENDIF   ! endif dimn == 2

! Verify that an outlet is not on a face that is connected to a periodic
! boundary condition.  If so, write error message and exit.
! No Xew outlet with X direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'XW' .OR. DES_MO_CLASS(BCV_I) == 'XE') &
          .AND. DES_PERIODIC_WALLS_X) THEN
         IF (DMP_LOG) THEN
            WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_X'
            WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_X'
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF
! No Ysn outlet with Y direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'YS' .OR. DES_MO_CLASS(BCV_I) == 'YN') &
          .AND. DES_PERIODIC_WALLS_Y) THEN
         IF (DMP_LOG) THEN
            WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_Y'
            WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_Y'
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF
! No Zbt outlet with Z direction periodic walls
      IF((DES_MO_CLASS(BCV_I) == 'ZB' .OR. DES_MO_CLASS(BCV_I) == 'ZT') &
          .AND. DES_PERIODIC_WALLS_Z) THEN
         IF (DMP_LOG) THEN
            WRITE (UNIT_LOG, 1500) BCV, 'DES_PERIODIC_WALLS_Z'
            WRITE (*, 1500) BCV, 'DES_PERIODIC_WALLS_Z'
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

 1500 FORMAT(/1X,70('*')//,' From: DES_MO_CLASSIFY -',/10X,&
         'DEM outlets can not be placed on a periodic boundary.',/10X,&
         'Check DEM boundary condtion ',I2,' and ',A,'.',/1X,70('*')/)

      RETURN
      END SUBROUTINE DES_MO_CLASSIFY




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_DES_LE_BC                                        !
!                                                                      !
!  Purpose: Check/set parameters for DES Lees Edeards BC.              !
!                                                                      !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Comments: *** DES Lees Edwards BC funcionality has been lost. ***   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_LE_BC

      use discretelement
      use mpi_utility

      IMPLICIT NONE

! Lees Edwards BC functionality has been lost in current DEM code
      IF(DES_LE_BC) THEN
         IF (DES_CONTINUUM_COUPLED) THEN
            WRITE(UNIT_LOG, 1064)
             CALL MFIX_EXIT(myPE)
         ENDIF 
         IF (DES_NEIGHBOR_SEARCH .NE. 4) THEN
            WRITE(UNIT_LOG, 1060)
            CALL MFIX_EXIT(myPE)
         ENDIF
! not all possible shear directions are fully coded         
         IF (DIMN .EQ. 2) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY' .AND. &
               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX') THEN
               WRITE(UNIT_LOG, 1061)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF(DIMN.EQ.3) THEN
            IF(TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDY') THEN ! .AND. & 
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DUDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DVDZ' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDX' .AND. &
!               TRIM(DES_LE_SHEAR_DIR) .NE. 'DWDY') THEN
               WRITE(UNIT_LOG, 1062)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
         IF (DES_PERIODIC_WALLS) THEN
            DES_PERIODIC_WALLS = .FALSE.
            DES_PERIODIC_WALLS_X = .FALSE.
            DES_PERIODIC_WALLS_Y = .FALSE.
            DES_PERIODIC_WALLS_Z = .FALSE.            
            WRITE(UNIT_LOG, 1063)
            WRITE(*,1063)
         ENDIF
      ENDIF

      RETURN

 1060 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Only the grid based search option is allowed when using',&
         'using',/10X,'Lees & Edwards BC.',/1X,70('*')/)

 1061 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=2 shear options are DUDY or DVDX',/1X,70('*')/)

 1062 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'Invalid option for DES_LE_SHEAR_DIR in mfix.dat. When',/10X,&
         'DIMN=3 shear options are DUDY, DUDZ, DVDX, DVDZ, DWDX or',&
         'DWDY.',/1X,70('*')/)

 1063 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'WARNING: DES_PERIODIC_WALLS set to false when DES_LE_BC.',&
         /10X,'DES_LE_BC implies periodic walls, however, the ',&
         'periodicity is handled',/10X, 'independently of ',&
         'DES_PERIODIC_WALLS option and so it is shut off.',&
         /1X,70('*')/)

 1064 FORMAT(/1X,70('*')//' From: CHECK_DES_DATA',/' Message: ',&
         'DES_CONTINUUM_COUPLED cannot be true when using ',&
         'DES_LE_BC.',/1X,70('*')/)

      END SUBROUTINE CHECK_DES_LE_BC
