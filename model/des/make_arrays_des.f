!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays 
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MAKE_ARRAYS_DES

      USE param1
      USE funits
      USE run
      USE compar
      USE discretelement
      USE cutcell
      use desmpi
      use mpi_utility
      USE geometry
      USE des_rxns
      USE des_thermo
      USE des_stl_functions

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, IJK, PC, SM_CELL
      INTEGER :: lface, lcurpar, lpip_all(0:numpes-1), lglobal_id  , lparcnt

      INTEGER :: FACTOR

! MPPIC related quantities
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'

! cfassign and des_init_bc called before reading the particle info
      CALL CFASSIGN

! parallelization: desmpi_init needs to be called after des_init_bc
! since it relies on setting/checking of des_mio
      call desgrid_init
      call desmpi_init

! Set the initial particle data.
      IF(RUN_TYPE == 'NEW' .and. PARTICLES /= 0) THEN
         IF(GENER_PART_CONFIG) THEN
            CALL COPY_PARTICLE_CONFIG_FROMLISTS
         ELSE
            CALL READ_PAR_INPUT
         ENDIF

! Set the global ID for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0
         lpip_all(mype) = pip
         call global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         imax_global_id = 0
         do lcurpar  = 1,pip
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id
            imax_global_id = iglobal_id(pip)
         end do
         call global_all_max(imax_global_id)

! Initialize old values
         omega_old(:,:)   = zero
         omega_new(:,:)   = zero
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)

! Read the restart file.
      ELSEIF(RUN_TYPE == 'RESTART_1') THEN

         CALL READ_RES0_DES
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)

! Initizlie the old values.
         omega_old(:,:)   = omega_new(:,:)
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)
         IF(ENERGY_EQ) DES_T_s_OLD(:) = DES_T_s_NEW(:)

      ELSE

         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1100 FORMAT('Error 1100: Unsupported RUN_TYPE for DES.')

      ENDIF

! setting the global id for walls. this is required to handle
! particle-wall contact
      DO lface = 1, merge(4,6,DO_K)
         iglobal_id(max_pip+lface) = -lface
      ENDDO

! setting additional particle properties now that the particles
! have been identified
      DO L = 1, MAX_PIP
! Skip 'empty' locations when populating the particle property arrays.
         IF(.NOT.PEA(L,1)) CYCLE
         IF(PEA(L,4)) CYCLE
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
! the following is used aid visualization of mixing but can be employed
! for other purposes if desired
         MARK_PART(L) = 1
         IF(DES_POS_NEW(2,L).LE.YLENGTH/2.d0) MARK_PART(L) = 0
      ENDDO

      CALL SET_PHASE_INDEX
      CALL INIT_PARTICLES_IN_CELL

! do_nsearch should be set before calling particle in cell
      DO_NSEARCH =.TRUE.
      CALL PARTICLES_IN_CELL

      IF(DEM_SOLIDS) THEN
         CALL NEIGHBOUR
         CALL INIT_SETTLING_DEM
      ENDIF


      IF(MPPIC) CALL CALC_DTPIC


      RETURN
      END SUBROUTINE MAKE_ARRAYS_DES
