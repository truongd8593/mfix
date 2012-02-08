!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MAKE_ARRAYS_DES                                        C
!  Purpose: DES - allocating DES arrays                                
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Added some calls that are necessary if INTERPOLATION IS ON C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE MAKE_ARRAYS_DES

!-----------------------------------------------
! Modules
!-----------------------------------------------            
      USE param1
      USE funits
      USE run
      USE compar      
      USE discretelement
      USE cutcell 
      use desmpi 
      use mpi_utility
      USE geometry 
      USE des_ic
      USE des_rxns

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER :: I, J, K, L, IJK, PC, SM_CELL 
      INTEGER :: lface,lcurpar,lpip_all(0:numpes-1),lglobal_id 
      
! MPPIC related quantities      
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------      
      INCLUDE 'function.inc'
!-----------------------------------------------


      CALL DES_ALLOCATE_ARRAYS
      CALL DES_INIT_ARRAYS

! cfassign and des_init_bc called before reading the particle info 
      CALL CFASSIGN

! Make the necessary calculations for the mass inflow/outflow boundary
! conditions.  DTSOLID is needed so call is made after cfassign.f
      CALL DES_INIT_BC

! parallelization: desmpi_init needs to be called after des_init_bc
! since it relies on setting/checking of des_mio
      call desgrid_init 
      call desmpi_init       

      CALL DES_WALLBC_PREPROCSSING 

      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '---------- START MAKE_ARRAYS_DES ---------->'

! If no particles are in the system then there is no need to read 
! particle_input.dat or call generate_particle_config. Note, if no 
! particles are in the system and no dem inlet is specified, then 
! the run will have already been aborted from checks conducted in 
! check_des_bc                     
      IF(RUN_TYPE == 'NEW' .and. particles /= 0) THEN ! Fresh run

         IF(.NOT.GENER_PART_CONFIG) THEN 
            CALL READ_PAR_INPUT
         ELSE
            CALL GENERATE_PARTICLE_CONFIG
         ENDIF
         
! Further initialization now the particles have been specified         
! for run type new set the global id for the particles and set the ghost cnt
         ighost_cnt = 0
         lpip_all = 0 
         lpip_all(mype) = pip 
         call global_all_sum(lpip_all)
         lglobal_id = sum(lpip_all(0:mype-1))
         do lcurpar  = 1,pip 
            lglobal_id = lglobal_id + 1
            iglobal_id(lcurpar) = lglobal_id 
         end do 
         imax_global_id = iglobal_id(pip)
         call global_all_max(imax_global_id)

! setting the old values 
         omega_old(:,:)   = zero
         omega_new(:,:)   = zero  
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)
         DES_VEL_OOLD(:,:) = DES_VEL_NEW(:,:)
! setting an initial radius for reacting particles
         IF(ANY_DES_SPECIES_EQ) CORE_RAD(:) = DES_RADIUS(:)

! overriding initial particle velocities with velocities assigned from a
! Gaussian distribution based on usr specified standard deviation and
! mean
         IF(PVEL_StDev > ZERO) CALL INIT_PARTICLES_JN         

      ELSEIF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
          
         call des_read_restart 
     
         IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A,G17.8)') &
            'DES_RES file read at Time= ', TIME
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)

! setting the old values : what about des_vel_oold?
         omega_old(:,:)   = omega_new(:,:)
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)
      ELSEIF (RUN_TYPE == 'RESTART_2') THEN 
         IF(DMP_LOG) WRITE(UNIT_LOG,'(3X,A)') &
            'Restart 2 is not implemented with DES'
         CALL MFIX_EXIT(myPE)
      ENDIF

! setting the global id for walls. this is required to handle 
! particle-wall contact
      DO lface = 1,dimn*2
         iglobal_id(max_pip+lface) = -lface 
      ENDDO

! setting additional particle properties now that the particles
! have been identified
      DO L = 1, MAX_PIP
! If RESTART_1 is being used with DEM inlets/outlets, then it is possible
! that the particle arrays have indices without data (without particles).
! Skip 'empty' locations when populating the particle property arrays.      
         IF(.NOT.PEA(L,1)) CYCLE  
         IF(PEA(L,4)) CYCLE  
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
! the following is used aid visualization of mixing but can be employed
! for other purposes if desired
         MARK_PART(L) = 1
         IF(DES_POS_NEW(L,2).LE.YLENGTH/2.d0) MARK_PART(L) = 0         
      ENDDO



! This call will delete the particles outside the domain. It will then
! re-arrange the arrays such that the active particles are in a block.
! This works only for MPPIC as the cell information is added to 
! particles in generate_particle_config_mppic itself. For DEM particles,
! the cell information is not added in generate_particle-config but is
! done only during first call to particles_in_cell. 
! Call this before particles_in_cell as EP_G is computed which 
! might become less than zero if the particles outside the domain 
! are not removed first       
      IF(CARTESIAN_GRID) CALL CG_DEL_OUTOFDOMAIN_PARTS

! do_nsearch should be set before calling particle in cell  
      DO_NSEARCH =.TRUE.
      CALL PARTICLES_IN_CELL

! Set initial conditions obtained from mfix.dat file. (ENERGY/SPECIES)
      IF(RUN_TYPE == 'NEW' .AND. DES_IC_EXIST) CALL DES_SET_IC

! If cut-cell then remove the particles that are outside of the cut-cell
! faces.
! Call this after particles_in_cell so that the particles are already
! assigned grid id's.  I thought mppic sets grid id's in 
      IF(TRIM(RUN_TYPE) == 'RESTART_1'.AND. CARTESIAN_GRID) THEN 
         open(1000, file='parts_out.dat', form="formatted")
         
         DO L = 1, PIP 
            SM_CELL = 0.d0
            IF(.NOT.FLUID_AT(PIJK(L,4))) THEN 
               IF(SMALL_CELL_AT(PIJK(L,4))) SM_CELL = 1.d0
               WRITE(1000, '(10(2x,g17.8))')&
                  SM_CELL, (DES_POS_NEW(L, I), I=1,DIMN),&
                  (DES_VEL_NEW(L, I), I=1,DIMN), REAL(L)
            ENDIF
         ENDDO
         close(1000, status = 'keep')
      ENDIF

      IF (TRIM(RUN_TYPE) == 'NEW') THEN
         CALL WRITEIC
      ENDIF

      IF(MPPIC) THEN 
         DTPIC_CFL = LARGE_NUMBER 
         PC = 1 
         DO L = 1, MAX_PIP
            IF(PC.GT.PIP) EXIT
            IF(.NOT.PEA(L,1)) CYCLE 
            PC = PC+1
            IF(PEA(L,4)) CYCLE 

            DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/&
               (ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
            DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/&
               (ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
            DTPIC_TMPZ = LARGE_NUMBER 
            IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/&
               (ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)

            DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         ENDDO
         CALL global_all_max(DTPIC_CFL)

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)
         
         DTSOLID = DTPIC_MAX
         
         IF(DMP_LOG) THEN 
            WRITE(*,'(A40, 2x, 2(2x,g17.8))') &
               'DTPIC BASED ON CFL AND TAUP = ', dtpic_cfl, dtpic_taup
            WRITE(*,'(A40, 2x, 2(2x,g17.8))') 'DTSOLID SET TO ', DTSOLID
         ENDIF

         CNP_ARRAY = ZERO 
      ENDIF


      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(1X,A)')&
         '<---------- END MAKE_ARRAYS_DES ----------'

      RETURN

      END SUBROUTINE MAKE_ARRAYS_DES 


!------------------------------------------------------------------------
! Subroutine       : read_par_input 
! Purpose          : reads the particle input and broadcasts the particle  
!                    to respective processors
!------------------------------------------------------------------------

      subroutine read_par_input
            
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      use funits 
      use compar
      use desmpi 
      use cdist 
      use mpi_utility

      implicit none 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices      
      integer :: i,j,k
! index of particle
      INTEGER :: lcurpar
! local unit      
      integer :: lunit 
! local filename      
      character(30) lfilename 
!-----------------------------------------------



      IF(DMP_LOG) WRITE(UNIT_LOG,'(3x,a,/,5x,a)') &
       'Reading particle configuration from the supplied input file'

      lunit = 10

! based on distributed or single IO open the input file
      IF (bdist_io) THEN
         write(lfilename,'("particle_input_",I4.4,".dat")')mype
         open(unit=lunit, file=lfilename, form="formatted")
      ELSE
         IF (mype .eq. pe_io) THEN
            lfilename= "particle_input.dat"
            open(unit=lunit, file=lfilename, form="formatted")
         ENDIF 
      ENDIF  

! Read the file
! in the distribute IO first line of the file will be number of particles in that processor 
      IF (bdist_io) then 
         read(lunit,*) pip
         DO lcurpar = 1,pip
            pea(lcurpar,1) = .true.
            read (lunit,*) (des_pos_new(lcurpar,k),k=1,dimn),&
               des_radius(lcurpar), ro_sol(lcurpar),&
               (des_vel_new(lcurpar,k),k=1,dimn)
         ENDDO

! Read into temporary variable and scatter 
      ELSE
         IF (mype .eq. pe_io) THEN
! temporary variable                 
            ALLOCATE (dpar_pos(particles,dimn))
            ALLOCATE (dpar_vel(particles,dimn))
            ALLOCATE (dpar_rad(particles))
            ALLOCATE (dpar_den(particles))
! Initialize            
            dpar_pos=0.0
            dpar_vel=0.0
            dpar_rad=0.0
            dpar_den = 0.0
            DO lcurpar = 1, particles
               read (lunit,*) (dpar_pos(lcurpar,k),k=1,dimn),&
                            dpar_rad(lcurpar), & 
                            dpar_den(lcurpar),(dpar_vel(lcurpar,k),k=1,dimn)
            ENDDO
         ENDIF 
         call des_scatter_particle
         IF(mype.eq.pe_io) deallocate (dpar_pos,dpar_vel,dpar_rad,dpar_den)
      ENDIF 

      IF(bdist_io .or. mype.eq.pe_io) close(lunit)

      RETURN

  999 IF(dmp_log)write(unit_log,"(/1X,70('*')//,A,/10X,A,/1X,70('*')/)")&
         ' From: read_par_input -',&
         ' particle_input.dat file is missing.  Terminating run.'
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE READ_PAR_INPUT
      
