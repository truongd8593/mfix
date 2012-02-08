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
      INTEGER I, J, K, L, IJK, PC, SM_CELL 
      integer lface,lcurpar,lpip_all(0:numpes-1),lglobal_id 
      
      double precision :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ 
      
      INCLUDE 'function.inc'
!-----------------------------------------------
! pradeep moved from get_Data to here.
      CALL DES_ALLOCATE_ARRAYS
      CALL DES_INIT_ARRAYS
!pradeep moving cfassign and des_init_bc before reading the particle info 
      CALL CFASSIGN
! J.Musser
! Make the necessary calculations for the mass inflow/outflow boundary
! conditions.  DTSOLID is needed so call is made after cfassign.f
      CALL DES_INIT_BC
      call desgrid_init 
      call desmpi_init       

      CALL DES_WALLBC_PREPROCSSING 

      if(dmp_log.and.debug_des)write(unit_log,'(1X,A)')&
         '---------- START MAKE_ARRAYS_DES ---------->'
! J.Musser 
! If no particles are in the system then there is no need to read 
! particle_input.dat or call generate_particle_config. Note, if no 
! particles are in the system and no dem inlet is specified, then 
! the run will have already been aborted from checks conducted in 
! check_des_bc                     
      IF(RUN_TYPE == 'NEW' .and. particles /= 0) THEN ! Fresh run
         IF(.NOT.GENER_PART_CONFIG) THEN 
           call read_par_input
         ELSE
           call generate_particle_config
         ENDIF
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
! set the old values 
         omega_old(:,:)   = zero
         omega_new(:,:)   = zero  
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)
! Set an initial radius for reacting particles
         IF(ANY_DES_SPECIES_EQ) CORE_RAD(:) = DES_RADIUS(:)
      ELSEIF(RUN_TYPE == 'RESTART_1') THEN !  Read Restart
          
         call des_read_restart 
     
         if(dmp_log)write(unit_log,'(3X,A,G17.8)') 'DES_RES file read at Time= ', TIME
         imax_global_id = maxval(iglobal_id(1:pip))
         call global_all_max(imax_global_id)
! set the old values 
         omega_old(:,:)   = omega_new(:,:)
         des_pos_old(:,:) = des_pos_new(:,:)
         des_vel_old(:,:) = des_vel_new(:,:)
      ELSEIF (RUN_TYPE == 'RESTART_2') THEN 
         if(dmp_log)write(unit_log,'(3X,A)') 'Restart 2 is not implemented with DES'
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF(USE_COHESION .AND. VAN_DER_WAALS) THEN ! calculate these for all run types.
         SURFACE_ENERGY = HAMAKER_CONSTANT/&
                 (24d0*PI*VDW_INNER_CUTOFF*VDW_INNER_CUTOFF)
         WALL_SURFACE_ENERGY = WALL_HAMAKER_CONSTANT/&
                 (24d0*PI*WALL_VDW_INNER_CUTOFF*WALL_VDW_INNER_CUTOFF)
      END IF 
!set the global id for walls. this is required to handle particle-wall contact
      do lface = 1,dimn*2
         iglobal_id(max_pip+lface) = -lface 
      end do 
!set particle properties - this code moved from cfassign   
      DO L = 1, MAX_PIP
         IF(.NOT.PEA(L,1)) CYCLE  
         IF(PEA(L,4)) CYCLE  
         PVOL(L) = (4.0D0/3.0D0)*PI*DES_RADIUS(L)**3
         PMASS(L) = PVOL(L)*RO_SOL(L)
         OMOI(L) = 2.5D0/(PMASS(L)*DES_RADIUS(L)**2) !ONE OVER MOI
      ENDDO
 
      IF(CARTESIAN_GRID) CALL CG_DEL_OUTOFDOMAIN_PARTS
!the above call will delete the particles outside the domain.
!it will then re-arrange the arrays such that the active particles
!are in a block. Finally it will call particles_in_cell for
!other consistency requirements
!the above will work only for MPPIC as the cell information 
!is added to particles in generarate_particle_config itself. 
!For DEM particles, the cell information is not added and is done
!only during first call to particles_in_cell. Don't want to call 
!particles_in_cell here as that will also compute EP_G which might 
!become less that zero if the particles outside the domain are not 
!removed first 

!Pradeep: do_nsearch should be set before calling particle in cell  
      do_nsearch = .true.
      CALL PARTICLES_IN_CELL

! J.Musser
! Set initial conditions obtained from mfix.dat file. (ENERGY/SPECIES)
      IF(RUN_TYPE == 'NEW' .AND. DES_IC_EXIST) CALL DES_SET_IC

!rahul:
!If cut-cell then remove the particles that are outside of the 
!cut-cell faces. Do this after particles_in_cell so that the particles are 
!are already assigned grid id's
      IF(TRIM(RUN_TYPE) == 'RESTART_1'.AND. CARTESIAN_GRID) THEN 
         open(1000, file='parts_out.dat', form="formatted")
         
         DO L = 1, PIP 
            SM_CELL = 0.d0
            IF(.NOT.FLUID_AT(PIJK(L,4))) THEN 
               IF(SMALL_CELL_AT(PIJK(L,4))) SM_CELL = 1.d0
               WRITE(1000, '(10(2x,g17.8))') SM_CELL, (DES_POS_NEW(L, I), I=1,DIMN), (DES_VEL_NEW(L, I), I=1,DIMN), REAL(L)
            ENDIF
         ENDDO
         close(1000, status = 'keep')
      ENDIF


! Overrides initial particle velocity with velocities assigned from a
! Gaussian distribution based on usr specified standard deviation and
! mean
      IF(TRIM(RUN_TYPE) == 'NEW') THEN 
         IF(PVEL_StDev.GT.ZERO) CALL INIT_PARTICLES_JN
         !Rahul: I just want to visualize the initial configuration
         !through tecplot. so playing with this variable here. 
         !the following 3 lines shud be removed in the final version.
         !DES_OUTPUT_TYPE = 'TECPLOT'
         !CALL write_des_data 
         !DES_OUTPUT_TYPE = UNDEFINED_C

         !only leave this line to write out the initial distribution
        
         ! Pradeep: writeic is modified to MPI, so it can be called when required 
         !CALL write_des_data 
         CALL WRITEIC
         !writeic is not up2date for MPI verison. it is better
         !to write the initial configuration through write_des_data
         !DES_VEL_NEW(:,1) = 1.
         !DES_VEL_NEW(:,2) = 2.
         !CALL PARTICLES_IN_CELL
         
         !CALL WRITE_DES_DATA
         !CALL WRITE_VTK_FILE
         
         !stop
      ENDIF

      IF(MPPIC) THEN 
         DTPIC_CFL = LARGE_NUMBER 
         PC = 1 
         DO L = 1, MAX_PIP
            IF(PC.GT.PIP) EXIT
            IF(.NOT.PEA(L,1)) CYCLE 
            PC = PC+1
! skip ghost particles
            IF(PEA(L,4)) CYCLE 

            DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
            DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
            DTPIC_TMPZ = LARGE_NUMBER 
            IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)

            DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         ENDDO
         CALL global_all_max(DTPIC_CFL)

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)
         
         DTSOLID = DTPIC_MAX
         
         IF(myPE.eq.pe_IO) THEN 
            WRITE(*,'(A40, 2x, 2(2x,g17.8))') 'DTPIC BASED ON CFL AND TAUP = ', dtpic_cfl, dtpic_taup
            WRITE(*,'(A40, 2x, 2(2x,g17.8))') 'DTSOLID SET TO ', DTSOLID
         ENDIF

         !read(*,*) 
         CNP_ARRAY = ZERO 
      end IF

!      Stop 

      if(dmp_log.and.debug_des)write(unit_log,'(1X,A)')&
         '<---------- END MAKE_ARRAYS_DES ----------'

      RETURN

! Flag that file particle_input.dat is missing and exit         
 999  if(dmp_log)write(unit_log,"(/1X,70('*')//,A,/10X,A,/1X,70('*')/)")&
         ' From: MAKE_ARRAYS_DES -',&
         ' particle_input.dat file is missing.  Terminating run.'
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE MAKE_ARRAYS_DES 


!------------------------------------------------------------------------
! Subroutine       : read_par_input 
! Purpose          : reads the particle input and bcasts the particle to 
!                    respective processors
!------------------------------------------------------------------------
      subroutine read_par_input
            
      USE discretelement
      use funits 
      use compar
      use desmpi 
      use cdist 
      use mpi_utility

      implicit none 
      integer i,j,k,lcurpar
      integer lunit 
      character(30) lfilename 
            
! based on distributed or single IO open the input file      
      if(dmp_log)write(unit_log,'(3x,a,/,5x,a)') &
       'Reading particle configuration from the supplied input file'

      lunit = 10
      if (bdist_io) then 
         write(lfilename,'("particle_input_",I4.4,".dat")')mype
         open(unit=lunit, file=lfilename, form="formatted")
      else 
         if (mype .eq. pe_io) then  
            lfilename= "particle_input.dat"
            open(unit=lunit, file=lfilename, form="formatted")
         end if 
      end if  

! Read the file
! in the distribute IO first line of the file will be number of particles in that processor 
      if (bdist_io) then 
         read(lunit,*) pip
         do lcurpar = 1,pip
            pea(lcurpar,1) = .true.
            read (lunit,*) (des_pos_new(lcurpar,k),k=1,dimn),des_radius(lcurpar), & 
                            ro_sol(lcurpar),(des_vel_new(lcurpar,k),k=1,dimn)
         enddo
      else 
! Read into temporary variable and scatter 
         if (mype .eq. pe_io) then 
         allocate (dpar_pos(particles,dimn),dpar_vel(particles,dimn),dpar_rad(particles),dpar_den(particles))
         dpar_pos=0.0;dpar_vel=0.0;dpar_rad = 0.0;dpar_den = 0.0
         do lcurpar = 1, particles
            read (lunit,*) (dpar_pos(lcurpar,k),k=1,dimn),dpar_rad(lcurpar), & 
                            dpar_den(lcurpar),(dpar_vel(lcurpar,k),k=1,dimn)
         end do
         end if 
         call des_scatter_particle
         if(mype.eq.pe_io) deallocate (dpar_pos,dpar_vel,dpar_rad,dpar_den)
      end if 
      if(bdist_io .or. mype.eq.pe_io) close(lunit)
      return 

  999 if(dmp_log)write(unit_log,"(/1X,70('*')//,A,/10X,A,/1X,70('*')/)")&
         ' From: read_par_input -',&
         ' particle_input.dat file is missing.  Terminating run.'
      CALL MFIX_EXIT(myPE)
      end subroutine 
      
