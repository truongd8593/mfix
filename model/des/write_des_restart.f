!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!                                       
!  module name: des_write_restart        
!  purpose: writing des data for restart  
!                         
!  Author : Pradeep G      
!  Purpose : Reads either single restart file or multiple restart files     
!            (based on bdist_io) flag 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      subroutine des_write_restart

      use param1      
      use compar
      use discretelement
      use run
      use des_bc
      use des_rxns
      use des_thermo
      use desmpi 
      use machine 
      use cdist 
      use mpi_utility

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lboundnum 
      integer li,lj,lres_unit,lnext_rec,ltot_pip,ltor_dimn,lsize
      integer lproc,lparcnt,lglocnt,lgathercnts(0:numpes-1)
      logical lassoc
      character(30) lfilename

!-----------------------------------------------
! remove *********************
      ltor_dimn = 1 + (dimn-2)*2 

! Create the restart file name and open. A restart can either be
! distributed, where each process writes a restart file. Or, a restart
! can be serial, where a single restart file is written. 
      if (bdist_io) then 
         lres_unit = 901 
         write(lfilename,'(A,I4.4,A)') trim(run_name)//'_DES',mype,'.res'
         open (unit=lres_unit,file=lfilename,form='unformatted',status='unknown', &
               access='direct',recl=open_n1)
      else 
         if(mype.eq.pe_io) then 
            lres_unit = 901
            write(lfilename,'(A)') trim(run_name)//'_DES.res'
            open (unit=lres_unit,file=lfilename,form='unformatted',status='unknown', &
                  access='direct',recl=open_n1)
         end if 
      end if 

! Distributed restart files:
!---------------------------------------------------------------------//
      if (bdist_io) then

! Operational conditions:
!--------------------------------------------------------------------->>
! pip - Partiles on processor, including ghost particles
! ighost_cnt - Number of ghost particles on processor
! vtp_findex - Output file count for *.vtp type files
! tecplot_findex - Output file count for tecplot type files.
! dtsolid - Value of solids time step based on particle properties.
         write(lres_unit,rec=1) pip, ighost_cnt, vtp_findex, &
            tecplot_findex, dtsolid 
         lnext_rec = 6

! Particle properties:
!--------------------------------------------------------------------->>
! Particle positions: All x-coordinates, then all y-coordinates and for
! 3D simulations, all z-coordinates.
         do li = 1,dimn
            call des_writepar(lres_unit,des_pos_new(:,li),pip,lnext_rec)
         end do 
! Global id of all particles for parallel processing.
         call des_writepar(lres_unit,iglobal_id,pip,lnext_rec)
! Dynamic particle count elements:
! 2 - Classifies particles as 'new'
! 4 - Classifies particles as ghost particles
         do li = 2,4
            call des_writepar(lres_unit,pea(:,li),pip,lnext_rec)
         end do
! Particle velocities: All x-components, all y-components, and for 3D
! simulations, all z-components
         do li = 1,dimn
            call des_writepar(lres_unit,des_vel_new(:,li),pip,lnext_rec)
         end do
! Rotational velocities: 1 degree of fredom for 2D, 3 DOF for 3D
         do li = 1,ltor_dimn
            call des_writepar(lres_unit,omega_new(:,li),pip,lnext_rec)
         end do 
! Particle radi
         call des_writepar(lres_unit,des_radius,pip,lnext_rec)
! Particle densities
         call des_writepar(lres_unit,ro_sol,pip,lnext_rec)
! MPPIC Only: Statistical weight of each particle.
         IF(MPPIC) call des_writepar(lres_unit,des_stat_wt(:),pip, lnext_rec)

! Logicals indicating the restart is written containing energy equation
! and reaction data.
         write (lres_unit,rec=lnext_rec) DES_ENERGY_EQ, &
            ANY_DES_SPECIES_EQ, MAX_DES_NMAX
         lnext_rec=lnext_rec+1
! Write out particle temperatures.
         IF(DES_ENERGY_EQ) &
            CALL des_writepar(lres_unit,DES_T_s_NEW,pip,lnext_rec)
         IF(ANY_DES_SPECIES_EQ) THEN
! Solids phase index: This is needed for reacting flows as the diameter
! (shrinking_particle) or density (variable_density) reaction models
! prevent identifying particles from initial phase specifications.
            CALL des_writepar(lres_unit,PIJK(:,5),pip,lnext_rec)
! Write out species mass fractions.
            DO li = 1, MAX_DES_NMAX
               CALL des_writepar(lres_unit,DES_X_s(:,li),pip,lnext_rec)
            ENDDO
         ENDIF

! Neighbor/collision data:
!--------------------------------------------------------------------->>
! Number of neighbors
         call des_writepar(lres_unit,neighbours(:,1),pip,lnext_rec)
! Number of particles in contact with each particle.
         call des_writepar(lres_unit,pn(:,1),pip,lnext_rec)
! Actual neighbors, contact pairs, ...
         do li = 2,maxneighbors
            call des_writepar(lres_unit,neighbours(:,li),pip,lnext_rec)
            call des_writepar(lres_unit,pn(:,li),pip,lnext_rec)
            call des_writepar(lres_unit,pv(:,li),pip,lnext_rec)
! Accumulated tangential displacement that occured during collision
            do lj = 1,dimn  
               call des_writepar(lres_unit,pft(:,li,lj),pip,lnext_rec)
            end do 
         end do 

! Mass inlet/outlet parameters:
!--------------------------------------------------------------------->>
         if(des_mi)then

!  *******************************************************************
!  * Warning: DES_MI_TIME is an array allocated to des_bcmi.         *
!  *          Could this be a problem?                               *
!  *******************************************************************
! Time to seed next particle.
            write (lres_unit,rec=lnext_rec) des_mi_time
            lnext_rec=lnext_rec+1

! For each mass inlet, store:
!   mi_factor - number of inlet grid cells containing particles
!   mi_window - grid cell length (> diameter)
!   particle_plcmnt - indicates seeding approach (ordered vs. random)
            do lboundnum =1, des_bcmi
               write (lres_unit,rec=lnext_rec) mi_factor(lboundnum), &
                  mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 

               if(particle_plcmnt(lboundnum) == 'ORDR')then
! For ordered inlets, store:
!   mi_order pointer is 'associated', and if so, store:
!     number of grid cells
!     mi_order - order to populate grid cells with incoming particles
                  lassoc = associated(mi_order(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(mi_order(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
!     number of cells in i-direction
!     number of cells in i-direction
!     i indices of mass inlet grid
                  lassoc = associated(i_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(i_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
!   j_of_mi is 'associated', and if so, store:
!     number of cells in j-direction
!     j indices of mass inlet grid
                  lassoc = associated(j_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(j_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif

! Serial restart file:
!---------------------------------------------------------------------//
       else
! Root process is informed of the total number of (real) particles.
         lglocnt = 10
         lparcnt = pip - ighost_cnt
         call global_sum(lparcnt,lglocnt)

! Allocate variables for gather.
! double precision arryas:
         allocate (dprocbuf(lparcnt)) ! local process particle count
         allocate (drootbuf(lglocnt)) ! global particle count (root)
! integer arrays:
         allocate (iprocbuf(lparcnt)) ! local process particle count
         allocate (irootbuf(lglocnt)) ! global particle count (root)

! Construct an array for the Root process that states the number of
! (real) particles on each process.
         igath_sendcnt = lparcnt 
         lgathercnts = 0
         lgathercnts(mype) = lparcnt
         call global_sum(lgathercnts,igathercnts)

! Calculate the displacements for each process in the global array.
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 

! Store data in restart file:

! Operational conditions:
!--------------------------------------------------------------------->>
! 1) total number of real particles
! 2) Output file count for *.vtp type files
! 3) Output file count for tecplot type files.
! 4) Value of solids time step based on particle properties.
         if(mype.eq.pe_io) write(lres_unit,rec=1) lglocnt, vtp_findex, &
            tecplot_findex,dtsolid 
         lnext_rec = 6

! Position data (one coordinate at a time).
         do li = 1,dimn
            call des_gatherwrite(lres_unit,des_pos_new(:,li),lglocnt,lnext_rec)
         end do
! Global particle index numbers.
         call des_gatherwrite(lres_unit,iglobal_id,lglocnt,lnext_rec)
! Dynamic particle count elements:
! 2 - Classifies particles as 'new'
! 4 - Classifies particles as ghost particles
         do li = 2,4
            call des_gatherwrite(lres_unit,pea(:,li),lglocnt,lnext_rec)
         end do 
! Translational velocities (one coordinate at a time).
         do li = 1,dimn
            call des_gatherwrite(lres_unit,des_vel_new(:,li),lglocnt,lnext_rec)
         end do 
! Rotational velocity.
         do li = 1,ltor_dimn
            call des_gatherwrite(lres_unit,omega_new(:,li),lglocnt,lnext_rec)
         end do 
! Particle radi.
         call des_gatherwrite(lres_unit,des_radius,lglocnt,lnext_rec)
! Particle densities.
         call des_gatherwrite(lres_unit,ro_sol,lglocnt,lnext_rec)

! MPPIC Only: Statistical weight of each particle.
         IF(MPPIC) call des_gatherwrite(lres_unit,des_stat_wt,lglocnt,lnext_rec)

! Write out particle temperatures.
         IF(DES_ENERGY_EQ) call des_gatherwrite(lres_unit, &
            DES_T_s_NEW,lglocnt,lnext_rec)
! Solids phase index: This is needed for reacting flows as the diameter
! (shrinking_particle) or density (variable_density) reaction models
! prevent identifying particles from initial phase specifications.
         IF(ANY_DES_SPECIES_EQ) THEN
            call des_gatherwrite(lres_unit, &
               PIJK(:,5),lglocnt,lnext_rec)
! Write out species mass fractions.
            DO li = 1, MAX_DES_NMAX
               call des_gatherwrite(lres_unit, &
                  DES_X_s(:,li),lglocnt,lnext_rec)
            ENDDO
         ENDIF

! Neighbor/collision data:
!--------------------------------------------------------------------->>
! Number of neighbors
         call des_gatherwrite(lres_unit,neighbours(:,1),lglocnt,lnext_rec)
! Number of particles in contact with each particle.
         call des_gatherwrite(lres_unit,pn(:,1),lglocnt,lnext_rec)

! *** Use global particle IDs instead of local IDs for neighbours.
         do li = 2,maxneighbors
! Global IDs of neighbor particles.
            call des_gatherwrite(lres_unit,neighbours(:,li),lglocnt,lnext_rec,ploc2glb=.true.)
! Neighbor contact data by global ID.
            call des_gatherwrite(lres_unit,pn(:,li),lglocnt,lnext_rec,ploc2glb=.true.)
            call des_gatherwrite(lres_unit,pv(:,li),lglocnt,lnext_rec)
            do lj = 1,dimn
! Accumulated tangential displacement that occured during collision
               call des_gatherwrite(lres_unit,pft(:,li,lj),lglocnt,lnext_rec)
            end do 
         end do

         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf)


! Mass inlet/outlet parameters:
!--------------------------------------------------------------------->>
         if(des_mi .and. mype .eq. pe_io)then

!  *******************************************************************
!  * Warning: DES_MI_TIME is an array allocated to des_bcmi.         *
!  *          Could this be a problem?                               *
!  *******************************************************************
! Time to seed next particle.
            write (lres_unit,rec=lnext_rec) des_mi_time
            lnext_rec=lnext_rec+1

! For each mass inlet, store:
!   mi_factor - number of inlet grid cells containing particles
!   mi_window - grid cell length (> diameter)
!   particle_plcmnt - indicates seeding approach (ordered vs. random)
            do lboundnum =1, des_bcmi
               write (lres_unit,rec=lnext_rec) mi_factor(lboundnum), &
                  mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 

               if(particle_plcmnt(lboundnum) == 'ORDR')then
! For ordered inlets, store:
!   mi_order pointer is 'associated', and if so, store:
!     number of grid cells
!     mi_order - order to populate grid cells with incoming particles
                  lassoc = associated(mi_order(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec)lassoc
                  lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(mi_order(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize
                     lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit, &
                        mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
!   i_of_mi is 'associated', and if so, store:
!     number of cells in i-direction
!     i indices of mass inlet grid
                  lassoc = associated(i_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc
                  lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(i_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize
                     lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit, &
                        i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
!   j_of_mi is 'associated', and if so, store:
!     number of cells in j-direction
!     j indices of mass inlet grid
                  lassoc = associated(j_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc
                  lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(j_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize
                     lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit, &
                        j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif
      end if 

      if(bdist_io .or.mype .eq. pe_io) close(lres_unit)
       
      end subroutine des_write_restart 

      
