!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!                                                                    
!  module name: des_read_restart                                      
!  purpose: writing des data for restart                               
!                                                                      
!  Author  : Pradeep G                                                  
!  Purpose : Reads either single restart file or multiple restart files     
!            (based on bdist_io) flag 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      subroutine des_read_restart
!-----------------------------------------------
! Modules
!-----------------------------------------------
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
      integer lboundnum,lsize     
      logical lassoc   
      integer li,lj,lres_unit,lnext_rec,ltot_pip,ltor_dimn
      integer lproc,lparcnt,lglocnt,lscattercnts(0:numpes-1)
      integer,dimension(:),allocatable::ltemparr
      character(30) lfilename
      integer lmax_pip 
!-----------------------------------------------

      ltor_dimn = 1 + (dimn-2)*2 
 
! Create the restart file name and open. A restart can either be
! distributed, where each process writes a restart file. Or, a restart
! can be serial, where a single restart file is written. 
      lres_unit = 901 
      if (bdist_IO) then 
         write(lfilename,'(A,I4.4,A)') trim(run_name)//'_DES',mype,'.res'
         open (unit=lres_unit,file=lfilename,form='unformatted',status='old', &
               access='direct',recl=open_n1)
      else 
         if(mype.eq.pe_io) then 
            write(lfilename,'(A)') trim(run_name)//'_DES.res'
            open (unit=lres_unit,file=lfilename,form='unformatted',status='old', &
                  access='direct',recl=open_n1)
         end if 
      end if 


! Distributed restart files:
!---------------------------------------------------------------------//
      if (bdist_IO) then 

! Operational conditions:
!--------------------------------------------------------------------->>
! pip - Partiles on processor, including ghost particles
! ighost_cnt - Number of ghost particles on processor
! vtp_findex - Output file count for *.vtp type files
! tecplot_findex - Output file count for tecplot type files.
! dtsolid - Value of solids time step based on particle properties.
         read(lres_unit,rec=1) pip, ighost_cnt, vtp_findex, &
            tecplot_findex, dtsolid 

         pea(1:pip,1) = .true.
         if (pip .gt. max_pip) then 
            write(*,*) "From des_read_restart:"
            write(*,*) "Error: The pip is grater than current max_pip"
            write(*,*) "pip=" ,pip,"; max_pip =", max_pip 
            call mfix_exit(mype)
         end if 

         lnext_rec = 6

! Particle properties:
!--------------------------------------------------------------------->>
! Particle positions: All x-coordinates, then all y-coordinates and for
! 3D simulations, all z-coordinates.
         do li = 1,dimn
            call in_bin_512(lres_unit,des_pos_new(:,li),pip,lnext_rec)
         end do 
! Global id of all particles for parallel processing.
         call in_bin_512i(lres_unit,iglobal_id,pip,lnext_rec)
         allocate(ltemparr(pip))

! Dynamic particle count elements:
! 2 - Classifies particles as 'new'
! 4 - Classifies particles as ghost particles
         do li = 2,4
            call in_bin_512i(lres_unit,ltemparr,pip,lnext_rec)
            do lj =1,pip
               if(ltemparr(lj).eq.1) then 
                  pea(lj,li)= .true.
               else 
                  pea(lj,li)= .false.
               end if 
            end do 
         end do 
         deallocate(ltemparr)

! Particle velocities: All x-components, all y-components, and for 3D
! simulations, all z-components
         do li = 1,dimn
            call in_bin_512(lres_unit,des_vel_new(:,li),pip,lnext_rec)
         end do 

! Rotational velocities: 1 degree of fredom for 2D, 3 DOF for 3D
         do li = 1,ltor_dimn
            call in_bin_512(lres_unit,omega_new(:,li),pip,lnext_rec)
         end do
! Particle radi
         call in_bin_512(lres_unit,des_radius,pip,lnext_rec)
! Particle densities
         call in_bin_512(lres_unit,ro_sol,pip,lnext_rec)
! MPPIC Only: Statistical weight of each particle.
         IF(MPPIC) call in_bin_512(lres_unit,des_stat_wt,pip,lnext_rec)

! Particle temperatures.
         IF(DES_ENERGY_EQ) &
            CALL in_bin_512(lres_unit,DES_T_s_NEW,pip,lnext_rec)
         IF(ANY_DES_SPECIES_EQ) THEN
! Solids phase index: This is needed for reacting flows as the diameter
! (shrinking_particle) or density (variable_density) reaction models
! prevent identifying particles from initial phase specifications.
            call in_bin_512i(lres_unit,PIJK(:,5),pip,lnext_rec)
! Species mass fractions.
            DO li = 1, MAX_DES_NMAX
               CALL in_bin_512(lres_unit,DES_X_s(:,li),pip,lnext_rec)
            ENDDO
         ENDIF

! Neighbor/collision data:
!--------------------------------------------------------------------->>
! Number of neighbors
         call in_bin_512i(lres_unit,neighbours(:,1),pip,lnext_rec)
! Number of particles in contact with each particle.
         call in_bin_512i(lres_unit,pn(:,1),pip,lnext_rec)
! Actual neighbors, contact pairs, ...
         do li = 2,maxneighbors
            call in_bin_512i(lres_unit,neighbours(:,li),pip,lnext_rec)
            call in_bin_512i(lres_unit,pn(:,li),pip,lnext_rec)
            call in_bin_512i(lres_unit,pv(:,li),pip,lnext_rec)
! Accumulated tangential displacement that occured during collision
            do lj = 1,dimn  
               call in_bin_512(lres_unit,pft(:,li,lj),pip,lnext_rec)
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
            read (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
! For each mass inlet, store:
            do lboundnum =1, des_bcmi
!   mi_factor - number of inlet grid cells containing particles
!   mi_window - grid cell length (> diameter)
!   particle_plcmnt - indicates seeding approach (ordered vs. random)
               read (lres_unit,rec=lnext_rec) mi_factor(lboundnum),mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 
! For ordered inlets:
               if(particle_plcmnt(lboundnum) == 'ORDR')then
!   mi_order pointer is 'associated', and if so, store:
!     number of grid cells
!     mi_order - order to populate grid cells with incoming particles
                  read (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(mi_order(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
!   i_of_mi is associated
!     number of cells in i-direction
!     i indices of mass inlet grid
                  read (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(i_of_mi(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
!   j_of_mi is associated
!     number of cells in j-direction
!     j indices of mass inlet grid
                  read (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(j_of_mi(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif

! Serial restart file:
!---------------------------------------------------------------------//
      else   ! else branch if not bdist_IO

         lglocnt = 10

! Operational conditions:
!--------------------------------------------------------------------->>
! 1) total number of real particles
! 2) Output file count for *.vtp type files
! 3) Output file count for tecplot type files.
! 4) Value of solids time step based on particle properties.
         if(mype.eq.pe_io) read(lres_unit,rec=1) lglocnt,vtp_findex,tecplot_findex,dtsolid 
         call bcast(vtp_findex)
         call bcast(tecplot_findex)
         call bcast(dtsolid)

         allocate(dpar_pos(lglocnt,dimn))
         allocate(irestartmap(lglocnt))

! for single file read the particle info first and then read the global
! id and position and devlop the mapping between the reading order 
! and local id and proc number 

         if(mype.eq.pe_io) then 
            lnext_rec = 6
! Position data (one coordinate at a time).
            do li = 1,dimn
               call in_bin_512(lres_unit,dpar_pos(:,li),lglocnt,lnext_rec)
            end do
         end if 
! The call to des_restart_map allocates and sets the restartmap. It
! sets pip, and pea(pip,1) to true
         call des_restart_map(lglocnt) 
         deallocate(dpar_pos)

! At this point, each process knows how many and which particles 
! it contains. Since only the root process reads the restart file, 
! parameters for scattering the data are set.

! Allocate variables for scatter.
! double precision arryas:
         allocate(dprocbuf(pip))
         allocate(drootbuf(lglocnt))
! integer arrays:
         allocate(iprocbuf(pip))
         allocate(irootbuf(lglocnt))

! Construct an array for the Root process that states the number of
! (real) particles on each process.
         lscattercnts(:) = 0 
         lscattercnts(mype)=pip 
         iscr_recvcnt = pip
         call global_sum(lscattercnts,iscattercnts)

! Calculate the displacements for each process in the global array.
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)  
         end do 

! Read in particle properties:
!--------------------------------------------------------------------->>
! Global particle index numbers.
         call des_readscatter(lres_unit,iglobal_id,lglocnt,lnext_rec)
! Dynamic particle count elements:
! 2 - Classifies particles as 'new'
! 4 - Classifies particles as ghost particles
         do li = 2,4
            call des_readscatter(lres_unit,pea(:,li),lglocnt,lnext_rec)
         end do 
! Translational velocities (one coordinate at a time).
         do li = 1,dimn
            call des_readscatter(lres_unit,des_vel_new(:,li),lglocnt,lnext_rec)
         end do 
! Rotational velocity.
         do li = 1,ltor_dimn
            call des_readscatter(lres_unit,omega_new(:,li),lglocnt,lnext_rec)
         end do 
! Particle radi.
         call des_readscatter(lres_unit,des_radius,lglocnt,lnext_rec)
! Particle densities.
         call des_readscatter(lres_unit,ro_sol,lglocnt,lnext_rec)
! MPPIC Only: Statistical weight of each particle.
         IF(MPPIC) call des_readscatter(lres_unit,des_stat_wt,lglocnt,lnext_rec)

! Particle temperatures.
         IF(DES_ENERGY_EQ)  call des_readscatter(lres_unit, &
            DES_T_s_NEW, lglocnt, lnext_rec)
         IF(ANY_DES_SPECIES_EQ) THEN
! Solids phase index: This is needed for reacting flows as the diameter
! (shrinking_particle) or density (variable_density) reaction models
! prevent identifying particles from initial phase specifications.
            call des_readscatter(lres_unit,PIJK(:,5),lglocnt,lnext_rec)
! Species mass fractions.
            DO li = 1, MAX_DES_NMAX
               call des_readscatter(lres_unit, &
                  DES_X_s(:,li),lglocnt,lnext_rec)
            ENDDO
         ENDIF

! Neighbor/collision data:
!--------------------------------------------------------------------->>
! Number of neighbors
         call des_readscatter(lres_unit,neighbours(:,1),lglocnt,lnext_rec)
! Number of particles in contact with each particle.
         call des_readscatter(lres_unit,pn(:,1),lglocnt,lnext_rec)

! *** Use global particle IDs for neighbors were stored. After the 
! *** restart is completely read, the IDs read from the restart are
! *** converted from global to local.
         do li = 2,maxneighbors
! Neighbor particles.
            call des_readscatter(lres_unit,neighbours(:,li),lglocnt,lnext_rec)
! Neighbor contact data.
            call des_readscatter(lres_unit,pn(:,li),lglocnt,lnext_rec)
            call des_readscatter(lres_unit,pv(:,li),lglocnt,lnext_rec)
            do lj = 1,dimn  
! Accumulated tangential displacement that occured during collision
               call des_readscatter(lres_unit,pft(:,li,lj),lglocnt,lnext_rec)
            end do 
         end do 

! Mass inlet/outlet parameters:
!--------------------------------------------------------------------->>
         if(des_mi) then 
!  *******************************************************************
!  * Warning: DES_MI_TIME is an array allocated to des_bcmi.         *
!  *          Could this be a problem?                               *
!  *******************************************************************
! Time to seed next particle.
            if(mype .eq. pe_io)then
               read (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
            end if 
            call bcast(des_mi_time)

! For each mass inlet, store:
!   mi_factor - number of inlet grid cells containing particles
!   mi_window - grid cell length (> diameter)
!   particle_plcmnt - indicates seeding approach (ordered vs. random)
            do lboundnum =1, des_bcmi
               if(mype.eq.pe_io) then 
                  read (lres_unit,rec=lnext_rec) mi_factor(lboundnum), &
                     mi_window(lboundnum),particle_plcmnt(lboundnum)
                  lnext_rec = lnext_rec + 1 
               end if  
               call bcast (mi_factor)
               call bcast (mi_window)
               call bcast (particle_plcmnt)
               if(particle_plcmnt(lboundnum) == 'ORDR')then
! For ordered inlets, store:
!   mi_order pointer is 'associated', and if so, store:
!     number of grid cells
!     mi_order - order to populate grid cells with incoming particles
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc
                     lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize
                        lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(mi_order(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit, &
                        mi_order(lboundnum)%value,lsize,lnext_rec)
                     call bcast(mi_order(lboundnum)%value)
                  endif
!   i_of_mi is 'associated', and if so, store:
!     number of cells in i-direction
!     i indices of mass inlet grid
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc
                     lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize
                        lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(i_of_mi(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit, &
                        i_of_mi(lboundnum)%value,lsize,lnext_rec)
                     call bcast(i_of_mi(lboundnum)%value)
                  endif

!   j_of_mi is 'associated', and if so, store:
!     number of cells in j-direction
!     j indices of mass inlet grid
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc
                     lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize
                        lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(j_of_mi(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit, &
                        j_of_mi(lboundnum)%value,lsize,lnext_rec)
                     call bcast(j_of_mi(lboundnum)%value)
                  endif
               endif  
            enddo
         endif

! Clean up arrays:
         deallocate (dprocbuf)
         deallocate (drootbuf)
         deallocate (iprocbuf)
         deallocate (irootbuf)
         deallocate (irestartmap)

! Convert the  neighbors and contact data from global IDs to local IDs.
         call des_restart_neigh

      end if   ! end if/else bdist_io
!-----------------------------------------------------------------<<<

! Copy values from "NEW" arrays to "OLD" arrays.
! Higher order integration (Adams-Bashforth) will default to Euler for
! the first time step since 'OLD' arrays values were not stored.
      omega_old(:,:) = omega_new(:,:) 
      des_pos_old(:,:) = des_pos_new(:,:)
      des_vel_old(:,:) = des_vel_new(:,:)
      DES_T_s_OLD(:) = DES_T_s_NEW(:)

! Close the restart file.
      if(bdist_io .or.mype .eq. pe_io) close(lres_unit)

      end subroutine des_read_restart 
