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
 
! open the restart file  
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


!----------------------------------------------------------------->>>
      if (bdist_IO) then 
         read(lres_unit,rec=1) pip,ighost_cnt,vtp_findex,tecplot_findex,dtsolid 
         pea(1:pip,1) = .true.
         if (pip .gt. max_pip) then 
            write(*,*) "From des_read_restart:"
            write(*,*) "Error: The pip is grater than current max_pip"
            write(*,*) "pip=" ,pip,"; max_pip =", max_pip 
            call mfix_exit(mype)
         end if 

         lnext_rec = 6
         do li = 1,dimn
            call in_bin_512(lres_unit,des_pos_new(:,li),pip,lnext_rec)
         end do 
         call in_bin_512i(lres_unit,iglobal_id,pip,lnext_rec)
         allocate(ltemparr(pip))
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
         do li = 1,dimn
            call in_bin_512(lres_unit,des_vel_new(:,li),pip,lnext_rec)
         end do 
         do li = 1,ltor_dimn
            call in_bin_512(lres_unit,omega_new(:,li),pip,lnext_rec)
         end do 
         call in_bin_512(lres_unit,des_radius,pip,lnext_rec)
         call in_bin_512(lres_unit,ro_sol,pip,lnext_rec)
         IF(MPPIC) call in_bin_512(lres_unit,des_stat_wt,pip,lnext_rec)
         
         call in_bin_512i(lres_unit,neighbours(:,1),pip,lnext_rec)
         call in_bin_512i(lres_unit,pn(:,1),pip,lnext_rec)
         do li = 2,maxneighbors
            call in_bin_512i(lres_unit,neighbours(:,li),pip,lnext_rec)
            call in_bin_512i(lres_unit,pn(:,li),pip,lnext_rec)
            call in_bin_512i(lres_unit,pv(:,li),pip,lnext_rec)
            do lj = 1,dimn  
               call in_bin_512(lres_unit,pft(:,li,lj),pip,lnext_rec)
            end do 
         end do 

! read details for mass inlet and outlet parameters
         if(des_mi)then
            read (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
            do lboundnum =1, des_bcmi
               read (lres_unit,rec=lnext_rec) mi_factor(lboundnum),mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 
               if(particle_plcmnt(lboundnum) == 'ORDR')then
                  read (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(mi_order(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
                  read (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(i_of_mi(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
                  read (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     allocate(j_of_mi(lboundnum)%value(lsize))
                     call in_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif
!-----------------------------------------------------------------<<<
      else   ! else branch if not bdist_IO
!----------------------------------------------------------------->>>

! for single file read the particle info first and then read the global
! id and position and devlop the mapping between the reading order 
! and local id and proc number 
         lglocnt = 10
         if(mype.eq.pe_io) read(lres_unit,rec=1) lglocnt,vtp_findex,tecplot_findex,dtsolid 
         call bcast(vtp_findex)
         call bcast(tecplot_findex)
         call bcast(dtsolid)
         allocate(dpar_pos(lglocnt,dimn),irestartmap(lglocnt))
         if(mype.eq.pe_io) then 
            lnext_rec = 6
            do li = 1,dimn
               call in_bin_512(lres_unit,dpar_pos(:,li),lglocnt,lnext_rec)
            end do
         end if 
! this call to des_restart_map allocates and sets the restartmap,
! sets pip, and sets pea(pip,1) to true          
         call des_restart_map(lglocnt) 
         deallocate(dpar_pos)
 
! set the parameters for scatter
         allocate(dprocbuf(pip),drootbuf(lglocnt),iprocbuf(pip),irootbuf(lglocnt))
         lscattercnts(:) = 0 
         lscattercnts(mype)=pip 
         iscr_recvcnt = pip
         call global_sum(lscattercnts,iscattercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)  
         end do 

! read the variables in the order it is written 
         call des_readscatter(lres_unit,iglobal_id,lglocnt,lnext_rec)
         do li = 2,4
            call des_readscatter(lres_unit,pea(:,li),lglocnt,lnext_rec)
         end do 
         do li = 1,dimn
            call des_readscatter(lres_unit,des_vel_new(:,li),lglocnt,lnext_rec)
         end do 
         do li = 1,ltor_dimn
            call des_readscatter(lres_unit,omega_new(:,li),lglocnt,lnext_rec)
         end do 
         call des_readscatter(lres_unit,des_radius,lglocnt,lnext_rec)
         call des_readscatter(lres_unit,ro_sol,lglocnt,lnext_rec)
         IF(MPPIC) call des_readscatter(lres_unit,des_stat_wt,lglocnt,lnext_rec)

! read neighbour and contact details 
         call des_readscatter(lres_unit,neighbours(:,1),lglocnt,lnext_rec)
         call des_readscatter(lres_unit,pn(:,1),lglocnt,lnext_rec)
         do li = 2,maxneighbors
            call des_readscatter(lres_unit,neighbours(:,li),lglocnt,lnext_rec)
            call des_readscatter(lres_unit,pn(:,li),lglocnt,lnext_rec)
            call des_readscatter(lres_unit,pv(:,li),lglocnt,lnext_rec)
            do lj = 1,dimn  
               call des_readscatter(lres_unit,pft(:,li,lj),lglocnt,lnext_rec)
            end do 
         end do 

! read details for mass inlet and outlet parameters
         if(des_mi) then 
            if(mype .eq. pe_io)then
               read (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
            end if 
            call bcast(des_mi_time)
            do lboundnum =1, des_bcmi
               if(mype.eq.pe_io) then 
                  read (lres_unit,rec=lnext_rec) mi_factor(lboundnum),mi_window(lboundnum),particle_plcmnt(lboundnum)
                  lnext_rec = lnext_rec + 1 
               end if  
               call bcast (mi_factor)
               call bcast (mi_window)
               call bcast (particle_plcmnt)
               if(particle_plcmnt(lboundnum) == 'ORDR')then
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(mi_order(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                     call bcast(mi_order(lboundnum)%value)
                  endif
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(i_of_mi(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                     call bcast(i_of_mi(lboundnum)%value)
                  endif
                  if (mype.eq.pe_io) then 
                     read (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                     if(lassoc)then
                        read (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     end if
                  end if 
                  call bcast(lassoc)
                  call bcast (lsize)
                  if (lassoc) then 
                     allocate(j_of_mi(lboundnum)%value(lsize))
                     if(mype.eq.pe_io) call in_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                     call bcast(j_of_mi(lboundnum)%value)
                  endif
               endif  
            enddo
         endif

         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf,irestartmap)
! Change neighbours and contact details from global id to local 
! particle number 
         call des_restart_neigh

      end if   ! end if/else bdist_io
!-----------------------------------------------------------------<<<

! setting old positions to new positions 
! will higher order integration be incorrect since des_vel_oold will not
! be the old value?
      omega_old(:,:) = omega_new(:,:) 
      des_pos_old(:,:) = des_pos_new(:,:)
      des_vel_old(:,:) = des_vel_new(:,:)
      if(bdist_io .or.mype .eq. pe_io) close(lres_unit)

      end subroutine des_read_restart 
