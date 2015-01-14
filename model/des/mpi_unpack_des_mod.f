!------------------------------------------------------------------------
! Module           : desmpi
! Purpose          : Contains wrapper class for mpi communications- send,recv
!
! Author           : Pradeep.G
!
! Purpose          : Module contains subroutines and variables related to
!                    des mpi communication.
!
! Comments         : do_nsearch flag should be set to true before calling
!                    des_par_exchange; when do_nsearch is true ghost particles of the
!                    system will be updated, which will be later used to generate
!                    neighbour list.
!------------------------------------------------------------------------
      module mpi_unpack_des

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use parallel_mpi
      use mpi_utility
      use discretelement
      use desgrid
      use compar
      use physprop
      use sendrecv
      use des_bc
      use desmpi_wrapper
      use sendrecvnode
      use mfix_pic
      use des_thermo
      use run, only: ENERGY_EQ,ANY_SPECIES_EQ
      use param, only: DIMENSION_N_s
      use des_rxns
      use desmpi

      use mpi_comm_des, only: desmpi_sendrecv_init
      use mpi_comm_des, only: desmpi_sendrecv_wait

      
      interface unpack_dbuf
         module procedure unpack_db0, unpack_db1,unpack_i0,unpack_i1
      end interface unpack_dbuf


      contains


!----------------------------------------------------------------------!
!Unpack subroutine for single real variables                           !
!----------------------------------------------------------------------!
      subroutine unpack_db0(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata

      idata = drecvbuf(lbuf,pface)
      lbuf = lbuf + 1

      return
      end subroutine unpack_db0

!----------------------------------------------------------------------!
!Unpack subroutine for real arrays                                     !
!----------------------------------------------------------------------!
      subroutine unpack_db1(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(lbuf:lbuf+lsize-1,pface)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_db1


!----------------------------------------------------------------------!
!Unpack subroutine for single integer variables                        !
!----------------------------------------------------------------------!
      subroutine unpack_i0(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata

      idata = drecvbuf(lbuf,pface)
      lbuf = lbuf + 1

      return
      end subroutine unpack_i0

!----------------------------------------------------------------------!
!Unpack subroutine for integer arrays                                  !
!----------------------------------------------------------------------!
      subroutine unpack_i1(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(lbuf:lbuf+lsize-1,pface)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_i1


!------------------------------------------------------------------------
! Subroutine       : desmpi_unpack_ghostpar
! Purpose          : unpacks the ghost particle from the recv buffer
! Parameter        : pface - value from 1 to 6 represents faces
!
!------------------------------------------------------------------------
      subroutine desmpi_unpack_ghostpar(pface)
      use constant, only: PI
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: pface
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparid,lprvijk,lijk,lparijk,lparcnt,ltot_ind
      integer :: lpacketsize,lbuf,lindx,llocpar,lnewcnt,lpicloc
      logical,dimension(:),allocatable :: lfound
      integer,dimension(:),allocatable :: lnewspot,lnewpic
!-----------------------------------------------

! unpack the particles:
! if it already exists update the position
! if not and do_nsearch is true then add to the particle array


      lpacketsize = 2*dimn + 3+ 5
      lparcnt = drecvbuf(1,pface)
      lnewcnt = lparcnt
      allocate (lfound(lparcnt),lnewspot(lparcnt),lnewpic(dg_ijksize2))
      lfound(:) = .false.
      lnewspot(:) =0
      lnewpic = 0

      do lcurpar = 1,lparcnt
         lbuf = (lcurpar-1)*lpacketsize+ibufoffset

! 1) Global ID
         call unpack_dbuf(lbuf,lparid,pface)
! 2) DES Grid IJK
         call unpack_dbuf(lbuf,lparijk,pface)
! 3) DES Grid IJK - Previous
         call unpack_dbuf(lbuf,lprvijk,pface)

! Determine if this particle already exists on this process as a
! ghost particle. If so, (lfound), then the current infomration is
! updated on the current process. Otherwise (.NOT.lfound) a new
! ghost particle is created on this process.
         lfound(lcurpar) = locate_par(lparid,lprvijk,llocpar)
         if (lparijk .ne. lprvijk .and. .not.lfound(lcurpar)) then
            lfound(lcurpar) = locate_par(lparid,lparijk,llocpar)
         endif

         if(lfound(lcurpar)) then
! Store the local variables
            dg_pijk(llocpar) = lparijk
            dg_pijkprv(llocpar) = lprvijk

! 4) Radious
            call unpack_dbuf(lbuf,des_radius(llocpar),pface)
! 5) Phase index
            call unpack_dbuf(lbuf,pijk(llocpar,5),pface)
! 6) Position
            call unpack_dbuf(lbuf,des_pos_new(1:dimn,llocpar),pface)
! 7) Translational Velocity
            call unpack_dbuf(lbuf,des_vel_new(1:dimn,llocpar),pface)
! 8) Rotational Velocity
            call unpack_dbuf(lbuf,omega_new(1:3,llocpar),pface)
! 9) Temperature
            if(ENERGY_EQ)then
               call unpack_dbuf(lbuf,des_t_s_new(llocpar),pface)
            endif
! 10) Species Composition
            if(ANY_SPECIES_EQ)then
               call unpack_dbuf(lbuf,des_x_s(llocpar,1:dimension_n_s),pface)
            endif

! 11) User Variables
            call unpack_dbuf(lbuf,des_usr_var(1:3,llocpar),pface)

! Calculate the volume of the ghost particle.
            PVOL(llocpar) = (4.0D0/3.0D0)*PI*DES_RADIUS(llocpar)**3
! Flag that the ghost particle was updated.
            ighost_updated(llocpar) = .true.
            lnewcnt = lnewcnt-1

! Copy the current value to the previous value if needed.
            IF (DO_OLD) THEN
               des_pos_old(:,llocpar)= des_pos_new(:,llocpar)
               des_vel_old(:,llocpar)= des_vel_new(:,llocpar)
               if(ENERGY_EQ)des_t_s_old(llocpar)= des_t_s_new(llocpar)
               omega_old(:,llocpar)= omega_new(:,llocpar)
            ENDIF

         else
            lnewpic(lparijk) = lnewpic(lparijk) + 1
         endif
      enddo

! if do_nsearch is on then add new particles and clean up ghost particles
      if (do_nsearch) then
         if((max_pip-pip).lt.lnewcnt) call redim_par(pip+lnewcnt)
         ighost_cnt = ighost_cnt + lnewcnt
         pip = pip + lnewcnt
         do lcurpar = 1,lparcnt
            if(lfound(lcurpar)) cycle
            lbuf = (lcurpar-1)*lpacketsize+ibufoffset

!  1) Global particle ID
            call unpack_dbuf(lbuf,lparid,pface)
!  2) DES grid IJK
            call unpack_dbuf(lbuf,lparijk,pface)
!  3) DES grid IJK - Previous
            call unpack_dbuf(lbuf,lprvijk,pface)
! Locate the first open space in the particle array.
            do while(pea(ispot,1))
               ispot = ispot + 1
            enddo
! Set the flags for the ghost particle and store the local variables.
            pea(ispot,1) = .true.
            pea(ispot,2) = .false.
            pea(ispot,3) = .false.
            pea(ispot,4) = .true.
            iglobal_id(ispot)  = lparid
            dg_pijk(ispot) = lparijk
            dg_pijkprv(ispot) = lprvijk
!  4) Particle radius
            call unpack_dbuf(lbuf,des_radius(ispot),pface)
!  5) Particle phase index
            call unpack_dbuf(lbuf,pijk(ispot,5),pface) 
!  6) Particle position
            call unpack_dbuf(lbuf,des_pos_new(1:dimn,ispot),pface)
!  7) Particle velocity
            call unpack_dbuf(lbuf,des_vel_new(1:dimn,ispot),pface)
!  8) Particle rotational velocity
            call unpack_dbuf(lbuf,omega_new(1:dimn,ispot),pface)
!  9) Particle temperature.
            if(ENERGY_EQ)then
               call unpack_dbuf(lbuf,des_t_s_new(ispot),pface)
            endif
! 10) Particle species composition
            if(ANY_SPECIES_EQ)then
               call unpack_dbuf(lbuf,des_x_s(ispot,1:dimension_n_s),pface)
            endif
! 11) User varaible
            call unpack_dbuf(lbuf,des_usr_var(1:3,ispot),pface)

            ighost_updated(ispot) = .true.
            lnewspot(lcurpar) = ispot

            PVOL(ispot) = (4.0D0/3.0D0)*PI*DES_RADIUS(ispot)**3

            IF (DO_OLD) THEN
               des_pos_old(1:dimn,ispot) = des_pos_new(1:dimn,ispot)
               des_vel_old(1:dimn,ispot) = des_vel_new(1:dimn,ispot)
               omega_old(1:3,ispot) = omega_new(1:3,ispot)
               if(ENERGY_EQ) des_t_s_old(ispot) = des_t_s_new(ispot)
            ENDIF
         enddo
      endif

!deallocate temporary variablies
      deallocate (lfound,lnewspot,lnewpic)

      end subroutine desmpi_unpack_ghostpar


!------------------------------------------------------------------------
! Subroutine       : desmpi_unpack_parcross
! Purpose          : pack the particle crossing the boundary
! Parameter        : pface - value from 1 to 6 represents faces
!
!------------------------------------------------------------------------
      subroutine desmpi_unpack_parcross(pface)

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: pface
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lijk,lcurpar,lparcnt,llocpar,lparid,lparijk,lprvijk
      integer :: lneighindx,lneigh,lcontactindx,lcontactid,lcontact,&
                 lneighid,lneighijk,lneighprvijk
      logical :: lfound
      integer :: lpacketsize,lbuf,ltmpbuf,lcount
      logical :: lcontactfound,lneighfound
      integer :: cc,ii,kk,num_pairs_sent
!-----------------------------------------------

! loop through particles and locate them and make changes
      lpacketsize = 9*dimn + 3*4 + 15
      lparcnt = drecvbuf(1,pface)

! if mppic make sure enough space available
      if(mppic .and. (max_pip-pip).lt.lparcnt) call redim_par(pip+lparcnt)

      do lcurpar =1,lparcnt
         lfound = .false.
         lbuf = (lcurpar-1)*lpacketsize + ibufoffset
         call unpack_dbuf(lbuf,lparid,pface)
         call unpack_dbuf(lbuf,lparijk,pface)
         call unpack_dbuf(lbuf,lprvijk,pface)

! if mppic add the particles to free spots else locate the particles
         if (mppic) then
            do while(pea(ispot,1))
               ispot = ispot + 1
            enddo
            llocpar = ispot
         else
            lfound  = locate_par(lparid,lprvijk,llocpar)
            if (.not. lfound) then
               WRITE(*,700) ineighproc(pface), mype
               call des_mpi_stop
            endif
            ighost_cnt = ighost_cnt - 1
         endif

! convert the local particle from ghost to existing and update its position
         pea(llocpar,1) = .true.
         pea(llocpar,4) = .false.
         dg_pijk(llocpar) = lparijk
         dg_pijkprv(llocpar) = lprvijk
         call unpack_dbuf(lbuf,des_radius(llocpar),pface)
         call unpack_dbuf(lbuf,pijk(llocpar,1:5),pface)
!         pea(llocpar,2:3)     = drecvbuf(lbuf:lbuf+1,pface) ; lbuf=lbuf+2
         pea(llocpar,2:3) = .false.
         if (drecvbuf(lbuf,pface).eq.1) pea(llocpar,2) = .true. ; lbuf = lbuf + 1
         if (drecvbuf(lbuf,pface).eq.1) pea(llocpar,3) = .true. ; lbuf = lbuf + 1
         call unpack_dbuf(lbuf,ro_sol(llocpar),pface)
         call unpack_dbuf(lbuf,pvol(llocpar),pface)
         call unpack_dbuf(lbuf,pmass(llocpar),pface)
         call unpack_dbuf(lbuf,omoi(llocpar),pface)
         call unpack_dbuf(lbuf,des_pos_new(1:dimn,llocpar),pface)
         call unpack_dbuf(lbuf,des_vel_new(1:dimn,llocpar),pface)

         if(ENERGY_EQ) then
            call unpack_dbuf(lbuf,des_t_s_old(llocpar),pface)
            call unpack_dbuf(lbuf,des_t_s_new(llocpar),pface)
         endif

         if(ANY_SPECIES_EQ)then
            call unpack_dbuf(lbuf,des_x_s(llocpar,1:dimension_n_s),pface)
         endif

         call unpack_dbuf(lbuf,des_usr_var(1:3,llocpar),pface)

         call unpack_dbuf(lbuf,omega_new(1:3,llocpar),pface)
         IF (DO_OLD) THEN
            call unpack_dbuf(lbuf,des_pos_old(1:dimn,llocpar),pface)
            call unpack_dbuf(lbuf,des_vel_old(1:dimn,llocpar),pface)
            call unpack_dbuf(lbuf,omega_old(1:3,llocpar),pface)
            call unpack_dbuf(lbuf,des_acc_old(1:dimn,llocpar),pface)
            call unpack_dbuf(lbuf,rot_acc_old(1:3,llocpar),pface)
         ENDIF
         call unpack_dbuf(lbuf,fc(:,llocpar),pface)
         call unpack_dbuf(lbuf,tow(1:3,llocpar),pface)

      end do

      lbuf = lparcnt*lpacketsize + ibufoffset

      call unpack_dbuf(lbuf,num_pairs_sent,pface)

      do cc = 1, num_pairs_sent

         call unpack_dbuf(lbuf,lparid,pface)

         call unpack_dbuf(lbuf,lparijk,pface)

         if (.not. locate_par(lparid,lparijk,llocpar)) then
            print *,"at buffer location",lbuf," pface = ",pface
            print *,"COULD NOT FIND PARTICLE ",lparid," IN IJK ",lparijk
            call des_mpi_stop
         endif

         call unpack_dbuf(lbuf,lneighid,pface)

         call unpack_dbuf(lbuf,lneighijk,pface)

         if (.not. locate_par(lneighid,lneighijk,lneigh)) then
            if (.not. exten_locate_par(lneighid,lparijk,lneigh)) then
               print *,"at buffer location",lbuf," pface = ",pface
               print *,"COULD NOT FIND NEIGHBOR ",lneighid," IN IJK ",lneighijk
               call des_mpi_stop
            endif
         endif

         call add_pair(llocpar,lneigh)

         pv_pair(pair_num)=merge(.true.,.false.,0.5 < drecvbuf(lbuf,pface))
         lbuf = lbuf+1 

         do ii=1,DIMN
            call unpack_dbuf(lbuf,pfn_pair(ii,pair_num),pface)
            call unpack_dbuf(lbuf,pft_pair(ii,pair_num),pface)
         enddo
      enddo

 700 FORMAT(/2X,'From: DESMPI_UNPACK_PARCROSS: ',/2X,&
         'ERROR: Unable to locate particles moving from ',I4.4,&
         ' to ', I4.4)
 701 FORMAT(/2X,'From: DESMPI_UNPACK_PARCROSS: ',/2X,&
         'WARNING: Unable to locate neighbor for particles ',&
         'crossing boundary')
 702 FORMAT(/2X,'From: DESMPI_UNPACK_PARCROSS: ',/2X,&
         'WARNING: Unable to locate neighbor for particles ',&
         'crossing boundary.'/2X,'Contact particle ID =',I10)

      END SUBROUTINE desmpi_unpack_parcross


!------------------------------------------------------------------------
! Function         : locate_par
! Purpose          : locates particle in ijk and returns true if found
! Parameter        : pglobalid - global id of the particle (input)
!                    pijk - ijk of the cell (input)
!                    plocalno - local particle number (output)
!------------------------------------------------------------------------
      function locate_par(pglobalid,pijk,plocalno)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      logical :: locate_par
      integer :: pglobalid,pijk,plocalno
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lpicloc,lcurpar
!-----------------------------------------------

      locate_par = .false.
      if (pijk .lt. dg_ijkstart2 .or. pijk .gt. dg_ijkend2) then
         return
      endif

      do lpicloc = 1,dg_pic(pijk)%isize
         lcurpar = dg_pic(pijk)%p(lpicloc)
         if (iglobal_id(lcurpar) .eq. pglobalid) then
            plocalno = lcurpar
            locate_par = .true.
            return
         endif
      enddo

      return
      end function locate_par

!------------------------------------------------------------------------
! Function         : exten_locate_par
! Purpose          : locates particles extensively by searching all neighbouring
!                    cells similar to grid based search
! Parameter        : pglobalid - global id of the particle (input)
!                    pijk - ijk of the center cell (input)
!                    plocalno - localparticle number (output)
!------------------------------------------------------------------------
      function exten_locate_par(pglobalid,pijk,plocalno)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      logical :: exten_locate_par
      integer :: pglobalid,pijk,plocalno
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lpicloc,lcurpar
      integer :: lijk,li,lj,lk,lic,ljc,lkc,lkoffset
!-----------------------------------------------
      exten_locate_par = .false.
      lic = dg_iof_lo(pijk)
      ljc = dg_jof_lo(pijk)
      lkc = dg_kof_lo(pijk)
      lkoffset = merge(0, 1, NO_K)
      do  lk = lkc-lkoffset,lkc+lkoffset
      do  lj = ljc-1,ljc+1
      do  li = lic-1,lic+1
         lijk = dg_funijk(li,lj,lk)
         if (lijk .lt. dg_ijkstart2 .or. lijk .gt. dg_ijkend2) cycle
         do lpicloc = 1, dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)
            if (iglobal_id(lcurpar) .eq. pglobalid) then
               plocalno = lcurpar
               exten_locate_par = .true.
               return
            end if
         end do
      end do
      end do
      end do
      return
      end function exten_locate_par

      end module mpi_unpack_des
