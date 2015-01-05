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
      module mpi_funs_des

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

      contains

!------------------------------------------------------------------------
! Subroutine       : des_par_exchange
! Comments         : main subroutine controls entire exchange of particles
!                    between the processors
! Steps
!  --> bin the particles
!  --> check if the send and recv buffer size is enough
!  --> particle crossing boundary - will take place in following order
!      a. top-bottom interface
!      b. north-south interface
!      c. east-west interface
!      active particles(pea(par,1))in the ghost cell will be packed
!      and exchange to respective processor.
!      This order will also take care of particles crossing corners.
!      (for example particle crossing directly into northwest block)
!  --> bin the particles (if required)
!  --> ghost cell particles - exchange will take place in following order
!      a. east-west interface
!      b. north-south interface
!      c. top-bottom interface
!------------------------------------------------------------------------
      subroutine des_par_exchange()

      use mpi_pack_des, only: desmpi_pack_parcross
      use mpi_unpack_des, only: desmpi_unpack_parcross

      use mpi_pack_des, only: desmpi_pack_ghostpar
      use mpi_unpack_des, only: desmpi_unpack_ghostpar
!-----------------------------------------------
      implicit none

!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer linter,lface
      integer, save :: lcheckbuf = 0
!-----------------------------------------------

! bin the particles and check the bufsize
      call desgrid_pic(plocate=.true.)

! mpi_all_reduce is expensive so avoiding the check for buffer size
! further high factors are used for buffer and hence check at low
! frequency is enough
      if (mod(lcheckbuf,100) .eq. 0) then
         call desmpi_check_sendrecvbuf
         lcheckbuf = 0
      end if
      lcheckbuf = lcheckbuf + 1

! call particle crossing the boundary exchange in T-B,N-S,E-W order
      dsendbuf(1,:) = 0; drecvbuf(1,:) =0
      ispot = 1
      do linter = merge(2,3,NO_K),1,-1
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface))cycle
            call desmpi_pack_parcross(lface)
            call desmpi_sendrecv_init(lface)
         end do
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface)) cycle
            call desmpi_sendrecv_wait(lface)
            call desmpi_unpack_parcross(lface)
         end do
! update pic this is required for particles crossing corner cells
         do lface = linter*2-1,linter*2
            if(dsendbuf(1,lface).gt.0.or.drecvbuf(1,lface).gt.0) then
               call desgrid_pic(plocate=.false.)
               exit
            end if
         end do
      end do
      call des_mpi_barrier

!      call des_dbgmpi(5)


      IF(.NOT.MPPIC) THEN
! call ghost particle exchange in E-W, N-S, T-B order
         dsendbuf(1,:) = 0; drecvbuf(1,:) =0
         ighost_updated(:) = .false.
         ispot = 1
         do linter = 1,merge(2,3,NO_K)
            do lface = linter*2-1,linter*2
               if(.not.iexchflag(lface))cycle
               call desmpi_pack_ghostpar(lface)
               call desmpi_sendrecv_init(lface)
            end do
            do lface = linter*2-1,linter*2
               if(.not.iexchflag(lface)) cycle
               call desmpi_sendrecv_wait(lface)
               call desmpi_unpack_ghostpar(lface)
            end do
! update pic required as particles in ghost cell can move between ghost cells
            do lface = linter*2-1,linter*2
               if(dsendbuf(1,lface).gt.0.or.drecvbuf(1,lface).gt.0) then
                  call desgrid_pic(plocate=.false.)
                  exit
               end if
            end do
         end do
         if(do_nsearch) call desmpi_cleanup
         call des_mpi_barrier
      ENDIF   ! end if(.not.mppic)

!      call des_dbgmpi(2)
!      call des_dbgmpi(3)
!      call des_dbgmpi(4)
!      call des_dbgmpi(6)
!      call des_dbgmpi(7)

      end subroutine des_par_exchange


!------------------------------------------------------------------------
! Subroutine       : desmpi_check_sendrecvbuf
! Purpose          : checks if the sendrecvbuf size is enough if not redefine
!                    the arrays and sets recvcnts
!
!------------------------------------------------------------------------
      subroutine desmpi_check_sendrecvbuf
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!local variables
!-----------------------------------------------
      integer:: lijk,ltot_ind,lparcnt,lmaxcnt,lface,ltordimn,lpacketsize,lindx
      real :: lfactor = 1.5
!-----------------------------------------------

      lmaxcnt = 0
      ltordimn = merge(1,3,NO_K)
      lpacketsize = 2*dimn + ltordimn+ 5
      do lface = 1,2*dimn
         ltot_ind = isendindices(1,lface)
         lparcnt = 0
         do lindx = 2,ltot_ind+1
            lijk = isendindices(lindx,lface)
            lparcnt = lparcnt + dg_pic(lijk)%isize
         enddo
         if(lparcnt.gt.lmaxcnt) lmaxcnt = lparcnt
      enddo

      call global_all_max(lmaxcnt)
      if (imaxbuf .lt. lmaxcnt*lpacketsize+ibufoffset) then
         imaxbuf = lmaxcnt*lpacketsize*lfactor
         if(allocated(dsendbuf)) deallocate(dsendbuf,drecvbuf)
         allocate(dsendbuf(imaxbuf,2*dimn),drecvbuf(imaxbuf,2*dimn))
      endif

      end subroutine desmpi_check_sendrecvbuf

!------------------------------------------------------------------------
! Subroutine       : desmpi_cleanup
! Purpose          : Cleans the ghost particles
! Parameter        :
!------------------------------------------------------------------------
      subroutine desmpi_cleanup
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer ltot_ind,lface,lindx,lijk,lcurpar,lpicloc
!-----------------------------------------------

      do lface = 1,dimn*2
         if(.not.iexchflag(lface))cycle
         ltot_ind = irecvindices(1,lface)
         do lindx = 2,ltot_ind+1
            lijk = irecvindices(lindx,lface)
            do lpicloc =1,dg_pic(lijk)%isize
               lcurpar = dg_pic(lijk)%p(lpicloc)
               if(ighost_updated(lcurpar)) cycle
               pip = pip - 1
               ighost_cnt = ighost_cnt-1
               pea(lcurpar,1:4) = .false.
               fc(:,lcurpar) = 0.0
               des_pos_new(:,lcurpar)=0
               IF (DO_OLD) THEN
                  des_pos_old(:,lcurpar)=0
                  des_vel_old(:,lcurpar)=0
               ENDIF
               des_vel_new(:,lcurpar)=0
               omega_new(:,lcurpar)=0

               if(ENERGY_EQ) then
                  des_t_s_new(lcurpar)=0
                  des_t_s_old(lcurpar)=0
               endif

               if(ANY_SPECIES_EQ)then
                  des_x_s(lcurpar,1:dimension_n_s)= 0
               endif

               des_usr_var(1:3,lcurpar)= 0

            end do
         end do
      end do
      end subroutine desmpi_cleanup



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



      end module mpi_funs_des
