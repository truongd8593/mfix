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

! Contains following subroutines:
!    des_par_exchange
!    desmpi_init, desmpi_setcomm, desmpi_sendrecv_init,
!    desmpi_sendrecv_wait, desmpi_gatherv, desmpi_scatterv
!    des_scatter_particle, des_restart_map, desmpi_check_sendrecvbuf,
!    desmpi_pack_ghostpar, desmpi_unpack_ghostpar, desmpi_cleanup,
!    desmpi_pack_parcross, desmpi_unpack_parcross,
!    des_addnodevalues, des_addnodevalues2,
!    des_gather_d,l,i, des_gatherwrite_d,l,i, des_writepar_d,l,i
!    des_scatter_d,l,i, des_readscatter_d,l,i,
!    des_restart_neigh, redim_par, des_dbgmpi
! Contains following functions:
!    locate_par, exten_locate_par
!------------------------------------------------------------------------
      module desmpi

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
!-----------------------------------------------

! flags and constants for interfaces
      integer,dimension(:),allocatable   :: ineighproc
      logical,dimension(:),allocatable   :: iexchflag

! offset for periodic boundaries
      double precision,dimension(:,:),allocatable   :: dcycl_offset

! following variables used for sendrecv ghost particles and particle exchange
      double precision, dimension(:,:), allocatable :: dsendbuf,drecvbuf
      integer,dimension(:),allocatable:: isendcnt
      integer,dimension(:),allocatable:: isendreq,irecvreq
      integer,parameter :: ibufoffset = 2
      integer :: imaxbuf, ispot

! following variables are used for gather and scatter
      double precision, dimension(:), allocatable :: drootbuf,dprocbuf
      integer, dimension(:), allocatable :: irootbuf,iprocbuf
      integer,dimension(:), allocatable:: iscattercnts,idispls,igathercnts
      integer :: iscr_recvcnt,igath_sendcnt

! following variables are used to identify the cell number for ghost cells
      integer,dimension(:,:),allocatable :: isendindices,irecvindices

! variable to clean the ghost cells
      logical,dimension(:),allocatable :: ighost_updated

! variables used to read initial particle properties
      double precision, dimension(:,:),allocatable:: dpar_pos,dpar_vel
      double precision, dimension(:),allocatable::dpar_den,dpar_rad
! Variables used for reading restart file
      integer,dimension(:),allocatable:: irestartmap
      integer,dimension(:),allocatable  :: itempglobal_id

! generic interface definition
      interface des_gather
         module procedure des_gather_l,des_gather_i,des_gather_d
      end interface
      interface des_readscatter
         module procedure des_readscatter_l,des_readscatter_i,des_readscatter_d
      end interface
      interface des_scatter
         module procedure des_scatter_l,des_scatter_i,des_scatter_d
      end interface


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
! Subroutine       : desmpi_init
! Purpose          : allocates and initializes the variables related to send and recv
!                    sets the flag related to periodic and processor interfaces
!------------------------------------------------------------------------
      subroutine desmpi_init
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lpacketsize,ltordimn,lfaces,lfactor=4
      integer :: lmaxlen1,lmaxlen2,lmaxarea,lmaxghostpar
!-----------------------------------------------

! determine initial size of send and recv buffer based on max_pip,
! total cells max number of boundary cells, and the packet size
      ltordimn = merge(1,3,NO_K)
      lpacketsize = 2*dimn + ltordimn+ 5
      lfaces = dimn*2
      lmaxlen1 = dg_iend2-dg_istart2+1
      lmaxlen2 = dg_jend2-dg_jstart2+1
      if (do_K) then
         lmaxlen1 = max(lmaxlen1,dg_kend2 -dg_kstart2+1)
         lmaxlen2 = max(lmaxlen2,dg_kend2 -dg_kstart2+1)
      else
         lmaxlen1 = max(lmaxlen1,lmaxlen2)
         lmaxlen2 = 1
      end if
      lmaxarea = lmaxlen1*lmaxlen2 + 10 !10 is added for buffer and is required for send and recv indices
      lmaxghostpar = (max_pip/dg_ijksize2)* lfactor
      if(lmaxghostpar.lt.100) lmaxghostpar = 100
      imaxbuf = lmaxghostpar*lmaxarea*lpacketsize

      allocate (isendindices(lmaxarea,lfaces),irecvindices(lmaxarea,lfaces))
      isendindices =0
      irecvindices=0

      allocate (dsendbuf(imaxbuf,lfaces),drecvbuf(imaxbuf,lfaces), &
                isendreq(lfaces),irecvreq(lfaces),isendcnt(lfaces))
      dsendbuf=0.0
      drecvbuf=0.0
      isendreq =0
      irecvreq=0
      isendcnt=0

      allocate (dcycl_offset(lfaces,dimn),ineighproc(lfaces),iexchflag(lfaces))
      ineighproc=0
      iexchflag=.false.
      dcycl_offset=0.0

! allocate variables related to scattergather
      allocate(iscattercnts(0:numpes-1),idispls(0:numpes-1),igathercnts(0:numpes-1))
      iscattercnts=0
      idispls=0
      igathercnts=0

! allocate variables related to ghost particles
      allocate(ighost_updated(max_pip))

! call node exchange init in case
! this could be needed if des_interp_on is true (i.e., drag is interpolated)
! or DES_INTERP_MEAN_FIELDS is true (i.e., mean fields are interpolated)
      IF(DES_INTERP_ON.OR.DES_INTERP_MEAN_FIELDS) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,'(/,5x,A,/,5x,A,/)') 'In desmpi_mod, &
         &setting the node indices &
         &for MPI communication', 'of nodal information needed for backward &
         &interpolation'
         IF(myPE.eq.pe_IO) WRITE(*, '(/,5x,A,/,5x,A,/)') 'In desmpi_mod, &
         &setting the node indices &
         &for MPI communication', 'of nodal information needed for backward &
         &interpolation'
         call des_setnodeindices
         ENDIF


! set the communication flags
      call desmpi_setcomm

!      call des_dbgmpi(1)

      end subroutine desmpi_init


!------------------------------------------------------------------------
! Subroutine       : desmpi_setcomm
! Purpose          : sets the flags required for interprocessor communication
!------------------------------------------------------------------------
      subroutine desmpi_setcomm()
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lijkproc,liproc,ljproc,lkproc
      integer li,lj,lk,lis,lie,ljs,lje,lks,lke,lcount,lface,litmp,ljtmp,lktmp
      integer listart1,liend1,ljstart1,ljend1,lkstart1,lkend1
      integer listart2,liend2,ljstart2,ljend2,lkstart2,lkend2
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
!-----------------------------------------------

! set flags for interprocessor boundaries and set the corresponding to proc
      lijkproc = mype
      liproc = iofproc(lijkproc)
      ljproc = jofproc(lijkproc)
      lkproc = kofproc(lijkproc)
      iexchflag(:) = .false.
      ineighproc(:) = 0
      if(liproc.gt.0) then
         iexchflag(2) = .true.
         ineighproc(2)=procijk(liproc-1,ljproc,lkproc)
      end if
      if(liproc.lt.nodesi-1) then
         iexchflag(1) = .true.
         ineighproc(1)=procijk(liproc+1,ljproc,lkproc)
      end if
      if(ljproc.gt.0) then
         iexchflag(4)= .true.
         ineighproc(4)=procijk(liproc,ljproc-1,lkproc)
      end if
      if(ljproc.lt.nodesj-1) then
         iexchflag(3) = .true.
         ineighproc(3)=procijk(liproc,ljproc+1,lkproc)
      end if
      if(lkproc.gt.0) then
         iexchflag(6)=.true.
         ineighproc(6)=procijk(liproc,ljproc,lkproc-1)
      end if
      if(lkproc.lt.nodesk-1) then
         iexchflag(5) =.true.
         ineighproc(5)=procijk(liproc,ljproc,lkproc+1)
      end if

!set flags for cyclic boundary conditions and corresponding to proc
      dcycl_offset(:,:) = 0
      if (des_periodic_walls_x) then
         if(liproc.eq.0) then
            iexchflag(2)=.true.
            ineighproc(2)= procijk(nodesi-1,ljproc,lkproc)
            dcycl_offset(2,1)= xlength
         end if
         if(liproc.eq.nodesi-1) then
            iexchflag(1)=.true.
            ineighproc(1)= procijk(0,ljproc,lkproc)
            dcycl_offset(1,1)=-xlength
         end if
      end if
      if (des_periodic_walls_y) then
         if(ljproc.eq.0) then
            iexchflag(4)=.true.
            ineighproc(4)= procijk(liproc,nodesj-1,lkproc)
            dcycl_offset(4,2)= ylength
         end if
         if(ljproc.eq.nodesj-1) then
            iexchflag(3)=.true.
            ineighproc(3)= procijk(liproc,0,lkproc)
            dcycl_offset(3,2)=-ylength
         end if
      end if
      if (des_periodic_walls_z) then
         if(lkproc.eq.0) then
            iexchflag(6)=.true.
            ineighproc(6)= procijk(liproc,ljproc,nodesk-1)
            dcycl_offset(6,3)=zlength
         end if
         if(lkproc.eq.nodesk-1) then
            iexchflag(5)=.true.
            ineighproc(5)= procijk(liproc,ljproc,0)
            dcycl_offset(5,3)=-zlength
         end if
      end if

! For mass inlet and outlet, the cells where the particle injected are
! considered as part of the domain; This avoids flagging newly injected particles
! and outgoing particles as ghost particles PEA(:,4)
      listart1=dg_istart1;liend1=dg_iend1;listart2=dg_istart2;liend2=dg_iend2
      ljstart1=dg_jstart1;ljend1=dg_jend1;ljstart2=dg_jstart2;ljend2=dg_jend2
      lkstart1=dg_kstart1;lkend1=dg_kend1;lkstart2=dg_kstart2;lkend2=dg_kend2
!      if (dem_mio) then
!         if (dem_mi_x .or. dem_mo_x) then
!            if(listart1.eq.dg_imin1) listart1 = dg_imin1-1
!            if(liend1.eq.dg_imax1) liend1 = dg_imax1+1
!         end if
!         if (dem_mi_y .or. dem_mo_y) then
!            if(ljstart1.eq.dg_jmin1) ljstart1 = dg_jmin1-1
!            if(ljend1.eq.dg_jmax1) ljend1 = dg_jmax1+1
!         end if
!         if (dem_mi_z .or. dem_mo_z) then
!            if(lkstart1.eq.dg_kmin1) lkstart1 = dg_kmin1-1
!            if(lkend1.eq.dg_kmax1) lkend1 = dg_kmax1+1
!         end if
!      end if

! set the ghost cell indices for e-w, n-s and t-b
! for east and west faces
      lks = lkstart1
      lke = lkend1
      ljs = ljstart1
      lje = ljend1

!east face
      lface = 1
      li = liend1
      lcount = 1
      do lk = lks,lke
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li+1,lj,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

!west face
      lface = 2
      li = listart1
      lcount = 1
      do lk = lks,lke
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li-1,lj,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

! for north and south faces
      lks = lkstart1
      lke = lkend1
      lis = listart2
      lie = liend2

!north face
      lface = 3
      lj = ljend1
      lcount = 1
      do lk = lks,lke
         do li = lis,lie
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj+1,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1


!south face
      lface = 4
      lj = ljstart1
      lcount = 1
      do lk = lks,lke
         do li = lis,lie
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj-1,lk)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

! for top and bottom
      if (no_k) return
      lis = listart2
      lie = liend2
      ljs = ljstart2
      lje = ljend2

!top face
      lface = 5
      lk = lkend1
      lcount = 1
      do li = lis,lie
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj,lk+1)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1


!bottom face
      lface = 6
      lk = lkstart1
      lcount = 1
      do li = lis,lie
         do lj = ljs,lje
            lcount = lcount + 1
            isendindices(lcount,lface) = dg_funijk(li,lj,lk)
            irecvindices(lcount,lface) = dg_funijk(li,lj,lk-1)
         end do
      end do
      isendindices(1,lface) = lcount - 1
      irecvindices(1,lface) = lcount - 1

      return
      end subroutine desmpi_setcomm

!------------------------------------------------------------------------
! Subroutine       : desmpi_sendrecv_init
! Purpose          : posts asynchronous send and recv and updates the request id
!
! Parameters       : pface - face number (1to6)
!                    debug - for printing debug statments
!------------------------------------------------------------------------
      subroutine desmpi_sendrecv_init(pface,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer,intent(in) :: pface
      integer,intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'desmpi_sendrecv_init'
      integer :: ldebug,ltag,lerr,lsource,ldest,lrecvface
      integer :: message_tag
!-----------------------------------------------

      message_tag(lsource,ldest,lrecvface) = lsource+numpes*ldest+numpes*numpes*lrecvface+100

! set the debug flag
      ldebug = 0
      if (present(pdebug)) then
        ldebug = pdebug
      endif

!direct copy in case of single processor
      lrecvface = pface+mod(pface,2)-mod(pface+1,2)
      if (ineighproc(pface).eq.mype) then
         drecvbuf(1:isendcnt(pface),lrecvface)= dsendbuf(1:isendcnt(pface),pface)
      else
         ltag = message_tag(ineighproc(pface),mype,pface)
         call des_mpi_irecv(drecvbuf(:,pface),imaxbuf, &
                            ineighproc(pface),ltag,irecvreq(pface),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )

         ltag = message_tag(mype,ineighproc(pface),lrecvface)
         call des_mpi_isend(dsendbuf(:,pface),isendcnt(pface), &
                        ineighproc(pface),ltag,isendreq(pface),lerr)
         call mpi_check( name //':mpi_isend ', lerr )
      end if
      return
      end subroutine desmpi_sendrecv_init

!------------------------------------------------------------------------
! Subroutine       : desmpi_sendrecv_wait
! Purpose          : waits for the communication for the specified interface
!
! Parameters       : pface - face number (1to6)
!                    debug - for printing debug statments
!------------------------------------------------------------------------
      subroutine desmpi_sendrecv_wait(pface,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer,intent(in) :: pface
      integer,intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'desmpi_sendrecv_wait'
      integer :: ldebug,ltag,lerr,lsource,ldest
!-----------------------------------------------

! set the debug flag
      ldebug = 0
      if (present(pdebug)) then
        ldebug = pdebug
      endif

! wait for both send and recv request completes
      if (ineighproc(pface).ne.mype) then
         call des_mpi_wait(isendreq(pface),lerr)
         call mpi_check( name //':mpi_wait-send', lerr )
         call des_mpi_wait(irecvreq(pface),lerr)
         call mpi_check( name //':mpi_wait-recv', lerr )
      end if
      return
      end subroutine desmpi_sendrecv_wait


!------------------------------------------------------------------------
! Subroutine       : desmpi_scatterv
! Purpose          : scatters the particle from PE_IO
! Parameters       : ptype - flag for datatype integer (1) or double precision (2)
!                    pdebug - optional flag for debugging
!------------------------------------------------------------------------
      subroutine desmpi_scatterv(ptype,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: ptype
      integer, intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lroot,lidebug,lerr
      character(len=80), parameter :: name = 'desmpi_scatterv'
!-----------------------------------------------

      lroot = pe_io
      if (.not. present(pdebug)) then
         lidebug = 0
      else
         lidebug = pdebug
      endif

      if (ptype .eq. 1) then
         call des_MPI_Scatterv(irootbuf,iscattercnts,idispls, &
                               iprocbuf,iscr_recvcnt,lroot,lerr)
      else
         call des_MPI_Scatterv(drootbuf,iscattercnts,idispls, &
                               dprocbuf,iscr_recvcnt,lroot,lerr)
      end if
      call MPI_Check( name //':MPI_Scatterv', lerr )

      return
      end subroutine desmpi_scatterv


!------------------------------------------------------------------------
! Subroutine       : desmpi_gatherv
! Purpose          : gathers the particle from local proc to root proc
! Parameters       : ptype - flag for datatype integer (1) or double precision (2)
!                    pdebug - optional flag for debugging
!------------------------------------------------------------------------
      subroutine desmpi_gatherv(ptype,pdebug)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: ptype
      integer, intent(in),optional :: pdebug
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lroot,lidebug,lerr
      character(len=80), parameter :: name = 'des_gather'
!-----------------------------------------------

      lroot = pe_io
      if (.not. present(pdebug)) then
         lidebug = 0
      else
         lidebug = pdebug
      endif
      if(ptype.eq.1) then
         call des_MPI_Gatherv(iprocbuf,igath_sendcnt,irootbuf, &
                              igathercnts,idispls,lroot,lerr)
      else
         call des_MPI_Gatherv(dprocbuf,igath_sendcnt,drootbuf, &
                              igathercnts,idispls,lroot,lerr)
      end if
      call MPI_Check( name //':MPI_Gatherv', lerr )

      return
      end subroutine desmpi_gatherv


!------------------------------------------------------------------------
! Subroutine       : des_scatter_particle
! Purpose          : Scatters the particles read from input file or
!                    generated based on specified volume fraction
! Comments         : In the main thread (pe_io) particles will be seperated
!                    based on the location and it will be packed in drootbuf
!                    Each proc receives its particles along with particle count
!------------------------------------------------------------------------
      subroutine des_scatter_particle
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lcurpar,lproc,lbuf,lpacketsize,lface
      integer lproc_parcnt(0:numpes-1),lpar_proc(particles)
!-----------------------------------------------

      integer :: rdimn

      rdimn = merge(2,3, NO_K)

! set the packet size for transfer
      lpacketsize = 2*rdimn + 2

! build the send buffer in PE_IO proc
! first pass to get the count of particles
      lpar_proc(:) =-1
      lproc_parcnt(:) = 0
      if(myPE.eq.pe_io) then
         if (no_k) then
            do lcurpar = 1,particles
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.xe(istart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,1).lt.xe(iend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,2).ge.yn(jstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,2).lt.yn(jend1_all(lproc))) then
                     lpar_proc(lcurpar) = lproc
                     lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                     exit
                  endif
               enddo
               if (lpar_proc(lcurpar).eq.-1) then
                  WRITE(*,500) lcurpar
                  call des_mpi_stop
               endif
            enddo
         else
            do lcurpar = 1,particles
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.xe(istart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,1).lt.xe(iend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,2).ge.yn(jstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,2).lt.yn(jend1_all(lproc))     &
                   .and. dpar_pos(lcurpar,3).ge.zt(kstart1_all(lproc)-1) &
                   .and. dpar_pos(lcurpar,3).lt.zt(kend1_all(lproc))) then
                     lpar_proc(lcurpar) = lproc
                     lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
                     exit
                  end if
               end do
               if (lpar_proc(lcurpar).eq.-1) then
                  WRITE(*,501) lcurpar
                  call des_mpi_stop
               endif
            enddo
         endif  ! if (no_k)
      endif ! if (my_pe.eq.pe_io)
      call bcast(lproc_parcnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lproc_parcnt(mype)
      if (pip .gt. max_pip) then
         WRITE(*,502) pip, max_pip
         call des_mpi_stop
      endif
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(particles*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lproc_parcnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lproc_parcnt(lproc)*lpacketsize
         end do
         lproc_parcnt(:) = 0
         do lcurpar = 1,particles
            lproc = lpar_proc(lcurpar)
            lbuf = idispls(lproc)+lproc_parcnt(lproc)*lpacketsize+1
            drootbuf(lbuf) = dpar_rad(lcurpar); lbuf = lbuf + 1
            drootbuf(lbuf) = dpar_den(lcurpar); lbuf = lbuf + 1
            drootbuf(lbuf:lbuf+rdimn-1) = dpar_pos(lcurpar,1:rdimn); lbuf = lbuf + rdimn
            drootbuf(lbuf:lbuf+rdimn-1) = dpar_vel(lcurpar,1:rdimn); lbuf = lbuf + rdimn
            lproc_parcnt(lproc) = lproc_parcnt(lproc) + 1
         enddo
      endif
      call desmpi_scatterv(ptype=2)

! unpack the particles in each processor and set the pip
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize+1
         des_radius(lcurpar) = dprocbuf(lbuf); lbuf = lbuf+1
         ro_sol(lcurpar) = dprocbuf(lbuf); lbuf = lbuf+1
         des_pos_new(1:rdimn,lcurpar) = dprocbuf(lbuf:lbuf+rdimn-1); lbuf = lbuf+rdimn
         des_vel_new(1:rdimn,lcurpar) = dprocbuf(lbuf:lbuf+rdimn-1); lbuf = lbuf+rdimn
         pea(lcurpar,1) = .true.
      enddo
      deallocate (dprocbuf,drootbuf)

 500  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 501  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 502  FORMAT(/2X,'From: DES_SCATTER_PARTICLE: ',/2X,&
         'ERROR: Particles in the processor ',I10,&
         'exceeds MAX_PIP', I10)

      RETURN
      end subroutine des_scatter_particle


!------------------------------------------------------------------------
! Subroutine       : des_restart_map
! Purpose          : Scatters the particles read from restart file and
!                    generates mapping which will be used in des_readscatter
!                    routines
! Parameters       : pglocnt - total particles read from restart file
!------------------------------------------------------------------------
      subroutine des_restart_map (pglocnt)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer pglocnt
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lcurpar,lproc,lbuf,lpacketsize,lface
      integer lpar_cnt(0:numpes-1)
      double precision,dimension(0:numpes-1):: lxmin,lxmax,lymin,lymax,lzmin,lzmax
!-----------------------------------------------

! set the domain range for each processor
      do lproc= 0,numpes-1
         lxmin(lproc) = xe(istart1_all(lproc)-1)
         lxmax(lproc) = xe(iend1_all(lproc))
         lymin(lproc) = yn(jstart1_all(lproc)-1)
         lymax(lproc) = yn(jend1_all(lproc))
         lzmin(lproc) = zt(kstart1_all(lproc)-1)
         lzmax(lproc) = zt(kend1_all(lproc))
      end do
! modify the range for mass inlet and outlet, as particles injected
! can lie outside the domain and not ghost particles
      if (dem_mio) then
         do lproc= 0,numpes-1
!            if(des_mi_x .or. des_mo_x) then
               if(istart1_all(lproc).eq.imin1) then
                  lxmin(lproc) = xe(istart1_all(lproc)-2)
               endif
               if(iend1_all(lproc).eq.imax1) then
                  lxmax(lproc) = xe(iend1_all(lproc)+1)
               endif
!            end if
!            if(des_mi_y .or. des_mo_y) then
               if(jstart1_all(lproc).eq.jmin1) then
                  lymin(lproc) = yn(jstart1_all(lproc)-2)
               end if
               if(jend1_all(lproc).eq.jmax1) then
                  lymax(lproc) = yn(jend1_all(lproc)+1)
               end if
!            end if
!            if(des_mi_z .or. des_mo_z) then
               if(kstart1_all(lproc).eq.kmin1) then
                  lzmin(lproc) = zt(kstart1_all(lproc)-2)
               end if
               if(kend1_all(lproc).eq.kmax1) then
                  lzmax(lproc) = zt(kend1_all(lproc)+1)
               end if
!            end if
         end do
      end if



! set the packet size for transfer
      lpacketsize = dimn

! build the send buffer in PE_IO proc
! first pass to get the count of particles
      irestartmap(:) = -1
      lpar_cnt(:) = 0
      if(myPE.eq.pe_io) then
         if (no_k) then
            do lcurpar = 1,pglocnt
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.lxmin(lproc) &
                   .and. dpar_pos(lcurpar,1).lt.lxmax(lproc) &
                   .and. dpar_pos(lcurpar,2).ge.lymin(lproc) &
                   .and. dpar_pos(lcurpar,2).lt.lymax(lproc)) then
                     lpar_cnt(lproc) = lpar_cnt(lproc) + 1
                     irestartmap(lcurpar) = lproc
                     exit
                  end if
               end do
               if (irestartmap(lcurpar).eq.-1) then
                  WRITE(*,600) lcurpar
                  call des_mpi_stop
               endif
            enddo
         else
            do lcurpar = 1,pglocnt
               do lproc= 0,numpes-1
                  if (   dpar_pos(lcurpar,1).ge.lxmin(lproc) &
                   .and. dpar_pos(lcurpar,1).lt.lxmax(lproc) &
                   .and. dpar_pos(lcurpar,2).ge.lymin(lproc) &
                   .and. dpar_pos(lcurpar,2).lt.lymax(lproc) &
                   .and. dpar_pos(lcurpar,3).ge.lzmin(lproc) &
                   .and. dpar_pos(lcurpar,3).lt.lzmax(lproc)) then
                     lpar_cnt(lproc) = lpar_cnt(lproc) + 1
                     irestartmap(lcurpar) = lproc
                     exit
                  end if
               end do
               if (irestartmap(lcurpar).eq.-1) then
                  WRITE(*,611) lcurpar
                  irestartmap(lcurpar) = PE_IO

!                  call des_mpi_stop
               endif
            enddo
         endif  ! if (no_k)
      endif ! if (my_pe.eq.pe_io)
      call bcast(lpar_cnt(0:numpes-1),pe_io)

! second pass: set and allocate scatter related variables
      pip = lpar_cnt(mype)
      if (pip .gt. max_pip) then
         WRITE(*,602) pip, max_pip
         call des_mpi_stop
      endif
      iscr_recvcnt = pip*lpacketsize
      allocate (dprocbuf(iscr_recvcnt))
      if (mype.eq.pe_io) then
         allocate (drootbuf(pglocnt*lpacketsize))
      else
         allocate (drootbuf(10))
      endif

! in the IO processor build the drootbuffer and idispls required
! for mpi communication
      if(mype.eq.pe_io) then
         idispls(0) = 0
         iscattercnts(0) = lpar_cnt(0)*lpacketsize
         do lproc = 1,numpes-1
            idispls(lproc) = idispls(lproc-1) + iscattercnts(lproc-1)
            iscattercnts(lproc) = lpar_cnt(lproc)*lpacketsize
         enddo
         lpar_cnt(:) = 0
         do lcurpar = 1,pglocnt
            lproc = irestartmap(lcurpar)
            lbuf = idispls(lproc)+lpar_cnt(lproc)*lpacketsize+1
            drootbuf(lbuf:lbuf+dimn-1) = dpar_pos(lcurpar,1:dimn); lbuf = lbuf + dimn
            lpar_cnt(lproc) = lpar_cnt(lproc) + 1
         enddo
      endif
      call desmpi_scatterv(ptype=2)

! unpack the particles in each processor and set the pip
      do lcurpar = 1,pip
         lbuf = (lcurpar-1)*lpacketsize+1
         des_pos_new(1:dimn,lcurpar) = dprocbuf(lbuf:lbuf+dimn-1); lbuf = lbuf+dimn
         pea(lcurpar,1) = .true.
      enddo
      deallocate (drootbuf,dprocbuf)

 600  FORMAT(/2X,'From: DES_RESTART_MAP: (0)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')
 601  FORMAT(/2X,'From: DES_RESTART_MAP: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain')

 611  FORMAT(/2X,'From: DES_RESTART_MAP: (1)',/2X,&
         'ERROR: Unable to locate the particle (no. ',I10,&
         ') inside the domain'/,2x,'... DELETE IT!')



 602  FORMAT(/2X,'From: DES_RESTART_MAP: ',/2X,&
         'ERROR: Particles in the processor ',I10,&
         'exceeds MAX_PIP', I10)

      RETURN
      end subroutine des_restart_map


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
! Subroutine       : desmpi_pack_ghostpar
! Purpose          : packs the ghost particle in the buffer based on the flag
! Parameter        : face - value from 1 to 6 represents faces
!------------------------------------------------------------------------
      subroutine desmpi_pack_ghostpar(pface)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: pface
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lijk,lindx,ltot_ind,lpicloc,lpar_cnt,lcurpar
      integer :: ltordimn,lpacketsize,lbuf
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
!-----------------------------------------------
      ltordimn = merge(1,3,NO_K)
      lpacketsize = 2*dimn + ltordimn+ 5
      lpar_cnt = 0
      ltot_ind = isendindices(1,pface)
      do lindx = 2,ltot_ind+1
         lijk = isendindices(lindx,pface)
         do lpicloc =1,dg_pic(lijk)%isize
            lbuf = lpar_cnt*lpacketsize+ibufoffset
            lcurpar = dg_pic(lijk)%p(lpicloc)
            if(pea(lcurpar,4) .and. .not.ighost_updated(lcurpar) ) cycle
            dsendbuf(lbuf,pface) = iglobal_id(lcurpar);lbuf = lbuf +1
            dsendbuf(lbuf,pface) = dg_ijkconv(lijk,pface,ineighproc(pface))
            lbuf = lbuf +1
            dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface))
            lbuf = lbuf +1
            dsendbuf(lbuf,pface) = des_radius(lcurpar); lbuf = lbuf + 1
            dsendbuf(lbuf,pface) = pijk(lcurpar,5); lbuf = lbuf + 1
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
            lbuf = lbuf + dimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_new(1:dimn,lcurpar)
            lbuf = lbuf + dimn
            dsendbuf(lbuf:lbuf+ltordimn-1,pface) = omega_new(1:ltordimn,lcurpar)
            lbuf = lbuf + ltordimn
            lpar_cnt = lpar_cnt + 1
         end do
      end do
      dsendbuf(1,pface)=lpar_cnt
      isendcnt(pface) = lpar_cnt*lpacketsize+ibufoffset

      end subroutine desmpi_pack_ghostpar


!------------------------------------------------------------------------
! Subroutine       : desmpi_unpack_ghostpar
! Purpose          : unpacks the ghost particle from the recv buffer
! Parameter        : pface - value from 1 to 6 represents faces
!
!------------------------------------------------------------------------
      subroutine desmpi_unpack_ghostpar(pface)
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
      integer :: ltordimn,lpacketsize,lbuf,lindx,llocpar,lnewcnt,lpicloc
      logical,dimension(:),allocatable :: lfound
      integer,dimension(:),allocatable :: lnewspot,lnewpic
!-----------------------------------------------

! unpack the particles:
! if it already exists update the position
! if not and do_nsearch is true then add to the particle array

      ltordimn = merge(1,3,NO_K)

      lpacketsize = 2*dimn + ltordimn+ 5
      lparcnt = drecvbuf(1,pface)
      lnewcnt = lparcnt
      allocate (lfound(lparcnt),lnewspot(lparcnt),lnewpic(dg_ijksize2))
      lfound(:) = .false.
      lnewspot(:) =0
      lnewpic = 0

      do lcurpar = 1,lparcnt
         lbuf = (lcurpar-1)*lpacketsize+ibufoffset
         lparid  = drecvbuf(lbuf,pface); lbuf = lbuf + 1
         lparijk = drecvbuf(lbuf,pface); lbuf = lbuf + 1
         lprvijk = drecvbuf(lbuf,pface); lbuf = lbuf + 1
! locate the particles first based on previous ijk and then based on current ijk
         lfound(lcurpar) = locate_par(lparid,lprvijk,llocpar)
         if (lparijk .ne. lprvijk .and. .not.lfound(lcurpar)) then
            lfound(lcurpar) = locate_par(lparid,lparijk,llocpar)
         endif
         if(lfound(lcurpar)) then
            dg_pijk(llocpar) = lparijk
            dg_pijkprv(llocpar) = lprvijk
            des_radius(llocpar) = drecvbuf(lbuf,pface)
            lbuf = lbuf + 1
            pijk(llocpar,5) = drecvbuf(lbuf,pface)
            lbuf = lbuf + 1
            des_pos_old(:,llocpar)= des_pos_new(:,llocpar)
            des_vel_old(:,llocpar)= des_vel_new(:,llocpar)
            omega_old(:,llocpar)= omega_new(:,llocpar)
            des_pos_new(1:dimn,llocpar)= drecvbuf(lbuf:lbuf+dimn-1,pface)
            lbuf = lbuf + dimn
            des_vel_new(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
            lbuf = lbuf + dimn
            omega_new(1:ltordimn,llocpar) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
            lbuf = lbuf + ltordimn
            ighost_updated(llocpar) = .true.
            lnewcnt = lnewcnt-1
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
            lparid  = drecvbuf(lbuf,pface); lbuf = lbuf + 1
            lparijk = drecvbuf(lbuf,pface); lbuf = lbuf + 1
            lprvijk = drecvbuf(lbuf,pface); lbuf = lbuf + 1
            do while(pea(ispot,1))
               ispot = ispot + 1
            enddo
            pea(ispot,1) = .true.
            pea(ispot,2) = .false.
            pea(ispot,3) = .false.
            pea(ispot,4) = .true.
            iglobal_id(ispot)  = lparid
            dg_pijk(ispot) = lparijk
            dg_pijkprv(ispot) = lprvijk
            des_radius(ispot) = drecvbuf(lbuf,pface) ; lbuf = lbuf + 1
            pijk(ispot,5) = drecvbuf(lbuf,pface) ; lbuf = lbuf + 1
            des_pos_new(1:dimn,ispot)= drecvbuf(lbuf:lbuf+dimn-1,pface)
            lbuf = lbuf + dimn
            des_vel_new(1:dimn,ispot) = drecvbuf(lbuf:lbuf+dimn-1,pface)
            lbuf = lbuf + dimn
            omega_new(1:ltordimn,ispot) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
            lbuf = lbuf + ltordimn
            ighost_updated(ispot) = .true.
            lnewspot(lcurpar) = ispot
            des_pos_old(1:dimn,ispot) = des_pos_new(1:dimn,ispot)
            des_vel_old(1:dimn,ispot) = des_vel_new(1:dimn,ispot)
            omega_old(1:ltordimn,ispot) = omega_new(1:ltordimn,ispot)
         enddo
      endif

!deallocate temporary variablies
      deallocate (lfound,lnewspot,lnewpic)

      end subroutine desmpi_unpack_ghostpar

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
               pn(:,lcurpar) = 0 ; pv(:,lcurpar) = .false.
               pft(lcurpar,:,:) = 0
               des_pos_new(:,lcurpar)=0
               des_pos_old(:,lcurpar)=0
               des_vel_new(:,lcurpar)=0
               des_vel_old(:,lcurpar)=0
               omega_new(:,lcurpar)=0
               neighbours(lcurpar,:)=0
            end do
         end do
      end do
      end subroutine desmpi_cleanup


!------------------------------------------------------------------------
! Subroutine       : desmpi_pack_parcross
! Purpose          : packs the particle crossing the boundary
! Parameter        : pface - value from 1 to 6 represents faces
!------------------------------------------------------------------------
      subroutine desmpi_pack_parcross(pface)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: pface
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: ltot_ind,lindx,ijk
      integer :: lneighindx,lcontactindx,lneigh,lcontact,lijk,&
                 lpicloc,lparcnt,lcurpar
      integer :: lpacketsize,lbuf,ltordimn,ltmpbuf
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
      INCLUDE '../function.inc'
!-----------------------------------------------

! pack the particle crossing the boundary
      ltordimn = merge(1,3,NO_K)
      lpacketsize = 9*dimn + ltordimn*4 + maxneighbors * (dimn+5) + 15
      ltot_ind = irecvindices(1,pface)
      lparcnt = 0
      do lindx = 2,ltot_ind+1
         lijk = irecvindices(lindx,pface)
         do lpicloc = 1,dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)
            if (pea(lcurpar,4)) cycle ! if ghost particle then cycle
            lbuf = lparcnt*lpacketsize + ibufoffset
            dsendbuf(lbuf,pface) = iglobal_id(lcurpar)
            lbuf = lbuf+1
            dsendbuf(lbuf,pface) = dg_ijkconv(lijk,pface,ineighproc(pface))
            lbuf = lbuf+1
            dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface))
            lbuf = lbuf+1
            dsendbuf(lbuf,pface) = des_radius(lcurpar)  ;lbuf = lbuf+1
            li = pijk(lcurpar,1) + icycoffset(pface,1)
            lj = pijk(lcurpar,2) + icycoffset(pface,2)
            lk = pijk(lcurpar,3) + icycoffset(pface,3)
            dsendbuf(lbuf,pface) = li ; lbuf = lbuf+1
            dsendbuf(lbuf,pface) = lj ; lbuf = lbuf+1
            dsendbuf(lbuf,pface) = lk ; lbuf = lbuf+1
            dsendbuf(lbuf,pface) = funijk_proc(li,lj,lk,ineighproc(pface))
            lbuf = lbuf+1
            dsendbuf(lbuf,pface) = pijk(lcurpar,5) ;lbuf = lbuf+1
!            dsendbuf(lbuf:lbuf+1,pface) = pea(lcurpar,2:3);lbuf=lbuf+2
            dsendbuf(lbuf:lbuf+1,pface) = 0
            if (pea(lcurpar,2)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            if (pea(lcurpar,3)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            dsendbuf(lbuf,pface) = ro_sol(lcurpar)      ;lbuf = lbuf+1
            dsendbuf(lbuf,pface) = pvol(lcurpar)        ;lbuf = lbuf+1
            dsendbuf(lbuf,pface) = pmass(lcurpar)       ;lbuf = lbuf+1
            dsendbuf(lbuf,pface) = omoi(lcurpar)        ;lbuf = lbuf+1
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_old(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_old(1:dimn,lcurpar)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_new(1:dimn,lcurpar)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+ltordimn-1,pface) = omega_old(1:ltordimn,lcurpar)
            lbuf = lbuf+ltordimn
            dsendbuf(lbuf:lbuf+ltordimn-1,pface) = omega_new(1:ltordimn,lcurpar)
            lbuf = lbuf+ltordimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_acc_old(1:dimn,lcurpar)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+ltordimn-1,pface) = rot_acc_old(1:ltordimn,lcurpar)
            lbuf = lbuf+ltordimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = fc(:,lcurpar)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+ltordimn-1,pface) = tow(1:ltordimn,lcurpar)
            lbuf = lbuf+ltordimn

! build the neighbour with global number and current and previous pijk
            dsendbuf(lbuf,pface) = neighbours(lcurpar,1)
            ltmpbuf = lbuf+1
            do lneighindx = 2,neighbours(lcurpar,1)+1
               lneigh = neighbours(lcurpar,lneighindx)
               dsendbuf(ltmpbuf,pface) = iglobal_id(lneigh)
               ltmpbuf = ltmpbuf+1
               dsendbuf(ltmpbuf,pface) = dg_ijkconv(dg_pijk(lneigh),pface,ineighproc(pface))
               ltmpbuf = ltmpbuf+1
               dsendbuf(ltmpbuf,pface) = dg_ijkconv(dg_pijkprv(lneigh),pface,ineighproc(pface))
               ltmpbuf = ltmpbuf+1
            enddo

            lbuf = lbuf+3*maxneighbors
! build contact list with global number
            dsendbuf(lbuf,pface) = pn(1,lcurpar);ltmpbuf=lbuf+1
            do lcontactindx = 2,pn(1,lcurpar)+1
               lcontact = pn(lcontactindx,lcurpar)
               dsendbuf(ltmpbuf,pface) = iglobal_id(lcontact)
               ltmpbuf=ltmpbuf+1
               dsendbuf(ltmpbuf,pface) = merge(1,0,pv(lcontactindx,lcurpar))
               ltmpbuf=ltmpbuf+1
               dsendbuf(ltmpbuf:ltmpbuf+dimn-1,pface) = pft(lcurpar,lcontactindx,1:dimn)
               ltmpbuf=ltmpbuf+dimn
            enddo
            lbuf = lbuf+(2+dimn)*maxneighbors
! Pradeep remove ********************
!             print * ,"--------------------------------------------------------------------"
!             print * , "inside pack par cross" ,mype
!             print *, "particle moved is", iglobal_id(lcurpar),des_pos_new(:,lcurpar)
!             print *, "number of neighbours",neighbours(lcurpar,1)
!             do lneighindx = 2,neighbours(lcurpar,1)+1
!             print *, "neighbourid" , iglobal_id(neighbours(lcurpar,lneighindx))
!             end do
!             print *, "number of contacts ", pn(1,lcurpar)
!             do lcontactindx = 2,pn(1,lcurpar)+1
!             print *, "contactid" , iglobal_id(pn(lcontactindx,lcurpar))
!             print *, "contact status", pv(lcontactindx,lcurpar)
!             print *, "accumulated force ",pft(lcurpar,lcontactindx,1),pft(lcurpar,lcontactindx,2)
!             end do
!             print * ,"--------------------------------------------------------------------"

! In case of mppic remove the particles else
! Convert the particle as ghost and set the forces zero
            if (mppic) then
               pea(lcurpar,1:4) = .false.
            else
               pea(lcurpar,4) = .true.
               ighost_cnt = ighost_cnt + 1
            end if
            fc(:,lcurpar) = 0.
            neighbours(lcurpar,:)=0
            pn(:,lcurpar) = 0
            pv(:,lcurpar) = .false.
            pft(lcurpar,:,:) = 0

            lparcnt = lparcnt + 1
         end do
      end do
      dsendbuf(1,pface) = lparcnt
      isendcnt(pface) = lparcnt*lpacketsize+ibufoffset

! following unused variables are not sent across the processor
! well_depth
! is_linked
! links
! part_grid

      end subroutine desmpi_pack_parcross


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
      integer :: lpacketsize,lbuf,ltordimn,ltmpbuf,lcount
      logical :: lcontactfound,lneighfound
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
!-----------------------------------------------

! loop through particles and locate them and make changes
      ltordimn = merge(1,3,NO_K)
      lpacketsize = 9*dimn + ltordimn*4 + maxneighbors * (dimn+5) + 15
      lparcnt = drecvbuf(1,pface)

! if mppic make sure enough space available
      if(mppic .and. (max_pip-pip).lt.lparcnt) call redim_par(pip+lparcnt)

      do lcurpar =1,lparcnt
         lfound = .false.
         lbuf = (lcurpar-1)*lpacketsize + ibufoffset
         lparid  = drecvbuf(lbuf,pface)
         lbuf = lbuf+1
         lparijk = drecvbuf(lbuf,pface)
         lbuf = lbuf+1
         lprvijk = drecvbuf(lbuf,pface)
         lbuf = lbuf+1

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
         des_radius(llocpar)  = drecvbuf(lbuf,pface)
         lbuf = lbuf + 1
         pijk(llocpar,1:5)    = drecvbuf(lbuf:lbuf+4,pface)
         lbuf = lbuf+5
!         pea(llocpar,2:3)     = drecvbuf(lbuf:lbuf+1,pface) ; lbuf=lbuf+2
         pea(llocpar,2:3) = .false.
         if (drecvbuf(lbuf,pface).eq.1) pea(llocpar,2) = .true. ; lbuf = lbuf + 1
         if (drecvbuf(lbuf,pface).eq.1) pea(llocpar,3) = .true. ; lbuf = lbuf + 1
         ro_sol(llocpar)      = drecvbuf(lbuf,pface)
         lbuf = lbuf + 1
         pvol(llocpar)        = drecvbuf(lbuf,pface)
         lbuf = lbuf + 1
         pmass(llocpar)       = drecvbuf(lbuf,pface)
         lbuf = lbuf + 1
         omoi(llocpar)        = drecvbuf(lbuf,pface)
         lbuf = lbuf + 1
         des_pos_old(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         des_pos_new(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         des_vel_old(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         des_vel_new(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         omega_old(1:ltordimn,llocpar) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
         lbuf = lbuf + ltordimn
         omega_new(1:ltordimn,llocpar) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
         lbuf = lbuf + ltordimn
         des_acc_old(1:dimn,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         rot_acc_old(1:ltordimn,llocpar) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
         lbuf = lbuf + ltordimn
         fc(:,llocpar) = drecvbuf(lbuf:lbuf+dimn-1,pface)
         lbuf = lbuf + dimn
         tow(1:ltordimn,llocpar) = drecvbuf(lbuf:lbuf+ltordimn-1,pface)
         lbuf = lbuf + ltordimn

! get the neighbour id based on its global number, current and previous pijk and extensive search
         neighbours(llocpar,1) = drecvbuf(lbuf,pface)
         ltmpbuf=lbuf +1
         lcount = 0
         do lneighindx = 2,neighbours(llocpar,1)+1
            lneighfound = .false.
            lneighid = drecvbuf(ltmpbuf,pface)
            ltmpbuf = ltmpbuf+1
            lneighijk = drecvbuf(ltmpbuf,pface)
            ltmpbuf = ltmpbuf+1
            lneighprvijk = drecvbuf(ltmpbuf,pface)
            ltmpbuf = ltmpbuf+1
            lneighfound = locate_par(lneighid,lneighprvijk,lneigh)
            if (.not.lneighfound) lneighfound = locate_par(lneighid,lneighijk,lneigh)
            if (.not.lneighfound) lneighfound = exten_locate_par(lneighid,lparijk,lneigh)
            if (.not.lneighfound) then
               WRITE(*,701)
               cycle
            endif
            lcount = lcount+1
            neighbours(llocpar,lcount+1) = lneigh
         enddo
         neighbours(llocpar,1) = lcount
         lbuf = lbuf+3*maxneighbors

! loop through contact list and find local particle number using neighbor list
         pn(1,llocpar) = drecvbuf(lbuf,pface);ltmpbuf=lbuf+1
         pv(1,llocpar) = .false.
         lcount = 0
         do lcontactindx = 2,pn(1,llocpar)+1
            lcontactfound = .false.
            lcontactid = drecvbuf(ltmpbuf,pface)
            ltmpbuf=ltmpbuf+1
            do lneighindx = 2,neighbours(llocpar,1)+1
               if (iglobal_id(neighbours(llocpar,lneighindx)).eq.lcontactid) then
                  lcontact = neighbours(llocpar,lneighindx)
                  lcontactfound = .true.
                  exit
               endif
            enddo
            if (.not.lcontactfound) then
!check for wall contact and if not print warning message
               if(lcontactid .lt. 0) then
                  lcontact = max_pip + (-1) * lcontactid
               else
                  WRITE(*,702) lcontactid
                  ltmpbuf = ltmpbuf + 1 + dimn ! necessary as pv and pft not yet read for this particle
                  cycle
               endif
            endif
            lcount = lcount+1
            pn(lcount+1,llocpar) = lcontact
            pv(lcount+1,llocpar) = merge(.true.,.false.,drecvbuf(ltmpbuf,pface).gt.0.5)
            ltmpbuf=ltmpbuf+1
            pft(llocpar,lcount+1,1:dimn) = drecvbuf(ltmpbuf:ltmpbuf+dimn-1,pface)
            ltmpbuf=ltmpbuf+dimn
         enddo
         pn(1,llocpar)=lcount
         lbuf = lbuf+(2+dimn)*maxneighbors

!debuging print PRADEEP remove
!          print * ,"--------------------------------------------------------------------"
!          print * , "inside unpack par cross"
!          print *, "particle moved is", iglobal_id(llocpar)
!          print *, "number of neighbours",neighbours(llocpar,1)
!          do lneighindx = 2,neighbours(llocpar,1)+1
!          print *, "neighbourid" , iglobal_id(neighbours(llocpar,lneighindx))
!          end do
!          print *, "number of contacts ", pn(1,llocpar)
!          do lcontactindx = 2,pn(1,llocpar)+1
!          print *, "contactid" , iglobal_id(pn(lcontactindx,llocpar))
!          print *, "contact status", pv(lcontactindx,llocpar)
!          print *, "accumulated force ",pft(llocpar,lcontactindx,1),pft(llocpar,lcontactindx,2)
!          end do
!          print * ,"--------------------------------------------------------------------"

      end do

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
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
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
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
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


!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues_mean_fields
! Purpose          : This routine is specially used for computing mean
!                    fields by backward interpolation.
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues_mean_fields()
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lm,ijk,lface,lijkmin,lijkmax
      integer :: linode,ljnode,lknode,lijknode
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! fill the temporary buffer
      DO LM = 1,DES_MMAX
         CALL DES_EXCHANGENODE(DES_ROPS_NODE(:,LM),PADD=.TRUE.)
         DO LI =1,DIMN
            CALL DES_EXCHANGENODE(DES_VEL_NODE(:,LI,LM),PADD=.TRUE.)
         END DO
      END DO

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            des_rops_node(lijkmin,:)  = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_vel_node(lijkmin,:,:) = des_vel_node(lijkmin,:,:)+des_vel_node(lijkmax,:,:)
            des_rops_node(lijkmax,:)  = des_rops_node(lijkmin,:)
            des_vel_node(lijkmax,:,:) = des_vel_node(lijkmin,:,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues_mean_fields



!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues
! Purpose          : This routine is specially used for des_drag_gs
!                    The backward interpolation in des_drag_gs computes
!                    the grid node values of drag_am and drag_bm
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues()
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lm,ijk,lface,lijkmin,lijkmax
      integer :: linode,ljnode,lknode,lijknode
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! fill the temporary buffer
      do lm = 1,DES_MMAX
         call des_exchangenode(drag_am(:,lm),padd=.true.)
         do li =1,dimn
            call des_exchangenode(drag_bm(:,li,lm),padd=.true.)
         end do
      end do

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            drag_am(lijkmin,:) = drag_am(lijkmin,:)+drag_am(lijkmax,:)
            drag_bm(lijkmin,:,:) = drag_bm(lijkmin,:,:)+drag_bm(lijkmax,:,:)
            drag_am(lijkmax,:) = drag_am(lijkmin,:)
            drag_bm(lijkmax,:,:) = drag_bm(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            drag_am(lijkmin,:) = drag_am(lijkmin,:)+drag_am(lijkmax,:)
            drag_bm(lijkmin,:,:) = drag_bm(lijkmin,:,:)+drag_bm(lijkmax,:,:)
            drag_am(lijkmax,:) = drag_am(lijkmin,:)
            drag_bm(lijkmax,:,:) = drag_bm(lijkmin,:,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            drag_am(lijkmin,:) = drag_am(lijkmin,:)+drag_am(lijkmax,:)
            drag_bm(lijkmin,:,:) = drag_bm(lijkmin,:,:)+drag_bm(lijkmax,:,:)
            drag_am(lijkmax,:) = drag_am(lijkmin,:)
            drag_bm(lijkmax,:,:) = drag_bm(lijkmin,:,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues


!------------------------------------------------------------------------
! Subroutine       : des_addnodevalues2
! Purpose          : This routine is specially used for calc_des_rop_s
!                    The backward interpolation in calc_des_rop_s computes
!                    the grid node values of des_rops_node
!                    node values are from istart2 to iend1;
!                    hence a separate module is created to exchange
!                    node values
!
! Parameters       : None
!------------------------------------------------------------------------
      subroutine des_addnodevalues2()
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lm,ijk,lface,lijkmin,lijkmax
      integer :: linode,ljnode,lknode,lijknode
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! fill the temporary buffer
      do lm = 1,DES_MMAX
         call des_exchangenode(des_rops_node(:,lm),padd=.true.)
      end do

! adjust for periodic boundaries with no domain decomposition
      if (des_periodic_walls_x .and. nodesi.eq.1) then
         do lk = kstart2,kend2
         do lj = jstart2,jend2
            lijkmin = funijk(1,lj,lk)
            lijkmax = funijk(imax1,lj,lk)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_y .and. nodesj.eq.1) then
         do lk = kstart2,kend2
         do li = istart2,iend2
            lijkmin = funijk(li,1,lk)
            lijkmax = funijk(li,jmax1,lk)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if
      if (des_periodic_walls_z .and. nodesk.eq.1 .and. do_K) then
         do li = istart2,iend2
         do lj = jstart2,jend2
            lijkmin = funijk(li,lj,1)
            lijkmax = funijk(li,lj,kmax1)
            des_rops_node(lijkmin,:) = des_rops_node(lijkmin,:)+des_rops_node(lijkmax,:)
            des_rops_node(lijkmax,:) = des_rops_node(lijkmin,:)
         end do
         end do
      end if

      return

      end subroutine des_addnodevalues2


!------------------------------------------------------------------------
! Subroutine       : des_gather_d
! Purpose          : gathers double precision array from local to root
! Parameters       :
!                    parray - array to be writen
!------------------------------------------------------------------------
      subroutine des_gather_d(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
!-----------------------------------------------

! pack the variables in case of
      lparcount = 1
      lcount = 0
      do lcurpar = 1, max_pip
         if (lparcount.gt.pip) exit
         if (.not. pea(lcurpar,1)) cycle
         lparcount = lparcount +1
         if (pea(lcurpar,4)) cycle
         lcount = lcount + 1
         dprocbuf(lcount) = parray(lcurpar)
      end do
      call desmpi_gatherv(ptype=2)
      end subroutine des_gather_d


!------------------------------------------------------------------------
! Subroutine       : des_gather_l
! Purpose          : gathers logical array from local to root
! Parameters       :
!                    parray - array to be writen
!------------------------------------------------------------------------
      subroutine des_gather_l(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      logical, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
!-----------------------------------------------

! pack the variables in proc buffer
      lparcount = 1
      lcount = 0
      do lcurpar = 1, max_pip
         if (lparcount.gt.pip) exit
         if (.not. pea(lcurpar,1)) cycle
         lparcount = lparcount +1
         if (pea(lcurpar,4)) cycle
         lcount = lcount + 1
         if(parray(lcurpar)) then
            iprocbuf(lcount) = 1
         else
            iprocbuf(lcount) = 0
         end if
      end do
      call desmpi_gatherv(ptype=1)

      end subroutine des_gather_l


!------------------------------------------------------------------------
! Subroutine       : des_gather_i
! Purpose          : gathers integer array from local to root
! Parameters       :
!                    parray - array to be writen
!                    ploc2glb - this flag is used to conver local particle
!                    number into global particle number (used for history
!                    and neighbour terms)
!------------------------------------------------------------------------
      subroutine des_gather_i(parray,ploc2glb)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, dimension(:) :: parray
      logical,optional :: ploc2glb
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar,lparcount,lcount
      logical :: lloc2glb
!-----------------------------------------------

      if (present(ploc2glb)) then
         lloc2glb = ploc2glb
      else
         lloc2glb = .false.
      end if
! pack the variables in proc buffer
      lparcount = 1
      lcount = 0
      if (lloc2glb) then
         do lcurpar = 1, max_pip
            if (lparcount.gt.pip) exit
            if (.not. pea(lcurpar,1)) cycle
            lparcount = lparcount +1
            if (pea(lcurpar,4)) cycle
            lcount = lcount + 1
            if(parray(lcurpar).gt.0) then
               iprocbuf(lcount) = iglobal_id(parray(lcurpar))
            else
               iprocbuf(lcount) = 0
            end if
         end do
      else
         do lcurpar = 1, max_pip
            if (lparcount.gt.pip) exit
            if (.not. pea(lcurpar,1)) cycle
            lparcount = lparcount +1
            if (pea(lcurpar,4)) cycle
            lcount = lcount + 1
            iprocbuf(lcount) = parray(lcurpar)
         end do
      end if
      call desmpi_gatherv(ptype=1)

      end subroutine des_gather_i


!------------------------------------------------------------------------
! Subroutine       : des_scatter_i
! Purpose          : scatters integer array into parray
!                    note: used only for restart read
! Parameters       : parray - destination array
!------------------------------------------------------------------------
      subroutine des_scatter_i(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar
!-----------------------------------------------

! scatter the variables and map into local particles array
      call desmpi_scatterv(ptype=1)
      do lcurpar = 1,pip
         parray(lcurpar) = iprocbuf(lcurpar)
      end do
      end subroutine des_scatter_i


!------------------------------------------------------------------------
! Subroutine       : des_scatter_l
! Purpose          : scatters integer array into logical parray
!                    note: used only for restart read
! Parameters       : parray - destination array
!------------------------------------------------------------------------
      subroutine des_scatter_l(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      logical, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar
!-----------------------------------------------

! scatter the variables and map into local particles array
      call desmpi_scatterv(ptype=1)
      do lcurpar = 1,pip
         if(iprocbuf(lcurpar).eq.1) then
            parray(lcurpar) = .true.
         else
            parray(lcurpar) = .false.
         end if
      end do
      end subroutine des_scatter_l

!------------------------------------------------------------------------
! Subroutine       : des_scatter_d
! Purpose          : scatters double precision array into parray
!                    note: used only for restart read
! Parameters       :
!                    parray - destination array
!------------------------------------------------------------------------
      subroutine des_scatter_d(parray)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lcurpar
!-----------------------------------------------

! scatter the variables and map into local particles array
      call desmpi_scatterv(ptype=2)
      do lcurpar = 1,pip
         parray(lcurpar) = dprocbuf(lcurpar)
      end do
      end subroutine des_scatter_d


!------------------------------------------------------------------------
! Subroutine       : des_readscatter_i
! Purpose          : reads the file into temporary array and using irestartmap
!                    build irootbuf and call scatter
! Parameters       : punit - unit number
!                    parray - destination array
!                    ptotsize - total array size at rootproc to be read
!                    pnext_rec - next rec pointer for file
!------------------------------------------------------------------------
      subroutine des_readscatter_i(punit,parray,ptotsize,pnext_rec)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer:: punit,pnext_rec,ptotsize
      integer, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer, dimension(ptotsize) :: larray
      integer,dimension(0:numpes-1) :: lproc_parcnt
      integer :: lcurpar,lproc
!-----------------------------------------------
      if(mype.eq.pe_io) then
         call in_bin_512i(punit,larray,ptotsize,pnext_rec)
         lproc_parcnt(:)= 0
         do lcurpar = 1,ptotsize
            lproc = irestartmap(lcurpar)
            lproc_parcnt(lproc)= lproc_parcnt(lproc)+1
            irootbuf(idispls(lproc)+lproc_parcnt(lproc))=larray(lcurpar)
         end do
      end if
      call des_scatter(parray)
      end subroutine des_readscatter_i


!------------------------------------------------------------------------
! Subroutine       : des_readscatter_l
! Purpose          : reads the file into temporary array and using irestartmap
!                    build irootbuf and call scatter
! Parameters       : punit - unit number
!                    parray - destination array
!                    ptotsize - total array size at rootproc to be read
!                    pnext_rec - next rec pointer for file
!------------------------------------------------------------------------
      subroutine des_readscatter_l(punit,parray,ptotsize,pnext_rec)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer:: punit,pnext_rec,ptotsize
      logical, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer, dimension(ptotsize) :: larray
      integer,dimension(0:numpes-1) :: lproc_parcnt
      integer :: lcurpar,lproc
!-----------------------------------------------
      if(mype.eq.pe_io) then
         call in_bin_512i(punit,larray,ptotsize,pnext_rec)
         lproc_parcnt(:)= 0
         do lcurpar = 1,ptotsize
            lproc = irestartmap(lcurpar)
            lproc_parcnt(lproc)= lproc_parcnt(lproc)+1
            irootbuf(idispls(lproc)+lproc_parcnt(lproc))=larray(lcurpar)
         end do
      end if
      call des_scatter(parray)
      end subroutine des_readscatter_l


!------------------------------------------------------------------------
! Subroutine       : des_readscatter_d
! Purpose          : reads the file into temporary array and using irestartmap
!                    build drootbuf and calls scatter
! Parameters       : punit - unit number
!                    parray - destination array
!                    ptotsize - total array size at rootproc to be read
!                    pnext_rec - next rec pointer for file
!------------------------------------------------------------------------
      subroutine des_readscatter_d(punit,parray,ptotsize,pnext_rec)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer:: punit,pnext_rec,ptotsize
      double precision, dimension(:) :: parray
!-----------------------------------------------
! local variables
!-----------------------------------------------
      double precision, dimension(ptotsize) :: larray
      integer,dimension(0:numpes-1) :: lproc_parcnt
      integer :: lcurpar,lproc
!-----------------------------------------------
      if(mype.eq.pe_io) then
         call in_bin_512(punit,larray,ptotsize,pnext_rec)
         lproc_parcnt(:)= 0
         do lcurpar = 1,ptotsize
            lproc = irestartmap(lcurpar)
            lproc_parcnt(lproc)= lproc_parcnt(lproc)+1
            drootbuf(idispls(lproc)+lproc_parcnt(lproc))=larray(lcurpar)
         end do
      end if
      call des_scatter(parray)
      end subroutine des_readscatter_d


!------------------------------------------------------------------------
! Subroutine       : des_restart_neigh
! Purpose          : restart file contains neighbour information in terms
!                    global id. This routine converts the global id into
!                    local particle number
!                    steps
!                    1. Exchange the ghost particles (does not involve any neighbour info)
!                    2. loop through particles neighbour and contact list
!                       and convert global numbers to local numbers
! Parameters       : none
!------------------------------------------------------------------------

      subroutine des_restart_neigh
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer linter,lface
      integer lparcnt,lcurpar,lneigh,lneighid,lneighindx,lcontact,&
              lcontactid,lcontactindx
      integer lcurijk,lcount
      logical lneighfound,lcontactfound
!-----------------------------------------------
! set do_nsearch true so that the ghost cell will be updated
      do_nsearch = .true.
      call desgrid_pic(plocate=.true.)
      call desmpi_check_sendrecvbuf

!call ghost particle exchange in E-W, N-S, T-B order
      dsendbuf(1,:) = 0; drecvbuf(1,:) =0
      ighost_updated(:) = .false.
      ispot = 1
      do linter = 1,dimn
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface))cycle
            call desmpi_pack_ghostpar(lface)
            call desmpi_sendrecv_init(lface)
         enddo
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface)) cycle
            call desmpi_sendrecv_wait(lface)
            call desmpi_unpack_ghostpar(lface)
         enddo
! update pic required as particles in ghost cell can move between ghost cells
         do lface = linter*2-1,linter*2
            if(dsendbuf(1,lface).gt.0.or.drecvbuf(1,lface).gt.0) then
               call desgrid_pic(plocate=.false.)
               exit
            endif
         enddo
      enddo
      call des_mpi_barrier

! loop through particles neighbour and contact list and find the local particles number
      lparcnt = 1
      do lcurpar =1,max_pip
! pradeep skip ghost particles
         if(lparcnt.gt.pip) exit
         if(.not.pea(lcurpar,1)) cycle
         lparcnt = lparcnt+1
         if(pea(lcurpar,4)) cycle

         lcurijk = dg_pijk(lcurpar)
         lcount = 0
         do lneighindx = 2,neighbours(lcurpar,1)+1
            lneighfound = .false.
            lneighid = neighbours(lcurpar,lneighindx)
            lneighfound = locate_par(lneighid,lcurijk,lneigh)
            if (.not.lneighfound) lneighfound = &
               exten_locate_par(lneighid,lcurijk,lneigh)
            if (.not.lneighfound) then
               WRITE(*,800)
               cycle
            endif
            lcount = lcount + 1
            neighbours(lcurpar,lcount+1) = lneigh
         enddo
         neighbours(lcurpar,1) = lcount

! loop through contact list and find local particle number using neighbor list
         lcount = 0
         do lcontactindx = 2,pn(1,lcurpar)+1
            lcontactfound = .false.
            lcontactid = pn(lcontactindx,lcurpar)
            do lneighindx = 2,neighbours(lcurpar,1)+1
               if (iglobal_id(neighbours(lcurpar,lneighindx)).eq.lcontactid) then
                  lcontact = neighbours(lcurpar,lneighindx)
                  lcontactfound = .true.
                  exit
               endif
            enddo
            if (.not.lcontactfound) then
! check for wall contact and if not print warning message
               if(lcontactid .lt. 0) then
                  lcontact = max_pip + (-1) * lcontactid
               else
                  WRITE(*,801)
                  cycle
               endif
            endif
            lcount = lcount+1
            pn(lcount+1,lcurpar) = lcontact
         enddo
         pn(1,lcurpar) = lcount
      enddo

 800  FORMAT(/2X,'From: DES_RESTART_NEIGH: ',/2X,&
         'WARNING: Unable to locate neighbor during restart (0)',/)
 801  FORMAT(/2X,'From: DES_RESTART_NEIGH: ',/2X,&
         'WARNING: Unable to locate neighbor during restart (1)',/)

      end subroutine des_restart_neigh


!------------------------------------------------------------------------
! Subroutine       : redim_par
! Author           : Pradeep G
! Purpose          : subroutine to redimension the particle when the
!                    array size is not sufficient
! Parameter        : pmaxpip - particle size required
!
!------------------------------------------------------------------------
      subroutine redim_par(pmaxpip)
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer :: pmaxpip
!-----------------------------------------------

      WRITE(*,900) pmaxpip, max_pip
      call des_mpi_stop

 900  FORMAT(/2X,'From: REDIM_PAR: ',/2X,&
         'ERROR: Number of particles ',I10,/2X,&
         'exceeds allowable particles (MAX_PIP)', I10,/2X,&
         'Suggestion: increase PARTICLES_FACTOR in mfix.dat',/2X,&
         'Comment: error may be the result of too many ',&
         'particles moving',/2X,'across processors and/or ',&
         'result of periodic treatment')

      end  subroutine redim_par

!------------------------------------------------------------------------
! subroutine       : des_dbgmpi
! Purpose          : For printing the flags and values set for interface
!                    communication
! Parameters       : ptype - based on this following info is printed to
!                    the file
!                    1 - interface flags
!                    2 - send buffer for ghost particles
!                    3 - recv buffer for ghost particles
!                    4 - particle information
!                    5 - send buffer for particles exchanging processor
!                    6 - particles info
!                    7 - neighinfo
!------------------------------------------------------------------------
      subroutine des_dbgmpi(ptype)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer ptype
!-----------------------------------------------
! local varaiables
!-----------------------------------------------
      character (30) filename
      integer lcurpar,lpacketsize,lface,lparcnt,lbuf,lindx,ltordimn
      integer lcurijk
      integer lneighcnt,lneigh,lneighindx,lcontactcnt,lcontact,lcontactindx
      integer lstart,lsize
      double precision xpos,ypos
      integer li,lj,lparcount
!-----------------------------------------------
! include statement functions
!-----------------------------------------------
      INCLUDE 'desgrid_functions.inc'
!-----------------------------------------------

      write(filename,'("dbg_desmpi",I4.4,".dat")') mype
      open(44,file=filename)
      select case(ptype)
      case (1)
         write(44,*)&
            "------------------------------------------------------"
         write(44,*) "Flag Information"
         do lface =1,dimn*2
            write(44,*) "details for face =" , lface
            write(44,*) "Exchflag, cyclfac, neighproc" ,iexchflag(lface),ineighproc(lface)
         end do
         write(44,*) &
            "------------------------------------------------------"
      case (2)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 2*dimn + ltordimn+ 5
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = dsendbuf(1,lface)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "ghost send buffer for face", lface
               write(44,*) "Number of particles in sendbuf",lparcnt
               write(44,*) "particle number global_id ijk prvijk ",&
                  "radius material new_pos new_vel omega_new"
               write(44,*) &
                "-----------------------------------------------------"
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) lcurpar,(dsendbuf(lindx,lface),lindx=lbuf,lbuf+lpacketsize-1)
               end do
            end if
         end do
      case (3)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 2*dimn + ltordimn+ 5
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = drecvbuf(1,lface)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "ghost recv buffer for face", lface
               write(44,*) "Number of particles in recvbuf",lparcnt
               write(44,*) "particle number global_id ijk prvijk ",&
                  "radius material new_pos new_vel omega_new"
               write(44,*) &
                 "-----------------------------------------------------"
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) lcurpar,(drecvbuf(lindx,lface),lindx=lbuf,lbuf+lpacketsize-1)
               end do
            end if
         end do
      case (4)
          write(44,*) &
             "---------------------------------------------------------"
          write(44,*) "Particle info"
          write(44,*) "max_pip,pip =" , max_pip,pip
          write(44,*) "ghost position                        ",&
             "i       j     k    ijk"
          write(44,*) &
             "---------------------------------------------------------"
          lparcount = 1
          do lcurpar=1,max_pip
             if (lparcount.gt.pip) exit
             if (.not.pea(lcurpar,1))cycle
             lparcount=lparcount + 1
             xpos = des_pos_new(1,lcurpar)
             ypos = des_pos_new(2,lcurpar)
             li=iofpos(xpos);lj=jofpos(ypos)
             write(44,*)pea(lcurpar,4),xpos,ypos,li,lj,dg_funijk(li,lj,1)
          end do
      case (5)
         ltordimn = merge(1,3,NO_K)
         lpacketsize = 9*dimn + ltordimn*4 + maxneighbors * (dimn+5) + 13
         do lface =1,dimn*2
            if (.not.iexchflag(lface))cycle
            lparcnt = dsendbuf(1,lface)
            if (lparcnt .gt. 0) then
               write(44,*) &
                "------------------------------------------------------"
               write(44,*) "particle crossing info send buffer", lface
               write(44,*) "Number of particles in sendbuf",lparcnt
               do lcurpar = 1,lparcnt
                  lbuf = (lcurpar-1) * lpacketsize + ibufoffset
                  write(44,*) "global_id  ijk prvijk radius  i,j,k, ijk"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 8
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "phase density vol mass omoi pos_old"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 5+dimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "pos_new     vel_old   vel_new"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 3*dimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "omega_old     omega_new"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = ltordimn*2
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "acc_old     rot_acc_old   fc "
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 2*dimn + ltordimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

                  write(44,*) "fn ft tow"
                  write(44,*) &
                 "-----------------------------------------------------"
                  lsize = 2*dimn + ltordimn
                  write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                  lbuf = lbuf + lsize

! print neighbour information
                  lneighcnt =dsendbuf(lbuf,lface);lbuf = lbuf + 1
                  write(44,*) "total neighbour=",lneighcnt
                  write(44,*) "neighbou",lneighcnt
                  do lneighindx = 1, lneighcnt
                     lsize = 3
                     write(44,'(5(2x,f8.4))') (dsendbuf(lindx,lface),lindx=lbuf,lbuf+lsize-1)
                     lbuf = lbuf + lsize
                  enddo
               enddo
            endif
         enddo
      case (6)
         write(44,*) "-----------------------------------------------"
         write(44,*) "at Time =",s_time
         write(44,*) "Total paticles =",pip
         write(44,*) "Total ghost paticles =",ighost_cnt
         write(44,*) "do_nsearch =",do_nsearch
         lparcnt = 1
         do lcurpar = 1,max_pip
            if(lparcnt.gt.pip) exit
            lparcnt = lparcnt + 1
            write(44,*) "particle position =",des_pos_new(1:dimn,lcurpar)
         end do
         write(44,*) "-----------------------------------------------"
      case (7)
         write(44,*) "-----------------------------------------------"
         write(44,*) "pip and max_pip" , pip, max_pip,pea(1,1)
         write(44,*) s_time
         lparcnt = 1
         do lcurpar =1,max_pip
            if(lparcnt.gt.pip) exit
            if(.not.pea(lcurpar,1)) cycle
            lparcnt = lparcnt+1
            if(pea(lcurpar,4)) cycle
            write(44,*) "Info for particle", iglobal_id(lcurpar)
            write(44,*) "position new ", des_pos_new(:,lcurpar)
            lcurijk = dg_pijk(lcurpar)
            write(44,*) "Total Neighbours", neighbours(lcurpar,1)
            do lneighindx = 2,neighbours(lcurpar,1)+1
               write(44,*)"Neghibour ",neighbours(lcurpar,lneighindx)
               write(44,*)"Neighbour par position",des_pos_new(:,neighbours(lcurpar,lneighindx))
            end do
            write(44,*) "Total contacts", pn(1,lcurpar)
            do lcontactindx = 2,pn(1,lcurpar)+1
               write(44,*)"contact ", pn(lcontactindx,lcurpar)
               write(44,*)"incontact ", pv(lcontactindx,lcurpar)
               write(44,*)"contact par position",des_pos_new(:,pn(lcontactindx,lcurpar))
            end do
         end do
         write(44,*) "-----------------------------------------------"
      end select
      close(44)
      end subroutine des_dbgmpi


      end module
