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
!    desmpi_sendrecv_wait, desmpi_gatherv, desmpi_scatterv
!    desmpi_check_sendrecvbuf,
!    desmpi_pack_ghostpar, desmpi_unpack_ghostpar, desmpi_cleanup,
!    desmpi_pack_parcross, desmpi_unpack_parcross,
!    des_addnodevalues, des_addnodevalues2,
!    des_gather_d,l,i

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
      use des_thermo
      use run, only: ENERGY_EQ,ANY_SPECIES_EQ
      use param, only: DIMENSION_N_s
      use des_rxns



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


      integer,dimension(:),allocatable  :: itempglobal_id

      contains

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
         lpacketsize = 9*dimn + ltordimn*4 + 13
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
         end do
         write(44,*) "-----------------------------------------------"
      end select
      close(44)
      end subroutine des_dbgmpi


      end module
