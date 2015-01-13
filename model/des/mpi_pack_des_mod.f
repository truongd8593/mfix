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
      module mpi_pack_des

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

      interface pack_dbuf
         module procedure pack_db0, pack_db1
      end interface pack_dbuf

      contains


!----------------------------------------------------------------------!
!Pack subroutine for single variables                                                                      !
!----------------------------------------------------------------------!
      subroutine pack_db0(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata

      dsendbuf(lbuf,pface) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_db0

!----------------------------------------------------------------------!
!Pack subroutine for arrays                                                                      !
!----------------------------------------------------------------------!
      subroutine pack_db1(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(lbuf:lbuf+lsize-1,pface) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_db1


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
      integer :: lpacketsize,lbuf
!-----------------------------------------------
      lpacketsize = 2*dimn + 3+ 5
      lpar_cnt = 0
      ltot_ind = isendindices(1,pface)
      do lindx = 2,ltot_ind+1
         lijk = isendindices(lindx,pface)
         do lpicloc =1,dg_pic(lijk)%isize
            lbuf = lpar_cnt*lpacketsize+ibufoffset
            lcurpar = dg_pic(lijk)%p(lpicloc)
            if(pea(lcurpar,4) .and. .not.ighost_updated(lcurpar) ) cycle

! 1) Global ID
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(dg_ijkconv(lbuf,dg_pijkprv(lcurpar),pface,ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 6) Position
            call pack_dbuf(lbuf,des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn),pface)
! 7) Translational Velocity
            call pack_dbuf(lbuf,des_vel_new(1:dimn,lcurpar),pface)
! 8) Rotational Velocity
            call pack_dbuf(lbuf,omega_new(1:dimn,lcurpar),pface)

! 9) Temperature
            if(ENERGY_EQ)then
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
            endif
! 10) Species Composition
            if(ANY_SPECIES_EQ)then
               call pack_dbuf(lbuf,des_x_s(lcurpar,1:dimension_n_s),pface)
            endif

! 11) User Variable
            call pack_dbuf(lbuf,des_usr_var(1:3,lcurpar),pface)

            lpar_cnt = lpar_cnt + 1
         end do
      end do
      dsendbuf(1,pface)=lpar_cnt
      isendcnt(pface) = lpar_cnt*lpacketsize+ibufoffset

      end subroutine desmpi_pack_ghostpar


!------------------------------------------------------------------------
! Subroutine       : desmpi_pack_parcross
! Purpose          : packs the particle crossing the boundary
! Parameter        : pface - value from 1 to 6 represents faces
!------------------------------------------------------------------------
      subroutine desmpi_pack_parcross(pface)

      use functions

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      integer, intent(in) :: pface
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: li, lj, lk
      integer :: ltot_ind,lindx,ijk,cc,ii,ll,kk
      integer :: lneighindx,lcontactindx,lneigh,lcontact,lijk,&
                 lpicloc,lparcnt,lcurpar
      integer :: lpacketsize,lbuf,ltmpbuf,num_pairs_to_send,lpairsize

      logical, allocatable, dimension(:) :: going_to_send

!-----------------------------------------------

! pack the particle crossing the boundary
      lpacketsize = 9*dimn + 3*4 + 15
      ltot_ind = irecvindices(1,pface)
      lparcnt = 0

      allocate(going_to_send(max_pip))
      going_to_send(:) = .false.

      do lindx = 2,ltot_ind+1
         lijk = irecvindices(lindx,pface)
         do lpicloc = 1,dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)

            if (pea(lcurpar,4)) cycle ! if ghost particle then cycle

            going_to_send(lcurpar) = .true.

       IF(iGLOBAL_ID(lcurpar) == 1641 .OR. iGLOBAL_ID(lcurpar) == 777) THEN

          write(*,*) 'I am packing particle: ',iGLOBAL_ID(lcurpar),'on',myPE
       ENDIF

            lbuf = lparcnt*lpacketsize + ibufoffset
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,ineighproc(pface)),pface)
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface)),pface)
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
            li = pijk(lcurpar,1) + icycoffset(pface,1)
            lj = pijk(lcurpar,2) + icycoffset(pface,2)
            lk = pijk(lcurpar,3) + icycoffset(pface,3)
            call pack_dbuf(lbuf,li,pface)
            call pack_dbuf(lbuf,lj,pface)
            call pack_dbuf(lbuf,lk,pface)
            call pack_dbuf(lbuf,funijk_proc(li,lj,lk,ineighproc(pface)),pface)
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
!            dsendbuf(lbuf:lbuf+1,pface) = pea(lcurpar,2:3);lbuf=lbuf+2
            dsendbuf(lbuf:lbuf+1,pface) = 0
            if (pea(lcurpar,2)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            if (pea(lcurpar,3)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            call pack_dbuf(lbuf,ro_sol(lcurpar),pface)
            call pack_dbuf(lbuf,pvol(lcurpar),pface)
            call pack_dbuf(lbuf,pmass(lcurpar),pface)
            call pack_dbuf(lbuf,omoi(lcurpar),pface)
            call pack_dbuf(lbuf,des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn),pface)
            call pack_dbuf(lbuf,des_vel_new(1:dimn,lcurpar),pface)

            if(ENERGY_EQ) then
               call pack_dbuf(lbuf,des_t_s_old(lcurpar),pface)
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
            endif

            if(ANY_SPECIES_EQ)then
               call pack_dbuf(lbuf,des_x_s(lcurpar,1:dimension_n_s),pface)
            endif

            call pack_dbuf(lbuf, des_usr_var(1:3,lcurpar),pface)

            call pack_dbuf(lbuf,omega_new(1:3,lcurpar),pface)
            IF (DO_OLD) THEN
               call pack_dbuf(lbuf,des_pos_old(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn),pface)
               call pack_dbuf(lbuf,des_vel_old(1:dimn,lcurpar),pface)
               call pack_dbuf(lbuf,omega_old(1:3,lcurpar),pface)
               call pack_dbuf(lbuf,des_acc_old(1:dimn,lcurpar),pface)
               call pack_dbuf(lbuf,rot_acc_old(1:3,lcurpar),pface)
            ENDIF
            call pack_dbuf(lbuf,fc(:,lcurpar),pface)
            call pack_dbuf(lbuf,tow(1:3,lcurpar),pface)

! In case of mppic remove the particles else
! Convert the particle as ghost and set the forces zero
            if (mppic) then
               pea(lcurpar,1:4) = .false.
            else
               pea(lcurpar,4) = .true.
               ighost_cnt = ighost_cnt + 1
            end if
            fc(:,lcurpar) = 0.

            lparcnt = lparcnt + 1
         end do
      end do

      lbuf = lparcnt*lpacketsize + ibufoffset

      num_pairs_to_send = 0
      do cc = 1, pair_num
         LL = PAIRS(1,CC)
         if (going_to_send(LL)) then
            num_pairs_to_send = num_pairs_to_send + 1
         endif
      enddo

      call pack_dbuf(lbuf,num_pairs_to_send,pface)

      do cc = 1, pair_num
         lcurpar = PAIRS(1,CC)
         if (.not. going_to_send(lcurpar)) cycle

         call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
!         dsendbuf(lbuf,pface) = dg_ijkconv(lijk,pface,ineighproc(pface))
!         lbuf = lbuf+1
         call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface)),pface)

         lneigh = PAIRS(2,CC)

         call pack_dbuf(lbuf,iglobal_id(lneigh),pface)
!         dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijk(lneigh),pface,ineighproc(pface))
!         lbuf = lbuf+1
         call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lneigh),pface,ineighproc(pface)),pface)

         call pack_dbuf(lbuf,merge(1,0,pv_PAIR(CC)),pface)
         do ii=1,DIMN
            call pack_dbuf(lbuf,PFN_PAIR(II,CC),pface)
            call pack_dbuf(lbuf,PFT_PAIR(II,CC),pface)
         enddo
      enddo

      deallocate(going_to_send)

      lpairsize = 6 + 2*DIMN

      dsendbuf(1,pface) = lparcnt
      isendcnt(pface) = lparcnt*lpacketsize+num_pairs_to_send*lpairsize+ibufoffset + 3

! following unused variables are not sent across the processor
! well_depth
! is_linked
! links
! part_grid

      end subroutine desmpi_pack_parcross

      end module mpi_pack_des
