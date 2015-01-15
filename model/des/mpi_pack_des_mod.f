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

      contains


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

! Do not send particle data for a ghost particle whose owner has not yet
! updated the particle's data on this processor.
            if(pea(lcurpar,4) .and. .not.ighost_updated(lcurpar) ) cycle

! 1) Global ID
            dsendbuf(lbuf,pface) = iglobal_id(lcurpar)
            lbuf = lbuf +1
! 2) DES grid IJK
            dsendbuf(lbuf,pface) = dg_ijkconv(lijk,pface,ineighproc(pface))
            lbuf = lbuf +1
! 3) DES grid IJK - previous
            dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface))
            lbuf = lbuf +1
! 4) Radius
            dsendbuf(lbuf,pface) = des_radius(lcurpar)
            lbuf = lbuf + 1
! 5) Phase index
            dsendbuf(lbuf,pface) = pijk(lcurpar,5)
            lbuf = lbuf + 1
! 6) Position
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
            lbuf = lbuf + dimn
! 7) Translational Velocity
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_new(1:dimn,lcurpar)
            lbuf = lbuf + dimn
! 8) Rotational Velocity
            dsendbuf(lbuf:lbuf+dimn-1,pface) = omega_new(1:dimn,lcurpar)
            lbuf = lbuf + dimn
! 9) Exiting particle flag
            dsendbuf(lbuf,pface) = merge(1,0,pea(lcurpar,3))
            lbuf = lbuf + 1
! 10) Temperature
            if(ENERGY_EQ)then
               dsendbuf(lbuf,pface) = des_t_s_new(lcurpar)
              lbuf = lbuf +1
            endif
! 11) Species Composition
            if(ANY_SPECIES_EQ)then
               dsendbuf(lbuf:lbuf+dimension_n_s-1,pface) = &
                  des_x_s(lcurpar,1:dimension_n_s)
               lbuf = lbuf+dimension_n_s
            endif

! 12) User Variable
            dsendbuf(lbuf:lbuf+3-1,pface) = des_usr_var(1:3,lcurpar)
            lbuf = lbuf+3

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

! Location in the buffer where the number of pair data is specified.
      integer :: num_pairs_send_buf_loc


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
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_new(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_new(1:dimn,lcurpar)
            lbuf = lbuf+dimn

            if(ENERGY_EQ) then
               dsendbuf(lbuf,pface) = des_t_s_old(lcurpar)
               lbuf = lbuf+1
               dsendbuf(lbuf,pface) = des_t_s_new(lcurpar)
               lbuf = lbuf+1
            endif

            if(ANY_SPECIES_EQ)then
               dsendbuf(lbuf:lbuf+dimension_n_s-1,pface) = &
                  des_x_s(lcurpar,1:dimension_n_s)
               lbuf = lbuf + dimension_n_s
            endif

            dsendbuf(lbuf:lbuf+3-1,pface) = des_usr_var(1:3,lcurpar)
            lbuf = lbuf+3

            dsendbuf(lbuf:lbuf+3-1,pface) = omega_new(1:3,lcurpar)
            lbuf = lbuf+3
            IF (DO_OLD) THEN
               dsendbuf(lbuf:lbuf+dimn-1,pface) = des_pos_old(1:dimn,lcurpar)+dcycl_offset(pface,1:dimn)
               lbuf = lbuf+dimn
               dsendbuf(lbuf:lbuf+dimn-1,pface) = des_vel_old(1:dimn,lcurpar)
               lbuf = lbuf+dimn
               dsendbuf(lbuf:lbuf+3-1,pface) = omega_old(1:3,lcurpar)
               lbuf = lbuf+3
               dsendbuf(lbuf:lbuf+dimn-1,pface) = des_acc_old(1:dimn,lcurpar)
               lbuf = lbuf+dimn
               dsendbuf(lbuf:lbuf+3-1,pface) = rot_acc_old(1:3,lcurpar)
               lbuf = lbuf+3
            ENDIF
            dsendbuf(lbuf:lbuf+dimn-1,pface) = fc(:,lcurpar)
            lbuf = lbuf+dimn
            dsendbuf(lbuf:lbuf+3-1,pface) = tow(1:3,lcurpar)
            lbuf = lbuf+3

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


! Calculate the location in buffer where the number of pair data is
! stored and skip specifying the entry. After all the pair data is
! packed, then this value is set.
      lbuf = lparcnt*lpacketsize + ibufoffset
      num_pairs_send_buf_loc = lbuf
      lbuf = lbuf+1

      num_pairs_to_send = 0
      do cc = 1, pair_num
         lcurpar = PAIRS(1,CC)
! Only packup pairing data for particles being transfered.
         if (.not. going_to_send(lcurpar)) cycle
! Do not send pairing data if the pair no longer exists or if the
! particle is exiting as it may be locatable during unpacking.
         lneigh = PAIRS(2,CC)
         if(.not.PEA(lneigh,1)) cycle
         if(PEA(lneigh,3)) cycle
! Global ID of particle bing packed.
         dsendbuf(lbuf,pface) = iglobal_id(lcurpar)
         lbuf = lbuf+1
! DES grid IJK of cell receiving the particle.
         dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijkprv(lcurpar),pface,ineighproc(pface))
         lbuf = lbuf+1
! Global ID of particle pair. This particle may or may not live on the
! the current or destination processor.
         dsendbuf(lbuf,pface) = iglobal_id(lneigh)
         lbuf = lbuf+1
! DES grid IJK of cell containing the particle-pair.
         dsendbuf(lbuf,pface) = dg_ijkconv(dg_pijkprv(lneigh),pface,ineighproc(pface))
         lbuf = lbuf+1
! Convernt the logical flag to an integer.
         dsendbuf(lbuf,pface) = merge(1,0,pv_PAIR(CC))
         lbuf = lbuf+1
! Pack the normal and tangential collision histories.
         do ii=1,DIMN
            dsendbuf(lbuf,pface) = PFN_PAIR(II,CC)
            lbuf = lbuf+1
            dsendbuf(lbuf,pface) = PFT_PAIR(II,CC)
            lbuf = lbuf+1
         enddo
! Increment the number of pairs being sent.
         num_pairs_to_send = num_pairs_to_send + 1
      enddo

! Store the number of pair datasets being sent. This information is
! stored before the pairing data so the receiving process knows the
! amount of data to 'unpack.'
      dsendbuf(num_pairs_send_buf_loc,pface) = num_pairs_to_send

      lpairsize = 6 + 2*DIMN

      dsendbuf(1,pface) = lparcnt
      isendcnt(pface) = lparcnt*lpacketsize+num_pairs_to_send*lpairsize+ibufoffset + 3

      deallocate(going_to_send)

      return
      end subroutine desmpi_pack_parcross

      end module mpi_pack_des
