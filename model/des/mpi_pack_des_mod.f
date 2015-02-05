!----------------------------------------------------------------------!
!  Module: MPI_PACK_DES                                                !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for packing real and ghost particles     !
!     into the MPI send buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_PACK_DES


      PRIVATE
      PUBLIC :: DESMPI_PACK_PARCROSS, DESMPI_PACK_GHOSTPAR

      interface pack_dbuf
        module procedure pack_db0
        module procedure pack_db1
        module procedure pack_i0
        module procedure pack_i1
      end interface pack_dbuf

      CONTAINS

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_GHOSTPAR                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs ghost particle in the send buffer.                    !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_GHOSTPAR(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! Flag indicating that the ghost particle was updated
      use desmpi, only: iGHOST_UPDATED
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SPECIES_EQ
! Runtime flag for MPPIC solids
      use mfix_pic, only: MPPIC
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! Dimenions of DES grid
      use desgrid, only: DG_IJKSIZE2
! DES grid cell containing each particle: current/previous
      use desgrid, only: DG_PIJK, DG_PIJKPRV
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW, OMEGA_OLD
! Particle species composition
      use des_rxns, only: DES_X_s
! Particle tempertures. current/previous
      use des_thermo, only: DES_T_s_NEW, DES_T_s_OLD
! Particle radius, volume
      use discretelement, only: DES_RADIUS, PVOL
! Flags indicate the state of the particle
      use discretelement, only: PEA
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Flag to send/recv old (previous) values
      use discretelement, only: DO_OLD
! Flag to conduct a new neighbor search.
      use discretelement, only: DO_NSEARCH
! Number of particles on the process (max particle array size)
      use discretelement, only: PIP, MAX_PIP
! Number of ghost particles on the current process
      use discretelement, only: iGHOST_CNT
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR
! Function to convert DES grid IJK to new proc value.
      use desgrid, only: dg_ijkconv
! Size of the send buffer
      use desmpi, only: isendcnt
! Offset for particles with cycle BCs (otherwise zero)
      use desmpi, only: dcycl_offset
! Map of particles to DES grid
      use desgrid, only: DG_PIC
! Cell number of ghost particles
      use desmpi, only: iSENDINDICES

! Global Constants:
!---------------------------------------------------------------------//
      use constant, only: PI
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN


      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lijk,lindx,ltot_ind,lpicloc,lpar_cnt,lcurpar
      integer :: lpacketsize,lbuf
!......................................................................!

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
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 6) Position
            call pack_dbuf(lbuf,des_pos_new(:,lcurpar)+                &
               dcycl_offset(pface,:),pface)
! 7) Translational Velocity
            call pack_dbuf(lbuf,des_vel_new(:,lcurpar),pface)
! 8) Rotational Velocity
            call pack_dbuf(lbuf,omega_new(:,lcurpar),pface)
! 9) Exiting particle flag
            call pack_dbuf(lbuf,merge(1,0,pea(lcurpar,3)),pface)
! 10) Temperature
            if(ENERGY_EQ)then
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
            endif
! 11) Species Composition
            if(ANY_SPECIES_EQ)then
               call pack_dbuf(lbuf,des_x_s(lcurpar,:),pface)
            endif
! 12) User Variable
            call pack_dbuf(lbuf,des_usr_var(1:3,lcurpar),pface)

            lpar_cnt = lpar_cnt + 1
         end do
      end do
      dsendbuf(1,pface)=lpar_cnt
      isendcnt(pface) = lpar_cnt*lpacketsize+ibufoffset

      end subroutine desmpi_pack_ghostpar


!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_PARCROSS                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs real particle in the send buffer.                     !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_PARCROSS(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! Flag indicating that the ghost particle was updated
      use desmpi, only: iGHOST_UPDATED
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SPECIES_EQ
! Runtime flag for MPPIC solids
      use mfix_pic, only: MPPIC
! Dimenions of DES grid
      use desgrid, only: DG_IJKSIZE2
! DES grid cell containing each particle: current/previous
      use desgrid, only: DG_PIJK, DG_PIJKPRV
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! The statistical weight of each particle.
      use mfix_pic, only: DES_STAT_WT
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW, OMEGA_OLD
! Particle radius, volume, density, mass
      use discretelement, only: DES_RADIUS, PVOL, RO_SOL, PMASS
! Previous value for particle acceleration (tangential/rotational)
      use discretelement, only: DES_ACC_OLD, ROT_ACC_OLD
! Particle species composition
      use des_rxns, only: DES_X_s
! Particle tempertures. current/previous
      use des_thermo, only: DES_T_s_NEW, DES_T_s_OLD
! Force arrays acting on the particle
      use discretelement, only: FC, TOW
! One of the moment of inertia
      use discretelement, only: OMOI
! Flags indicate the state of the particle
      use discretelement, only: PEA
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Flag to send/recv old (previous) values
      use discretelement, only: DO_OLD
! Flag to conduct a new neighbor search.
      use discretelement, only: DO_NSEARCH
! Number of particles on the process (max particle array size)
      use discretelement, only: PIP, MAX_PIP
! Number of ghost particles on the current process
      use discretelement, only: iGHOST_CNT
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR
! Particle pair (neighborhood) arrays:
      use discretelement, only: PAIR_NUM, PAIRS
! Pair collision history information
      use discretelement, only: PV_PAIR, PFN_PAIR, PFT_PAIR
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN
! The ID of the current process
      use compar, only: myPE



      use desgrid, only: dg_ijkconv, icycoffset
      use desmpi, only: dcycl_offset, isendcnt
      use desgrid, only: DG_PIC
      use desmpi, only: iSENDINDICES
      use desmpi, only: irecvindices

      use functions

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: li, lj, lk
      integer :: ltot_ind,lindx,cc,ii
      integer :: lneigh,lijk,&
                 lpicloc,lparcnt,lcurpar
      integer :: lpacketsize,lbuf,num_pairs_to_send,lpairsize

      logical, allocatable, dimension(:) :: going_to_send

! Location in the buffer where the number of pair data is specified.
      integer :: num_pairs_send_buf_loc
!......................................................................!

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
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)

            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
            li = pijk(lcurpar,1) + icycoffset(pface,1)
            lj = pijk(lcurpar,2) + icycoffset(pface,2)
            lk = pijk(lcurpar,3) + icycoffset(pface,3)
            call pack_dbuf(lbuf,li,pface)
            call pack_dbuf(lbuf,lj,pface)
            call pack_dbuf(lbuf,lk,pface)
            call pack_dbuf(lbuf,funijk_proc(li,lj,lk,                  &
               ineighproc(pface)),pface)
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
!            dsendbuf(lbuf:lbuf+1,pface) = pea(lcurpar,2:3);lbuf=lbuf+2
            dsendbuf(lbuf:lbuf+1,pface) = 0
            if (pea(lcurpar,2)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            if (pea(lcurpar,3)) dsendbuf(lbuf,pface) = 1 ; lbuf = lbuf+1
            call pack_dbuf(lbuf,ro_sol(lcurpar),pface)
            call pack_dbuf(lbuf,pvol(lcurpar),pface)
            call pack_dbuf(lbuf,pmass(lcurpar),pface)
            call pack_dbuf(lbuf,omoi(lcurpar),pface)
            call pack_dbuf(lbuf,des_pos_new(:,lcurpar) +          &
               dcycl_offset(pface,:),pface)
            call pack_dbuf(lbuf,des_vel_new(:,lcurpar),pface)

            if(ENERGY_EQ) then
               call pack_dbuf(lbuf,des_t_s_old(lcurpar),pface)
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
            endif

            if(ANY_SPECIES_EQ)then
               call pack_dbuf(lbuf,des_x_s(lcurpar,:),pface)
            endif

            call pack_dbuf(lbuf, des_usr_var(1:3,lcurpar),pface)

            call pack_dbuf(lbuf,omega_new(1:3,lcurpar),pface)
            IF (DO_OLD) THEN
               call pack_dbuf(lbuf,des_pos_old(:,lcurpar) +            &
                  dcycl_offset(pface,:),pface)
               call pack_dbuf(lbuf,des_vel_old(:,lcurpar),pface)
               call pack_dbuf(lbuf,omega_old(:,lcurpar),pface)
               call pack_dbuf(lbuf,des_acc_old(:,lcurpar),pface)
               call pack_dbuf(lbuf,rot_acc_old(:,lcurpar),pface)
            ENDIF
            call pack_dbuf(lbuf,fc(:,lcurpar),pface)
            call pack_dbuf(lbuf,tow(:,lcurpar),pface)

! PIC particles are removed and the number of particles on the processor
! is decremented.
            IF (MPPIC) THEN
               call pack_dbuf(lbuf,des_stat_wt(lcurpar),pface)
               pea(lcurpar,1:4) = .false.
               pip = pip - 1

! DEM particles are converted to ghost particles. This action does not
! change the number of particles.
            ELSE
               pea(lcurpar,4) = .true.
               ighost_cnt = ighost_cnt + 1
            END IF

! Clear out the force array.
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
         call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! DES grid IJK of cell receiving the particle.
         call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,     &
            ineighproc(pface)),pface)
! Global ID of particle pair. This particle may or may not live on the
! the current or destination processor.
         call pack_dbuf(lbuf,iglobal_id(lneigh),pface)
! DES grid IJK of cell containing the particle-pair.
         call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lneigh),pface,      &
            ineighproc(pface)),pface)
! Convernt the logical flag to an integer.
         call pack_dbuf(lbuf,merge(1,0,pv_PAIR(CC)),pface)
! Pack the normal and tangential collision histories.

         do ii=1,DIMN
            call pack_dbuf(lbuf,PFN_PAIR(II,CC),pface)
            call pack_dbuf(lbuf,PFT_PAIR(II,CC),pface)
         enddo
! Increment the number of pairs being sent.
         num_pairs_to_send = num_pairs_to_send + 1
      enddo

! Store the number of pair datasets being sent. This information is
! stored before the pairing data so the receiving process knows the
! amount of data to 'unpack.'
      lbuf = num_pairs_send_buf_loc
      call pack_dbuf(lbuf,num_pairs_to_send,pface)

      lpairsize = 6 + 2*DIMN

      dsendbuf(1,pface) = lparcnt
      isendcnt(pface) = lparcnt*lpacketsize +                          &
         num_pairs_to_send*lpairsize+ibufoffset + 3

      deallocate(going_to_send)


      END SUBROUTINE DESMPI_PACK_PARCROSS

!----------------------------------------------------------------------!
! PACK SUBROUTINE FOR SINGLE REAL VARIABLES                            !
!----------------------------------------------------------------------!
      subroutine pack_db0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata

      dsendbuf(lbuf,pface) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_db0

!----------------------------------------------------------------------!
! Pack subroutine for real arrays                                      !
!----------------------------------------------------------------------!
      subroutine pack_db1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(lbuf:lbuf+lsize-1,pface) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_db1


!----------------------------------------------------------------------!
! Pack subroutine for single integer variables                         !
!----------------------------------------------------------------------!
      subroutine pack_i0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata

      dsendbuf(lbuf,pface) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_i0

!----------------------------------------------------------------------!
! Pack subroutine for integer arrays                                   !
!----------------------------------------------------------------------!
      subroutine pack_i1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(lbuf:lbuf+lsize-1,pface) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_i1

      end module mpi_pack_des
