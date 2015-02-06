!----------------------------------------------------------------------!
!  Module: MPI_UNPACK_DES                                              !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for unpacking real and ghost particles   !
!     from the MPI recv buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_UNPACK_DES


      PRIVATE
      PUBLIC :: DESMPI_UNPACK_PARCROSS, DESMPI_UNPACK_GHOSTPAR

      interface unpack_dbuf
         module procedure unpack_db0 ! real scalars
         module procedure unpack_db1 ! real arrays
         module procedure unpack_i0  ! integer scalars
         module procedure unpack_i1  ! integer arrays
         module procedure unpack_l0  ! logical scalars
      end interface unpack_dbuf


      CONTAINS


!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_GHOSTPAR                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks ghost particle from the recv buffer.                !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_GHOSTPAR(pface)


! Global Variables:
!---------------------------------------------------------------------//
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! Flag indicating that the ghost particle was updated
      use desmpi, only: iGHOST_UPDATED
! The MPI receive buffer
      use desmpi, only: dRECVBUF
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
      use discretelement, only: DG_PIJK, DG_PIJKPRV
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

      use des_allocate

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
      integer :: lcurpar,lparid,lprvijk,lijk,lparijk,lparcnt,ltot_ind
      integer :: lpacketsize,lbuf,lindx,llocpar,lnewcnt,lpicloc
      logical,dimension(:),allocatable :: lfound
      integer,dimension(:),allocatable :: lnewspot,lnewpic
!......................................................................!

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
! 9) Exiting particle flag
            call unpack_dbuf(lbuf,pea(llocpar,3),pface)
! 10) Temperature
            if(ENERGY_EQ)then
               call unpack_dbuf(lbuf,des_t_s_new(llocpar),pface)
            endif
! 11) Species Composition
            if(ANY_SPECIES_EQ)then
               call unpack_dbuf(lbuf,des_x_s(llocpar,:),pface)
            endif

! 12) User Variables
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

! iAdd new particles and clean up ghost particles if DO_NSEARCH is set.
      if (do_nsearch) then
         call PARTICLE_GROW(pip+lnewcnt)
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
!  9) Exiting particle flag
!           pea(ispot,3) = (drecvbuf(lbuf,pface) > 0.5)
            call unpack_dbuf(lbuf,pea(ispot,3),pface) ! (Need to check the logic)
! 10) Particle temperature.
            if(ENERGY_EQ)then
               call unpack_dbuf(lbuf,des_t_s_new(ispot),pface)
            endif
! 11) Particle species composition
            if(ANY_SPECIES_EQ)then
               call unpack_dbuf(lbuf,des_x_s(ispot,:),pface)
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


!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_PARCROSS                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks real particle from the recv buffer.                 !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_PARCROSS(pface)

! Global Variables:
!---------------------------------------------------------------------//
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! Flag indicating that the ghost particle was updated
      use desmpi, only: iGHOST_UPDATED
! The MPI receive buffer
      use desmpi, only: dRECVBUF
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
      use discretelement, only: DG_PIJK, DG_PIJKPRV
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

! Module Procedures:
!---------------------------------------------------------------------//
      use des_allocate
      use desmpi_wrapper, only: DES_MPI_STOP

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lijk,lcurpar,lparcnt,llocpar,lparid,lparijk,lprvijk
      integer :: lneighindx,lneigh,lcontactindx,lcontactid,lcontact,&
                 lneighid,lneighijk,lneighprvijk
      logical :: lfound
      integer :: lpacketsize,lbuf,ltmpbuf,lcount
      logical :: lcontactfound,lneighfound
      integer :: cc,ii,kk,num_pairs_sent

      integer :: pair_match
      logical :: do_add_pair
!......................................................................!

! loop through particles and locate them and make changes
      lpacketsize = 9*dimn + 3*4 + 15
      lparcnt = drecvbuf(1,pface)

! if mppic make sure enough space available
      call PARTICLE_GROW(pip+lparcnt)

      do lcurpar =1,lparcnt
         lfound = .false.
         lbuf = (lcurpar-1)*lpacketsize + ibufoffset
         call unpack_dbuf(lbuf,lparid,pface)
         call unpack_dbuf(lbuf,lparijk,pface)
         call unpack_dbuf(lbuf,lprvijk,pface)

! PIC particles are always 'new' to the receiving process. Find the
! first available array position and store the global ID. Increment
! the PIP counter to include the new particle.
         IF(MPPIC) THEN
            DO WHILE(PEA(ISPOT,1))
               ISPOT = ISPOT + 1
            ENDDO
            lLOCPAR = iSPOT
            iGLOBAL_ID(lLOCPAR) = lPARID
            PIP = PIP + 1

! A DEM particle should already exist on the current processor as a
! ghost particle. Match the sent particle to the local ghost particle
! by matching the global IDs. Decrement the iGHOST_CNT counter to
! account for the switch from ghost to real particle.
         else
            lFOUND  = LOCATE_PAR(lPARID,lPRVIJK,lLOCPAR)
            IF (.NOT. LFOUND) THEN
               WRITE(*,700) iNEIGHPROC(PFACE), MYPE
               CALL DES_MPI_STOP
            ENDIF
            iGHOST_CNT = iGHOST_CNT - 1
         ENDIF

! convert the local particle from ghost to existing and update its position
         pea(llocpar,1) = .true.
         pea(llocpar,4) = .false.
         dg_pijk(llocpar) = lparijk
         dg_pijkprv(llocpar) = lprvijk
         call unpack_dbuf(lbuf,des_radius(llocpar),pface)
         call unpack_dbuf(lbuf,pijk(llocpar,1:5),pface)
!         pea(llocpar,2:3)     = drecvbuf(lbuf:lbuf+1,pface) ; lbuf=lbuf+2
         pea(llocpar,2:3) = .false.
         if (drecvbuf(lbuf,pface).eq.1) &
            pea(llocpar,2) = .true. ; lbuf = lbuf + 1
         if (drecvbuf(lbuf,pface).eq.1) &
            pea(llocpar,3) = .true. ; lbuf = lbuf + 1
         call unpack_dbuf(lbuf,ro_sol(llocpar),pface)
         call unpack_dbuf(lbuf,pvol(llocpar),pface)
         call unpack_dbuf(lbuf,pmass(llocpar),pface)
         call unpack_dbuf(lbuf,omoi(llocpar),pface)
         call unpack_dbuf(lbuf,des_pos_new(:,llocpar),pface)
         call unpack_dbuf(lbuf,des_vel_new(:,llocpar),pface)

         if(ENERGY_EQ) then
            call unpack_dbuf(lbuf,des_t_s_old(llocpar),pface)
            call unpack_dbuf(lbuf,des_t_s_new(llocpar),pface)
         endif

         if(ANY_SPECIES_EQ)then
            call unpack_dbuf(lbuf,des_x_s(llocpar,:),pface)
         endif

         call unpack_dbuf(lbuf,des_usr_var(:,llocpar),pface)

         call unpack_dbuf(lbuf,omega_new(:,llocpar),pface)
         IF (DO_OLD) THEN
            call unpack_dbuf(lbuf,des_pos_old(:,llocpar),pface)
            call unpack_dbuf(lbuf,des_vel_old(:,llocpar),pface)
            call unpack_dbuf(lbuf,omega_old(:,llocpar),pface)
            call unpack_dbuf(lbuf,des_acc_old(:,llocpar),pface)
            call unpack_dbuf(lbuf,rot_acc_old(:,llocpar),pface)
         ENDIF
         call unpack_dbuf(lbuf,fc(:,llocpar),pface)
         call unpack_dbuf(lbuf,tow(:,llocpar),pface)

         IF(MPPIC) call unpack_dbuf(lbuf,des_stat_wt(llocpar),pface)

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

               print *,"  "
               print *,"  "
               print *," fail on  ", myPE
               print *,"at buffer location",lbuf," pface = ",pface
               print *,"COULD NOT FIND NEIGHBOR ",lneighid," IN IJK ",lneighijk
               call des_mpi_stop
            endif
         endif

! If the neighbor particle is a 'real' particle on this processor, then
! the pair data may already exist. Check before addeding it.
         do_add_pair = .TRUE.
         if(pea(lneigh,1) .and. .not.pea(lneigh,4)) then
            do ii=1,pair_num
               if(PAIRS(1,II) == lneigh) then
                  if(PAIRS(2,II) == llocpar) then
                     do_add_pair = .FALSE.
                     pair_match = II
                     exit
                  endif
               endif
            enddo
         endif

         if(do_add_pair) then
            call add_pair(llocpar,lneigh)
            pair_match = pair_num
         endif

         call unpack_dbuf(lbuf,pv_pair(pair_num),pface)

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

!----------------------------------------------------------------------!
! Function: LOCATE_PAR                                                 !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION LOCATE_PAR(pGLOBALID, pIJK, pLOCALNO)

      use discretelement, only: iGLOBAL_ID
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use discretelement, only: dg_pic

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Global ID of the particle
      INTEGER, INTENT(IN) :: pGlobalID
! IJK of DES grid cell containing the particle
      INTEGER, INTENT(IN) :: pIJK
! Local ID for the particle.
      INTEGER, INTENT(OUT) :: pLocalNO

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: lpicloc, lcurpar

! Initialize the result.
      locate_par = .false.

! Verify that the passied IJK value is within a valid range.
      if(pIJK < dg_ijkstart2 .or. pIJK > dg_ijkend2)  RETURN

! Loop the the particles in DES grid cell pIJK. Return to the calling
! routine if the passed global ID matches the global ID of one of
! the local particles.
      DO lpicloc = 1,dg_pic(pijk)%isize
         lcurpar = dg_pic(pijk)%p(lpicloc)
         IF(iGLOBAL_ID(lcurpar) == pGlobalID) THEN
            plocalno = lcurpar
            locate_par = .true.
            RETURN
         ENDIF
      ENDDO

      RETURN
      end function locate_par

!----------------------------------------------------------------------!
! Function: EXTEN_LOCATE_PAR                                           !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!------------------------------------------------------------------------
      LOGICAL FUNCTION EXTEN_LOCATE_PAR(pGlobalID, pIJK, pLocalNO)

      use discretelement, only: iGLOBAL_ID, dg_pic
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use desgrid, only: dg_Iof_LO, DG_Jof_LO, DG_Kof_LO
      use geometry, only: NO_K

      use desgrid, only: dg_funijk

      implicit none

! Dummy variables:
!---------------------------------------------------------------------//
! The global ID of the particle to be matched locally
      INTEGER, INTENT(IN) :: pGlobalId
! The DES grid cell index expected to contain the particle.
      INTEGER, INTENT(IN) :: pIJK
! The local ID of the matching particle.
      INTEGER, INTENT(OUT) :: pLocalNo

! Local variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: lijk, li, lj, lk, lic, ljc, lkc, lkoffset
      INTEGER :: lpicloc,lcurpar

      exten_locate_par = .false.

      lic = dg_iof_lo(pijk)
      ljc = dg_jof_lo(pijk)
      lkc = dg_kof_lo(pijk)
      lkoffset = merge(0, 1, NO_K)
      DO  lk = lkc-lkoffset,lkc+lkoffset
      DO  lj = ljc-1,ljc+1
      DO  li = lic-1,lic+1
         lijk = dg_funijk(li,lj,lk)
         IF (lijk .lt. dg_ijkstart2 .or. lijk .gt. dg_ijkend2) CYCLE
         DO lpicloc = 1, dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)
            IF (iglobal_id(lcurpar) .eq. pglobalid) THEN
               plocalno = lcurpar
               exten_locate_par = .true.
               RETURN
            END IF
         END DO
      END DO
      END DO
      END DO

      RETURN
      END FUNCTION EXTEN_LOCATE_PAR


!----------------------------------------------------------------------!
! Unpack subroutine for single real variables                          !
!----------------------------------------------------------------------!
      subroutine unpack_db0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata

      idata = drecvbuf(lbuf,pface)
      lbuf = lbuf + 1

      return
      end subroutine unpack_db0

!----------------------------------------------------------------------!
! Unpack subroutine for real arrays                                    !
!----------------------------------------------------------------------!
      subroutine unpack_db1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
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
! Unpack subroutine for single integer variables                       !
!----------------------------------------------------------------------!
      subroutine unpack_i0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata

      idata = drecvbuf(lbuf,pface)
      lbuf = lbuf + 1

      return
      end subroutine unpack_i0

!----------------------------------------------------------------------!
! Unpack subroutine for integer arrays                                 !
!----------------------------------------------------------------------!
      subroutine unpack_i1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(lbuf:lbuf+lsize-1,pface)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_i1

!----------------------------------------------------------------------!
! Unpack subroutine for logical variables                              !
!----------------------------------------------------------------------!
      subroutine unpack_l0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      logical, intent(inout) :: idata

      idata = merge(.true.,.false.,0.5<drecvbuf(lbuf,pface))
      lbuf = lbuf + 1

      return
      end subroutine unpack_l0


      end module mpi_unpack_des
