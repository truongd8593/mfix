!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      DOUBLE PRECISION :: D(DIMN), DIST, &
                          NEIGHBOR_SEARCH_DIST
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
!-----------------------------------------------


      IF(MPPIC) THEN
         IF(MPPIC_SOLID_STRESS_SNIDER) THEN
            CALL CFNEWVALUES_MPPIC_SNIDER
         ELSE
! call the coloring function like approach
            CALL CFNEWVALUES_MPPIC
         ENDIF
         RETURN
      ENDIF


!$omp parallel do if(max_pip .ge. 10000) default(shared)        &
!$omp private(l,d,dist,neighbor_search_dist)                    &
!$omp reduction(.or.:do_nsearch) schedule (auto)

      DO L = 1, MAX_PIP
! only process particles that exist
         IF(.NOT.PEA(L,1)) CYCLE
! skip ghost particles
         IF(PEA(L,4)) CYCLE

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.PEA(L,2))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
            IF(USE_COHESION .AND. VAN_DER_WAALS) &
               FC(:,L) = FC(:,L) + Fcohesive(L,:)/PMASS(L)
         ELSE
            FC(:,L) = ZERO
            TOW(:,L) = ZERO
         ENDIF


! Advance particle position, velocity
        IF (INTG_EULER) THEN
! first-order method
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + FC(:,L)*DTSOLID
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
               DES_VEL_NEW(L,:)*DTSOLID
! following is equivalent to x=xold + vold*dt + 1/2acc*dt^2
!         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
!             (DES_VEL_NEW(L,:)+DES_VEL_OLD(L,:))*DTSOLID
            OMEGA_NEW(L,:)   = OMEGA_OLD(L,:) + TOW(:,L)*OMOI(L)*DTSOLID

            ELSEIF (INTG_ADAMS_BASHFORTH) THEN
! T.Li:  second-order Adams-Bashforth scheme
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
               ( 3.d0*DES_VEL_OLD(L,:)-DES_VEL_OOLD(L,:) )*DTSOLID
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + 0.5d0*&
               ( 3.d0*FC(:,L)-DES_ACC_OLD(L,:) )*DTSOLID
            OMEGA_NEW(L,:)   =  OMEGA_OLD(L,:) + 0.5d0*&
               ( 3.d0*TOW(:,L)*OMOI(L)-ROT_ACC_OLD(L,:) )*DTSOLID
            DES_ACC_OLD(L,:) = FC(:,L)
            ROT_ACC_OLD(L,:) = TOW(:,L)*OMOI(L)
         ENDIF


! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
            DIST = SQRT(DES_DOTPRDCT(D,D))
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
         ENDIF


! Check if the particle has moved a distance greater than or equal to
! its radius during one solids time step. if so, call stop
         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         DIST = SQRT(DES_DOTPRDCT(D,D))
         IF(DIST.GE.DES_RADIUS(L)) THEN
            WRITE(*,1002) L, DIST, DES_RADIUS(L)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'old particle pos = ', DES_POS_OLD(L,:)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'new particle pos = ', DES_POS_NEW(L,:)
            WRITE(*,'(5X,A,3(ES17.9))')&
               'new particle vel = ', DES_VEL_NEW(L,:)
            WRITE(*,1003)
            STOP
         ENDIF

! Reset total contact force and torque
         FC(:,L) = ZERO
         TOW(:,L) = ZERO

      ENDDO
!$omp end parallel do


 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

      RETURN
      END SUBROUTINE CFNEWVALUES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES_MPPIC_SNIDER                               C
!
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES_MPPIC_SNIDER

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mppic_wallbc
      USE fldvar
      USE cutcell
      USE mfix_pic
      USE randomno
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM
      INTEGER I, J, K, IJK, IJK_OLD, IJK2, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, DP_BAR, COEFF_EN, MEANVEL(DIMN), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), MEAN_FREE_PATH, PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, MEANUS(DIMN, MMAX)
      DOUBLE PRECISION :: DPS_DXE, DPS_DXW, DPS_DYN, DPS_DYS, DPS_DZT, DPS_DZB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, epg_min2, velf_part(dimn)

      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT, count_bc

      DOUBLE PRECISION MEANUS_e(DIMN, MMAX), MEANUS_w(DIMN, MMAX),MEANUS_n(DIMN, MMAX),MEANUS_s(DIMN, MMAX),MEANUS_t(DIMN, MMAX), MEANUS_b(DIMN, MMAX)
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5), epg_min_loc(1)


      double precision  sig_u, mean_u,ymid
      double precision, allocatable, dimension(:,:) ::  rand_vel
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      PC = 1
      FOCUS_PARTICLE = -1
      DTPIC_CFL = LARGE_NUMBER

      if(dimn.eq.2) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(dimn.eq.3) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)

      allocate(rand_vel(MAX_PIP, DIMN))
      do idim = 1, dimn
         mean_u = zero
         sig_u = 1.d0
         CALL NOR_RNO(RAND_VEL(1:MAX_PIP, IDIM), MEAN_U, SIG_U)
      enddo

      DO L = 1, MAX_PIP
         DELETE_PART = .false.
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle
         pc = pc+1
         if(pea(l,4)) cycle

         DES_LOC_DEBUG = .FALSE.

         !if(mppic) then
         !   IJK = PIJK(L, 4)

         !   COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
         !   IF(COUNT_BC.GE.1) CALL MPPIC_ADD_FRIC_FORCE(L)
         !ENDIF
! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.PEA(L,2))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
         ENDIF

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
         !By comparing MFIX equations and equations in those papers,
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO

         M = PIJK(L,5)
         IJK = PIJK(L,4)

         IJK_OLD = IJK

         PIJK_OLD(:) = PIJK(L,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         COEFF_EN =  MPPIC_COEFF_EN1
         UPRIMETAU(:) = ZERO

         VEL_ORIG(:) = DES_VEL_NEW(L,:)

         DES_VEL_NEW(L,:) = (DES_VEL_OLD(L,:) + &
         & FC(:,L)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         do idim = 1, dimn
            SIG_U = 0.05D0
            rand_vel(L, idim)  = sig_u*DES_VEL_NEW(L, IDIM)*rand_vel(L, idim)
         enddo


         IF(L.EQ.FOCUS_PARTICLE) THEN

            WRITE(*,'(A20,2x,3(2x,i4))') 'CELL ID = ', PIJK(L,1:3)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'EPS = ', 1.d0 - EP_g(PIJK(L,4))
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', DES_VEL_NEW(L,:)

            WRITE(*,'(A20,2x,3(2x,g17.8))') 'FC = ', FC(:,L)
         ENDIF

         !MEANVEL(1) = DES_U_S(IJK_OLD,M)
         !MEANVEL(2) = DES_V_S(IJK_OLD,M)
         !IF(DIMN.EQ.3) MEANVEL(3) = DES_W_S(IJK_OLD,M)

         PS_FORCE(:) = PS_GRAD(L, :)
         DELUP(:) = -( DTSOLID*PS_FORCE(:))/((1.d0+DP_BAR*DTSOLID))
         DELUP(:) = DELUP(:)/ROP_S(IJK_OLD,M)

         MEANVEL(:) = AVGSOLVEL_P(L,:)

         DO IDIM = 1, DIMN
            IF(PS_FORCE(IDIM).LE.ZERO) THEN
               UPRIMETAU(IDIM) = MIN(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(L,IDIM)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MAX(UPRIMETAU(IDIM), ZERO)
            ELSE
               UPRIMETAU(IDIM) = MAX(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(L,IDIM)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MIN(UPRIMETAU(IDIM), ZERO)
            END IF

         ENDDO

         DES_VEL_NEW(L,:) = DES_VEL_NEW(L,:) + UPRIMETAU(:)
         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
         DES_VEL_NEW(L,:)*DTSOLID


         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(L,1:DIMN), DES_VEL_NEW(L, 1:DIMN)))

         RAD_EFF = DES_RADIUS(L)
               !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
         MEAN_FREE_PATH = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2

         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then
         !   DES_VEL_NEW(L,:) = (DES_VEL_NEW(L,:)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

          D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)

          DIST = SQRT(DES_DOTPRDCT(D,D))


         CALL MPPIC_FIND_NEW_CELL(L)

         IJK = PIJK(L,4)

         INSIDE_DOMAIN = .true.

         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD) then

            IF(CUT_CELL_AT(IJK)) then
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + rand_vel(L, :)*dtsolid
               DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
            ELSE
               !IF(IJK.NE.IJK_OLD) THEN
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + rand_vel(L, :)*dtsolid
               DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
               !ENDIF
            ENDIF
         ENDIF

         PIJK(L,:) = PIJK_OLD(:)
         DIST = SQRT(DES_DOTPRDCT(D,D))

          IF(DIST.GT.MEAN_FREE_PATH) THEN
          !WRITE(*,*) 'UPRIME OLD = ', UPRIMETAU(:)
          !WRITE(*,*) 'DIST GT MEAN FREE PATH= ', DIST, MEAN_FREE_PATH

             DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
             DES_VEL_NEW(L,:)*DTSOLID*MEAN_FREE_PATH/DIST

             !D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)

                                !DIST = SQRT(DES_DOTPRDCT(D,D))
             !WRITE(*,*) 'new moved distance  = ', dist,1. -  ep_g(ijk)
          ENDIF


         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         D_GRIDUNITS(1) = ABS(D(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(D(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DIMN.EQ.3) D_GRIDUNITS(3) = ABS(D(3))/DZ(PIJK(L,3))

         DIST = SQRT(DES_DOTPRDCT(D,D))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit

         DO IDIM = 1, DIMN
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN

                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  DELETE_PART = .true.

               ENDIF
               !CALL mfix_exit(myPE)
            END IF

         END DO
         IF(.not.DELETE_PART) DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         FC(:,L) = ZERO

         IF(DELETE_PART) THEN
            PEA(L,1) = .false.
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      IF(MPPIC) THEN
         CALL global_all_max(DTPIC_CFL)
         PIP = PIP - PIP_DEL_COUNT

         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
            IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'

            WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'
            !DO IPROC = 0, NUMPES-1
            !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            !ENDDO

         ENDIF

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN
            !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN

            !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
         ELSE
            !WRITE(*,'(A40,2x,g17.8)') 'DT
            !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         END IF


      ENDIF
 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

 2003 FORMAT(/10x, &
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)


      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC_SNIDER


      SUBROUTINE CFNEWVALUES_MPPIC

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mppic_wallbc
      USE mfix_pic
      USE randomno
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM
      INTEGER I, J, K, IJK, IJK_OLD, IJK2, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, DP_BAR, COEFF_EN, MEANVEL(DIMN), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), MEAN_FREE_PATH, PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, MEANUS(DIMN, MMAX), POS_Z
      DOUBLE PRECISION :: DPS_DXE, DPS_DXW, DPS_DYN, DPS_DYS, DPS_DZT, DPS_DZB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, epg_min2, velf_part(dimn), VELP_INT(DIMN)

      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT, count_bc

      DOUBLE PRECISION MEANUS_e(DIMN, MMAX), MEANUS_w(DIMN, MMAX),MEANUS_n(DIMN, MMAX),MEANUS_s(DIMN, MMAX),MEANUS_t(DIMN, MMAX), MEANUS_b(DIMN, MMAX)
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5), epg_min_loc(1)

      double precision  sig_u, mean_u,ymid, part_taup,  DTPIC_MIN_X,  DTPIC_MIN_Y,  DTPIC_MIN_Z
      double precision, allocatable, dimension(:,:) ::  rand_vel

      double precision :: norm1, norm2, norm3
      Logical :: OUTER_STABILITY_COND, DES_FIXED_BED
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      OUTER_STABILITY_COND = .false.
      DES_FIXED_BED = .false.
      PC = 1
      FOCUS_PARTICLE = -1
      DTPIC_CFL = LARGE_NUMBER
      DTPIC_MIN_X =  LARGE_NUMBER
      DTPIC_MIN_Y =  LARGE_NUMBER
      DTPIC_MIN_Z =  LARGE_NUMBER

      if(dimn.eq.2) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(dimn.eq.3) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)
      allocate(rand_vel(MAX_PIP, DIMN))
      do idim = 1, dimn
         mean_u = zero
         sig_u = 1.d0
         CALL NOR_RNO(RAND_VEL(1:MAX_PIP, IDIM), MEAN_U, SIG_U)
      enddo
      !WRITE(*, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))
      !WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))

      DO L = 1, MAX_PIP
         DELETE_PART = .false.
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle
         pc = pc+1
         if(pea(l,4)) cycle

         DES_LOC_DEBUG = .FALSE.

         IF(.NOT.PEA(L,2))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
         ENDIF

         IF(DES_FIXED_BED) FC(:,L) = ZERO

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
         !By comparing the MFIX and equations in those papers,
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         IF(DES_ONEWAY_COUPLED) F_gp(L) = ZERO
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO


         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO
         M = PIJK(L,5)
         IJK = PIJK(L,4)
         IJK_OLD = IJK
         PIJK_OLD(:) = PIJK(L,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         VEL_ORIG(:) = DES_VEL_NEW(L,:)


         DES_VEL_NEW(L,:) = (DES_VEL_OLD(L,:) + &
         & FC(:,L)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         IF(DES_FIXED_BED) DES_VEL_NEW(L, :) = ZERO

         !MPPIC_VPTAU(L,:) = DES_VEL_NEW(L,:)

         VELP_INT(:) = DES_VEL_NEW(L,:)

         MEANVEL(1) = DES_U_S(IJK_OLD,M)
         MEANVEL(2) = DES_V_S(IJK_OLD,M)
         IF(DIMN.EQ.3) MEANVEL(3) = DES_W_S(IJK_OLD,M)

         RAD_EFF = DES_RADIUS(L)
         !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
         MEAN_FREE_PATH  = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2

         DO IDIM = 1, DIMN
            !SIG_U = 0.05D0*MEANVEL(IDIM)
            !DES_VEL_NEW(L, IDIM) = DES_VEL_NEW(L, IDIM) + SIG_U*RAND_VEL(L, IDIM )
            !PART_TAUP = RO_Sol(L)*((2.d0*DES_RADIUS(L))**2.d0)/(18.d0* MU_G(IJK))
            SIG_U = 0.005D0      !*MEANVEL(IDIM)
            RAND_VEL(L, IDIM)  = SIG_U*RAND_VEL(L, IDIM)*DES_VEL_NEW(L, IDIM)
            IF(DES_FIXED_BED) RAND_VEL(L,IDIM) = ZERO
            !rand_vel(L, idim)  = sig_u*rand_vel(L, idim)/part_taup
            !rand_vel(L, idim)  = sig_u* mean_free_path*rand_vel(L, idim)/part_taup
            DES_VEL_NEW(L, idim) = DES_VEL_NEW(L, idim) + rand_vel(L, idim)
         enddo

         !DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + VELP_INT(:)*DTSOLID

                                !if(mod(L,400).eq.0) write(*,'(A,2x,10(2x,g17.8))') 'rand vel = ', rand_vel(1,:), sig_u


         IF(.not.DES_ONEWAY_COUPLED.and.(.not.des_fixed_bed)) CALL MPPIC_APPLY_PS_GRAD_PART(L)


         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(L,1:DIMN), DES_VEL_NEW(L, 1:DIMN)))

         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then
         !   DES_VEL_NEW(L,:) = (DES_VEL_NEW(L,:)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
         DES_VEL_NEW(L,:)*DTSOLID !+ rand_vel(L, :)*dtsolid

         IF(DES_FIXED_BED) DES_POS_NEW(L,:) = DES_POS_OLD(L,:)

         CALL MPPIC_FIND_NEW_CELL(L)

         IJK = PIJK(L,4)

         !IF((EP_G(IJK).LT.0.35.and.fluid_at(ijk)).or.(ep_g(ijk_old).lt.0.35)) then !.and.(ijk.ne.ijk_old)) then
         !IF((EP_G(IJK).LT.EP_STAR.and.fluid_at(ijk)).and.(ijk.ne.ijk_old)) then
         INSIDE_DOMAIN = .true.
         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(CUT_CELL_AT(IJK)) THEN
            POS_Z = zero
            IF(DIMN.eq.3) POS_Z = DES_POS_NEW(L,3)
            CALL GET_DEL_H_DES(IJK,'SCALAR', &
            & DES_POS_NEW(L,1),  DES_POS_NEW(L,2), &
            & POS_Z, &
            & DIST, NORM1, NORM2, NORM3, .true.)

            IF(DIST.LE.ZERO) INSIDE_DOMAIN = .false.
         ENDIF

         !IF((EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN)) then
         !IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN) then
         IF(1.d0 - EP_G(IJK).GT. 1.3d0*(1.d0 - EP_STAR).and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD.and.OUTER_STABILITY_COND) then

            IF(CUT_CELL_AT(IJK)) then
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:) !+ rand_vel(L, :)*dtsolid
               DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
               !DES_VEL_NEW(L,:) = VELP_INT(:)
            ELSE
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               !DES_POS_NEW(L,:) = DES_POS_OLD(L,:) !+ rand_vel(L, :)*dtsolid
               DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
               !DES_VEL_NEW(L,:) = VELP_INT(:)
            ENDIF
         ENDIF

         PIJK(L,:) = PIJK_OLD(:)

         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)

         DIST = SQRT(DES_DOTPRDCT(D,D))

         !IF(DIST.GT.MEAN_FREE_PATH) THEN
          !WRITE(*,*) 'UPRIME OLD = ', UPRIMETAU(:)
          !WRITE(*,*) 'DIST GT MEAN FREE PATH= ', DIST, MEAN_FREE_PATH

         !    DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
         !    DES_VEL_NEW(L, :)*DTSOLID*MEAN_FREE_PATH/DIST

             !DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
             !VELP_INT(:)*DTSOLID*MEAN_FREE_PATH/DIST

                                !D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)

                                !DIST = SQRT(DES_DOTPRDCT(D,D))
             !WRITE(*,*) 'new moved distance  = ', dist,1. -  ep_g(ijk)
        !  ENDIF



         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         D_GRIDUNITS(1) = ABS(D(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(D(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DIMN.EQ.3) D_GRIDUNITS(3) = ABS(D(3))/DZ(PIJK(L,3))

         DIST = SQRT(DES_DOTPRDCT(D,D))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit

         DO IDIM = 1, DIMN
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN

                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(L,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(l,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(L,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC          = ', FC(:,L)

                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)

                  WRITE(*, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(L,:)
                  WRITE(*, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(l,:)
                  WRITE(*, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(L,:)
                  WRITE(*, '(A,2x,3(g17.8))') 'FC          = ', FC(:,L)
                  read(*,*)
                  DELETE_PART = .true.

               ENDIF
               !CALL mfix_exit(myPE)
            END IF

         END DO

         IF(.not.DELETE_PART) then
            DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
            DTPIC_MIN_X = MIN(DTPIC_MIN_X, DTPIC_TMPX)
            DTPIC_MIN_Y = MIN(DTPIC_MIN_Y, DTPIC_TMPY)
            DTPIC_MIN_Z = MIN(DTPIC_MIN_Z, DTPIC_TMPZ)
         ENDIF
         FC(:,L) = ZERO

         IF(DELETE_PART) THEN
            PEA(L,1) = .false.
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      DEALLOCATE(RAND_VEL)

      CALL global_all_max(DTPIC_CFL)
      PIP = PIP - PIP_DEL_COUNT

      LPIP_DEL_COUNT_ALL(:) = 0
      LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
      CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
      IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
         IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
         & FREQUENTLY: MONITOR THIS MESSAGE'

         WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
         & FREQUENTLY: MONITOR THIS MESSAGE'
         !DO IPROC = 0, NUMPES-1
         !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
         !ENDDO

      ENDIF

      DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

      IF(DTSOLID.GT.DTPIC_MAX) THEN
         !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
         IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
      ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN

         !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
         IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
      ELSE
        !WRITE(*,'(A40,2x,g17.8)') 'DT
        !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
      END IF

      WRITE(UNIT_LOG, '(10x,A,2x,3(g17.8))') 'DTPIC MINS IN EACH DIRECTION = ', DTPIC_MIN_X, DTPIC_MIN_Y, DTPIC_MIN_Z
      WRITE(*, '(10x,A,2x,3(g17.8))') 'DTPIC MINS IN EACH DIRECTION = ', DTPIC_MIN_X, DTPIC_MIN_Y, DTPIC_MIN_Z


 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

 2003 FORMAT(/10x, &
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)


      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC


!------------------------------------------------------------------------
! subroutine       : des_dbgpic
! Author           : Pradeep G.
! Purpose          : For debugging the pic values
! Parameters       : pstart - start indices of the particle
!                    pend - end indices of the particle
!                    pfreq - optional frequency (when the local count matches the
!                    frequency the filw will be written)
!                    if not send then it prints the file
!------------------------------------------------------------------------

      subroutine des_dbgpic (pstart,pend,pfreq)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement
      USE fldvar
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: pstart,pend
      INTEGER, INTENT(IN), OPTIONAL :: pfreq
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0
      character(30) :: filename
!-----------------------------------------------
      if (present(pfreq)) then
         lfreq = lfreq+1
         if (lfreq .ne. pfreq) return
         lfreq =0
      end if
      lfcount = lfcount + 1
      write(filename,'("debug",I3.3)') lfcount
      open (unit=100,file=filename)
      do lp = pstart,pend
         if (pea(lp,1) .and. .not.pea(lp,4)) then
            lijk = pijk(lp,4)
            write(100,*)"positon =",lijk,pijk(lp,1),pijk(lp,2), &
               pijk(lp,3),ep_g(lijk),DES_U_s(lijk,1)
            write(100,*)"forces =", FC(2,lp),tow(1,lp)
         end if
      end do
      close (100)

      RETURN
      END SUBROUTINE des_dbgpic


!------------------------------------------------------------------------
! subroutine       : des_dbgtecplot
! Author           : Pradeep G.
! Purpose          : prints the tecplot file for particle location
! Parameters       : pstart - start indices of the particle
!                    pend - end indices of the particle
!                    pfreq - optional frequency (when the local count matches the
!                    frequency the filw will be written)
!                    if not send then it prints the file
!------------------------------------------------------------------------

      subroutine des_dbgtecplot (pstart,pend,pfreq)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement
      USE fldvar
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: pstart,pend
      INTEGER, INTENT(IN), OPTIONAL :: pfreq
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0
      character(30) :: filename
!-----------------------------------------------

      if (present(pfreq)) then
         lfreq = lfreq+1
         if (lfreq .ne. pfreq) return
         lfreq =0
      end if
      lfcount = lfcount + 1
      write(filename,'("new_tec",I3.3,".dat")') lfcount
      open (unit=100,file=filename)
      write(100,'(9(A,3X),A)') &
         'VARIABLES = ', '"ijk"', '"x"', '"y"', '"vx"', '"vy"', &
         '"ep_g"', '"FCX"' ,'"FCY"', '"TOW"'
      write(100,'(A,F14.7,A)') 'zone T = "' , s_time , '"'
      do lp = pstart,pend
         if (pea(lp,1)) then
            lijk = pijk(lp,4)
            write(100,*)lijk,des_pos_new(lp,1),des_pos_new(lp,2), &
               des_vel_new(lp,1),des_vel_new(lp,2),ep_g(lijk),&
               fc(1,lp),fc(2,lp),tow(lp,1)
         endif
      enddo
      close (100)
      RETURN
      END SUBROUTINE DES_DBGTECPLOT



