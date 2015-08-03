!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CFNEWVALUES                                            !
!                                                                      !
!  Purpose: DES - Calculate the new values of particle velocity,       !
!           position, angular velocity etc                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTEGRATE_TIME_PIC

      USE discretelement
      USE mfix_pic

      IMPLICIT NONE

! call the coloring function like approach
      IF(MPPIC_SOLID_STRESS_SNIDER) THEN
         CALL INTEGRATE_TIME_PIC_SNIDER
      ELSE
         CALL INTEGRATE_TIME_PIC_GARG
      ENDIF

      RETURN
      END SUBROUTINE INTEGRATE_TIME_PIC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INTEGRATE_TIME_PIC_SNIDER                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INTEGRATE_TIME_PIC_SNIDER

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE fldvar
      USE cutcell
      USE mfix_pic
      USE randomno
      USE geometry, only: DO_K, NO_K
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: NP, LC
      INTEGER :: M, IDIM
      INTEGER I, J, K, IJK

      DOUBLE PRECISION :: DP_BAR, ENp1, SLIPVEL(3), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(3), UPRIMETAU(3), MFP
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION :: THREEINTOSQRT2
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ

      LOGICAL :: DELETE_PART
      INTEGER :: PIP_DEL_COUNT

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)

      DOUBLE PRECISION :: DIST(3), DIST_MAG, VEL(3)

      DOUBLE PRECISION :: VELOLD1(3), POSOLD(3)
      DOUBLE PRECISION :: VELOLD2(3)
      DOUBLE PRECISION :: VELOLD3(3)
      DOUBLE PRECISION :: VELOLD4(3)
      DOUBLE PRECISION :: VELOLD5(3)
      DOUBLE PRECISION :: VELOLD6(3)

!-----------------------------------------------

      PC = 1
      DTPIC_CFL = LARGE_NUMBER

      if(NO_K) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(DO_K) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)

      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      DES_LOC_DEBUG = .FALSE.

      ENp1 = MPPIC_COEFF_EN1 + 1.0d0

      DO NP = 1, MAX_PIP
         DELETE_PART = .FALSE.

! pradeep skip ghost particles
         IF(PC.GT.PIP) EXIT
         IF(IS_NONEXISTENT(NP)) CYCLE

         PC = PC+1


! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(IS_ENTERING(NP))THEN
            FC(:,NP) = ZERO
         ELSE
            FC(:,NP) = FC(:,NP)/PMASS(NP) + GRAV(:)
         ENDIF

! DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
! By comparing MFIX equations and equations in those papers,
! D_p = Beta/(EP_S*RHOP)
! F_gp in drag_fgs.f  = Beta*PVOL/EP_S
! Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(NP)/(PMASS(NP))

         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO
         IF(.NOT.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO

         M = PIJK(NP,5)
         IJK = PIJK(NP,4)

!         IF(ep_g(IJK) < ep_star) then
!            write(*,*) 'Close pack: ',NP
!         endif

! Numerically integrated particle velocity without particle normal
! stress force :: Snider (Eq 38)
         VEL(:) = (DES_VEL_NEW(:,NP) + &
            FC(:,NP)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

! Estimated discrete particle velocity from continuum particle
! normal stress graident :: Snider (Eq 39)
         DELUP(:) = -(DTSOLID*PS_GRAD(NP,:)) / &
            ((1.0-EP_G(IJK))*RO_S0(1)*(1.d0+DP_BAR*DTSOLID))

! Slip velocity between parcel and bulk solids
         SLIPVEL = AVGSOLVEL_P(:,NP) - VEL

! Calculate particle velocity from particle normal stress.
         DO LC = 1, merge(2,3,NO_K)
! Snider (40)
            IF(PS_GRAD(NP,LC) <= ZERO) THEN
               UPRIMETAU(LC) = min(DELUP(LC), ENp1*SLIPVEL(LC))
               UPRIMETAU(LC) = max(UPRIMETAU(LC), ZERO)
            ELSE
! Snider (41)
               UPRIMETAU(LC) = max(DELUP(LC), ENp1*SLIPVEL(LC))
               UPRIMETAU(LC) = min(UPRIMETAU(LC), ZERO)
            ENDIF
         ENDDO

! Total particle velocity is the sum of the particle velocity from
! normal stress and velocity from all other foces :: Sinder (Eq 37)
         DES_VEL_NEW(:,NP) = VEL + UPRIMETAU

! Total distance traveled over the current time step
         DIST(:) = DES_VEL_NEW(:,NP)*DTSOLID
         DIST_MAG = dot_product(DIST, DIST)

! Calculate the mean free path :: Snider (43)
         MFP = DES_RADIUS(NP)/(THREEINTOSQRT2*(ONE-EP_G(IJK)))

! Limit parcel movement to mean free path if needed.
         IF(MFP**2 < DIST_MAG) THEN




!      if(mfp*2.0 < sqrt(dist_mag)) then
!      write(*,3000) NP, MFP, sqrt(dist_mag), &
!         DES_VEL_NEW(:,NP), DIST/DTSOLID
!
! 3000 FORMAT(2/,3x,'Limiting movement: ',I8,/&
!                5x,'Mean Free Path: ',f15.6,/&
!                5x,'DIST_MAG:       ',f15.6,/&
!                5x,'Velocity: ',3(2x,f15.4),/&
!                5x,'Dist/DT:  ',3(2x,f15.4))
!      endif

            DIST = MFP*DIST/sqrt(DIST_MAG)
            DES_VEL_NEW(:,NP) = DIST/DTSOLID
         ENDIF

! Update the parcel position.
         DES_POS_NEW(:,NP) = DES_POS_NEW(:,NP) + DIST(:)

! Bin the particle position to the fuild grid.
         CALL PIC_FIND_NEW_CELL(NP)


         D_GRIDUNITS(1) = ABS(DIST(1))/DX(PIJK(NP,1))
         D_GRIDUNITS(2) = ABS(DIST(2))/DY(PIJK(NP,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DO_K) D_GRIDUNITS(3) = ABS(DIST(3))/DZ(PIJK(NP,3))

         DTPIC_TMPX = (CFL_PIC*DX(PIJK(NP,1)))/(ABS(DES_VEL_NEW(1,NP))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(NP,2)))/(ABS(DES_VEL_NEW(2,NP))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(NP,3)))/(ABS(DES_VEL_NEW(3,NP))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit
         DO LC = 1, merge(2,3,NO_K)
            IF(D_GRIDUNITS(LC).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN


                  write(*,1200)NP,  VELOLD1, VELOLD2, VELOLD3, VELOLD4,&
                     POSOLD,  DES_POS_NEW(:,NP), &
                     DIST, sqrt(DIST_MAG), sqrt(MFP) 

 1200 FORMAT(2/,2x,70('*')/2x,'Undesired parcel movement: ',I8,&
         /5x,'VEL1: ',3(3x,f12.4),/5x,'VEL2: ',3(3x,f9.3), &
         /5x,'VEL3: ',3(3x,f12.4),/5x,'VEL4: ',3(3x,f9.3), & 
         /5x,'POLD: ',3(3x,f12.4),/5x,'PNEW: ',3(3x,f9.3), & 
         /5x,'DD:   ',3(3x,f12.4),&
         /5x,'D_MAG:',1(3x,f12.4),& 
         /5x,'MFP:  ',1(3x,f12.4)) 
        
               ENDIF
            END IF

         END DO

         IF(.not.DELETE_PART) DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         FC(:,NP) = ZERO

         IF(DELETE_PART) THEN
            CALL SET_NONEXISTENT(NP)
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      CALL global_all_max(DTPIC_CFL)
      PIP = PIP - PIP_DEL_COUNT

      LPIP_DEL_COUNT_ALL(:) = 0
      LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
      CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
      IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
         IF(PRINT_DES_SCREEN) THEN
            WRITE(*,'(/,2x,A,2x,i10,/,A)') &
               'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                    'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
            ENDIF

            WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') &
                 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                 'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'

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
 1001 FORMAT(3X,'<---------- END INTEGRATE_TIME_PIC ----------')

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN INTEGRATE_TIME_PIC: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

      RETURN
      END SUBROUTINE INTEGRATE_TIME_PIC_SNIDER






























      SUBROUTINE INTEGRATE_TIME_PIC_GARG

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
      USE randomno
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER NP, M, LC
      INTEGER I, J, K, IJK, IJK_OLD

      DOUBLE PRECISION DD(3), DIST, &
                       DP_BAR, MEANVEL(3), D_GRIDUNITS(3)

      DOUBLE PRECISION MEAN_FREE_PATH
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION UPRIMEMOD

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, POS_Z
      DOUBLE PRECISION :: VELP_INT(3)

      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5)

      double precision  sig_u, mean_u, DTPIC_MIN_X,  DTPIC_MIN_Y,  DTPIC_MIN_Z
      double precision, allocatable, dimension(:,:) ::  rand_vel

      double precision :: norm1, norm2, norm3
      Logical :: OUTER_STABILITY_COND, DES_FIXED_BED
!-----------------------------------------------

      OUTER_STABILITY_COND = .false.
      DES_FIXED_BED = .false.
      PC = 1
      DTPIC_CFL = LARGE_NUMBER
      DTPIC_MIN_X =  LARGE_NUMBER
      DTPIC_MIN_Y =  LARGE_NUMBER
      DTPIC_MIN_Z =  LARGE_NUMBER

      if(NO_K) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(DO_K) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)
      allocate(rand_vel(3, MAX_PIP))
      do LC = 1, merge(2,3,NO_K)
         mean_u = zero
         sig_u = 1.d0
         CALL NOR_RNO(RAND_VEL(LC, 1:MAX_PIP), MEAN_U, SIG_U)
      enddo
      !WRITE(*, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))
      !WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))

      DO NP = 1, MAX_PIP
         DELETE_PART = .false.
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(is_nonexistent(NP)) cycle
         pc = pc+1
         if(is_ghost(NP) .or. is_entering_ghost(NP) .or. is_exiting_ghost(NP)) cycle

         DES_LOC_DEBUG = .FALSE.

         IF(.NOT.IS_ENTERING(NP))THEN
            FC(:,NP) = FC(:,NP)/PMASS(NP) + GRAV(:)
         ELSE
            FC(:,NP) = ZERO
         ENDIF

         IF(DES_FIXED_BED) FC(:,NP) = ZERO

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
         !By comparing the MFIX and equations in those papers,
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(NP)/(PMASS(NP))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO
         M = PIJK(NP,5)
         IJK = PIJK(NP,4)
         IJK_OLD = IJK
         PIJK_OLD(:) = PIJK(NP,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         DES_VEL_NEW(:,NP) = (DES_VEL_NEW(:,NP) + &
         & FC(:,NP)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         IF(DES_FIXED_BED) DES_VEL_NEW(:,NP) = ZERO

         !MPPIC_VPTAU(NP,:) = DES_VEL_NEW(:,NP)

         VELP_INT(:) = DES_VEL_NEW(:,NP)

         MEANVEL(1) = DES_U_S(IJK_OLD,M)
         MEANVEL(2) = DES_V_S(IJK_OLD,M)
         IF(DO_K) MEANVEL(3) = DES_W_S(IJK_OLD,M)

         RAD_EFF = DES_RADIUS(NP)
         !RAD_EFF = (DES_STAT_WT(NP)**(1.d0/3.d0))*DES_RADIUS(NP)
         IF(.not.DES_ONEWAY_COUPLED) then
            MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
            MEAN_FREE_PATH  = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2
         endif

         DO LC = 1, merge(2,3,NO_K)
            !SIG_U = 0.05D0*MEANVEL(LC)
            !DES_VEL_NEW(LC,NP) = DES_VEL_NEW(LC,NP) + SIG_U*RAND_VEL(LC, L )
            !PART_TAUP = RO_Sol(NP)*((2.d0*DES_RADIUS(NP))**2.d0)/(18.d0* MU_G(IJK))
            SIG_U = 0.005D0      !*MEANVEL(LC)
            RAND_VEL(LC, NP)  = SIG_U*RAND_VEL(LC, NP)*DES_VEL_NEW(LC,NP)
            IF(DES_FIXED_BED) RAND_VEL(LC,NP) = ZERO
            !rand_vel(LC, NP)  = sig_u*rand_vel(LC, NP)/part_taup
            !rand_vel(LC, NP)  = sig_u* mean_free_path*rand_vel(LC, NP)/part_taup
            DES_VEL_NEW(LC,NP) = DES_VEL_NEW(LC,NP) + rand_vel(LC, NP)
         enddo

         IF(.not.DES_ONEWAY_COUPLED.and.(.not.des_fixed_bed)) CALL MPPIC_APPLY_PS_GRAD_PART(NP)

         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(1:,NP), DES_VEL_NEW(1:,NP)))

         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then
         !   DES_VEL_NEW(:,NP) = (DES_VEL_NEW(:,NP)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

         IF(DES_FIXED_BED) THEN
            DD(:) = ZERO
         ELSE
            DD(:) = DES_VEL_NEW(:,NP)*DTSOLID !+ rand_vel(:, NP)*dtsolid
         ENDIF

         CALL PIC_FIND_NEW_CELL(NP)

         IJK = PIJK(NP,4)

         !IF((EP_G(IJK).LT.0.35.and.fluid_at(ijk)).or.(ep_g(ijk_old).lt.0.35)) then !.and.(ijk.ne.ijk_old)) then
         !IF((EP_G(IJK).LT.EP_STAR.and.fluid_at(ijk)).and.(ijk.ne.ijk_old)) then
         INSIDE_DOMAIN = .true.
         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(CUT_CELL_AT(IJK)) THEN
            POS_Z = zero
            IF(DO_K) POS_Z = DES_POS_NEW(3,NP)
            CALL GET_DEL_H_DES(IJK,'SCALAR', &
            & DES_POS_NEW(1,NP),  DES_POS_NEW(2,NP), &
            & POS_Z, &
            & DIST, NORM1, NORM2, NORM3, .true.)

            IF(DIST.LE.ZERO) INSIDE_DOMAIN = .false.
         ENDIF

         !IF((EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN)) then
         !IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN) then
         IF(1.d0 - EP_G(IJK).GT. 1.3d0*(1.d0 - EP_STAR).and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD.and.OUTER_STABILITY_COND) then
            DD(:) = ZERO
            DES_VEL_NEW(:,NP) = 0.8d0*DES_VEL_NEW(:,NP)
         ENDIF

         PIJK(NP,:) = PIJK_OLD(:)

         DES_POS_NEW(:,NP) = DES_POS_NEW(:,NP) + DD(:)

         DIST = SQRT(DOT_PRODUCT(DD,DD))

         D_GRIDUNITS(1) = ABS(DD(1))/DX(PIJK(NP,1))
         D_GRIDUNITS(2) = ABS(DD(2))/DY(PIJK(NP,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DO_K) D_GRIDUNITS(3) = ABS(DD(3))/DZ(PIJK(NP,3))

         DIST = SQRT(DOT_PRODUCT(DD,DD))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(NP,1)))/(ABS(DES_VEL_NEW(1,NP))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(NP,2)))/(ABS(DES_VEL_NEW(2,NP))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(NP,3)))/(ABS(DES_VEL_NEW(3,NP))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit

         DO LC = 1, merge(2,3,NO_K)
            IF(D_GRIDUNITS(LC).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN

                  WRITE(UNIT_LOG, 2001) NP, D_GRIDUNITS(:), DES_VEL_NEW(:,NP)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(:,NP)
                  IF (DO_OLD) WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(:,NP)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(:,NP)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC          = ', FC(:,NP)

                  WRITE(*, 2001) NP, D_GRIDUNITS(:), DES_VEL_NEW(:,NP)

                  WRITE(*, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(:,NP)
                  IF (DO_OLD) WRITE(*, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(:,NP)
                  WRITE(*, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(:,NP)
                  WRITE(*, '(A,2x,3(g17.8))') 'FC          = ', FC(:,NP)
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
         FC(:,NP) = ZERO

         IF(DELETE_PART) THEN
            CALL SET_NONEXISTENT(NP)
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
         IF(PRINT_DES_SCREEN) THEN
            WRITE(*,'(/,2x,A,2x,i10,/,A)') &
                 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                 'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
         ENDIF

         WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') &
              'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), &
              'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
         !DO IPROC = 0, NUMPES-1
         !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
         !     'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
         !ENDDO

      ENDIF

      DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

      RETURN
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

 1001 FORMAT(3X,'<---------- END INTEGRATE_TIME_PIC ----------')

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN INTEGRATE_TIME_PIC: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

      RETURN
      END SUBROUTINE INTEGRATE_TIME_PIC_GARG


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
      character(255) :: filename
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
         if (is_normal(lp) .or. is_entering(lp) .or. is_exiting(lp)) then
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
      character(255) :: filename
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
         if (.not.is_nonexistent(lp)) then
            lijk = pijk(lp,4)
            write(100,*)lijk,des_pos_new(1,lp),des_pos_new(2,lp), &
               des_vel_new(1,lp),des_vel_new(2,lp),ep_g(lijk),&
               fc(1,lp),fc(2,lp),tow(lp,1)
         endif
      enddo
      close (100)
      RETURN
      END SUBROUTINE DES_DBGTECPLOT



