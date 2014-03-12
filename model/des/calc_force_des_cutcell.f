!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES_CUTCELL                                 C
!
!  Purpose: DES calculations of force acting on a particle for cutcell C
!           treatment of boundaries.
!           This routine is still under development. Works well        C
!           for reasonably difficult geometries.                       C
!
!                                                                      C
!  Reviewer: Rahul Garg                             Date: 01-Jan-2012  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES_CUTCELL

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE softspring_funcs_cutcell
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I, J, LL, II, IW, IDIM, WALL_COUNT
      INTEGER NI, NLIM, N_NOCON, NEIGH_L
      INTEGER OVERLAP_MAXP

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP
      DOUBLE PRECISION FRAC_OVERLAP1, FRAC_OVERLAP2
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
                       V_REL_NORM_OLD, VRN_OLD(DIMN)
      DOUBLE PRECISION PFT_TMP(DIMN), FT_TMP(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
! tmp variables to calculate magnitude of forces
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION NORMAL(DIMN), NORMAL_OLD(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION DIST(DIMN), DIST_OLD(DIMN), DISTMOD, &
                       DISTMOD_OLD, R_LM
      DOUBLE PRECISION DTSOLID_TMP

! logic flag telling whether contact pair is old
      LOGICAL ALREADY_NEIGHBOURS
! logic flag for local debug warnings
      LOGICAL DES_LOC_DEBUG
      LOGICAL ALREADY_EXISTS
! index to track accounted for particles
      INTEGER PC
! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER PHASEI, PHASELL, IJK
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES,  ETAT_DES, KN_DES, KT_DES

      Double precision :: OVERLAP_FRAC, OVERLAP_FRAC_MAX_ALLOWED

      LOGICAL PARTICLE_SLIDE
!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------

      INCLUDE 'function.inc'

! Calculate new values
!---------------------------------------------------------------------
      OVERLAP_FRAC_MAX_ALLOWED = 0.05
      OVERLAP_MAXP = UNDEFINED_I
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
      !DES_LOC_DEBUG = .true. ;    DEBUG_DES = .true.
      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.

      FOCUS_PARTICLE = 87
      IF (S_TIME.LE.DTSOLID) THEN
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
      ENDIF

      IF(USE_STL) then
         CALL CALC_FORCE_WITH_WALL_CUTFACE_STL
      ELSE
         CALL CALC_FORCE_WITH_WALL_CUTFACE
      ENDIF
!     Calculate contact force and torque
!---------------------------------------------------------------------

      PC = 1
      DO LL = 1, MAX_PIP

         IF(LL.eq.focus_particle) then
            !debug_des = .true.
         else
            debug_des = .false.
         endif
         IF (PEA(LL,1)) PC = PC + 1
         IF(.NOT.PEA(LL,1) .OR. PEA(LL,4)) CYCLE



         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO

         PFT_TMP(:) = ZERO

         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact
! after the previous call of CALC_FORCE_DES.
            DO NI = 2, NLIM
               IF(PV(LL,NI-N_NOCON).EQ.0) THEN
! For each particle in PN(LL,2:MAXNEIGHBORS) that is no longer in
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(LL,1) by one.
                  PN(LL,(NI-N_NOCON):MAXNEIGHBORS-1) = &
                     PN(LL,(NI-N_NOCON+1):MAXNEIGHBORS)
                  PV(LL,(NI-N_NOCON):(MAXNEIGHBORS-1)) = &
                     PV(LL,(NI-N_NOCON+1):MAXNEIGHBORS)
                  PFT(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFT(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
               ENDIF
            ENDDO
         ENDIF

! Initializing rest of the neighbor list which is not in contact and
! clean up after the above array left shifts
         NLIM = MAX(2,PN(LL,1) + 2)
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
         IF (PN(LL,1) .GT. NEIGH_MAX) NEIGH_MAX = PN(LL,1)


! Initializing the neighbor list contact information when particles are
! not in contact; i.e. when particle LL has no neighbors
         IF (PN(LL,1).EQ.0) THEN
            PFT(LL,:,:) = ZERO
         ENDIF

! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV(LL,2:MAXNEIGHBORS) = 0



! Check particle LL neighbour contacts
!---------------------------------------------------------------------
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)

               !IF(I.GT.LL .AND. PEA(I,1)) THEN

               IF(PEA(I,1)) THEN

                  ALREADY_NEIGHBOURS=.FALSE.

                  IF(PN(LL,1).GT.0) THEN
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF


                  R_LM = DES_RADIUS(LL) + DES_RADIUS(I)
                  DIST(:) = DES_POS_NEW(I,:) - DES_POS_NEW(LL,:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))


                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN

                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,'(5X,A,I10)') 'NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF

                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX) THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
                        OVERLAP_MAXP = LL
                     ENDIF

                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_RADIUS(LL)
                     FRAC_OVERLAP2 = (R_LM-DISTMOD)/DES_RADIUS(I)
                     IF (FRAC_OVERLAP1 > flag_overlap .OR. &
                         FRAC_OVERLAP2 > flag_overlap) THEN
!                        WRITE(*,'(5X,A,A,ES15.7)') &
!                           'WARNING: excessive overlap detected ', &
!                           'at time ', S_TIME
!                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
!                           'between particles ', LL, 'and ',&
!                           I, 'with'
!                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
!                           'overlap = ', (R_LM-DISTMOD), &
!                           ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,'(5X,A,I10,I10)') &
                           'DISTMOD is zero between particle-pair ',&
                           LL, I
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair.
                     CALL CFRELVEL2(LL, I, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD, 0)

! Overlap calculation changed from history based to current position
                     OVERLAP_N = R_LM-DISTMOD

                     IF(ALREADY_NEIGHBOURS) THEN
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE
                        IF(DEBUG_DES) THEN
                           IF (.NOT.DES_LOC_DEBUG) THEN
                              DES_LOC_DEBUG = .TRUE.
                              WRITE(*,1000)
                           ENDIF
                           WRITE(*,'(5X,A,2(I10,X),A,ES15.7)') &
                              'Normal overlap for particle pair ',&
                              LL, I, ' : ', OVERLAP_N
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                           DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          DIST_OLD(:)=DES_POS_OLD(I,:)-DES_POS_OLD(LL,:)
                          DISTMOD_OLD=SQRT(DES_DOTPRDCT(DIST_OLD,DIST_OLD))
                          NORMAL_OLD(:)=DIST_OLD(:)/DISTMOD_OLD
                          VRN_OLD(:)=DES_VEL_OLD(LL,:)-DES_VEL_OLD(I,:)
                          V_REL_NORM_OLD=DES_DOTPRDCT(VRN_OLD,NORMAL_OLD)
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I10,2X,A)') &
                             'for first contact between particles', LL, &
                             'and ', I, 'with'
                          WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM, &
                             'and V_REL_NORM_OLD = ', V_REL_NORM_OLD
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF
                  ELSE
                     GOTO 300
                  ENDIF

                  phaseLL = PIJK(LL,5)
                  phaseI = PIJK(I,5)

! T.Li : Hertz vs linear spring-dashpot contact model
                  IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
                     sqrt_overlap = SQRT(OVERLAP_N)
                     KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
                     KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
                     sqrt_overlap = SQRT(sqrt_overlap)
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap
                  ELSE
                     KN_DES = KN
                     KT_DES = KT
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)
                  ENDIF

                  FNS1(:) = -KN_DES * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
                  FNORM(:) = FNS1(:) + FNS2(:)

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed
                  PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T * TANGENT(:)
                  PFT_TMP(:) = PFT(LL,NI,:)   ! for an easy pass to des_dotprdct
                  PFT_TMP(:) = PFT(LL,NI,:) - &
                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
                  FTS1(:) = -KT_DES * PFT_TMP(:)
                  FTS2(:) = -ETAT_DES * V_REL_TRANS_TANG * TANGENT(:)
                  FTAN(:) = FTS1(:) + FTS2(:)

! Check for Coulombs friction law and limit the maximum value of the tangential
! force on a particle in contact with another particle
                  CALL CFSLIDE2(TANGENT, PARTICLE_SLIDE)

! Calculate the total force FC and TOW on a particle in a particle-particle collision
                  CALL CFFCTOW2(LL, I, NORMAL, DISTMOD)

! Save the tangential displacement history with the correction of Coulomb's law
                  IF (PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = -( FTAN(:) - FTS2(:) ) / KT_DES
                  ELSE
                     PFT(LL,NI,:) = PFT_TMP(:)
                  ENDIF

                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000)
                     ENDIF

                     WRITE(*,*)'NORMAL =         ;', NORMAL
                     WRITE(*,*)'TANGNT =         ;', TANGENT
                     PRINT*, '     EtaN, EtaT =  ', ETAN_DES, ETAT_DES
                     PRINT*, '     Percent overlap = ', (R_LM - DISTMOD)*100.d0/R_LM
                     PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, '     FTS1 and FTS2 = ', FTS1(:), FTS2(:)
                     PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, '     PFT = ', PFT(LL,NI,:)
                     PRINT*, '     FORCEST = ', FTAN(:)
                     PRINT*, '     FORCESN = ', FNORM(:)
                     READ(*,*)
                  ENDIF

                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FNORM(1), 'FNy=',FNORM(2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FNORM(1), 'FNy=',FNORM(2)
                     ENDIF
                     CLOSE (1)
                     PRINT*, 'PN', PN(LL,:)
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF         ! IF (I>LL .AND. PEA(I,1))

 300           CONTINUE
            ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF               ! IF(NEIGHBOURS(LL,1).GT.0)

!---------------------------------------------------------------------
! End check particle LL neighbour contacts
      ENDDO   ! end loop over paticles LL


      !IF(ABS(FC(1,1)).GT.ZERO)
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties

! Calculate gas-solids drag force on particle
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL CALC_DES_DRAG_GS
      ENDIF

! Calculate solids-solids drag force on particle
      IF (DES_CONTINUUM_HYBRID) THEN
         CALL CALC_DES_DRAG_SS
      ENDIF

! Update the old values of particle position and velocity with the new
! values computed
      CALL CFUPDATEOLD

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------')

      RETURN
      END SUBROUTINE CALC_FORCE_DES_CUTCELL
