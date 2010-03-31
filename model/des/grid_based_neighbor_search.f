!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GRID_BASED_NEIGHBOR_SEARCH                             C
!  Purpose: Cell linked search 
!                                                                      C
!  Author: Rahul Garg                                 Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 24-OCT-08  C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      
      SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
      
      USE param1
      USE discretelement
      USE geometry
      USE des_bc

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      LOGICAL ALREADY_NEIGHBOURS

      INTEGER PC ! Index to track accounted for particles
      INTEGER I, II, IP1, IM1   ! X-coordinate loop indices
      INTEGER J, JJ, JP1, JM1   ! Y-coordinate loop indices
      INTEGER K, KK, KP1, KM1   ! Z-coordinate loop indices
      INTEGER PNO ! Temp. particle number variable
      INTEGER NPG ! Temp. cell particle count
      INTEGER LL, NP, NEIGH_L  ! Loop Counters
      INTEGER NLIM ! 

! min and max fluid cell index for grid based search mesh
      INTEGER DESGS_IMIN1, DESGS_IMAX1, DESGS_JMIN1, DESGS_JMAX1, &
              DESGS_KMIN1, DESGS_KMAX1
! min cell index for grid based search mesh
      INTEGER DESGS_IMIN2,DESGS_JMIN2, DESGS_KMIN2

      DOUBLE PRECISION DISTVEC(DIMN), DIST, R_LM ! Contact variables
      DOUBLE PRECISION LX, LY, LZ ! System dimensions
      DOUBLE PRECISION XPER_FAC, YPER_FAC, ZPER_FAC !Periodic BC info
      DOUBLE PRECISION CORD_PNO(DIMN) !Periodic BC neighbor check

!-----------------------------------------------
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------  

! Determine system dimensions for periodic BC neighbor search
! Look at cells imin2=1 through imax1=imax+1
      LX = XLENGTH ! XE(IMAX1) - XE(1)
      LY = YLENGTH ! YN(JMAX1) - YN(1)
      LZ = ZLENGTH ! ZT(KMAX1) - ZT(1)


      DESGS_IMAX1 = DESGRIDSEARCH_IMAX+1
      DESGS_IMIN1 = 2
      DESGS_IMIN2 = 1
      DESGS_JMAX1 = DESGRIDSEARCH_JMAX+1
      DESGS_JMIN1 = 2
      DESGS_JMIN2 = 1
      IF (DIMN .EQ. 2) THEN
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX
         DESGS_KMIN1 = DESGRIDSEARCH_KMAX
         DESGS_KMIN2 = DESGRIDSEARCH_KMAX
      ELSE
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX+1
         DESGS_KMIN1 = 2
         DESGS_KMIN2 = 1
      ENDIF

! Loop only over particles in the system
!   cycle if there is no particle in the LL seat
!   skip new particles entering the system and particles 
!   exiting the system
      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         II = DESGRIDSEARCH_PIJK(LL,1)
         JJ = DESGRIDSEARCH_PIJK(LL,2)
         KK = DESGRIDSEARCH_PIJK(LL,3)
! Initialize and set II, JJ, KK indices
         IP1=II+1; IM1=II-1
! Initialize and set JJ indices
         JP1=JJ+1; JM1=JJ-1
! Initialize and set KK indices
         KP1=KK;   KM1=KK
         IF(DIMN.EQ.3)THEN 
            KP1 = KK+1;   KM1 = KK-1
         ENDIF

! Check fluid cells neighboring PIJK(LL,4) for potential neighbors of LL
         DO KK = KM1, KP1
            DO JJ = JM1, JP1
               DO II = IM1, IP1

! Shift loop index to new variables for manipulation
                  I = II;   J = JJ;   K = KK

! Initizalize PER_FAC values 
                  XPER_FAC = 0; YPER_FAC = 0; ZPER_FAC = 0

! Correct indices for neighbor cells that are not in the domain and
! also set up periodic BC neighbor search, unless DES_MI is true
! DES_MI currently not setup with periodic conditions
                  IF (.NOT. DES_MI) THEN

                     IF(II.GT.DESGS_IMAX1) THEN 
                        IF(DES_PERIODIC_WALLS_X) THEN 
                           I = DESGS_IMIN1  ! min fluid cell no. 
                           XPER_FAC = one
                        ELSE
                           I = DESGS_IMAX1  ! max fluid cell no.
                        ENDIF
                     ENDIF
                     IF(II.LT.DESGS_IMIN1) THEN 
                        IF(DES_PERIODIC_WALLS_X) THEN 
                           I = DESGS_IMAX1
                           XPER_FAC = -one
                        ELSE 
                           I = DESGS_IMIN1
                        ENDIF
                     ENDIF
                     IF(JJ.GT.DESGS_JMAX1) THEN
                        IF(DES_PERIODIC_WALLS_Y) THEN 
                           J = DESGS_JMIN1
                           YPER_FAC = one
                        ELSE
                           J = DESGS_JMAX1
                        ENDIF
                     ENDIF
                     IF(JJ.LT.DESGS_JMIN1) THEN
                        IF(DES_PERIODIC_WALLS_Y) THEN 
                           J = DESGS_JMAX1
                           YPER_FAC = -one
                        ELSE
                           J = DESGS_JMIN1
                        ENDIF
                     ENDIF
                     IF(DIMN.EQ.3) THEN 
                        IF(KK.GT.DESGS_KMAX1) THEN 
                           IF(DES_PERIODIC_WALLS_Z) THEN
                              K = DESGS_KMIN1
                              ZPER_FAC = one
                           ELSE
                              K = DESGS_KMAX1
                           ENDIF
                        ENDIF
                        IF(KK.LT.DESGS_KMIN1) THEN 
                           IF(DES_PERIODIC_WALLS_Z) THEN
                              K = DESGS_KMAX1
                              ZPER_FAC = -one
                           ELSE
                              K = DESGS_KMIN1
                           ENDIF
                        ENDIF
                     ENDIF

                  ELSE  ! DES_MI == TRUE

                     IF(II.GT.DESGS_IMAX2) THEN 
                        I = DESGS_IMAX2  ! max cell no.
                     ENDIF
                     IF(II.LT.DESGS_IMIN2) THEN 
                        I = DESGS_IMIN2  ! min cell no.
                     ENDIF
                     IF(JJ.GT.DESGS_JMAX2) THEN
                        J = DESGS_JMAX2
                     ENDIF
                     IF(JJ.LT.1) THEN
                        J = DESGS_JMIN2
                     ENDIF
                     IF(DIMN.EQ.3) THEN 
                        IF(KK.GT.DESGS_KMAX2) THEN 
                           K = DESGS_KMAX2
                        ENDIF
                        IF(KK.LT.DESGS_KMIN2) THEN 
                           K = DESGS_KMIN2
                        ENDIF
                     ENDIF
                  ENDIF


! If cell IJK contains particles, store the amount in NPG
                  IF(ASSOCIATED(DESGRIDSEARCH_PIC(I,J,K)%P))THEN
                     NPG = SIZE(DESGRIDSEARCH_PIC(I,J,K)%P)
                  ELSE
                     NPG = 0
                  ENDIF

! Loop over the particles in IJK cell to determine if they are
! neighbors to particle LL
                  DO NP = 1,NPG

                     PNO = DESGRIDSEARCH_PIC(I,J,K)%P(NP)
                     IF(.NOT.PEA(PNO,1)) CYCLE 

                     IF(PNO.GT.LL)THEN 

! Calculate the distance between particles, accounting for neighbors 
! with a periodic boundary condition (wrap arround)
                        R_LM = DES_RADIUS(LL) + DES_RADIUS(PNO)
                        R_LM = FACTOR_RLM*R_LM
                        CORD_PNO(1) = DES_POS_NEW(PNO,1) + XPER_FAC*LX
                        CORD_PNO(2) = DES_POS_NEW(PNO,2) + YPER_FAC*LY
                        IF(DIMN.EQ.3) THEN 
                           CORD_PNO(3) = DES_POS_NEW(PNO,3) + ZPER_FAC*LZ
                        ENDIF 
                        DISTVEC(:) = CORD_PNO(:) - DES_POS_NEW(LL,:)
                        DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))

! Check if LL and PNO are already neighbors
                        ALREADY_NEIGHBOURS = .FALSE.
                        DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1
                           IF(PNO .EQ. NEIGHBOURS(LL,NEIGH_L)) ALREADY_NEIGHBOURS=.TRUE.
                        ENDDO

! Check the contact condition       
                        IF( (DIST .LE. R_LM) .AND. &
                        (.NOT.ALREADY_NEIGHBOURS)) THEN 

! Check that the array is not full and store PNO                           
                           NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
                           NLIM  = NEIGHBOURS(LL,1) + 1
                           IF(NLIM.GT.MN) THEN 
                              WRITE(*,1000) NLIM, MN, LL
                              STOP 
                           ELSE
                              NEIGHBOURS(LL,NLIM) = PNO
                           ENDIF

! Check that the array is not full and store LL
                           NEIGHBOURS(PNO,1) = NEIGHBOURS(PNO,1) + 1
                           NLIM  = NEIGHBOURS(PNO,1) + 1
                           IF(NLIM.GT.MN) THEN 
                              WRITE(*,1000) NLIM, MN, PNO
                              STOP 
                           ELSE
                              NEIGHBOURS(PNO,NLIM) = LL
                           ENDIF
                        ENDIF  !contact condition

                     ENDIF  !PNO.GT.LL

                  ENDDO  !NP

               ENDDO  ! II cell loop
            ENDDO  ! JJ cell loop
         ENDDO  ! KK cell loop
         PC = PC + 1

      ENDDO  ! Particles in system loop

 1000 FORMAT(/1X,70('*')//&
         ' From: GRID_BASED_NEIGHBOR_SEARCH -',/&
         ' Message: Neighbors(=', I4,') > MN(=', I4,&
         ') for particle: ', I,/&
         ' Either reduce the factor R_LM or increase ',&
         'MN in mfix.dat'/,&
         1X,70('*')/)         

      END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
