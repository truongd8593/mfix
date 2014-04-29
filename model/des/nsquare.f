!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NSQUARE                                                C
!>  Purpose: DES - N-Square neighbor search  
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NSQUARE

      USE param1
      USE discretelement
      USE geometry
      USE des_bc      
      Use des_thermo
      use geometry, only: DO_K

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I, J, K
      INTEGER L, LL, NEIGH_L, NEIGH_LL
      DOUBLE PRECISION DISTVEC(3), DIST, R_LM
! Temporary variables to adjust particle position in event of periodic
! boundaries
      DOUBLE PRECISION XPOS(2), YPOS(2), ZPOS(2), TMPPOS(3)
! Max loop limit in each coordinate direction
      INTEGER II, JJ, KK
! Index to track accounted for particles
      INTEGER PC
! Index to track unaccounted for particles
      INTEGER PNPC
!-----------------------------------------------
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------  

      PC=1
      DO L=1, MAX_PIP

         IF(PC .GE. PIP ) EXIT
         IF(.NOT.PEA(L,1)) CYCLE

         PNPC = PIP - PC        
         DO LL = L+1, MAX_PIP

            IF(PNPC .LE. 0) EXIT
            IF(.NOT.PEA(LL,1)) CYCLE

            R_LM = DES_RADIUS(L) + DES_RADIUS(LL)
            R_LM = FACTOR_RLM*R_LM

! the following section adjusts the neighbor check routine when
! any boundary is periodic
! ------------------------------ 
            IF (DES_PERIODIC_WALLS) THEN
               XPOS(:) = DES_POS_NEW(LL,1)
               YPOS(:) = DES_POS_NEW(LL,2)
               II = 1 
               JJ = 1
               KK = 1

               IF(DES_PERIODIC_WALLS_X) THEN
                  IF (DES_POS_NEW(L,1) + R_LM > XLENGTH) THEN 
                     II = 2
                     XPOS(II) = DES_POS_NEW(LL,1) + XLENGTH
                  ELSEIF (DES_POS_NEW(L,1) - R_LM < ZERO) THEN
                     II = 2
                     XPOS(II) = DES_POS_NEW(LL,1) - XLENGTH
                  ENDIF
               ENDIF                     
               IF(DES_PERIODIC_WALLS_Y) THEN 
                  IF (DES_POS_NEW(L,2) + R_LM > YLENGTH) THEN 
                     JJ = 2
                     YPOS(JJ) = DES_POS_NEW(LL,2) + YLENGTH
                  ELSEIF (DES_POS_NEW(L,2) - R_LM < YLENGTH) THEN
                     JJ = 2
                     YPOS(JJ) = DES_POS_NEW(LL,2) - YLENGTH
                  ENDIF
               ENDIF
               IF(DO_K) THEN 
                  ZPOS(:) = DES_POS_NEW(LL,3)
                  IF(DES_PERIODIC_WALLS_Z) THEN 
                     IF (DES_POS_NEW(L,3) + R_LM > ZLENGTH) THEN 
                        KK = 2
                        ZPOS(KK) = DES_POS_NEW(LL,3) + ZLENGTH
                     ELSEIF (DES_POS_NEW(L,3) - R_LM < ZERO) THEN
                        KK = 2
                        ZPOS(KK) = DES_POS_NEW(LL,3) - ZLENGTH
                     ENDIF
                  ENDIF
               ENDIF

! if particle L is within R_LM of a periodic boundary then check
! particles LL current position and its position shifted to the
! opposite boundary for neighbor contact
               OUTER: DO I = 1,II
                  DO J = 1,JJ
                     TMPPOS(1) = XPOS(I)
                     TMPPOS(2) = YPOS(J)
                     IF (DO_K) THEN
                        DO K = 1,KK
                           TMPPOS(3) = ZPOS(K)
                           DISTVEC(:) = TMPPOS(:) - DES_POS_NEW(L,:)
                           DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
                           IF (DIST.LE.R_LM) EXIT OUTER
                        ENDDO
                     ELSE
                        DISTVEC(:) = TMPPOS(:) - DES_POS_NEW(L,:)
                        DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
                        IF (DIST.LE.R_LM) EXIT OUTER
                     ENDIF
                  ENDDO
               ENDDO OUTER

            ELSE   ! if .not.des_periodic_walls
               DISTVEC(:) = DES_POS_NEW(LL,:) - DES_POS_NEW(L,:)
               DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC)) 
            ENDIF    ! endif des_periodic_walls
! ------------------------------ 

            IF (DIST.LE.R_LM) THEN

               NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
               NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1)+1
               NEIGH_L = NEIGHBOURS(L,1) 
               NEIGH_LL = NEIGHBOURS(LL,1)

               IF (NEIGH_L.LE.MN) THEN
                  NEIGHBOURS(L,NEIGH_L+1) = LL
               ELSE
                  WRITE(*,1000) 
                  PRINT *, L, ':', (NEIGHBOURS(L,II), II=1,MN+1) 
                  WRITE(*,1001)
                  STOP
               ENDIF

               IF (NEIGH_LL.LE.MN) THEN
                  NEIGHBOURS(LL,NEIGH_LL+1) = L
               ELSE 
                  WRITE(*,1000) 
                  PRINT *, LL, ':', (NEIGHBOURS(LL,II), II=1,MN+1) 
                  WRITE(*,1001)
                  STOP
               ENDIF

            ENDIF
            PNPC = PNPC - 1
         ENDDO   ! end loop over LL

         PC = PC + 1
      ENDDO   ! end loop over L

 1000 FORMAT(/1X,70('*')//&
         ' From: NSQUARE -',/&
         ' Message: Neighbors GT MN')
 1001 FORMAT(/,1X,70('*')/)

      RETURN
      END SUBROUTINE NSQUARE



