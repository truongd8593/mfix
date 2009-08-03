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
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER L, I, K, KK, II
      DOUBLE PRECISION DISTVEC(DIMN), DIST, R_LM

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
      DO L=1, MAX_PIS

         IF(PC .GE. PIS ) EXIT
         IF(.NOT.PEA(L)) CYCLE

         PNPC = PIS - PC        
         DO I = L+1, MAX_PIS

            IF(PNPC .LE. 0) EXIT
            IF(.NOT.PEA(I)) CYCLE

            DISTVEC(:) = DES_POS_NEW(I,:) - DES_POS_NEW(L,:)
            DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
            R_LM = DES_RADIUS(L) + DES_RADIUS(I)
            R_LM = FACTOR_RLM*R_LM

            IF (DIST.LE.R_LM) THEN

               NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
               NEIGHBOURS(I,1) = NEIGHBOURS(I,1)+1
               K = NEIGHBOURS(L,1) 
               KK = NEIGHBOURS(I,1)

               IF (K.LE.MN) THEN
                  NEIGHBOURS(L,K+1) = I
               ELSE 
                  WRITE(*,*) '---------- FROM NSQUARE ---------->'
                  WRITE(*,*) '     NEIGHBORS GT MN'
                  PRINT *, L,':',(NEIGHBOURS(L,II), II=1,MAXNEIGHBORS) 
                  STOP
               ENDIF

               IF (KK.LE.MN) THEN
                  NEIGHBOURS(I,KK+1) = L
               ELSE 
                  WRITE(*,*) '---------- FROM NSQUARE ---------->'
                  WRITE(*,*) '     NEIGHBORS GT MN'
                  PRINT *, I,':',(NEIGHBOURS(I,II), II=1,MAXNEIGHBORS) 
                  STOP
               ENDIF

            ENDIF

            PNPC = PNPC - 1
         ENDDO

         PC = PC + 1
      ENDDO      

      RETURN
      END SUBROUTINE NSQUARE


