!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QUADTREE(PARTS)                                        C
!  Purpose: To find particle neighbors in 2D using quadtree method     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE QUADTREE(PARTS)

      USE run
      USE param1
      USE constant
      USE discretelement
      IMPLICIT NONE

      INTEGER I, J, K, TC1, TC2, TCR, TCM, PARTS, PQ
      DOUBLE PRECISION CT

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(S_TIME.EQ.TIME) THEN
            INQC = INIT_QUAD_COUNT
         END IF
      END IF

	 IF(INQC.EQ.INIT_QUAD_COUNT) THEN

            CALL INIT_QUAD(PARTICLES)
            
      	    DO I = 1, PARTS
            CALL ADD_QUAD(I)
            END DO

            DO I = 1, PARTS
               J = 1
               CALL QUAD_NEIGHBOURS(I, J)
            END DO

            INQC = INQC - 1
            IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT 

         ELSE

            DO I = 1, PARTS
               PQ = PQUAD(I)
               IF(.NOT.((DES_POS_NEW(I,1).GT.CQUAD(PQ,1)).AND.&
                 (DES_POS_NEW(I,1).LE.CQUAD(PQ,2)).AND.& 
                 (DES_POS_NEW(I,2).GT.CQUAD(PQ,3)).AND.&
	         (DES_POS_NEW(I,2).LE.CQUAD(PQ,4)))) THEN

                 CALL MOVE_PARTICLE(I,PQ)
                 CALL ADD_QUAD(I)

               END IF
            END DO

            DO I = 1, PARTS
               J = 1
               CALL QUAD_NEIGHBOURS(I, J)
            END DO
	
	    INQC = INQC - 1		
            IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT

         END IF
    
      RETURN
      END SUBROUTINE QUADTREE 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_QUAD(PARTS)                                       C
!  Purpose: Initialize arrays for quadtree search                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE INIT_QUAD(PARTS)

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE

      INTEGER I, J, PARTS
      DOUBLE PRECISION A, OMEGA, OOMEGA2, ASINOMEGAT, BOXLIMIT

! Vibrated bottom wall 
      IF(DES_F.NE.ZERO) THEN
         A = ZERO 
         OMEGA = ZERO
         ASINOMEGAT = ZERO
         OMEGA = 2.0D0*Pi*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         A = DES_GAMMA*GRAV(2)*OOMEGA2
      END IF
! End - Bottom wall vibration calculations 
      
      
      DO I = 1, MAXQUADS
         DO J = 1, NMQD 
            LQUAD(I,J) = 0
         END DO
         DO J = 1, NWALLS
            CQUAD(I,J) = ZERO
         END DO
      END DO

      NQUAD = 1
      BOXLIMIT = RADIUS_EQ 
      CQUAD(1,1) = WX1 - BOXLIMIT
      CQUAD(1,2) = BY1 - BOXLIMIT - A !Adding minimum displacement to the bottom wall
      CQUAD(1,3) = EX2 + BOXLIMIT
      CQUAD(1,4) = TY2 + BOXLIMIT


      RETURN
      END SUBROUTINE INIT_QUAD 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MOVE_PARTICLE(II,PQD)                                     C
!  Purpose: To move a particle from its quad                           C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MOVE_PARTICLE(II,PQD)

      USE param1
      USE discretelement
      IMPLICIT NONE
      
      INTEGER II, PQD, J, K

       LQUAD(PQD,7) =  LQUAD(PQD,7) - 1

       DO J = 1, 4 
          IF(LQUAD(PQD,J).EQ.II) THEN
             IF(J.NE.4) THEN
                DO K = J, 3
                   LQUAD(PQD,K) = LQUAD(PQD,K+1) 
                END DO
             ELSE
                LQUAD(PQD,J) = 0
             END IF
          END IF
       END DO

      RETURN
      END SUBROUTINE MOVE_PARTICLE 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_QUAD(II)                                           C
!  Purpose: To find the 'quad' for the particle                        C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADD_QUAD(II)

      USE param1
      USE discretelement
      IMPLICIT NONE
      
      INTEGER II
      INTEGER J ,IQ, IQS, IFOUND, K

      IQS = 1
      IFOUND = 0
      DO WHILE(IFOUND.EQ.0)
         CALL FIND_QUAD(II, IQS)
         J = LQUAD(IQS,7)
         IF(J.EQ.4) THEN
	    CALL SPLIT_QUAD(IQS)
         ELSE IF(J.LT.0) THEN
            PRINT *,'ERROR IN SPLIT_QUAD'
         ELSE 
            LQUAD(IQS,7) = J+1
            LQUAD(IQS,J+1) = II
            IFOUND = 1
         END IF
      END DO            
      RETURN
      END SUBROUTINE ADD_QUAD 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE FIND_QUAD(I,Q)
  
        USE param1
	USE discretelement
        IMPLICIT NONE
        
	INTEGER I, Q, I_F, IC, IX, IY, NC


	I_F = 0
	DO WHILE(I_F.EQ.0) 
	   IF(LQUAD(Q,7).GE.0) THEN
 	     I_F = 1
           ELSE
	     IC = LQUAD(Q,1)
	     IF(DES_POS_NEW(I,1).GE.CQUAD(IC,3)) THEN
		IX = 1
	     ELSE 
		IX = 0
	     END IF
             IF(DES_POS_NEW(I,2).GE.CQUAD(IC,4)) THEN
		IY = 1
	     ELSE
                IY = 0
             END IF
	     NC = IX + 2*IY
	     Q = IC + NC
           END IF
	END DO
        PQUAD(I) = Q

	RETURN
	END SUBROUTINE FIND_QUAD

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C           

	SUBROUTINE SPLIT_QUAD(Q)	
        
        USE param1   
	USE discretelement
        IMPLICIT NONE
        
	INTEGER Q
	INTEGER I_FQ, J, K, ORPHANS(4), NC
	DOUBLE PRECISION XL, YL, XU, YU, XMID, YMID

	LQUAD(Q,7) = -1
	DO J = 1, 4
	   ORPHANS(J) = LQUAD(Q,J)
	   LQUAD(Q,J) = NQUAD + J
	   LQUAD(NQUAD+J,6) = Q
	   LQUAD(NQUAD+J,5) = J
	END DO
	XL = CQUAD(Q,1)
	YL = CQUAD(Q,2)
	XU = CQUAD(Q,3)
	YU = CQUAD(Q,4) 
	XMID = (XU-XL)/2
	YMID = (YU-YL)/2

        CQUAD(NQUAD+1,1) = XL
	CQUAD(NQUAD+1,2) = YL
	CQUAD(NQUAD+1,3) = XL+XMID
	CQUAD(NQUAD+1,4) = YL+YMID

        CQUAD(NQUAD+2,1) = XL+XMID	
	CQUAD(NQUAD+2,2) = YL
	CQUAD(NQUAD+2,3) = XU
	CQUAD(NQUAD+2,4) = YL+YMID

        CQUAD(NQUAD+3,1) = XL
	CQUAD(NQUAD+3,2) = YL+YMID
	CQUAD(NQUAD+3,3) = XL+XMID
	CQUAD(NQUAD+3,4) = YU

        CQUAD(NQUAD+4,1) = XL+XMID
	CQUAD(NQUAD+4,2) = YL+YMID
	CQUAD(NQUAD+4,3) = XU
	CQUAD(NQUAD+4,4) = YU
	NQUAD = NQUAD +4

	DO J = 1, 4
	   I_FQ = Q
	   K = ORPHANS(J)
	   CALL FIND_QUAD(K,I_FQ) 
	   NC = LQUAD(I_FQ,7)
	   LQUAD(I_FQ,7) = NC+1
	   LQUAD(I_FQ,NC+1) = K 
	END DO

!	IF(I_FQ.GT.MAXQUADS-100) THEN 
!           PRINT *,'I_FQ.GT.MAXQUADS. CALLING REARRANGE'
!           CALL REARRANGE_QUADTREE(I_FQ) 
!        END IF

        IF(I_FQ.GT.MAXQUADS-100) THEN
          PRINT *,'QUADTREE: NQUADS GREATER THAN MAXQUADS SPECIFIED', I_FQ
	  STOP
	END IF

	RETURN
	END SUBROUTINE SPLIT_QUAD

          
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DELETE_QUADS(II)                                       C
!  Purpose: To check and delete empty quads                            C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DELETE_QUADS(II)

      USE param1
      USE discretelement
      IMPLICIT NONE
      
      INTEGER II, J

! If all child quads of II are empty then remove the child quads
      IF((LQUAD(LQUAD(II,1),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,2),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,3),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,4),7).EQ.0)) THEN

          DO J = 1, 4
            LQUAD(LQUAD(II,J),6) = -2 
            LQUAD(II,J) = 0
          END DO
          LQUAD(II,7) = 0

      END IF      

      RETURN
      END SUBROUTINE DELETE_QUADS 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REARRANGE_QUADTREE(II)                                 C
!  Purpose: To rearrange the quadtree                                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE REARRANGE_QUADTREE(MNQD)

      USE param1
      USE discretelement
      IMPLICIT NONE
      
      INTEGER MNQD, I, J, K, LQ

      DO I = 1, MNQD
         IF(LQUAD(I,7).EQ.-1) CALL DELETE_QUADS(I)
      END DO

      I = 1
      DO WHILE(I.LE.MNQD)
         IF(LQUAD(I,6).EQ.-2) THEN
           DO J = I, MNQD-4
              LQ = LQUAD(J+4,6)
              DO K = 1,7
                 LQUAD(J,K) = LQUAD(J+4,K) ! moving the lquad array
                 IF(K.LE.4) THEN
                    CQUAD(J,K) = CQUAD(J+4,K) ! moving the cquad array
                    IF(LQUAD(LQ,7).LT.0) LQUAD(LQ,K) = LQUAD(LQ,K) - 4 ! moved quads: parent's children should be changed
                 END IF
              END DO
           END DO
           MNQD = MNQD - 4
           I = I + 4
         ELSE
           I = I + 1
         END IF
      END DO

      RETURN
      END SUBROUTINE REARRANGE_QUADTREE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QUAD_NEIGHBOURS(II, NQ)                                C
!  Purpose: Recursive subroutine to find neighbors for particles       C
!           based on neighbor quads                                    C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      RECURSIVE SUBROUTINE QUAD_NEIGHBOURS(II, NQ)
      
      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER II, I, N, LL, IC, J, JJ, K, POS, IXU, IXL, IYU, IYL, KK, NQ
      DOUBLE PRECISION XM, YM, DIST, R_LM, XL, XU, YL, YU 
      INTEGER DIA_RATIO, cc, cc1

      DIA_RATIO = INT(RADIUS_EQ/DES_RADIUS(II))

      IF (LQUAD(NQ,7).GT.0) THEN

         DO I = 1, LQUAD(NQ,7)
            N = LQUAD(NQ,I)
            IF(N.GT.II) THEN
               DIST = ZERO 
               DO K = 1, DIMN
                  DIST = DIST + (DES_POS_NEW(II,K)-DES_POS_NEW(N,K))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(II)+DES_RADIUS(N)
               IF(DIST.LE.R_LM) THEN
		  NEIGHBOURS(II,1) = NEIGHBOURS(II,1) + 1
		  NEIGHBOURS(N,1) = NEIGHBOURS(N,1) + 1
		  J = NEIGHBOURS(II,1)
		  JJ = NEIGHBOURS(N,1)
                  IF(J.LE.MN) THEN
                     NEIGHBOURS(II,J+1) = N
                     PN_DIST(II,J+1) = DIST
                     PN_RLM(II,J+1) = R_LM
	          ELSE 
                     PRINT *,'QUADTREE - NEIGHBORS GT MN'
                     PRINT *, II,':', (NEIGHBOURS(II,LL),LL=1,MAXNEIGHBORS)
                     STOP
	          END IF
                  IF(JJ.LE.MN) THEN
                     NEIGHBOURS(N,JJ+1) = II 
                     PN_DIST(N,JJ+1) = DIST
                     PN_RLM(N,JJ+1) = R_LM
	          ELSE 
                     PRINT *,'QUADTREE - NEIGHBORS GT MN'
                     PRINT *, II,':', (NEIGHBOURS(N,LL),LL=1,MAXNEIGHBORS)
                     STOP
	          END IF
               END IF
            END IF
         END DO

      ELSE IF(LQUAD(NQ,7).EQ.-1) THEN
         IC = LQUAD(NQ,1)
         XM = CQUAD(IC,3)
         YM = CQUAD(IC,4)

         XL = DES_POS_NEW(II,1) - (DIA_RATIO + 1) * DES_RADIUS(II)
         YL = DES_POS_NEW(II,2) - (DIA_RATIO + 1) * DES_RADIUS(II)
         XU = DES_POS_NEW(II,1) + (DIA_RATIO + 1) * DES_RADIUS(II)
         YU = DES_POS_NEW(II,2) + (DIA_RATIO + 1) * DES_RADIUS(II)

         IF(XL.LT.CQUAD(1,1)) XL=CQUAD(1,1)
         IF(XU.GT.CQUAD(1,3)) XU=CQUAD(1,3)
         IF(YL.LT.CQUAD(1,2)) YL=CQUAD(1,2)
         IF(YU.GT.CQUAD(1,4)) YU=CQUAD(1,4)

         IF(XL.GE.XM) THEN
            IXL = 1
         ELSE 
            IXL = 0
         END IF
         IF(YL.GE.YM) THEN
            IYL = 1
         ELSE 
            IYL = 0
         END IF
         IF(XU.GE.XM) THEN
            IXU = 1
         ELSE
            IXU = 0
         END IF
         IF(YU.GE.YM) THEN
            IYU = 1
         ELSE
            IYU = 0
         END IF
         POS = IXL + 2*IYL + 4*IXU + 8*IYU

         IF(POS.EQ.0) THEN 
            CALL QUAD_NEIGHBOURS(II, IC) 

         ELSE IF(POS.EQ.4) THEN
            CALL QUAD_NEIGHBOURS(II,IC)
            CALL QUAD_NEIGHBOURS(II, IC+1)
            
         ELSE IF(POS.EQ.8) THEN
            CALL QUAD_NEIGHBOURS(II, IC) 
            CALL QUAD_NEIGHBOURS(II, IC+2) 
            
         ELSE IF(POS.EQ.12) THEN
            CALL QUAD_NEIGHBOURS(II, IC)
            CALL QUAD_NEIGHBOURS(II, IC+1) 
            CALL QUAD_NEIGHBOURS(II, IC+2) 
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSE IF(POS.EQ.5) THEN
            CALL QUAD_NEIGHBOURS(II, IC+1)
            
         ELSE IF(POS.EQ.13) THEN
            CALL QUAD_NEIGHBOURS(II, IC+1)
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSE IF(POS.EQ.10) THEN
            CALL QUAD_NEIGHBOURS(II, IC+2)
            
         ELSE IF(POS.EQ.14) THEN
            CALL QUAD_NEIGHBOURS(II, IC+2)
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSE IF(POS.EQ.15) THEN
            CALL QUAD_NEIGHBOURS(II, IC+3) 
            
         END IF
      END IF

      RETURN
      END SUBROUTINE QUAD_NEIGHBOURS 

