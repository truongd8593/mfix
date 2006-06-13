!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OCTREE(PARTS)                                          C
!  Purpose: To find neighbors in 3D using octree search method         C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OCTREE(PARTS)

      USE run
      USE param1
      USE constant
      USE discretelement
      IMPLICIT NONE

      INTEGER I, J, TC1, TC2, TCR, TCM, PARTS
      INTEGER PQ, PQM(10), PQXL, PQXU, PQYL, PQYU, PQZL, PQZU

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(S_TIME.EQ.TIME) THEN
            INQC = INIT_QUAD_COUNT
         END IF
      END IF

         IF(INQC.EQ.INIT_QUAD_COUNT) THEN
     
                  CALL INIT_OCT(PARTICLES)
     
                  DO I = 1, PARTS
                     PQ = 1
                     CALL ADD_OCT(I,PQ)
                  END DO
     
                  DO I = 1, PARTS
                     J = 1
                     CALL OCT_NEIGHBOURS(I, J)
                  END DO
     
                  INQC = INQC - 1
                  IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT
     
               ELSE

                 DO I = 1, PARTS
                    J = 1
                    PQ = PQUAD(I)
! Removing a moved particle from its old quad
                     IF(.NOT.((DES_POS_NEW(I,1).GE.CQUAD(PQ,1)).AND.&
                       (DES_POS_NEW(I,1).LT.CQUAD(PQ,3)).AND.&
                       (DES_POS_NEW(I,2).GE.CQUAD(PQ,2)).AND.&
                       (DES_POS_NEW(I,2).LT.CQUAD(PQ,4)).AND.&
                       (DES_POS_NEW(I,3).GE.CQUAD(PQ,5)).AND.&
                       (DES_POS_NEW(I,3).LT.CQUAD(PQ,6)))) THEN
                     CALL REMOVE_PARTICLE_OCT(I,PQ)
      ! Placing the particle in a new quad
                  IF(QLM.LE.2) THEN
                     PQ = 1
                     GO TO 40
                  ELSE
                     PQM(1) = PQ
                     IF(PQM(1).EQ.1) THEN
                        GO TO 40
                     ELSE
                        DO J = 2, QLM
                           PQM(J) = LQUAD(PQM(J-1),10)
                           IF(PQM(J).EQ.1) THEN
                              PQ = 1
                              GO TO 40
                           ELSE
                              PQ = PQM(J)
                              IF((DES_POS_NEW(I,1).GE.CQUAD(PQ,1)).AND.&
                                 (DES_POS_NEW(I,1).LT.CQUAD(PQ,3)).AND.&
                                 (DES_POS_NEW(I,2).GE.CQUAD(PQ,2)).AND.&
                                 (DES_POS_NEW(I,2).GE.CQUAD(PQ,4)).AND.&
                                 (DES_POS_NEW(I,3).GE.CQUAD(PQ,5)).AND.&
                                 (DES_POS_NEW(I,3).LT.CQUAD(PQ,6))) THEN
                                 PQ = LQUAD(PQM(J),10)
                                 GO TO 40
                              END IF
                           END IF
                        END DO
                        IF(J.GT.QLM) PQ = 1
                     END IF
                  END IF
 40   CONTINUE
                     CALL ADD_OCT(I,PQ)
                     END IF
                  END DO

! Neighbor search
                  DO I = 1, PARTS
                     J = 1
                   IF(QLN.LE.2) THEN
                     PQ = 1
                     GO TO 50
                   ELSE
                     PQM(1) = PQUAD(I)
                     IF(PQM(1).EQ.1) THEN
                        GO TO 50
                    ELSE
                       DO J = 2, QLN
                          PQM(J) = LQUAD(PQM(J-1),6)
                          IF(PQM(J).EQ.1) THEN
                             PQ = 1
                             GO TO 50
                          ELSE
                             PQ = PQM(J)
                             PQXL = CQUAD(PQ,1) + 2.1D0 * DES_RADIUS(I)
                             PQXU = CQUAD(PQ,3) - 2.1D0 * DES_RADIUS(I)
                             PQYL = CQUAD(PQ,2) + 2.1D0 * DES_RADIUS(I)
                             PQYU = CQUAD(PQ,4) - 2.1D0 * DES_RADIUS(I)
                             PQZL = CQUAD(PQ,5) + 2.1D0 * DES_RADIUS(I)
                             PQZU = CQUAD(PQ,6) - 2.1D0 * DES_RADIUS(I)
    
                             IF((DES_POS_NEW(I,1).GE.PQXL).AND.&
                                (DES_POS_NEW(I,1).LT.PQXU).AND.&
                                (DES_POS_NEW(I,2).GE.PQYL).AND.&
                                (DES_POS_NEW(I,2).GE.PQYU).AND.&
                                (DES_POS_NEW(I,3).GE.PQZL).AND.&
                                (DES_POS_NEW(I,3).LT.PQZU)) THEN
                                GO TO 50
                             END IF
                          END IF
                       END DO
                       IF(J.GT.QLN) PQ = 1
                    END IF
                  END IF
 50   CONTINUE
                    CALL OCT_NEIGHBOURS(I, PQ)
                 END DO
    
                 INQC = INQC - 1
                 IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT
    
              END IF


      RETURN
      END SUBROUTINE OCTREE 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_OCT(PARTS)                                        C
!  Purpose: Initializing arrays for octree search                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE INIT_OCT(PARTS)

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

      A = ZERO
      OMEGA = ZERO
      IF(DES_F.NE.ZERO) THEN
         OMEGA = 2.0D0*Pi*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         IF(UNITS == "CGS") THEN
            A = DES_GAMMA*GRAV(2)*OOMEGA2
         ELSE
            A = DES_GAMMA*GRAV(2)*OOMEGA2
         END IF         
      END IF

      ASINOMEGAT = A*SIN(OMEGA*S_TIME)
      BOXLIMIT = RADIUS_EQ

      DO I = 1, MAXQUADS
         DO J = 1, NMQD 
            LQUAD(I,J) = 0
         END DO
         DO J = 1, NWALLS 
            CQUAD(I,J) = ZERO
         END DO
      END DO

      NQUAD = 1

      CQUAD(1,1) = WX1 - BOXLIMIT
      CQUAD(1,2) = BY1- BOXLIMIT + ASINOMEGAT 
      CQUAD(1,3) = EX2 + BOXLIMIT 
      CQUAD(1,4) = TY2 + BOXLIMIT 
      CQUAD(1,5) = SZ1 - BOXLIMIT
      CQUAD(1,6) = NZ2 + BOXLIMIT 

      RETURN
      END SUBROUTINE INIT_OCT 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REMOVE_PARTICLE_OCT(II,PQD)                                  C
!  Purpose: To move a particle from its quad                           C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
           SUBROUTINE REMOVE_PARTICLE_OCT(II,PQD)
    
           USE param1
           USE discretelement
           IMPLICIT NONE
    
           INTEGER II, PQD, J, K
    
           LQUAD(PQD,11) =  LQUAD(PQD,11) - 1
    
            DO J = 1, 8 
               IF(LQUAD(PQD,J).EQ.II) THEN
                  IF(J.NE.8) THEN
                     DO K = J, 7 
                        LQUAD(PQD,K) = LQUAD(PQD,K+1)
                     END DO
                  ELSE
                     LQUAD(PQD,J) = 0
                  END IF
               END IF
            END DO
    
           RETURN
           END SUBROUTINE REMOVE_PARTICLE_OCT
   

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_OCT(II, IQS)                                            C
!  Purpose: To find the 'oct' for the particle                         C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADD_OCT(II, IQS)
      USE param1
      USE discretelement

      IMPLICIT NONE
      
      INTEGER II, IQS
      INTEGER J ,IQ, IFOUND, K


      IFOUND = 0
 10   CONTINUE
         CALL FIND_OCT(II, IQS)
         J = LQUAD(IQS,11)
         IF(J.EQ.8) THEN
	    CALL SPLIT_OCT(IQS)
         ELSE IF(J.LT.0) THEN
            PRINT *,'ERROR IN SPLIT_OCT'
         ELSE 
            LQUAD(IQS,11) = J+1
            LQUAD(IQS,J+1) = II
            IFOUND = 1
         END IF
      IF(IFOUND.EQ.0) GO TO 10
      RETURN
      END SUBROUTINE ADD_OCT 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE FIND_OCT(I,Q)
        USE param1
	USE discretelement
	INTEGER I, Q, I_F, IC, IX, IY, IZ

	I_F = 0
  20     CONTINUE
	   IF(LQUAD(Q,11).GE.0) THEN
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
	     IF(DES_POS_NEW(I,3).GE.CQUAD(IC,6)) THEN
		IZ = 1
	     ELSE
		IZ = 0
	     END IF
	     NC = IX + 2*IY + 4*IZ
	     Q = IC + NC
           END IF
        IF(I_F.EQ.0) GO TO 20
        PQUAD(I) = Q

	RETURN
	END SUBROUTINE FIND_OCT

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE SPLIT_OCT(Q)	
        USE param1
	USE discretelement
	INTEGER Q
	INTEGER I_FQ, J, K, ORPHANS(8), NC
	DOUBLE PRECISION XL, YL, XU, YU, ZL, ZU, XMID, YMID, ZMID
	LQUAD(Q,11) = -1
	DO J = 1, 8 
	   ORPHANS(J) = LQUAD(Q,J)
	   LQUAD(Q,J) = NQUAD + J
	   LQUAD(NQUAD+J,10) = Q
	   LQUAD(NQUAD+J,9) = J
	END DO
	XL = CQUAD(Q,1)
	YL = CQUAD(Q,2)
	XU = CQUAD(Q,3)
	YU = CQUAD(Q,4) 
	ZL = CQUAD(Q,5)
	ZU = CQUAD(Q,6)
		
	XMID = (XU-XL)/2
	YMID = (YU-YL)/2
	ZMID = (ZU-ZL)/2

	CQUAD(NQUAD+1,1) = XL
	CQUAD(NQUAD+1,2) = YL
	CQUAD(NQUAD+1,3) = XL+XMID
	CQUAD(NQUAD+1,4) = YL+YMID
	CQUAD(NQUAD+1,5) = ZL
	CQUAD(NQUAD+1,6) = ZL+ZMID

	CQUAD(NQUAD+2,1) = XL+XMID	
	CQUAD(NQUAD+2,2) = YL
	CQUAD(NQUAD+2,3) = XU
	CQUAD(NQUAD+2,4) = YL+YMID
	CQUAD(NQUAD+2,5) = ZL
	CQUAD(NQUAD+2,6) = ZL+ZMID

	CQUAD(NQUAD+3,1) = XL
	CQUAD(NQUAD+3,2) = YL+YMID
	CQUAD(NQUAD+3,3) = XL+XMID
	CQUAD(NQUAD+3,4) = YU
	CQUAD(NQUAD+3,5) = ZL
	CQUAD(NQUAD+3,6) = ZL+ZMID

	CQUAD(NQUAD+4,1) = XL+XMID
	CQUAD(NQUAD+4,2) = YL+YMID
	CQUAD(NQUAD+4,3) = XU
	CQUAD(NQUAD+4,4) = YU
	CQUAD(NQUAD+4,5) = ZL
	CQUAD(NQUAD+4,6) = ZL+ZMID

	CQUAD(NQUAD+5,1) = XL
	CQUAD(NQUAD+5,2) = YL
	CQUAD(NQUAD+5,3) = XL+XMID
	CQUAD(NQUAD+5,4) = YL+YMID
	CQUAD(NQUAD+5,5) = ZL+ZMID
	CQUAD(NQUAD+5,6) = ZU

	CQUAD(NQUAD+6,1) = XL+XMID	
	CQUAD(NQUAD+6,2) = YL
	CQUAD(NQUAD+6,3) = XU
	CQUAD(NQUAD+6,4) = YL+YMID
	CQUAD(NQUAD+6,5) = ZL+ZMID
	CQUAD(NQUAD+6,6) = ZU 

	CQUAD(NQUAD+7,1) = XL
	CQUAD(NQUAD+7,2) = YL+YMID
	CQUAD(NQUAD+7,3) = XL+XMID
	CQUAD(NQUAD+7,4) = YU
	CQUAD(NQUAD+7,5) = ZL+ZMID
	CQUAD(NQUAD+7,6) = ZU

	CQUAD(NQUAD+8,1) = XL+XMID
	CQUAD(NQUAD+8,2) = YL+YMID
	CQUAD(NQUAD+8,3) = XU
	CQUAD(NQUAD+8,4) = YU
	CQUAD(NQUAD+8,5) = ZL+ZMID
	CQUAD(NQUAD+8,6) = ZU
	
	NQUAD = NQUAD + 8
	DO J = 1, 8
	   I_FQ = Q
	   K = ORPHANS(J)
	   CALL FIND_OCT(K,I_FQ) 
	   NC = LQUAD(I_FQ,11)
	   LQUAD(I_FQ,11) = NC+1
	   LQUAD(I_FQ,NC+1) = K 
	END DO

        IF(I_FQ.GT.MAXQUADS-100) THEN
           PRINT *,'I_FQ.GT.MAXQUADS. CALLING REARRANGE'
           CALL REARRANGE_OCTREE(I_FQ)
        END IF
   
        IF(I_FQ.GT.MAXQUADS-100) THEN
           PRINT *,'OCTREE: NQUADS GREATER THAN MAXQUADS SPECIFIED', I_FQ
           STOP
        END IF

	RETURN
	END SUBROUTINE SPLIT_OCT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DELETE_OCTS(II)                                       C
!  Purpose: To check and delete empty quads                            C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
           SUBROUTINE DELETE_OCTS(II)
    
           USE param1
           USE discretelement
           IMPLICIT NONE
    
           INTEGER II, J
    
! If all child quads of II are empty then remove the child quads
           IF((LQUAD(LQUAD(II,1),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,2),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,3),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,4),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,5),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,6),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,7),11).EQ.0).AND.&
              (LQUAD(LQUAD(II,8),11).EQ.0)) THEN
    
               DO J = 1,8 
                 LQUAD(LQUAD(II,J),10) = -2
                 LQUAD(II,J) = 0
               END DO
               LQUAD(II,11) = 0
    
           END IF
    
           RETURN
           END SUBROUTINE DELETE_OCTS


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REARRANGE_OCTREE(II)                                 C
!  Purpose: To rearrange the quadtree                                  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
           SUBROUTINE REARRANGE_OCTREE(MNQD)
    
           USE param1
           USE discretelement
           IMPLICIT NONE
    
           INTEGER MNQD, I, J, K, LQ
    
           DO I = 1, MNQD
              IF(LQUAD(I,11).EQ.-1) CALL DELETE_OCTS(I)
           END DO
    
           I = 1
      30   CONTINUE
              IF(LQUAD(I,10).EQ.-2) THEN
                DO J = I, MNQD-8
                   LQ = LQUAD(J+8,10)
                   DO K = 1,11
                      LQUAD(J,K) = LQUAD(J+8,K) ! moving the lquad array
                      IF(K.LE.8) THEN
                         CQUAD(J,K) = CQUAD(J+8,K) ! moving the cquad array
                         IF(LQUAD(LQ,11).LT.0) LQUAD(LQ,K) = LQUAD(LQ,K) - 8 
                      END IF
                   END DO
                END DO
                MNQD = MNQD - 8 
                I = I + 8 
              ELSE
                I = I + 1
              END IF
           IF(I.LE.MNQD) GO TO 30
    
           RETURN
           END SUBROUTINE REARRANGE_OCTREE
       
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OCT_NEIGHBOURS(II, NQ)                                 C
!  Purpose: Recursive sub-routine to find neighbors for particles      C
!           based on neighbor octs                                     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      RECURSIVE SUBROUTINE OCT_NEIGHBOURS(II, NQ)

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER II, I, N, IC, J, JJ, K, POS, IXU, IXL, IYU, IYL, IZL, IZU
      DOUBLE PRECISION XM, YM, ZM, DIST, R_LM, XL, XU, YL, YU, ZL, ZU 
      INTEGER LL, NQ  

      IF(LQUAD(NQ,11).GT.0) THEN

         DO I = 1, LQUAD(NQ,11)
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
                     PRINT *,'QUAD - NEIGHBORS GT MN'
                     PRINT *, II, (NEIGHBOURS(II,LL),LL=1,MAXNEIGHBORS)
                     STOP
                  END IF
                  IF(JJ.LE.MN) THEN
                     NEIGHBOURS(II,JJ+1) = II 
                     PN_DIST(N,JJ+1) = DIST
                     PN_RLM(N,JJ+1) = R_LM
                  ELSE
                     PRINT *,'QUAD - NEIGHBORS GT MN'
                     PRINT *, II, (NEIGHBOURS(N,LL),LL=1,MAXNEIGHBORS)
                     STOP
                  END IF
               END IF
	    END IF
         END DO

      ELSE IF(LQUAD(NQ,11).EQ.-1) THEN
         IC = LQUAD(NQ,1)
         XM = CQUAD(IC,3)
         YM = CQUAD(IC,4)
         ZM = CQUAD(IC,6)
         XL = DES_POS_NEW(II,1) - 2.1D0 * DES_RADIUS(II)
         XU = DES_POS_NEW(II,1) + 2.1D0 * DES_RADIUS(II)
         YL = DES_POS_NEW(II,2) - 2.1D0 * DES_RADIUS(II)
         YU = DES_POS_NEW(II,2) + 2.1D0 * DES_RADIUS(II)
         ZL = DES_POS_NEW(II,3) - 2.1D0 * DES_RADIUS(II)
         ZU = DES_POS_NEW(II,3) + 2.1D0 * DES_RADIUS(II)

         IF(XL.LT.CQUAD(1,1)) XL=CQUAD(1,1) 		
         IF(XU.GT.CQUAD(1,3)) XU=CQUAD(1,3)
         IF(YL.LT.CQUAD(1,2)) YL=CQUAD(1,2)
         IF(YU.GT.CQUAD(1,4)) YU=CQUAD(1,4)   
         IF(ZL.LT.CQUAD(1,5)) ZL=CQUAD(1,5)
         IF(ZU.GT.cquad(1,6)) ZU=CQUAD(1,6)  

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
         IF(ZL.GE.ZM) THEN
            IZL = 1
         ELSE 
            IZL = 0
         END IF
         IF(ZU.GE.ZM) THEN
            IZU = 1
         ELSE
            IZU =0
         END IF
         POS = IXL + 2*IYL + 4*IZL + 8*IXU + 16*IYU + 32*IZU

         IF(POS.EQ.0) THEN 
            CALL OCT_NEIGHBOURS(II, IC)
         ELSE IF(POS.EQ.9) THEN
            CALL OCT_NEIGHBOURS(II, IC+1)
         ELSE IF(POS.EQ.18) THEN
            CALL OCT_NEIGHBOURS(II, IC+2)
         ELSE IF(POS.EQ.27) THEN
            CALL OCT_NEIGHBOURS(II, IC+3) 
         ELSE IF(POS.EQ.36) THEN 
            CALL OCT_NEIGHBOURS(II, IC+4)
         ELSE IF(POS.EQ.45) THEN
            CALL OCT_NEIGHBOURS(II, IC+5)
         ELSE IF(POS.EQ.54) THEN
            CALL OCT_NEIGHBOURS(II, IC+6)
         ELSE IF(POS.EQ.63) THEN
            CALL OCT_NEIGHBOURS(II, IC+7) 

         ELSE IF(POS.EQ.8) THEN
            CALL OCT_NEIGHBOURS(II,IC)
            CALL OCT_NEIGHBOURS(II, IC+1)      
         ELSE IF(POS.EQ.16) THEN
            CALL OCT_NEIGHBOURS(II, IC) 
            CALL OCT_NEIGHBOURS(II, IC+2) 
         ELSE IF(POS.EQ.32) THEN
            CALL OCT_NEIGHBOURS(II, IC) 
            CALL OCT_NEIGHBOURS(II, IC+4)

         ELSE IF(POS.EQ.25) THEN
            CALL OCT_NEIGHBOURS(II, IC+1)
            CALL OCT_NEIGHBOURS(II, IC+3)
         ELSE IF(POS.EQ.41) THEN
            CALL OCT_NEIGHBOURS(II, IC+1)
            CALL OCT_NEIGHBOURS(II, IC+5)
         ELSE IF(POS.EQ.26) THEN
            CALL OCT_NEIGHBOURS(II, IC+2)
            CALL OCT_NEIGHBOURS(II, IC+3)
         ELSE IF(POS.EQ.50) THEN
            CALL OCT_NEIGHBOURS(II, IC+2)
            CALL OCT_NEIGHBOURS(II, IC+6)
         ELSE IF(POS.EQ.59) THEN
            CALL OCT_NEIGHBOURS(II, IC+3)
            CALL OCT_NEIGHBOURS(II, IC+7)
         ELSE IF(POS.EQ.52) THEN
            CALL OCT_NEIGHBOURS(II, IC+4)
            CALL OCT_NEIGHBOURS(II, IC+5)
         ELSE IF(POS.EQ.44) THEN
            CALL OCT_NEIGHBOURS(II, IC+4)
            CALL OCT_NEIGHBOURS(II, IC+6)
         ELSE IF(POS.EQ.61) THEN
            CALL OCT_NEIGHBOURS(II, IC+5)
            CALL OCT_NEIGHBOURS(II, IC+7)
         ELSE IF(POS.EQ.62) THEN
            CALL OCT_NEIGHBOURS(II, IC+6)
            CALL OCT_NEIGHBOURS(II, IC+7)

         ELSE IF(POS.EQ.24) THEN
            CALL OCT_NEIGHBOURS(II, IC)
            CALL OCT_NEIGHBOURS(II, IC+1) 
            CALL OCT_NEIGHBOURS(II, IC+2) 
            CALL OCT_NEIGHBOURS(II, IC+3)
         ELSE IF(POS.EQ.60) THEN
            CALL OCT_NEIGHBOURS(II, IC+4)
            CALL OCT_NEIGHBOURS(II, IC+5) 
            CALL OCT_NEIGHBOURS(II, IC+6) 
            CALL OCT_NEIGHBOURS(II, IC+7)
         ELSE IF(POS.EQ.48) THEN
            CALL OCT_NEIGHBOURS(II, IC)
            CALL OCT_NEIGHBOURS(II, IC+2) 
            CALL OCT_NEIGHBOURS(II, IC+4) 
            CALL OCT_NEIGHBOURS(II, IC+6)
         ELSE IF(POS.EQ.57) THEN
            CALL OCT_NEIGHBOURS(II, IC+1)
            CALL OCT_NEIGHBOURS(II, IC+3) 
            CALL OCT_NEIGHBOURS(II, IC+5) 
            CALL OCT_NEIGHBOURS(II, IC+7)
         ELSE IF(POS.EQ.58) THEN
            CALL OCT_NEIGHBOURS(II, IC+2)
            CALL OCT_NEIGHBOURS(II, IC+3) 
            CALL OCT_NEIGHBOURS(II, IC+6) 
            CALL OCT_NEIGHBOURS(II, IC+7)
         ELSE IF(POS.EQ.40) THEN
            CALL OCT_NEIGHBOURS(II, IC)
            CALL OCT_NEIGHBOURS(II, IC+1) 
            CALL OCT_NEIGHBOURS(II, IC+4) 
            CALL OCT_NEIGHBOURS(II, IC+5)
            
         ELSE IF(POS.EQ.56) THEN
            CALL OCT_NEIGHBOURS(II, IC)
            CALL OCT_NEIGHBOURS(II, IC+1) 
            CALL OCT_NEIGHBOURS(II, IC+2) 
            CALL OCT_NEIGHBOURS(II, IC+3)
            CALL OCT_NEIGHBOURS(II, IC+4)
            CALL OCT_NEIGHBOURS(II, IC+5) 
            CALL OCT_NEIGHBOURS(II, IC+6) 
            CALL OCT_NEIGHBOURS(II, IC+7)          
         END IF
      END IF

      RETURN
      END SUBROUTINE OCT_NEIGHBOURS 

