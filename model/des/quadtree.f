!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QUADTREE                                               C
!>  Purpose: To find particle neighbors in 2D using quadtree method     
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C  
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE QUADTREE

      USE run
      USE param1
      USE constant
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------  
      INTEGER I, J, K
      INTEGER PQ, PQM(10), PQXL, PQXU, PQYL, PQYU
      DOUBLE PRECISION CT
!----------------------------------------------- 

      IF(DES_CONTINUUM_COUPLED) THEN
         IF(S_TIME.EQ.TIME) THEN
            INQC = INIT_QUAD_COUNT
         ENDIF
      ENDIF

      IF(INQC.EQ.INIT_QUAD_COUNT) THEN

         CALL INIT_QUAD
         DO I = 1, PARTICLES
            PQ = 1
            CALL ADD_QUAD(I,PQ)
         ENDDO
         DO I = 1, PARTICLES
            J = 1
            CALL QUAD_NEIGHBOURS(I, J)
         ENDDO
         INQC = INQC - 1
         IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT 

      ELSE

         DO I = 1, PARTICLES
            J = 1
            PQ = PQUAD(I)
! Removing a moved particle from its old quad 
            IF(.NOT.((DES_POS_NEW(I,1).GE.CQUAD(PQ,1)).AND.&
            (DES_POS_NEW(I,1).LT.CQUAD(PQ,3)).AND.& 
            (DES_POS_NEW(I,2).GE.CQUAD(PQ,2)).AND.&
            (DES_POS_NEW(I,2).LT.CQUAD(PQ,4)))) THEN
               CALL REMOVE_PARTICLE(I,PQ)
! Placing the particle in a new quad
               IF(QLM.LE.2) THEN
                  PQ = 1
                  GOTO 40
               ELSE
                  PQM(1) = PQ
                  IF(PQM(1).EQ.1) THEN
                     GOTO 40
                  ELSE
                     DO J = 2, QLM 
                        PQM(J) = LQUAD(PQM(J-1),6)
                        IF(PQM(J).EQ.1) THEN
                           PQ = 1
                           GOTO 40
                        ELSE
                           PQ = PQM(J) 
                           IF((DES_POS_NEW(I,1).GE.CQUAD(PQ,1)).AND.&
                           (DES_POS_NEW(I,1).LT.CQUAD(PQ,3)).AND.&
                           (DES_POS_NEW(I,2).GE.CQUAD(PQ,2)).AND.&
                           (DES_POS_NEW(I,2).LT.CQUAD(PQ,4))) THEN
                              PQ = LQUAD(PQM(J),6)
                              GOTO 40    
                           ENDIF
                        ENDIF
                     ENDDO
                     IF(J.GT.QLM) PQ = 1
                  ENDIF
               ENDIF
 40    CONTINUE
               CALL ADD_QUAD(I,PQ)
            ENDIF
         ENDDO

! Neighbor search
         DO I = 1, PARTICLES
            J = 1
            IF(QLN.LE.2) THEN
               PQ = 1
               GO TO 50
            ELSE
               PQM(1) = PQUAD(I)
               IF(PQM(1).EQ.1) THEN
                  GOTO 50
               ELSE
                  DO J = 2, QLN  
                     PQM(J) = LQUAD(PQM(J-1),6)
                     IF(PQM(J).EQ.1) THEN
                        PQ = 1
                        GOTO 50
                     ELSE
                        PQ = PQM(J) 
                        PQXL = CQUAD(PQ,1) + 2.1D0 * DES_RADIUS(I)
                        PQXU = CQUAD(PQ,3) - 2.1D0 * DES_RADIUS(I)
                        PQYL = CQUAD(PQ,2) + 2.1D0 * DES_RADIUS(I)
                        PQYU = CQUAD(PQ,4) - 2.1D0 * DES_RADIUS(I)
   
                        IF((DES_POS_NEW(I,1).GE.PQXL).AND.&
                        (DES_POS_NEW(I,1).LT.PQXU).AND.&
                        (DES_POS_NEW(I,2).GE.PQYL).AND.&
                        (DES_POS_NEW(I,2).LT.PQYU)) THEN
                           GOTO 50    
                        ENDIF
                     ENDIF
                  ENDDO
                  IF(J.GT.QLN) PQ = 1
               ENDIF
            ENDIF
 50   CONTINUE
            CALL QUAD_NEIGHBOURS(I, PQ)
         ENDDO
      
         INQC = INQC - 1
         IF(INQC.EQ.0) INQC = INIT_QUAD_COUNT

      ENDIF
      
      RETURN
      END SUBROUTINE QUADTREE 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_QUAD                                              C
!  Purpose: Initialize arrays for quadtree search                      C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE INIT_QUAD

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
!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER I, J
      DOUBLE PRECISION A, OMEGA, OOMEGA2, BOXLIMIT
!-----------------------------------------------

!     Vibrated bottom wall 
      A = ZERO 
      IF(DES_F.NE.ZERO) THEN
         OMEGA = ZERO
         OMEGA = 2.0D0*Pi*DES_F
         OOMEGA2 = ONE/(OMEGA**2)
         A = DES_GAMMA*GRAV(2)*OOMEGA2
      ENDIF
!     End - Bottom wall vibration calculations 
      
      DO I = 1, MAXQUADS
         DO J = 1, NMQD 
            LQUAD(I,J) = 0
         ENDDO
         DO J = 1, NWALLS
            CQUAD(I,J) = ZERO
         ENDDO
      ENDDO

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
!  Module name: REMOVE_PARTICLE(II,PQD)                                C
!  Purpose: To move a particle from its quad                           C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE REMOVE_PARTICLE(II,PQD)

      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER II, PQD, J, K
!----------------------------------------------- 

      LQUAD(PQD,7) =  LQUAD(PQD,7) - 1

      DO J = 1, 4 
         IF(LQUAD(PQD,J).EQ.II) THEN
            IF(J.NE.4) THEN
               DO K = J, 3
                  LQUAD(PQD,K) = LQUAD(PQD,K+1) 
               END DO
            ELSE
               LQUAD(PQD,J) = 0
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE REMOVE_PARTICLE 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_QUAD(II,IQS)                                       C
!  Purpose: To find the 'quad' for the particle                        C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ADD_QUAD(II,IQS)

      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------   
      INTEGER II, IQS
      INTEGER J ,IQ, IFOUND, K
!-----------------------------------------------

      IFOUND = 0
 10   CONTINUE
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
      ENDIF
      IF(IFOUND.EQ.0) GOTO 10 

      RETURN
      END SUBROUTINE ADD_QUAD 


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE FIND_QUAD(I,Q)
      
      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------        
      INTEGER I, Q, I_F, IC, IX, IY, NC
!-----------------------------------------------  

      I_F = 0
 20   CONTINUE
      IF(LQUAD(Q,7).GE.0) THEN
         I_F = 1
      ELSE
         IC = LQUAD(Q,1)
         IF(DES_POS_NEW(I,1).GE.CQUAD(IC,3)) THEN
            IX = 1
         ELSE 
            IX = 0
         ENDIF
         IF(DES_POS_NEW(I,2).GE.CQUAD(IC,4)) THEN
            IY = 1
         ELSE
            IY = 0
         ENDIF
         NC = IX + 2*IY
         Q = IC + NC
      ENDIF
      IF(I_F.EQ.0) GOTO 20 
      PQUAD(I) = Q

      RETURN
      END SUBROUTINE FIND_QUAD


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C           

      SUBROUTINE SPLIT_QUAD(Q)
      
      USE param1   
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER Q
      INTEGER I_FQ, J, K, ORPHANS(4), NC
      DOUBLE PRECISION XL, YL, XU, YU, XMID, YMID
!----------------------------------------------- 

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
      ENDDO

      IF(I_FQ.GT.MAXQUADS-100) THEN 
         PRINT *,'I_FQ.GT.MAXQUADS. CALLING REARRANGE'
         CALL REARRANGE_QUADTREE(I_FQ) 
      ENDIF

      IF(I_FQ.GT.MAXQUADS-100) THEN
         PRINT *,'QUADTREE: NQUADS GREATER THAN MAXQUADS SPECIFIED', I_FQ
         STOP
      ENDIF

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
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER II, J
!----------------------------------------------- 

! If all child quads of II are empty then remove the child quads
      IF((LQUAD(LQUAD(II,1),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,2),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,3),7).EQ.0).AND.&
         (LQUAD(LQUAD(II,4),7).EQ.0)) THEN

         DO J = 1, 4
            LQUAD(LQUAD(II,J),6) = -2 
            LQUAD(II,J) = 0
         ENDDO
         LQUAD(II,7) = 0
   
      ENDIF      

      RETURN
      END SUBROUTINE DELETE_QUADS 


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
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER MNQD, I, J, K, LQ
!-----------------------------------------------

      DO I = 1, MNQD
         IF(LQUAD(I,7).EQ.-1) CALL DELETE_QUADS(I)
      ENDDO

      I = 1
 30   CONTINUE
      IF(LQUAD(I,6).EQ.-2) THEN
         DO J = I, MNQD-4
            LQ = LQUAD(J+4,6)
            DO K = 1,7
               LQUAD(J,K) = LQUAD(J+4,K) ! moving the lquad array
               IF(K.LE.4) THEN
                  CQUAD(J,K) = CQUAD(J+4,K) ! moving the cquad array
                  IF(LQUAD(LQ,7).LT.0) LQUAD(LQ,K) = LQUAD(LQ,K) - 4 
               ENDIF
            ENDDO
         ENDDO
         MNQD = MNQD - 4
         I = I + 4
      ELSE
         I = I + 1
      ENDIF
      IF(I.LE.MNQD) GO TO 30 

      RETURN
      END SUBROUTINE REARRANGE_QUADTREE


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
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER II, I, N, LL, IC, J, JJ, K, POS, IXU, IXL, IYU, IYL, KK, NQ
      DOUBLE PRECISION XM, YM, DIST, R_LM, XL, XU, YL, YU 
!-----------------------------------------------

      IF (LQUAD(NQ,7).GT.0) THEN

         DO I = 1, LQUAD(NQ,7)
            N = LQUAD(NQ,I)
            IF(N.GT.II) THEN
               DIST = ZERO 
               DO K = 1, DIMN
                  DIST = DIST + (DES_POS_NEW(II,K)-DES_POS_NEW(N,K))**2
               ENDDO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(II)+DES_RADIUS(N)
               IF(DIST.LE.R_LM) THEN
                  NEIGHBOURS(II,1) = NEIGHBOURS(II,1) + 1
                  NEIGHBOURS(N,1) = NEIGHBOURS(N,1) + 1
                  J = NEIGHBOURS(II,1)
                  JJ = NEIGHBOURS(N,1)
                  IF(J.LE.MN) THEN
                     NEIGHBOURS(II,J+1) = N
                  ELSE 
                     PRINT *,'QUADTREE - NEIGHBORS GT MN'
                     PRINT *, II,':', (NEIGHBOURS(II,LL),LL=1,MAXNEIGHBORS)
                     STOP
                  ENDIF
                  IF(JJ.LE.MN) THEN
                     NEIGHBOURS(N,JJ+1) = II 
                  ELSE 
                     PRINT *,'QUADTREE - NEIGHBORS GT MN'
                     PRINT *, II,':', (NEIGHBOURS(N,LL),LL=1,MAXNEIGHBORS)
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

      ELSEIF(LQUAD(NQ,7).EQ.-1) THEN
         IC = LQUAD(NQ,1)
         XM = CQUAD(IC,3)
         YM = CQUAD(IC,4)

         XL = DES_POS_NEW(II,1) - 2.1D0 * DES_RADIUS(II)
         YL = DES_POS_NEW(II,2) - 2.1D0 * DES_RADIUS(II)
         XU = DES_POS_NEW(II,1) + 2.1D0 * DES_RADIUS(II)
         YU = DES_POS_NEW(II,2) + 2.1D0 * DES_RADIUS(II)

         IF(XL.LT.CQUAD(1,1)) XL=CQUAD(1,1)
         IF(XU.GT.CQUAD(1,3)) XU=CQUAD(1,3)
         IF(YL.LT.CQUAD(1,2)) YL=CQUAD(1,2)
         IF(YU.GT.CQUAD(1,4)) YU=CQUAD(1,4)

         IF(XL.GE.XM) THEN
            IXL = 1
         ELSE 
            IXL = 0
         ENDIF
         IF(YL.GE.YM) THEN
            IYL = 1
         ELSE 
            IYL = 0
         ENDIF
         IF(XU.GE.XM) THEN
            IXU = 1
         ELSE
            IXU = 0
         ENDIF
         IF(YU.GE.YM) THEN
            IYU = 1
         ELSE
            IYU = 0
         ENDIF
         POS = IXL + 2*IYL + 4*IXU + 8*IYU

         IF(POS.EQ.0) THEN 
            CALL QUAD_NEIGHBOURS(II, IC) 

         ELSEIF(POS.EQ.4) THEN
            CALL QUAD_NEIGHBOURS(II,IC)
            CALL QUAD_NEIGHBOURS(II, IC+1)
            
         ELSEIF(POS.EQ.8) THEN
            CALL QUAD_NEIGHBOURS(II, IC) 
            CALL QUAD_NEIGHBOURS(II, IC+2) 
            
         ELSEIF(POS.EQ.12) THEN
            CALL QUAD_NEIGHBOURS(II, IC)
            CALL QUAD_NEIGHBOURS(II, IC+1) 
            CALL QUAD_NEIGHBOURS(II, IC+2) 
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSEIF(POS.EQ.5) THEN
            CALL QUAD_NEIGHBOURS(II, IC+1)
            
         ELSEIF(POS.EQ.13) THEN
            CALL QUAD_NEIGHBOURS(II, IC+1)
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSEIF(POS.EQ.10) THEN
            CALL QUAD_NEIGHBOURS(II, IC+2)
            
         ELSEIF(POS.EQ.14) THEN
            CALL QUAD_NEIGHBOURS(II, IC+2)
            CALL QUAD_NEIGHBOURS(II, IC+3)
            
         ELSEIF(POS.EQ.15) THEN
            CALL QUAD_NEIGHBOURS(II, IC+3) 
            
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE QUAD_NEIGHBOURS 

