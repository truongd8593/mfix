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

      USE discretelement
      IMPLICIT NONE

      INTEGER I, J, TC1, TC2, TCR, TCM, PARTS
      DOUBLE PRECISION CT

!     CALL SYSTEM_CLOCK(TC1, TCR, TCM)

      CALL INIT_QUAD(PARTICLES)
      
      DO I = 1, PARTS
         CALL ADD_QUAD(I)
      END DO

      DO I = 1, PARTS
         J = 1
         CALL QUAD_NEIGHBOURS(I, J)
      END DO

!     CALL SYSTEM_CLOCK(TC2, TCR, TCM)

!     CT =  TC2-TC1
!     IF(CT.LE.0) THEN
!     CT = TCM + TC2 -TC1
!     END IF
!     CT = CT/TCR
!     QUADCT = CT

!     PRINT *,'QUAD:- CPU TIME TAKEN:',CT

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
      DOUBLE PRECISION STIME, A, OMEGA, ASINOMEGAT, BOXLIMIT

      STIME = CALLED*DTSOLID

      A = 0D0
      OMEGA = 0D0
      IF(DES_F.NE.0.0) THEN
         OMEGA = 2*22*DES_F/7
         IF(UNITS == "CGS") THEN
            A = DES_GAMMA*980/(OMEGA*OMEGA)
         ELSE
            A = DES_GAMMA*9.81/(OMEGA*OMEGA)
         END IF
      END IF

      ASINOMEGAT = A*SIN(OMEGA*STIME)
      BOXLIMIT = RADIUS_EQ 

      DO I = 1, MAXQUADS
         DO J = 1, NMQD 
            LQUAD(J,I) = 0D0
         END DO
         DO J = 1, NWALLS
            CQUAD(J,I) = 0D0
         END DO
      END DO

      NQUAD = 1

      CQUAD(1,1) = WX1 - BOXLIMIT
      CQUAD(2,1) = BY1 - BOXLIMIT + ASINOMEGAT
      CQUAD(3,1) = EX2 + BOXLIMIT
      CQUAD(4,1) = TY2 + BOXLIMIT


      RETURN
      END SUBROUTINE INIT_QUAD 


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

      USE discretelement
      IMPLICIT NONE
      
      INTEGER II
      INTEGER J ,IQ, IQS, IFOUND, K

      IQS = 1
      IFOUND = 0
      DO WHILE(IFOUND.EQ.0)
         CALL FIND_QUAD(II, IQS)
         J = LQUAD(7,IQS)
         IF(J.EQ.4) THEN
	    CALL SPLIT_QUAD(IQS)
         ELSE IF(J.LT.0) THEN
            PRINT *,'ERROR IN SPLIT_QUAD'
         ELSE 
            LQUAD(7,IQS) = J+1
            LQUAD(J+1,IQS) = II
            IFOUND = 1
         END IF
      END DO            
      RETURN
      END SUBROUTINE ADD_QUAD 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE FIND_QUAD(I,Q)
      
	USE discretelement
        IMPLICIT NONE
        
	INTEGER I, Q, I_F, IC, IX, IY, NC


	I_F = 0
	DO WHILE(I_F.EQ.0) 
	   IF(LQUAD(7,Q).GE.0) THEN
 	     I_F = 1
           ELSE
	     IC = LQUAD(1,Q)
	     IF(DES_POS_NEW(1,I).GE.CQUAD(3,IC)) THEN
		IX = 1
	     ELSE 
		IX = 0
	     END IF
             IF(DES_POS_NEW(2,I).GE.CQUAD(4,IC)) THEN
		IY = 1
	     ELSE
                IY = 0
             END IF
	     NC = IX + 2*IY
	     Q = IC + NC
           END IF
	END DO
	RETURN
	END SUBROUTINE FIND_QUAD

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C           

	SUBROUTINE SPLIT_QUAD(Q)	
           
	USE discretelement
        IMPLICIT NONE
        
	INTEGER Q
	INTEGER I_FQ, J, K, ORPHANS(4), NC
	DOUBLE PRECISION XL, YL, XU, YU, XMID, YMID

	LQUAD(7,Q) = -1
	DO J = 1, 4
	   ORPHANS(J) = LQUAD(J,Q)
	   LQUAD(J,Q) = NQUAD + J
	   LQUAD(6,NQUAD+J) = Q
	   LQUAD(5,NQUAD+J) = J
	END DO
	XL = CQUAD(1,Q)
	YL = CQUAD(2,Q)
	XU = CQUAD(3,Q)
	YU = CQUAD(4,Q) 
	XMID = (XU-XL)/2
	YMID = (YU-YL)/2
	CQUAD(1,NQUAD+1) = XL
	CQUAD(2,NQUAD+1) = YL
	CQUAD(3,NQUAD+1) = XL+XMID
	CQUAD(4,NQUAD+1) = YL+YMID
	CQUAD(1,NQUAD+2) = XL+XMID	
	CQUAD(2,NQUAD+2) = YL
	CQUAD(3,NQUAD+2) = XU
	CQUAD(4,NQUAD+2) = YL+YMID
	CQUAD(1,NQUAD+3) = XL
	CQUAD(2,NQUAD+3) = YL+YMID
	CQUAD(3,NQUAD+3) = XL+XMID
	CQUAD(4,NQUAD+3) = YU
	CQUAD(1,NQUAD+4) = XL+XMID
	CQUAD(2,NQUAD+4) = YL+YMID
	CQUAD(3,NQUAD+4) = XU
	CQUAD(4,NQUAD+4) = YU
	NQUAD = NQUAD +4

	DO J = 1, 4
	   I_FQ = Q
	   K = ORPHANS(J)
	   CALL FIND_QUAD(K,I_FQ) 
	   NC = LQUAD(7,I_FQ)
	   LQUAD(7,I_FQ) = NC+1
	   LQUAD(NC+1,I_FQ) = K 
	END DO

	IF(I_FQ.GT.MAXQUADS) THEN
          PRINT *,'CHECK I_FQ', I_FQ
	  STOP
	END IF

	RETURN
	END SUBROUTINE SPLIT_QUAD

          
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

      USE discretelement
      IMPLICIT NONE

      INTEGER II, I, N, LL, IC, J, K, POS, IXU, IXL, IYU, IYL, KK, NQ
      DOUBLE PRECISION XM, YM, DIST, R_LM, XL, XU, YL, YU 
      INTEGER DIA_RATIO

      DIA_RATIO = INT(RADIUS_EQ/DES_RADIUS(II))

      IF (LQUAD(7,NQ).GT.0) THEN
         DO I = 1, LQUAD(7,NQ)
            N = LQUAD(I,NQ)
            IF(II.NE.N) THEN
               DIST = 0
               DO K = 1, DIMN
               DIST = DIST + (DES_POS_NEW(K,II)-DES_POS_NEW(K,N))**2
               END DO
               DIST = SQRT(DIST)
               R_LM = DES_RADIUS(II)+DES_RADIUS(N)
               IF(DIST.LE.R_LM) THEN
		  NEIGHBOURS(1,II) = NEIGHBOURS(1,II) + 1
		  J = NEIGHBOURS(1,II)
		  IF(J.LT.MN) THEN
                     NEIGHBOURS(J+1,II) = N
	          ELSE IF(J.GE.MN) THEN
                     PRINT *,'QUAD', II, J, (NEIGHBOURS(LL,II),LL=1,J)
                     STOP
	          END IF
               END IF
            END IF
         END DO
      ELSE IF(LQUAD(7,NQ).EQ.-1) THEN
         IC = LQUAD(1,NQ)
         XM = CQUAD(3,IC)
         YM = CQUAD(4,IC)

         XL = DES_POS_NEW(1,II) - (DIA_RATIO + 1) * DES_RADIUS(II)
         XU = DES_POS_NEW(1,II) + (DIA_RATIO + 1) * DES_RADIUS(II)
         YL = DES_POS_NEW(2,II) - (DIA_RATIO + 1) * DES_RADIUS(II)
         YU = DES_POS_NEW(2,II) + (DIA_RATIO + 1) * DES_RADIUS(II)

         IF(XL.LT.CQUAD(1,1)) XL=CQUAD(1,1)
         IF(XU.GT.CQUAD(3,1)) XU=CQUAD(3,1)
         IF(YL.LT.CQUAD(2,1)) YL=CQUAD(2,1)
         IF(YU.GT.CQUAD(4,1)) YU=CQUAD(4,1)

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

