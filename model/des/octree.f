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

      USE discretelement
      IMPLICIT NONE

      INTEGER I, J, TC1, TC2, TCR, TCM, PARTS
      DOUBLE PRECISION CT

!     CALL SYSTEM_CLOCK(TC1, TCR, TCM)

      CALL INIT_OCT(PARTICLES)

      DO I = 1, PARTS
         CALL ADD_OCT(I)
      END DO
      
      DO I = 1, PARTS
         J = 1
         CALL OCT_NEIGHBOURS(I, J)
      END DO

!     CALL SYSTEM_CLOCK(TC2, TCR, TCM)

!     CT =  TC2-TC1
!     IF(CT.LE.0) THEN
!     CT = TCM + TC2 -TC1
!     END IF
!     CT = CT/TCR
!     OCTCT = CT

!     PRINT *,'OCT:- CPU TIME TAKEN:',CT, TCM, TCR
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
      IMPLICIT NONE

      INTEGER I, J, PARTS
      DOUBLE PRECISION TIME, A, OMEGA, ASINOMEGAT, BOXLIMIT

      TIME = CALLED*DTSOLID

      A = 0.0
      OMEGA = 0.0
      IF(DES_F.NE.0.0) THEN
         OMEGA = 2*22*DES_F/7
         A = DES_GAMMA*9.81/(OMEGA*OMEGA)
      END IF

      ASINOMEGAT = A*SIN(OMEGA*TIME)
      BOXLIMIT = RADIUS_EQ

      DO I = 1, MAXQUADS
         DO J = 1, NMQD 
            LQUAD(J,I) = 0
         END DO
         DO J = 1, NWALLS 
            CQUAD(J,I) = 0
         END DO
      END DO

      NQUAD = 1

      CQUAD(1,1) = WX1 - BOXLIMIT
      CQUAD(2,1) = BY1- BOXLIMIT + ASINOMEGAT 
      CQUAD(3,1) = EX2 + BOXLIMIT 
      CQUAD(4,1) = TY2 + BOXLIMIT 
      CQUAD(5,1) = SZ1 - BOXLIMIT
      CQUAD(6,1) = NZ2 + BOXLIMIT 

      RETURN
      END SUBROUTINE INIT_OCT 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_OCT(II)                                            C
!  Purpose: To find the 'oct' for the particle                         C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADD_OCT(II)

      USE discretelement

      IMPLICIT NONE
      
      INTEGER II
      INTEGER J ,IQ, IQS, IFOUND, K


      IQS = 1
      IFOUND = 0
      DO WHILE(IFOUND.EQ.0)
         CALL FIND_OCT(II, IQS)
         J = LQUAD(11,IQS)
         IF(J.EQ.8) THEN
	    CALL SPLIT_OCT(IQS)
         ELSE IF(J.LT.0) THEN
            PRINT *,'ERROR IN SPLIT_OCT'
         ELSE 
            LQUAD(11,IQS) = J+1
            LQUAD(J+1,IQS) = II
            IFOUND = 1
         END IF
      END DO            
      RETURN
      END SUBROUTINE ADD_OCT 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE FIND_OCT(I,Q)
	USE discretelement
	INTEGER I, Q, I_F, IC, IX, IY, IZ

	I_F = 0
	DO WHILE(I_F.EQ.0) 
	   IF(LQUAD(11,Q).GE.0) THEN
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
	     IF(DES_POS_NEW(3,I).GE.CQUAD(6,IC)) THEN
		IZ = 1
	     ELSE
		IZ = 0
	     END IF
	     NC = IX + 2*IY + 4*IZ
	     Q = IC + NC
           END IF
	END DO
	RETURN
	END SUBROUTINE FIND_OCT

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE SPLIT_OCT(Q)	
	USE discretelement
	INTEGER Q
	INTEGER I_FQ, J, K, ORPHANS(8), NC
	DOUBLE PRECISION XL, YL, XU, YU, ZL, ZU, XMID, YMID, ZMID
	LQUAD(11,Q) = -1
	DO J = 1, 8 
	   ORPHANS(J) = LQUAD(J,Q)
	   LQUAD(J,Q) = NQUAD + J
	   LQUAD(10,NQUAD+J) = Q
	   LQUAD(9,NQUAD+J) = J
	END DO
	XL = CQUAD(1,Q)
	YL = CQUAD(2,Q)
	XU = CQUAD(3,Q)
	YU = CQUAD(4,Q) 
	ZL = CQUAD(5,Q)
	ZU = CQUAD(6,Q)
		
	XMID = (XU-XL)/2
	YMID = (YU-YL)/2
	ZMID = (ZU-ZL)/2

	CQUAD(1,NQUAD+1) = XL
	CQUAD(2,NQUAD+1) = YL
	CQUAD(3,NQUAD+1) = XL+XMID
	CQUAD(4,NQUAD+1) = YL+YMID
	CQUAD(5,NQUAD+1) = ZL
	CQUAD(6,NQUAD+1) = ZL+ZMID

	CQUAD(1,NQUAD+2) = XL+XMID	
	CQUAD(2,NQUAD+2) = YL
	CQUAD(3,NQUAD+2) = XU
	CQUAD(4,NQUAD+2) = YL+YMID
	CQUAD(5,NQUAD+2) = ZL
	CQUAD(6,NQUAD+2) = ZL+ZMID

	CQUAD(1,NQUAD+3) = XL
	CQUAD(2,NQUAD+3) = YL+YMID
	CQUAD(3,NQUAD+3) = XL+XMID
	CQUAD(4,NQUAD+3) = YU
	CQUAD(5,NQUAD+3) = ZL
	CQUAD(6,NQUAD+3) = ZL+ZMID

	CQUAD(1,NQUAD+4) = XL+XMID
	CQUAD(2,NQUAD+4) = YL+YMID
	CQUAD(3,NQUAD+4) = XU
	CQUAD(4,NQUAD+4) = YU
	CQUAD(5,NQUAD+4) = ZL
	CQUAD(6,NQUAD+4) = ZL+ZMID

	CQUAD(1,NQUAD+5) = XL
	CQUAD(2,NQUAD+5) = YL
	CQUAD(3,NQUAD+5) = XL+XMID
	CQUAD(4,NQUAD+5) = YL+YMID
	CQUAD(5,NQUAD+5) = ZL+ZMID
	CQUAD(6,NQUAD+5) = ZU

	CQUAD(1,NQUAD+6) = XL+XMID	
	CQUAD(2,NQUAD+6) = YL
	CQUAD(3,NQUAD+6) = XU
	CQUAD(4,NQUAD+6) = YL+YMID
	CQUAD(5,NQUAD+6) = ZL+ZMID
	CQUAD(6,NQUAD+6) = ZU 

	CQUAD(1,NQUAD+7) = XL
	CQUAD(2,NQUAD+7) = YL+YMID
	CQUAD(3,NQUAD+7) = XL+XMID
	CQUAD(4,NQUAD+7) = YU
	CQUAD(5,NQUAD+7) = ZL+ZMID
	CQUAD(6,NQUAD+7) = ZU

	CQUAD(1,NQUAD+8) = XL+XMID
	CQUAD(2,NQUAD+8) = YL+YMID
	CQUAD(3,NQUAD+8) = XU
	CQUAD(4,NQUAD+8) = YU
	CQUAD(5,NQUAD+8) = ZL+ZMID
	CQUAD(6,NQUAD+8) = ZU
	
	NQUAD = NQUAD + 8
	DO J = 1, 8
	   I_FQ = Q
	   K = ORPHANS(J)
	   CALL FIND_OCT(K,I_FQ) 
	   NC = LQUAD(11,I_FQ)
	   LQUAD(11,I_FQ) = NC+1
	   LQUAD(NC+1,I_FQ) = K 
	END DO

       IF(I_FQ.GT.MAXQUADS) THEN
         PRINT *,'CHECK I_FQ', I_FQ
         STOP
       END IF

	RETURN
	END SUBROUTINE SPLIT_OCT

       
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

      USE discretelement
      IMPLICIT NONE

      INTEGER II, I, N, IC, J, K, POS, IXU, IXL, IYU, IYL, IZL, IZU
      DOUBLE PRECISION XM, YM, ZM, DIST, R_LM, XL, XU, YL, YU, ZL, ZU 
      INTEGER DIA_RATIO, LL, NQ  

      DIA_RATIO = INT(RADIUS_EQ/DES_RADIUS(II))

      IF(LQUAD(11,NQ).GT.0) THEN
         DO I = 1, LQUAD(11,NQ)
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
      ELSE IF(LQUAD(11,NQ).EQ.-1) THEN
         IC = LQUAD(1,NQ)
         XM = CQUAD(3,IC)
         YM = CQUAD(4,IC)
         ZM = CQUAD(6,IC)
         XL = DES_POS_NEW(1,II) - (DIA_RATIO + 1) * DES_RADIUS(II)
         XU = DES_POS_NEW(1,II) + (DIA_RATIO + 1) * DES_RADIUS(II)
         YL = DES_POS_NEW(2,II) - (DIA_RATIO + 1) * DES_RADIUS(II)
         YU = DES_POS_NEW(2,II) + (DIA_RATIO + 1) * DES_RADIUS(II)
         ZL = DES_POS_NEW(3,II) - (DIA_RATIO + 1) * DES_RADIUS(II)
         ZU = DES_POS_NEW(3,II) + (DIA_RATIO + 1) * DES_RADIUS(II)

         IF(XL.LT.CQUAD(1,1)) XL=CQUAD(1,1) 		
         IF(XU.GT.CQUAD(3,1)) XU=CQUAD(3,1)
         IF(YL.LT.CQUAD(2,1)) YL=CQUAD(2,1)
         IF(YU.GT.CQUAD(4,1)) YU=CQUAD(4,1)   
         IF(ZL.LT.CQUAD(5,1)) ZL=CQUAD(5,1)
         IF(ZU.GT.cquad(6,1)) ZU=CQUAD(6,1)  

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

