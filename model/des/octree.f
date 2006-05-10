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

      USE param1
      USE constant
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
         CALL OCT_NEIGHBOURS(I,J)
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
!  Module name: ADD_OCT(II)                                            C
!  Purpose: To find the 'oct' for the particle                         C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADD_OCT(II)
      USE param1
      USE discretelement

      IMPLICIT NONE
      
      INTEGER II
      INTEGER J ,IQ, IQS, IFOUND, K


      IQS = 1
      IFOUND = 0
      DO WHILE(IFOUND.EQ.0)
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
      END DO            
      RETURN
      END SUBROUTINE ADD_OCT 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

	SUBROUTINE FIND_OCT(I,Q)
        USE param1
	USE discretelement
	INTEGER I, Q, I_F, IC, IX, IY, IZ

	I_F = 0
	DO WHILE(I_F.EQ.0) 
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
	END DO
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

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER II, I, N, IC, J, JJ, K, POS, IXU, IXL, IYU, IYL, IZL, IZU
      DOUBLE PRECISION XM, YM, ZM, DIST, R_LM, XL, XU, YL, YU, ZL, ZU 
      INTEGER DIA_RATIO, LL, NQ  

      DIA_RATIO = INT(RADIUS_EQ/DES_RADIUS(II))

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
                  J = NEIGHBOURS(II,1)
                  JJ = NEIGHBOURS(N,1)
                  IF(J.LE.MN) THEN
                     NEIGHBOURS(II,J+1) = N
                  ELSE
                     PRINT *,'QUAD - NEIGHBORS GT MN'
                     PRINT *, II, (NEIGHBOURS(II,LL),LL=1,MAXNEIGHBORS)
                     STOP
                  END IF
                  IF(JJ.LE.MN) THEN
                     NEIGHBOURS(II,JJ+1) = II 
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
         XL = DES_POS_NEW(II,1) - (DIA_RATIO + 1) * DES_RADIUS(II)
         XU = DES_POS_NEW(II,1) + (DIA_RATIO + 1) * DES_RADIUS(II)
         YL = DES_POS_NEW(II,2) - (DIA_RATIO + 1) * DES_RADIUS(II)
         YU = DES_POS_NEW(II,2) + (DIA_RATIO + 1) * DES_RADIUS(II)
         ZL = DES_POS_NEW(II,3) - (DIA_RATIO + 1) * DES_RADIUS(II)
         ZU = DES_POS_NEW(II,3) + (DIA_RATIO + 1) * DES_RADIUS(II)

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

