!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: DES - Calculte the drag force and pressure force           C
!           on particles exerted by the gas. Cell centered             C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DRAG_FGS(NPARTS)

      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      IMPLICIT NONE
      
      INTEGER IPJK, IJPK, IJKP, IMJK, IJMK, IJKM 
      INTEGER L, LL, I, J, K, KK, M, MM, IJK
      INTEGER NPARTS
      DOUBLE PRECISION P_FORCE(DIMENSION_3,DIMN)
      DOUBLE PRECISION TEMP1, TEMP2 
      DOUBLE PRECISION UGC, VGC, WGC
      DOUBLE PRECISION OEPS, OVOL
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

       DO IJK = IJKSTART3, IJKEND3
      
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)

         IF(PINC(IJK).GT.0) THEN

             UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
             VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
            
           DO M = 1, MMAX
             SOLID_DRAG(IJK,M,1) = -F_GS(IJK,M)*&
                                         (DES_U_S(IJK,M)-UGC)
             SOLID_DRAG(IJK,M,2) = -F_GS(IJK,M)*&
                                         (DES_V_S(IJK,M)-VGC)
             IF(DIMN.EQ.3) THEN
                WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
                SOLID_DRAG(IJK,M,3) = -F_GS(IJK,M)*&
                                         (DES_W_S(IJK,M)-WGC)
             END IF
             OEPS = ONE/EP_S(IJK,M)
             SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
           END DO     
                 
            IF(I.EQ.IMIN1) THEN
               TEMP2 = (P_G(IJK)+P_G(IPJK))/2
               TEMP1 = (P_G(IJK)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
               P_FORCE(IJK,1) = (TEMP1-TEMP2)*DY(J)*DZ(K)
            ELSE IF(I.EQ.IMAX1) THEN
               TEMP2 = (P_G(IMJK)+P_G(IJK))/2
               TEMP1 = (P_G(IJK)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
               P_FORCE(IJK,1) = (TEMP2 - TEMP1)*DY(J)*DZ(K)
            ELSE IF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
               TEMP2 = (P_G(IJK)+P_G(IPJK))/2
               TEMP1 = (P_G(IMJK)+P_G(IJK))/2
               P_FORCE(IJK,1) = (TEMP1 - TEMP2)*DY(J)*DZ(K)
            END IF

            IF(J.EQ.JMIN1) THEN
               TEMP2 = (P_G(IJK)+P_G(IJPK))/2
               TEMP1 = (P_G(IJK)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
               P_FORCE(IJK,2) = (TEMP1 - TEMP2)*DX(I)*DZ(K)
            ELSE IF(J.EQ.JMAX1) THEN
               TEMP2 = (P_G(IJMK)+P_G(IJK))/2
               TEMP1 = (P_G(IJK)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
               P_FORCE(IJK,2) = (TEMP2 - TEMP1)*DX(I)*DZ(K)
            ELSE IF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
               TEMP2 = (P_G(IJK)+P_G(IJPK))/2
               TEMP1 = (P_G(IJMK)+P_G(IJK))/2
               P_FORCE(IJK,2) = (TEMP1 - TEMP2)*DX(I)*DZ(K)
            END IF

            IF(DIMN.EQ.3) THEN
            IF(K.EQ.KMIN1) THEN
               TEMP2 = (P_G(IJK)+P_G(IJKP))/2
               TEMP1 = (P_G(IJK)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
               P_FORCE(IJK,3) = (TEMP1 - TEMP2)*DX(I)*DY(J)
            ELSE IF(K.EQ.KMAX1) THEN
               TEMP2 = (P_G(IJKM)+P_G(IJK))/2
               TEMP1 = (P_G(IJK)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
               P_FORCE(IJK,3) = (TEMP2 - TEMP1)*DX(I)*DY(J)
            ELSE IF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
               TEMP2 = (P_G(IJK)+P_G(IJKP))/2
               TEMP1 = (P_G(IJKM)+P_G(IJK))/2
               P_FORCE(IJK,3) = (TEMP1 - TEMP2)*DX(I)*DY(J)
            END IF
            END IF
         END IF
       END DO
      
       DO KK = 1, PARTICLES
          IJK = PIJK(KK,4)
            K = PIJK(KK,3)
            MM = PIJK(KK,5)
            IF(DIMN.EQ.2) THEN
               OVOL = ONE/(VOL(IJK)*DZ(K))
            ELSE
               OVOL = ONE/VOL(IJK)
            END IF
            FC(KK,:)=FC(KK,:)+ (SOLID_DRAG(IJK,MM,:) +&
                               P_FORCE(IJK,:)*OVOL)*PVOL(KK)
       END DO
      
      RETURN
      END SUBROUTINE DRAG_FGS


