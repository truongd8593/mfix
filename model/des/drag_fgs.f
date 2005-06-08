!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: DES - Calculte the drag force and pressure force           C
!           on particles exerted by the gas                            C
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
      INTEGER L, LL, I, J, K, KK, M, MM, IJK, IJKL
      INTEGER NPARTS
      DOUBLE PRECISION SOLID_DRAG(NDIM,IJKMAX2,MMAX)
      DOUBLE PRECISION P_FORCE(NDIM,IJKMAX2)
      DOUBLE PRECISION PVOL, TEMP1, TEMP2 

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      DO IJKL = 1, IJKMAX2
      
         I = I_OF(IJKL)
         J = J_OF(IJKL)
         K = K_OF(IJKL)
         IMJK = IM_OF(IJKL)
         IPJK = IP_OF(IJKL)
         IJMK = JM_OF(IJKL)
         IJPK = JP_OF(IJKL)
         IJKM = KM_OF(IJKL)
         IJKP = KP_OF(IJKL)

         IF(PINC(IJKL).GT.0) THEN

            DO MM= 1, MMAX
               SOLID_DRAG(1,IJKL,MM) = 0 - F_GS(IJKL,MM)*&
                                           (U_S(IJKL,MM)-U_G(IJKL))
               SOLID_DRAG(2,IJKL,MM) = 0 - F_GS(IJKL,MM)*&
                                           (V_S(IJKL,MM)-V_G(IJKL))
               IF(DIMN.EQ.3) THEN
                  SOLID_DRAG(3,IJKL,MM) = 0 - F_GS(IJKL,MM)*&
                                           (W_S(IJKL,MM)-W_G(IJKL))
               END IF
               IF(EP_S(IJKL,MM).GT.0) THEN
                  SOLID_DRAG(1,IJKL,MM) = SOLID_DRAG(1,IJKL,MM)*&
                                          (1/EP_S(IJKL,MM)) 
                  SOLID_DRAG(2,IJKL,MM) = SOLID_DRAG(2,IJKL,MM)*&
                                          (1/EP_S(IJKL,MM)) 
                  IF(DIMN.EQ.3) THEN
                    SOLID_DRAG(3,IJKL,MM) = SOLID_DRAG(3,IJKL,MM)*&
                                            (1/EP_S(IJKL,MM)) 
                  END IF
               END IF
            END DO
            
            IF(I.EQ.IMIN1) THEN
               TEMP2 = (P_G(IJKL)+P_G(IPJK))/2
               TEMP1 = (P_G(IJKL)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
               P_FORCE(1,IJKL) = (TEMP1-TEMP2)*DY(J)*DZ(K)
            ELSE IF(I.EQ.IMAX1) THEN
               TEMP2 = (P_G(IMJK)+P_G(IJKL))/2
               TEMP1 = (P_G(IJKL)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
               P_FORCE(1,IJKL) = (TEMP2 - TEMP1)*DY(J)*DZ(K)
            ELSE IF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
               TEMP2 = (P_G(IJKL)+P_G(IPJK))/2
               TEMP1 = (P_G(IMJK)+P_G(IJKL))/2
               P_FORCE(1,IJKL) = (TEMP1 - TEMP2)*DY(J)*DZ(K)
            END IF

            IF(J.EQ.JMIN1) THEN
               TEMP2 = (P_G(IJKL)+P_G(IJPK))/2
               TEMP1 = (P_G(IJKL)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
               P_FORCE(2,IJKL) = (TEMP1 - TEMP2)*DX(I)*DZ(K)
            ELSE IF(J.EQ.JMAX1) THEN
               TEMP2 = (P_G(IJMK)+P_G(IJKL))/2
               TEMP1 = (P_G(IJKL)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
               P_FORCE(2,IJKL) = (TEMP2 - TEMP1)*DX(I)*DZ(K)
            ELSE IF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
               TEMP2 = (P_G(IJKL)+P_G(IJPK))/2
               TEMP1 = (P_G(IJMK)+P_G(IJKL))/2
               P_FORCE(2,IJKL) = (TEMP1 - TEMP2)*DX(I)*DZ(K)
            END IF

            IF(DIMN.EQ.3) THEN
            IF(K.EQ.KMIN1) THEN
               TEMP2 = (P_G(IJKL)+P_G(IJKP))/2
               TEMP1 = (P_G(IJKL)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
               P_FORCE(3,IJKL) = (TEMP1 - TEMP2)*DX(I)*DY(J)
            ELSE IF(K.EQ.KMAX1) THEN
               TEMP2 = (P_G(IJKM)+P_G(IJKL))/2
               TEMP1 = (P_G(IJKL)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
               P_FORCE(3,IJKL) = (TEMP2 - TEMP1)*DX(I)*DY(J)
            ELSE IF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
               TEMP2 = (P_G(IJKL)+P_G(IJKP))/2
               TEMP1 = (P_G(IJKM)+P_G(IJKL))/2
               P_FORCE(3,IJKL) = (TEMP1 - TEMP2)*DX(I)*DY(J)
            END IF
            END IF

         END IF

      END DO
      
         DO KK = 1, PARTICLES
         
            PVOL = (4.0/3.0)*(22.0/7.0)*DES_RADIUS(KK)**3
            IJKL = NEIGHBOURS(MAXNEIGHBORS,KK)

            IF(IJKL.GT.0) THEN             
            K = K_OF(IJKL)
            
            IF(EQUIVALENT_RADIUS) THEN
               MM = 1
            ELSE
               DO M = 1, MMAX
                  IF((2*DES_RADIUS(KK)-D_P(IJK,M)).LE.1E-5) THEN
                     MM = M
                  END IF
               END DO
            END IF
            
            IF(KUIPERS_DRAG) THEN
               FC(1,KK)=FC(1,KK)+ (SOLID_DRAG(1,IJKL,MM)*EP_S(IJKL,MM)/&
               PINC(IJKL)) + (P_FORCE(1,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               FC(2,KK)=FC(2,KK)+ (SOLID_DRAG(2,IJKL,MM)*EP_S(IJKL,MM)/&
               PINC(IJKL)) + (P_FORCE(2,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               IF(DIMN.EQ.3) THEN
                 FC(3,KK)=FC(3,KK)+ (SOLID_DRAG(3,IJKL,MM)*EP_S(IJKL,MM)/& 
                   PINC(IJKL)) + (P_FORCE(3,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               END IF
            ELSE
               FC(1,KK)=FC(1,KK)+ (SOLID_DRAG(1,IJKL,MM)*PVOL) +&
                                  (P_FORCE(1,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               FC(2,KK)=FC(2,KK)+ (SOLID_DRAG(2,IJKL,MM)*PVOL) +&
                                  (P_FORCE(2,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               IF(DIMN.EQ.3) THEN
                 FC(3,KK)=FC(3,KK)+ (SOLID_DRAG(3,IJKL,MM)*PVOL) +&
                                    (P_FORCE(3,IJKL))*(PVOL/(VOL(IJKL)*DZ(K)))
               END IF
            END IF
            END IF
         END DO
      
      RETURN
      END SUBROUTINE DRAG_FGS


