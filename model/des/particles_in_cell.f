!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: DES - Finding the fluid computational cell in which        C
!           a particle lies, to calculte void fraction and also        C
!           the volume averaged solids velocity of the cell            C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTICLES_IN_CELL(PARTS)

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
      USE sendrecv
      IMPLICIT NONE

      INTEGER L, I, J, K, M, MM, PARTS, IJKL, IJK
      INTEGER IPJK, IJPK, IJKP
      DOUBLE PRECISION XE(IMAX2), YN(JMAX2), ZT(KMAX2), PVOL
      DOUBLE PRECISION SOLVOLINC(IJKMAX2,MMAX)
      DOUBLE PRECISION VOL_AVE_US(IJKMAX2,MMAX)
      DOUBLE PRECISION VOL_AVE_VS(IJKMAX2,MMAX)
      DOUBLE PRECISION VOL_AVE_WS(IJKMAX2,MMAX)

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      
      DO IJKL = 1, IJKMAX2
         PINC(IJKL) = 0
         DO M = 1, MMAX
            SOLVOLINC(IJKL,M) = 0.0
            VOL_AVE_US(IJKL,M) = 0.0
            VOL_AVE_VS(IJKL,M) = 0.0
            VOL_AVE_WS(IJKL,M) = 0.0
         END DO
      END DO

      XE(1) = 0.0
      YN(1) = 0.0

      DO I = IMIN1, IMAX2
         XE(I) = XE(I-1) + DX(I)
      END DO

      DO J  = JMIN1, JMAX2
         YN(J) = YN(J-1) + DY(J)
      END DO

      IF(DIMN.EQ.3) THEN
         ZT(1) = 0.0
         DO K = KMIN1, KMAX2
            ZT(K) = ZT(K-1) + DZ(K)
         END DO
      END IF
      
      DO L = 1, PARTICLES
         IF(EQUIVALENT_RADIUS) THEN
            MM = 1
         ELSE
            DO M = 1, MMAX
               IF((2*DES_RADIUS(L)-D_P(M)).LE.1E-5) THEN
                  MM = M
               END IF
            END DO
         END IF
         PVOL = 0.0
         PVOL = (4.0/3.0)*(22.0/7.0)*DES_RADIUS(L)**3
         DO IJKL = 1, IJKMAX2
            I = I_OF(IJKL)
            J = J_OF(IJKL)
            K = K_OF(IJKL)
            
           IF((I.GE.IMIN1).AND.(J.GE.JMIN1)) THEN   
            IF((DES_POS_NEW(1,L).GE.XE(I-1)).AND.(DES_POS_NEW(1,L).LT.XE(I))) THEN
               IF((DES_POS_NEW(2,L).GE.YN(J-1)).AND.(DES_POS_NEW(2,L).LT.YN(J))) THEN

                  IF(DIMN.EQ.3) THEN
                     IF((DES_POS_NEW(3,L).GT.ZT(K-1)).AND.(DES_POS_NEW(3,L).LE.ZT(K))) THEN
                        PINC(IJKL) = PINC(IJKL) + 1
                        NEIGHBOURS(MAXNEIGHBORS,L) = IJKL
                        SOLVOLINC(IJKL,MM) = SOLVOLINC(IJKL,MM) +  PVOL
                        VOL_AVE_US(IJKL,MM) = VOL_AVE_US(IJKL,MM) + PVOL*DES_VEL_NEW(1,L)
                        VOL_AVE_VS(IJKL,MM) = VOL_AVE_VS(IJKL,MM) + PVOL*DES_VEL_NEW(2,L)
                        VOL_AVE_WS(IJKL,MM) = VOL_AVE_WS(IJKL,MM) + PVOL*DES_VEL_NEW(3,L)
                     END IF
                  ELSE
                     PINC(IJKL) = PINC(IJKL) + 1
                     NEIGHBOURS(MAXNEIGHBORS,L) = IJKL
                     SOLVOLINC(IJKL,MM) = SOLVOLINC(IJKL,MM) + PVOL
                     VOL_AVE_US(IJKL,MM) = VOL_AVE_US(IJKL,MM) + PVOL*DES_VEL_NEW(1,L)
                     VOL_AVE_VS(IJKL,MM) = VOL_AVE_VS(IJKL,MM) + PVOL*DES_VEL_NEW(2,L)
                  END IF
               END IF
            END IF
          END IF  
         END DO
      END DO
      
      DO IJKL = 1, IJKMAX2

         K = K_OF(IJKL)

         IF(DIMN.EQ.2) THEN
            IF(COORDINATES == 'CARTESIAN') THEN
               DZ(K) = 2*RADIUS_EQ
            ELSE IF(COORDINATES == 'CYLINDRICAL') THEN
               DZ(K) = 44.0/7.0
            END IF
         END IF
         
         DO M = 1, MMAX
            IF(DIMN.EQ.3) THEN
               IF(PINC(IJKL).GT.0) THEN
                  VOL_AVE_US(IJKL,M) = VOL_AVE_US(IJKL,M)/SOLVOLINC(IJKL,M)
                  VOL_AVE_VS(IJKL,M) = VOL_AVE_VS(IJKL,M)/SOLVOLINC(IJKL,M)
                  VOL_AVE_WS(IJKL,M) = VOL_AVE_WS(IJKL,M)/SOLVOLINC(IJKL,M)
               END IF
            ELSE
               IF(PINC(IJKL).GT.0) THEN
                  VOL_AVE_US(IJKL,M) = VOL_AVE_US(IJKL,M)/SOLVOLINC(IJKL,M)
                  VOL_AVE_VS(IJKL,M) = VOL_AVE_VS(IJKL,M)/SOLVOLINC(IJKL,M)
               END IF
            END IF

            RO_S(M) = ROs
            EP_G(IJKL) = 1   
            IF(VOL(IJKL).GT.0) THEN
             ROP_S(IJKL,M) = RO_S(M)*SOLVOLINC(IJKL,M)/(VOL(IJKL)*DZ(K))
            END IF
            
            IF(PINC(IJKL).GT.0) THEN
              EP_G(IJKL) = EP_G(IJKL) - EP_S(IJKL,M)
            END IF
            
         END DO
      END DO


      DO IJKL = 1, IJKMAX2
         I = I_OF(IJKL)
         J = J_OF(IJKL)
         K = K_OF(IJKL)
         IPJK = IP_OF(IJKL)
         IJPK = JP_OF(IJKL)
         IJKP = KP_OF(IJKL)

         DO M = 1, MMAX
            IF(I.EQ.1) THEN
               U_S(IJKL,M) = 0D0
            ELSE IF(I.EQ.IMAX2) THEN
               U_S(IJKL,M) = 0D0
            ELSE 
               U_S(IJKL,M) = (VOL_AVE_US(IJKL,M)*DX(I) +&
                             VOL_AVE_US(IPJK,M)*DX(I+1))/(DX(I)+DX(I+1))
            END IF

            IF(J.EQ.1) THEN
               V_S(IJKL,M) = 0D0
            ELSE IF(J.EQ.JMAX2) THEN
               V_S(IJKL,M) = 0D0
            ELSE 
               V_S(IJKL,M) = (VOL_AVE_VS(IJKL,M)*DY(J) +&
                             VOL_AVE_VS(IJPK,M)*DY(J+1))/(DY(J)+DY(J+1))
            END IF
            
            IF(DIMN.EQ.3) THEN
               IF(K.EQ.1) THEN
                  W_S(IJKL,M) = 0D0
               ELSE IF(K.EQ.KMAX2) THEN
                  W_S(IJKL,M) = 0D0
               ELSE 
                  W_S(IJKL,M) = (VOL_AVE_WS(IJKL,M)*DZ(K) +&
                                VOL_AVE_WS(IJKP,M)*DZ(K+1))/(DZ(K)+DZ(K+1))
               END IF
            END IF
         END DO
      END DO
      
      RETURN
      END SUBROUTINE PARTICLES_IN_CELL

