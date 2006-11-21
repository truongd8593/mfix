!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  PARTICLES_IN_CELL                                     C
!  Purpose: DES - Finding the fluid computational cell in which        C
!           a particle lies, to calculte void fraction and also        C
!           the volume averaged solids velocity of the cell            C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C 
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTICLES_IN_CELL

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

      INTEGER L, I, J, K, M, MM 
      INTEGER IJK, IPJK, IJPK, IJKP
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,MMAX), OSOLVOL
!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS1 = .TRUE.

      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      PINC(:) = 0
      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

      IF(FIRST_PASS1) THEN
         XE(1) = ZERO
         YN(1) = ZERO
         DO I = IMIN1, IMAX2
            XE(I) = XE(I-1) + DX(I)
         END DO
         DO J  = JMIN1, JMAX2
            YN(J) = YN(J-1) + DY(J)
         END DO
         IF(DIMN.EQ.3) THEN
            ZT(1) = ZERO
            DO K = KMIN1, KMAX2
               ZT(K) = ZT(K-1) + DZ(K)
            END DO
         END IF
      ENDIF

      DO L = 1, PARTICLES

	 IF(S_TIME.LE.DTSOLID.OR.FIRST_PASS1) THEN ! Brute force technique to determine the particle locations in the Eulerian grid
            
                                ! currently the code does not distinguish different density particles with exactly same size - needs to be fixed
            DO M = 1, MMAX
               IF(ABS(2.0d0*DES_RADIUS(L)-D_P0(M)).LT.SMALL_NUMBER) THEN
                  PIJK(L,5) = M 
                  RO_S(M) = RO_Sol(L)
               END IF
            END DO
 
            IF(PIJK(L,5).EQ.0) THEN
              WRITE(*,*) 'Problem determining the solids association in PIC', L, DES_RADIUS(L), D_P0
            ENDIF


            DO I = IMIN1, IMAX3
               IF((DES_POS_NEW(L,1).GE.XE(I-1)).AND.(DES_POS_NEW(L,1).LT.XE(I))) THEN
                  PIJK(L,1) = I
                  GO TO 10
               END IF
            END DO

 10         CONTINUE
            DO J = JMIN1, JMAX3
               IF((DES_POS_NEW(L,2).GE.YN(J-1)).AND.(DES_POS_NEW(L,2).LT.YN(J))) THEN
                  PIJK(L,2) = J
                  GO TO 20
               END IF
            END DO

 20         CONTINUE
!           write(*,*) 'pijk', L, PIJK(L,1),  PIJK(L,2)
            IF(DIMN.EQ.2) THEN
               PIJK(L,3)  = 1
               GO TO 30
            ELSE
               DO K = KMIN1, KMAX3
                  IF((DES_POS_NEW(L,3).GT.ZT(K-1)).AND.(DES_POS_NEW(L,3).LE.ZT(K))) THEN 
                     PIJK(L,3) = K
                     GO TO 30
                  END IF
               END DO
            END IF
           

         ELSE                   ! Incremental approach to determine the new location of the particles

            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
!           write(*,*) 'pijk2', L, PIJK(L,1),  PIJK(L,2) , XE(I-1), XE(I)
            IF((DES_POS_NEW(L,1).GE.XE(I-1)).AND.(DES_POS_NEW(L,1).LT.XE(I)).OR.I.EQ.1) THEN
               GO TO 40 
            ELSE IF(DES_VEL_NEW(L,1).GT.ZERO) THEN
               IF((DES_POS_NEW(L,1).GE.XE(I)).AND.(DES_POS_NEW(L,1).LT.XE(I+1))) PIJK(L,1) = I+1
            ELSE IF(DES_VEL_NEW(L,1).LT.ZERO) THEN
               IF((DES_POS_NEW(L,1).GE.XE(I-2)).AND.(DES_POS_NEW(L,1).LT.XE(I-1))) PIJK(L,1) = I-1
            ELSE 
               PRINT *,'des/particles_in_cell.f : CHECK CELL I' , PIJK(L,1),  PIJK(L,2), L, DES_POS_NEW(L,1), DES_VEL_NEW(L,1)
               STOP
            END IF
 40         CONTINUE
!           write(*,*) 'pijk2', L, PIJK(L,1),  PIJK(L,2) 
            IF((DES_POS_NEW(L,2).GE.YN(J-1)).AND.(DES_POS_NEW(L,2).LT.YN(J)).OR.J.EQ.1) THEN
               GO TO 50 
            ELSE IF(DES_VEL_NEW(L,2).GT.ZERO) THEN
               IF((DES_POS_NEW(L,2).GE.YN(J)).AND.(DES_POS_NEW(L,2).LT.YN(J+1))) PIJK(L,2) = J+1
            ELSE IF(DES_VEL_NEW(L,2).LT.ZERO) THEN
               IF((DES_POS_NEW(L,2).GE.YN(J-2)).AND.(DES_POS_NEW(L,2).LT.YN(J-1))) PIJK(L,2) = J-1
            ELSE
               PRINT *,'des/particles_in_cell.f : CHECK CELL J' 
               STOP
            END IF
 50         CONTINUE
            IF(DIMN.EQ.2) THEN
               PIJK(L,3) = 1
            ELSE
               IF((DES_POS_NEW(L,3).GE.ZT(K-1)).AND.(DES_POS_NEW(L,3).LT.ZT(K)).OR.K.EQ.1) THEN
                  GO TO 30 
               ELSE IF(DES_VEL_NEW(L,3).GT.ZERO) THEN
                  IF((DES_POS_NEW(L,3).GE.ZT(K)).AND.(DES_POS_NEW(L,3).LT.ZT(K+1))) PIJK(L,3) = K+1
               ELSE IF(DES_VEL_NEW(L,3).LT.ZERO) THEN
                  IF((DES_POS_NEW(L,3).GE.ZT(K-2)).AND.(DES_POS_NEW(L,3).LT.ZT(K-1))) PIJK(L,3) = K-1
               ELSE 
                  PRINT *,'des/particles_in_cell.f : CHECK CELL K'
                  STOP
               END IF
            END IF

         END IF

 30      CONTINUE
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         PIJK(L,4) = IJK
         PINC(IJK) = PINC(IJK) + 1
         MM = PIJK(L,5)
         SOLVOLINC(IJK,MM) = SOLVOLINC(IJK,MM) +  PVOL(L)
         DES_U_s(IJK,MM) = DES_U_s(IJK,MM) + PVOL(L)*DES_VEL_NEW(L,1)
         DES_V_s(IJK,MM) = DES_V_s(IJK,MM) + PVOL(L)*DES_VEL_NEW(L,2)
         IF(DIMN.EQ.3) DES_W_s(IJK,MM) = DES_W_s(IJK,MM) + PVOL(L)*DES_VEL_NEW(L,3)

      END DO

      DO IJK = IJKSTART3, IJKEND3
         K = K_OF(IJK) 
         DO M = 1, MMAX
            IF(PINC(IJK).GT.0) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)   
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DIMN.EQ.3) THEN
                  DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
               END IF
            END IF
            EP_G(IJK) = ONE   
            IF(VOL(IJK).GT.0) THEN
               ROP_S(IJK,M)  = RO_S(M)*SOLVOLINC(IJK,M)/(VOL(IJK)*DZ(K))
            END IF
            IF(PINC(IJK).GT.0) THEN
               EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
            END IF
         END DO
      END DO

      FIRST_PASS1 = .False.

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL

