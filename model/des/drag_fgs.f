!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DRAG_FGS                                               C
!>
!!  Purpose: DES - Calculte the drag force and pressure force           
!!           on particles exerted by the gas. Cell centered            
!<
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number 3                                  Date: 2-July-07  C
!  Author: Rahul Garg                                                  C
!>
!!  Purpose: Now the drag_fgs routine is called from calc_drag in model 
!!  directory as well as by calc_forces_des. Calling arguments have     
!!  also changed. Depending on the choice, once can obtain drag force   
!!  based on local velocities or averaged velocities
!<
!   
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_FGS
      
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
      USE interpolation
      
      IMPLICIT NONE
!Logical to decide if to calculate darg forces on particles
      LOGICAL :: focus
!INTEGER, INTENT(IN) , OPTIONAL :: ERROR_INDEX
      
      INTEGER IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IPJPK, IPJKP, IJPKP&
      &, IPJPKP, IJK_EAST, IJK_NORTH, IJK_TOP
      INTEGER L, LL, I, J, K, KK, M, MM, IJK, IER, IJKFOCUS, FOCUS_PARTICLE2

      DOUBLE PRECISION P_FORCE(DIMENSION_3,DIMN), D_FORCE(DIMN)

      DOUBLE PRECISION TEMP_U_G(JMAX2, KMAX2), TEMP_V_G(IMAX2,KMAX2), TEMP_W_G(IMAX2, JMAX2)
      DOUBLE PRECISION TEMP1, TEMP2, AVG_FACTOR, AVG_FACTOR2
      DOUBLE PRECISION UGC, VGC, WGC
      DOUBLE PRECISION OEPS, OVOL, F_pf(DIMN), drag_bm_tmp(DIMN)

      INTEGER, DIMENSION(3):: PCELL
      INTEGER:: ONEW,PART_IJK
      INTEGER:: IB, IE, JB, JE, KB, KE, II,JJ
      DOUBLE PRECISION, DIMENSION(3) :: DRAG_P 
      DOUBLE PRECISION:: EPS_P, VCELL
      INTEGER:: NP
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      
      FOCUS_PARTICLE2 = 400000
!J = JMAX1
!K = 1
!DO I = 1, IMAX1
!   IJK = funijk(I,J,K)
!   IJK_EAST = funijk(I,JMAX2,K)
!   IF(CALLFROMDES) PRINT*,'UG = ', U_G(IJK), U_G(IJK_EAST), (U_G(IJK)+U_G(IJK_EAST))/2.d0
!ENDDO

      IF(.NOT.CALC_FC) GOTO 500 
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
               IF(EP_S(IJK,M).GT.ZERO) THEN
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
               END IF
            END DO
            
            IF(I.EQ.IMIN1) THEN
               TEMP2 = (P_G(IJK)+P_G(IPJK))/2
               TEMP1 = (P_G(IJK)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
               P_FORCE(IJK,1) = (TEMP1-TEMP2)*DY(J)*DZ(K )
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

 500  CONTINUE
      if(.not.DES_INTERP_ON) goto 200

      IF(INTX_PER) THEN 
         I = 1
         DO J = 1, JMAX2
            DO K = 1, KMAX2
               IJK = funijk(I,J,K)
               IJK_EAST = funijk(IMAX1, J, K)
               TEMP_U_G(J,K) = U_G(IJK) 
               U_G(IJK) = U_G(IJK_EAST)
            end DO
         ENDDO
      ENDIF
      
      IF(INTY_PER) THEN 
         J = 1
         DO I = 1, IMAX2
            DO K = 1, KMAX2
               IJK = funijk(I,J,K)
               IJK_NORTH = funijk(I, JMAX1, K)
               TEMP_V_G(I,K) = V_G(IJK) 
               V_G(IJK) = V_G(IJK_NORTH)
            ENDDO
         ENDDO
      ENDIF 

      IF(INTZ_PER.AND.DIMN.EQ.3) THEN 
         K = 1
         DO J = 1, JMAX2
            DO I = 1, IMAX2
               IJK = funijk(I,J,K)
               IJK_TOP = funijk(I, J, KMAX1)
               TEMP_W_G(I,J) = W_G(IJK) 
               W_G(IJK) = W_G(IJK_TOP)
            ENDDO
         ENDDO
      ENDIF 

      call set_interpolation_scheme(2)
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)
      AVG_FACTOR2 = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)
      
      DRAG_AM = ZERO
      DRAG_BM = ZERO
      wtbar = zero
      DO NP = 1, PARTICLES
         FOCUS = .FALSE.
! CALCUALTE THE DRAG FORCE ON EACH PARTICLE USING THE PARTICLE VELOCITY
         I = PIJK(NP, 1)
         J = PIJK(NP, 2)
         K = PIJK(NP, 3)

         IJK = PIJK(NP, 4)
         M = PIJK(NP,5)
         
         PART_IJK = PINC(IJK)
         
         OVOL = ONE/VOL(IJK)
         PCELL(1) = I-1
         PCELL(2) = J-1
         PCELL(3) = (3-DIMN)*1+(DIMN-2)*(K-1)

         
         CALL SET_INTERPOLATION_STENCIL(PCELL, IB, IE, JB, JE, KB,&
         &KE, INTERP_SCHEME, DIMN, ORDERNEW = ONEW) 
         IF(NP.EQ.FOCUS_PARTICLE2) THEN 
            FOCUS = .TRUE.
            PRINT*, 'CELL = ', PIJK(NP, 1)-1,PIJK(NP, 2)-1, PIJK(NP, 3)-1
            PRINT*,' IB IE1  = ', IB,IE, JB, JE, KB, KE
         ENDIF
         CALL INTERPOLATE_QUANTS( IB, IE, JB, JE, KB, KE, onew,&
         & DES_POS_NEW(NP, 1:DIMN), VEL_FP(NP, 1:DIMN), FOCUS)
         
         
!=========================================================
!NOW CALCULATE THE PARTICLE CENTERED DRAG COEFFICIENT
!=========================================================
         CALL DES_DRAG_GS(NP, vel_fp(NP,1:DIMN),des_vel_new(np,1:dimn))
         
         IF(CALC_FC) then
            OVOL = (ONE/VOL(PIJK(NP,4)))
            
            
            D_FORCE(1:dimn) = (f_gp(np)*(vel_fp(np,1:dimn)-des_vel_new(np,1:dimn)))
            
!PRINT*, 'D_FORCE = ', D_FORCE(2), FC(NP,2), f_gp(np)
            FC(NP,:) = FC(NP,:) + D_FORCE(:)
            FC(NP,:) = FC(NP,:) +  (P_FORCE(PIJK(NP,4),:)*OVOL)*PVOL(NP)
         ENDIF
         
         IF(CALLFROMDES) GOTO 300 
         
         
         DO k = 1, (3-dimn)*1+(dimn-2)*onew
            DO j = 1, onew
               DO i = 1, onew
                  II = IB+I-1
                  JJ = JB+J-1
                  KK = KB+K-1
                  
                  IF(II.LT.1.and.intx_per) II = IMAX1+II-1
                  IF(II.GT.IMAX1.and.intx_per) II = II-IMAX1+1
                  IF(JJ.LT.1.and.inty_per) JJ = JMAX1+JJ-1
                  IF(JJ.GT.JMAX1.and.inty_per) JJ = JJ-JMAX1+1
                  IF(KK.LT.1.and.intz_per) KK = KMAX1+KK-1
                  IF(KK.GT.KMAX1.and.intz_per) KK = KK-KMAX1+1
                  
                  IJK = funijk(II,JJ,KK)
                  VCELL = VOL(PIJK(NP,4))
                  
                  IF(II.EQ.1.or.II.EQ.IMAX1) VCELL = 0.5d0*VCELL
                  IF(JJ.EQ.1.or.JJ.EQ.JMAX1) VCELL = 0.5d0*VCELL
                  IF(DIMN.EQ.3)THEN
                     IF(KK.EQ.1.or.KK.EQ.KMAX1) VCELL = 0.5d0*VCELL
                  ENDIF

                  OVOL = ONE/VCELL
                  
                  drag_am(ii,jj,kk,M) = drag_am(ii,jj,kk,M) + F_gp(np)  &
                  & *(weightp(i,j,k))*OVOL
! drag_am(ii,jj,kk,M) = zero
!First remove the velocity component at this grid point from the vel_fp
                  drag_bm_tmp(1:dimn) = VEL_FP(NP, 1:DIMN) - WEIGHTP(I,J,K)*VSTENCIL(I,J,K, 1:dimn)
! Now find the remaning drag force
                  drag_bm_tmp(1:dimn) = des_vel_new(np,1:dimn) !- DRAG_BM_TMP(1:DIMN)
                  
                  drag_bm(ii,jj,kk, 1:DIMN,M) = drag_bm(ii,jj,kk,1:DIMN,M) + &
                  & F_gp(np)*(DRAG_BM_TMP(1:DIMN)) * weightp(i&
                  &,j,k)*OVOL
                  
!drag_bm(ii,jj,kk, 1:DIMN,M) = zero 

                  wtbar(ii,jj,kk,M) = wtbar(ii,jj,kk,M) + weightp(i,j,k)*RO_S(M)*OVOL*PVOL(NP)

               enddo
            enddo
         enddo
         
 300     continue

      end DO                    !NP = 1, PARTICLES
      
      IF(callfromdes) goto 400
      IF(INTX_PER) THEN
         I = 1
         DO K = 1, KMAX1
            DO J = 1, JMAX1
               DRAG_BM(i,j,k, 1,1:MMAX) =HALF*( DRAG_BM(i,j,k, 1,1:MMAX) + DRAG_BM(IMAX1, J&
               &,K,1, 1:MMAX))
               DRAG_BM(IMAX1, j, k, 1, 1:MMAX) = DRAG_BM(1, J, K, 1, 1:MMAX)

               WTBAR(I,J,K,:) = HALF*(WTBAR(I,J,K,:) + WTBAR(IMAX1, J, K, :))
               WTBAR(IMAX1, J,K,:) = WTBAR(1,J,K,:)

               DRAG_AM(I,J,K,:) = HALF*(DRAG_AM(I,J,K,:) + DRAG_AM(IMAX1, J, K, :))
               DRAG_AM(IMAX1, J,K,:) = DRAG_AM(1,J,K,:)
            end DO
         ENDDO
      ENDIF
      
      IF(INTY_PER) THEN
         J = 1
         DO K = 1, KMAX1
            DO I = 1, IMAX1
               

               DRAG_BM(i,j,k, 2,1:MMAX) = HALF*(DRAG_BM(i,j,k,  2,1:MMAX) + DRAG_BM(I, JMAX1&
               &,K, 2,1:MMAX))
               DRAG_BM(I, JMAX1, K,  2,1:MMAX) = DRAG_BM(I, 1, K,  2,1:MMAX)
               WTBAR(I,J,K,:) =HALF*( WTBAR(I,J,K,:) + WTBAR(I, JMAX1, K, :))
               WTBAR(I, JMAX1,K,:) = WTBAR(I,1,K,:)
               
               DRAG_AM(I,J,K,:) =HALF*( DRAG_AM(I,J,K,:) + DRAG_AM(I, JMAX1, K, :))
               DRAG_AM(I, JMAX1,K,:) = DRAG_AM(I,1,K,:)
            end DO
         ENDDO
      ENDIF

      IF(INTZ_PER.AND.DIMN.EQ.3) THEN
         K = 1
         DO J = 1, JMAX1
            DO I = 1, IMAX1
               DRAG_BM(i,j,k, 3,1:MMAX) = HALF*(DRAG_BM(i,j,k, 3,1:MMAX) + DRAG_BM(I, J&
               &,KMAX1, 3,1:MMAX))
               DRAG_BM(I, J, KMAX1,  3,1:MMAX) = DRAG_BM(I, J, 1,  3,1:MMAX)
               WTBAR(I,J,K,:) = HALF*(WTBAR(I,J,K,:) + WTBAR(I, J, KMAX1, :))
               WTBAR(I, J,KMAX1,:) = WTBAR(I,J,1,:)

               DRAG_AM(I,J,K,:) = HALF*(DRAG_AM(I,J,K,:) + DRAG_AM(I, J, KMAX1, :))
               DRAG_AM(I, J,KMAX1,:) = DRAG_AM(I,J,1,:)
            end DO
         ENDDO
      ENDIF


      DO M = 1, MMAX
         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN
               i = i_of(ijk)
               j = j_of(ijk)
               k = k_of(ijk)
               
               F_GS(IJK,M) = AVG_FACTOR2*(DRAG_AM(i,j,k,M) +&
               & DRAG_AM(I,J-1,K,M) + DRAG_AM(I-1,J-1,K,M)&
               &+DRAG_AM(I-1,J,K,M))
               
               ROP_S(IJK,M) = AVG_FACTOR2*(WTBAR(i,j,k,M) +&
               & WTBAR(I,J-1,K,M) + WTBAR(I-1,J-1,K,M)&
               &+WTBAR(I-1,J,K,M))

               
               if(dimn.eq.3) then 
                  F_GS(IJK,M) = F_GS(IJK,M) + AVG_FACTOR2*(DRAG_AM(i&
                  &,j,k-1,M) + DRAG_AM(I,J-1,K-1,M) + DRAG_AM(I-1,J&
                  &-1,K-1,M)+DRAG_AM(I-1,J,K-1,M))
                  
                  ROP_S(IJK,M) = ROP_S(IJK,M) + AVG_FACTOR2*(WTBAR(i&
                  &,j,k-1,M) + WTBAR(I,J-1,K-1,M) + WTBAR(I-1,J&
                  &-1,K-1,M)+WTBAR(I-1,J,K-1,M))
               end IF
               
               
            end IF
!     IF(J.EQ.JMAX/2) PRINT*, F_GS(IJK,M), EP_S(IJK,M),(18*MU_g0*EP_S(IJK,M))/((2.d0*DES_RADIUS(1))**2.d0)
         end DO
      end DO
 400  continue        

      IF(INTX_PER) THEN 
         I = 1
         DO J = 1, JMAX2
            DO K = 1, KMAX2
               IJK = funijk(I,J,K)
               U_G(IJK) = TEMP_U_G(J, K)
            end DO
         ENDDO
      ENDIF
      
      IF(INTY_PER) THEN 
         J = 1
         DO K = 1, KMAX2
            DO I = 1, IMAX2
               IJK = funijk(I,J,K)
               V_G(IJK) = TEMP_V_G(I,K)
            ENDDO
         ENDDO
      ENDIF 

      IF(INTZ_PER.AND.DIMN.EQ.3) THEN 
         K = 1
         DO J = 1, JMAX2
            DO I = 1, IMAX2
               IJK = funijk(I,J,K)
               W_G(IJK) = TEMP_W_G(I,J)
            ENDDO
         ENDDO
      ENDIF 

      RETURN

 200  continue
      if(calc_fc) then
         
         DO KK = 1, PARTICLES
            IJK = PIJK(KK,4)
            K = PIJK(KK,3)
            MM = PIJK(KK,5)
            OVOL = ONE/VOL(IJK)
            FC(KK,:)=FC(KK,:)  + (SOLID_DRAG(IJK,MM,:)*PVOL(KK)) 
            FC(KK,:) = FC(KK,:) + (P_FORCE(IJK,:)*OVOL)*PVOL(KK) 
            
         END DO
      end if

      RETURN

      END SUBROUTINE DRAG_FGS


      SUBROUTINE INTERPOLATE_QUANTS( IB, IE, JB, JE, KB, KE, onew, POSP, fvel, focus)
      
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
      USE interpolation
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IB, IE, JB, JE, KB, KE, onew
      DOUBLE PRECISION :: AVG_FACTOR
      DOUBLE PRECISION, DIMENSION(DIMN) , INTENT(OUT) :: FVEL
      LOGICAL :: FOCUS
      DOUBLE PRECISION, DIMENSION(DIMN) , INTENT(IN) :: POSP
      INTEGER :: I,J,K, II, JJ, KK , IJK, IJK_PER, IJKP_PER

      INTEGER IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IPJPK, IPJKP, IJPKP&
      &, IPJPKP
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)

      DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
         DO J = 1, ONEW
            DO I = 1, ONEW
               II = IB+I-1
               JJ = JB+J-1
               KK = KB+K-1
               GSTENCIL(I,J,K,1) = XE(II)
               GSTENCIL(I,J,K,2) = YN(JJ)
               GSTENCIL(I,J,K,3) = ZT(KK)*(DIMN-2) + DZ(1)*(3-DIMN)
               
               IF(II.LT.1.and.intx_per) II = IMAX1+II
               IF(II.GT.IMAX1.and.intx_per) II = II-IMAX1
               IF(JJ.LT.1.and.inty_per) JJ = JMAX1+JJ
               IF(JJ.GT.JMAX1.and.inty_per) JJ = JJ-JMAX1
               IF(KK.LT.1.and.intz_per) KK = KMAX1+KK
               IF(KK.GT.KMAX1.and.intz_per) KK = KK-KMAX1
               
               ijk = funijk(II,JJ,KK)
               
               ipjk = ip_of (ijk)    
               ijpk = jp_of (ijk)
               ijkp = kp_of (ijk)
               ijpkp = kp_of(ijpk)
               ipjkp = kp_of(ipjk)
               ipjpk = jp_of(ipjk)
               ipjpkp = kp_of(ipjpk)
               
               vstencil(i,j,k,1) = AVG_FACTOR*( u_g(ijk) + u_g(ijpk)&
               & + (u_g(ijkp) + u_g(ijpkp))*(dimn-2) )
               
               
               vstencil(i,j,k,2) = AVG_FACTOR*( v_g(ijk) + v_g(ipjk)&
               & + (v_g(ijkp) + v_g(ipjkp))*(dimn-2) )
               
               IF(FOCUS) THEN 
                  print*,'II JJ KK = ', ii, jj, kk
                  PRINT*,'v_g = ', v_g(ijk), v_g(ipjk), v_g(ijkp), v_g(ipjkp), vstencil(i,j,k,2)
               ENDIF
               
               if(dimn.eq.3) then 
                  vstencil(i,j,k,3) = AVG_FACTOR*(w_g(ijk) +&
                  & w_g(ijpk) + w_g(ipjk) + w_g(ipjpk) )
               else 
                  vstencil(i,j,k,3) = 0.d0 !doesn't matter what
! ever value is put here
               end if
               
            END DO
         END DO
      END DO
      
      if(dimn.eq.2) then 
         CALL interpolator(gstencil(1:onew,1:onew,1,1:dimn)&
         &,vstencil(1:onew,1:onew,1,1:2),POSP(1:2),FVEL(1:2),onew&  
         &,Interp_scheme,weightp)
         
      else 
         CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn)&
         &,vstencil(1:onew,1:onew,1:onew,1:dimn)&
         &,POSP( 1:3),FVEL(1:3),onew&
         &,interp_scheme,weightp)   
      end if
      IF(FOCUS) THEN 
         PRINT*, ' '
         PRINT*,'IN INTERP'
         
         PRINT*,' IB IE2  = ', IB,IE, JB, JE, KB, KE
         DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
            DO J = 1, ONEW
               DO I = 1, ONEW
                  
                  PRINT*,' V and Weight= ', VSTENCIL(I,J,K, 2), weightp(I,J,K)
               end DO
            end DO
         end DO

         print*, ' fvel = ', fvel(2)
         PRINT*, ' PPOS = ', POSP(2) !, gstencil(1:onew,1:onew,1,1)
!PRINT*, 'weight =', weightp(1:onew,1:onew,1:onew), SUM(weightp)
      ENDIF

      END SUBROUTINE INTERPOLATE_QUANTS

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_DRAG_GS (NP)                                       C
!  Purpose: Calculate the gas-particle drag coefficient for des        C
!           calculation such that in drag correlation, exact values of C
!            u_g and  v_sm are used                                    C
! Comments: No BVK drag model in this subroutine. BVK requires average C
!           diameter which needs to be defined for DEM case            C

!  Author: R. Garg and Jin Sun                        date 06/28/07    C 
!  Comments: It is not used in this current version pending some more  C
!            tests. 
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Variables referenced: EP_g, RO_g, MU_g, D_p                         C
!  Variables modified: DRAG_gs                                         C
!                                                                      C
!  Local variables: A, B, V_rm, Re                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     

      SUBROUTINE DES_DRAG_GS(KK, fvel, des_vel) 

!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar  
      USE drag  
      USE sendrecv 
      USE discretelement
      USE ur_facs 

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     particle number
      INTEGER , INTENT(IN) ::          KK


!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1) 
!     defined in the data file.
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0 
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0 
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0 
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0 
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: fvel, des_vel
!     
!     Indices 
      INTEGER          I,  IJK, IMJK, IJMK, IJKM, Im, M
!     
!     Cell center value of U_sm 
      DOUBLE PRECISION USCM 
!     
!     Cell center value of U_g 
      DOUBLE PRECISION UGC 
!     
!     Cell center value of V_sm 
      DOUBLE PRECISION VSCM 
!     
!     Cell center value of V_g 
      DOUBLE PRECISION VGC 
!     
!     Cell center value of W_sm 
      DOUBLE PRECISION WSCM 
!     
!     Cell center value of W_g 
      DOUBLE PRECISION WGC 
!     
!     Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION VREL 
!     
!     Reynolds number 
      DOUBLE PRECISION Re 
!     
!     Ratio of settling velocity of a multiparticle 
!     system to that of a single particle 
      DOUBLE PRECISION V_rm 
!     
!     Function of EP_g 
      DOUBLE PRECISION A 
!     
!     Function of EP_g 
      DOUBLE PRECISION B 
!     
!     Single sphere drag coefficient x Re 
      DOUBLE PRECISION C_DsxRe, C_DsxReT 
!     
!     single sphere drag coefficient 
      DOUBLE PRECISION C_d 
!     
!     drag coefficient 
      DOUBLE PRECISION DgA  

!     --- Gidaspow switch function variables [ceaf 2006-03-23]
      DOUBLE PRECISION Ergun
      DOUBLE PRECISION WenYu
      DOUBLE PRECISION PHI_gs
!     --- end Gidaspow switch function variables

!     
!     Gas Laminar viscosity redefined here to set
!     viscosity at pressure boundaries
      DOUBLE PRECISION Mu
!     
!     Gidaspow Reynolds number
      DOUBLE PRECISION Re_g
!     
!***********************************************************
!     Declaration of variables relevant to the Koch and Hill
!     drag correlation, sof
!***********************************************************
!     Stokes Drag Force
      DOUBLE PRECISION F_STOKES
!     
!     zero Re function for low Reynolds number
      DOUBLE PRECISION F_0
!     
!     inertial function for low Reynolds number
      DOUBLE PRECISION F_1
!     
!     zero Re function for high Reynolds number
      DOUBLE PRECISION F_2
!     
!     inertial function for high Reynolds number
      DOUBLE PRECISION F_3
!     
!     dimensionless drag force F
      DOUBLE PRECISION F
!     
!     transition Reynolds numbers
      DOUBLE PRECISION Re_Trans_1, Re_Trans_2
!     
!     solids volume fraction
      DOUBLE PRECISION phis
!     
!     weighting factor to compute F_0 and F_2
      DOUBLE PRECISION w, D_p_av, Y_i
!     
!     Hill and Koch Reynolds number
      DOUBLE PRECISION Re_kh
      
!     
!     End of Koch and Hill variables declaration, sof
!***********************************************************
!     
!     
!     Current value of F_gs (i.e., without underrelaxation)

      DOUBLE PRECISION F_gstmp

      DOUBLE PRECISION:: EPS, DIAMETER

      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER) ! Tsuji drag
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2 ! Dalla Valle (1948) 
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0 * Re**0.657D0)      ! Turton and
!     &          + 0.413D0 * Re**2.09D0 / (Re**1.09D0 + 16300.D0) ! Levenspiel (1986)
!     
!     
!     
!PRINT*,'IN CALC_DRAG, DRAG_TYPE = ', DRAG_TYPE
      IJK = PIJK(KK,4)
      M = PIJK(KK,5)
      EPS  = EP_S(IJK, M)
      DIAMETER = 2.D0*DES_RADIUS(KK)
      
      IF(DIMN == 2)THEN
         VREL = SQRT((FVEL(1)- DES_VEL(1))**2 +&
         (FVEL(2) - DES_VEL(2))**2)
      ELSE
         VREL = SQRT((FVEL(1) - DES_VEL(1))**2 +&
         (FVEL(2) - DES_VEL(2))**2 +&
         (FVEL(3) - DES_VEL (3))**2)
      ENDIF
      
      IF (P_OUTFLOW_AT(IJK)) THEN
         IF( FLUID_AT(EAST_OF(IJK) )) THEN
            Mu = MU_G(EAST_OF(IJK))
         ELSE IF ( FLUID_AT(WEST_OF(IJK)) ) THEN
            Mu = MU_G(WEST_OF(IJK))
         ELSE IF ( FLUID_AT(NORTH_OF(IJK)) ) THEN
            Mu = MU_G(NORTH_OF(IJK))
         ELSE IF ( FLUID_AT(SOUTH_OF(IJK)) ) THEN
            Mu = MU_G(SOUTH_OF(IJK))
         ELSE IF ( FLUID_AT(TOP_OF(IJK)) ) THEN
            Mu = MU_G(TOP_OF(IJK))
         ELSE IF ( FLUID_AT(BOTTOM_OF(IJK)) ) THEN
            Mu = MU_G(BOTTOM_OF(IJK))
         ENDIF
      ELSE
         Mu = MU_G(IJK)
      ENDIF

!     Reynolds number
      if(Mu > ZERO)then
         RE = diameter * VREL*RO_G(IJK)/Mu 
         
!     Note the presence of gas volume fraction in ROP_G
         RE_G = (diameter * VREL*ROP_G(IJK))/Mu
         
!     Note Reynolds' number for Hill and Koch has an additional factor of 1/2 & ep_g
         RE_kh = (0.5D0*diameter*VREL*ROP_G(IJK))/Mu
         
!Print*,'eps= ',eps
      else 
         RE = LARGE_NUMBER 
         RE_G = LARGE_NUMBER
         RE_kh = LARGE_NUMBER
      endif
!     f_gp() =  single particle drag excluding vector(v_g - v_p)
!     

!---------------Begin Syamlal and O'Brien ---------------------------
!     
!     Calculate V_rm
!     
      IF(TRIM(DRAG_TYPE).EQ.'SYAM_OBRIEN') then

         IF (EP_s(IJK,M) <= ZERO) THEN 
            F_gstmp = ZERO 
         ELSE IF (EP_G(IJK) == ZERO) THEN 
            F_gstmp = ZERO 
         ELSE 
            A = EP_G(IJK)**4.14D0 
            IF (EP_G(IJK) <= 0.85D0) THEN 
               B = drag_c1*EP_G(IJK)**1.28D0 
            ELSE 
               B = EP_G(IJK)**drag_d1
            ENDIF 
            V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) 
!------------------Begin cluster correction --------------------------
!     uncomment the following four lines and comment the above line V_RM=... 
!     V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) & 
!     * ( ONE + C(1) * exp( -a2*(Re - Re_c)**2 &
!     - a3*(EP_g(IJK)-ep_c)**2 &
!     )       * Re * (1. - EP_g(IJK))                )
!------------------End cluster correction ----------------------------
!     
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
!     
            IF(TSUJI_DRAG) THEN
               IF(EP_G(IJK).LE.0.8D0) THEN
                  F_gstmp = (Mu*PVOL(KK)/(DIAMETER**2))*&
                  (150D0*(EP_S(IJK,M)/EP_G(IJK)) + 1.75D0*RE)
               ELSE IF(EP_G(IJK).GT.0.8D0) THEN
                  IF(RE*EP_G(IJK).GT.1000D0) THEN
                     F_gstmp = 0.75D0*0.43D0*Mu*PVOL(KK)*RE/(DIAMETER**2 *&
                     EP_G(IJK)**1.7D0)
                  ELSE IF(RE*EP_G(IJK).LE.1000D0) THEN
                     F_gstmp = 0.75D0*C_DSXRET(RE*EP_G(IJK))*Mu*PVOL(KK)*&
                     RE/(DIAMETER**2 *EP_G(IJK)**1.7D0)
                  END IF
               END IF 
            ELSE IF(MODEL_B) THEN 
               F_gstmp = 0.75D0*Mu*(PVOL(KK))*C_DSXRE(RE/V_RM)/(&
               V_RM*DIAMETER*DIAMETER) 
            ELSE
               F_gstmp = 0.75D0*Mu*(PVOL(KK))*EP_G(IJK)*C_DSXRE(RE&
               /V_RM)/(V_RM*DIAMETER*DIAMETER) 
            ENDIF
         ENDIF
!---------------End Syamlal and O'Brien ---------------------------
!     
!--------------------------Begin Gidaspow --------------------------
      ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') then
         IF(EP_g(IJK) .LE. 0.8D0) THEN
            DgA = 150D0 * (ONE - EP_g(IJK)) * Mu &
            / ( EP_g(IJK) * diameter**2 ) &
            + 1.75D0 * RO_g(IJK) * VREL / diameter
         ELSE
            IF(Re_G .LE. 1000D0)THEN
               C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF
            DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
            /diameter
         ENDIF
         
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         
         IF(Model_B)THEN
            F_gstmp = DgA*(PVOL(KK))/EP_g(IJK)
!F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
         ELSE
            F_gstmp = DgA * (PVOL(KK)) !EP_s(IJK, M)
         ENDIF

         
!--------------------------End Gidaspow --------------------------
!     
!-----------------------Begin Gidaspow_blend ---------------------
      ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND') then
!     Dense phase - EP_g < 0.8
         Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu &
         / ( EP_g(IJK) * diameter**2 ) &
         + 1.75D0 * RO_g(IJK) * VREL / diameter
!     
!     Dilute phase - EP_g >= 0.8
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         WenYu = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
         /diameter
!     
!     Switch function
         PHI_gs = ATAN(150D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
!     Blend the models
         DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
         
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         
         IF(Model_B)THEN
            F_gstmp = DgA*(PVOL(KK))/EP_g(IJK)
!F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
         ELSE
            F_gstmp = DgA * (PVOL(KK)) !EP_s(IJK, M)
         ENDIF
         
!-----------------------End Gidaspow_blend -----------------------
!     
!--------------------------Begin WEN_YU --------------------------
      ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU') then
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
         /diameter
         
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            
            F_gstmp = DgA*(PVOL(KK))/EP_g(IJK)
!F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
         ELSE
            F_gstmp = DgA *(PVOL(KK))
         ENDIF
         
!--------------------------End WEN_YU ----------------------------

!--------------------Begin Koch & Hill (2001) --------------------
!     
!!!   Added by Clay Sutton (Lehigh University) 7-14-04
!!!   
!!!   MODIFICATIONS:
!!!   
!!!   1) Declared new variables F_STOKES, F_0, F_1, F_3
!!!   
!!!   2) Added new drag closure lines
!!!   
!!!   Clay's implementation was modified by Sof (01-21-2005)
!!!   for a report explaining these changes contact sof@fluent.com
!     
      ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') then
!     
         F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)*EP_g(IJK)/diameter**2
         
         phis = ONE-EP_G(IJK)   ! EP_s(IJK,M) for polydisperse systems (sof --> 03-27-2007)
         w = EXP(-10.0D0*(0.4D0-phis)/phis)
         
         IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
            F_0 = (1.0D0-w) *                                           &
            (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
            *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*      &
            phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *   &
            10.0D0*phis/(1.0D0-phis)**3
            
         ELSE IF(phis >= 0.4D0) THEN
            F_0 = 10.0D0*phis/(1.0D0-phis)**3
         ENDIF
         
         IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
            F_1 = dsqrt(2.0D0/phis) / 40.0D0
         ELSE IF(phis > 0.1D0) THEN
            F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
         ENDIF
         
         IF(phis < 0.4D0) THEN
            F_2 = (1.0D0-w) *                                           &
            (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis    &
            *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*      &
            phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *  &
            10.0D0*phis/(1.0D0-phis)**3
            
         ELSE
            F_2 = 10.0D0*phis/(1.0D0-phis)**3
         ENDIF
         
         IF(phis < 0.0953D0) THEN
            F_3 = 0.9351D0*phis + 0.03667D0
         ELSE
            F_3 = 0.0673D0 + 0.212D0*phis +0.0232D0/(1.0-phis)**5
         ENDIF
         
         Re_Trans_1 = (F_2 - 1.0D0)/(3.0D0/8.0D0 - F_3)
         Re_Trans_2 = (F_3 + dsqrt(F_3*F_3 - 4.0D0*F_1 &
         *(F_0-F_2))) / (2.0D0*F_1)
         
         IF(phis <= 0.01D0 .AND. Re_kh <= Re_Trans_1) THEN
            F = 1.0D0 + 3.0D0/8.0D0*Re_kh
            
         ELSE IF(phis > 0.01D0 .AND. Re_kh <= Re_Trans_2) THEN
            F = F_0 + F_1*Re_kh*Re_kh
            
            
         ELSE IF(phis <= 0.01D0 .AND. Re_kh > Re_Trans_1 .OR.         &
            phis >  0.01D0 .AND. Re_kh > Re_Trans_2) THEN
            F = F_2 + F_3*Re_kh
            
         ELSE
            F = zero
         ENDIF
         
!     This is a check for phis (or eps_(ijk,m)) to be within physical range
         IF(phis <= ZERO .OR. phis > ONE) F = zero
         
         DgA = F * F_STOKES
!!!   
!!!   Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         
         IF(Model_B)THEN
            F_gstmp = DgA*(PVOL(KK))/EP_g(IJK)
         ELSE
            F_gstmp = DgA *(PVOL(KK))
         ENDIF
!-------------------End Koch and Hill (2001)---------------------------------------
         
      ELSE IF((DRAG_TYPE).EQ.'DILUTE_CASE') then
         C_d =  C_DSXRET(RE)
!DgA = (0.75D0 * C_d * VREL * RO_g(IJK))/diameter
         DgA = (0.75D0 * C_d*VREL * RO_g(IJK))/(diameter)
!              IF(RE.GT.1.d0.AND.CALLFROMDES) THEN 
!          PRINT*,'FOR DILUTE CASE'
!          PRINT*, 'RE>1', RE, ' PARICLE ID = ', KK
!          PRINT*, 'PART VEL. = ', DES_VEL(:)
!          PRINT*, 'FLU. VEL. = ', FVEL(:)
!          PRINT*, 'PCELL: ', PIJK(KK, 1) -1 , PIJK(KK, 2) -1
!          PRINT*, 'MIN FLUID VELOCITY: ', MINVAL(U_G), MINVAL(V_G)
!          PRINT*, 'MAX FLUID VELOCITY: ', MAXVAL(U_G), MAXVAL(V_G)
!       ENDIF
         IF(Model_B)THEN
            F_gstmp = DgA*(PVOL(KK))/EP_g(IJK)
         ELSE
            F_gstmp = DgA *(PVOL(KK))
         ENDIF
!-------------------End Koch and Hill (2001)---------------------------------------
!--------------------------End Dilue case --------------------------
      ELSE
         CALL START_LOG 
!IF(.not.DMP_LOG)call open_pe_log(ier)
         if(mype == pe_io) WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         CALL END_LOG 
         call mfix_exit(myPE)  
      ENDIF
!PRINT*,'re= ', RE, CALC_FC
      F_gp(kk) = (ONE - UR_F_gs) * F_gp(KK) + UR_F_gs * F_gstmp
! F_gs(IJK, M) = F_gs(IJK,M) + F_gp(KK)/VOL(IJK)
      end SUBROUTINE DES_DRAG_GS

    
