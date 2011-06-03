!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DRAG_FGS                                               C
!  Purpose: DES - Calculte the drag force and pressure force           
!           on particles exerted by the gas. Cell centered            
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number 3                                  Date: 2-July-07  C
!  Author: Rahul Garg                                                  C
!
!  Purpose: Now the drag_fgs routine is called from calc_drag in model 
!  directory as well as by calc_forces_des. Calling arguments have     
!  also changed. Depending on the choice, once can obtain drag force   
!  based on local velocities or averaged velocities
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_FGS
      
      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------         
! local variable used for debugging
      LOGICAL FOCUS 

      DOUBLE PRECISION P_FORCE(DIMENSION_3,DIMN)
      DOUBLE PRECISION D_FORCE(DIMN)
      DOUBLE PRECISION drag_bm_tmp(DIMN)

! temporary variables used to calculate pressure at scalar cell edge      
      DOUBLE PRECISION TEMP1, TEMP2

! average fluid velocity in x, y, z direction at scalar cell center      
      DOUBLE PRECISION UGC, VGC, WGC

! local variables to temporarily store gas velocities near periodic boundaries      
      DOUBLE PRECISION TEMP_U_G_X(JMAX2, KMAX2),&
         TEMP_V_G_X(JMAX2,KMAX2), TEMP_W_G_X(JMAX2, KMAX2)
      DOUBLE PRECISION TEMP_U_G_Y(IMAX2, KMAX2),&
         TEMP_V_G_Y(IMAX2,KMAX2), TEMP_W_G_Y(IMAX2, KMAX2)
      DOUBLE PRECISION TEMP_U_G_Z(IMAX2, JMAX2),&
          TEMP_V_G_Z(IMAX2,JMAX2), TEMP_W_G_Z(IMAX2, JMAX2)
          
! indices used with periodic boundaries
      INTEGER IJK_IMAX1, IJK_JMAX1, IJK_KMAX1

! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
              IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK

! i,j,k indices of the fluid cell the particle resides in minus 1 
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL

! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER IW, IE, JS, JN, KB, KTP

! order of interpolation set in the call to set_interpolation_scheme unless it
! is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW

!mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION MPG_CYCLIC(DIMN)

! constant whose value depends on dimension of system      
      DOUBLE PRECISION AVG_FACTOR     

! index of solid phase that particle NP belongs to      
      INTEGER M

! particle number index, used for looping      
      INTEGER NP

! index to track accounted for particles 
      INTEGER PC 

! volume of fluid cell 
      DOUBLE PRECISION VCELL

! solids volume fraction of phase M in fluid cell
      DOUBLE PRECISION EP_SM

! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION OEPS, OVOL 

! for error messages      
      INTEGER IER
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'


! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
      AVG_FACTOR = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)   
      MPG_CYCLIC(1:DIMN) = ZERO 
      IF(CYCLIC_X_PD) MPG_CYCLIC(1) = DELP_X/XLENGTH
      IF(CYCLIC_Y_PD) MPG_CYCLIC(2) = DELP_Y/YLENGTH
      IF(CYCLIC_Z_PD.AND.DIMN.EQ.3) MPG_CYCLIC(3) = DELP_Z/ZLENGTH

! if calc_fc is true, then the contact forces (FC) will be updated 
! to include gas-solids drag and gas pressure force on the particle
! in this section the gas pressure force on the particle is computed.
! if the user did not specify using an interpolated drag force (i.e.,
! des_interp_on=t), then the gas solid drag force on the particle is
! also computed based on cell average quantities
!-----------------------------------------------      
      IF(CALC_FC) THEN

! initialization              
         SOLID_DRAG(:,:,:) = ZERO
         P_FORCE(:,:) = ZERO

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

! todo: fix interpolation of p_force points for consistency with mfix
! and to allow non-uniform grids            
            IF(PINC(IJK).GT.0) THEN
               IF(IMIN1.EQ.IMAX1) THEN
                  P_FORCE(IJK,1) = MPG_CYCLIC(1)*DX(I)*DY(J)*DZ(K ) 
               ELSEIF(I.EQ.IMIN1) THEN
                  TEMP2 = (P_G(IJK)+P_G(IPJK))/2
                  TEMP1 = (P_G(IJK)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
                  P_FORCE(IJK,1) = (TEMP1-TEMP2)*DY(J)*DZ(K )  +  MPG_CYCLIC(1)*DX(I)*DY(J)*DZ(K )
               ELSEIF(I.EQ.IMAX1) THEN
                  TEMP2 = (P_G(IMJK)+P_G(IJK))/2
                  TEMP1 = (P_G(IJK)*DX(I) - TEMP2*DX(I)/2)/(DX(I)/2)
                  P_FORCE(IJK,1) = (TEMP2 - TEMP1)*DY(J)*DZ(K) +  MPG_CYCLIC(1)*DX(I)*DY(J)*DZ(K )
               ELSEIF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
                  TEMP2 = (P_G(IJK)+P_G(IPJK))/2
                  TEMP1 = (P_G(IMJK)+P_G(IJK))/2
                  P_FORCE(IJK,1) = (TEMP1 - TEMP2)*DY(J)*DZ(K) + MPG_CYCLIC(1)*DX(I)*DY(J)*DZ(K )
               ENDIF
               IF(JMIN1.EQ.JMAX1) THEN
                  P_FORCE(IJK,2) = MPG_CYCLIC(2)*DX(I)*DY(J)*DZ(K )
               ELSEIF(J.EQ.JMIN1) THEN
                  TEMP2 = (P_G(IJK)+P_G(IJPK))/2
                  TEMP1 = (P_G(IJK)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
                  P_FORCE(IJK,2) = (TEMP1 - TEMP2)*DX(I)*DZ(K) + MPG_CYCLIC(2)*DX(I)*DY(J)*DZ(K )
               ELSEIF(J.EQ.JMAX1) THEN
                  TEMP2 = (P_G(IJMK)+P_G(IJK))/2
                  TEMP1 = (P_G(IJK)*DY(J) - TEMP2*DY(J)/2)/(DY(J)/2)
                  P_FORCE(IJK,2) = (TEMP2 - TEMP1)*DX(I)*DZ(K) + MPG_CYCLIC(2)*DX(I)*DY(J)*DZ(K )
               ELSEIF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
                  TEMP2 = (P_G(IJK)+P_G(IJPK))/2
                  TEMP1 = (P_G(IJMK)+P_G(IJK))/2
                  P_FORCE(IJK,2) = (TEMP1 - TEMP2)*DX(I)*DZ(K) +  MPG_CYCLIC(2)*DX(I)*DY(J)*DZ(K )
               ENDIF
               IF(DIMN.EQ.3) THEN
                  IF(KMIN1.EQ.KMAX1) THEN
                     P_FORCE(IJK,3) = MPG_CYCLIC(3)*DX(I)*DY(J)*DZ(K )
                  ELSEIF(K.EQ.KMIN1) THEN
                     TEMP2 = (P_G(IJK)+P_G(IJKP))/2
                     TEMP1 = (P_G(IJK)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
                     P_FORCE(IJK,3) = (TEMP1 - TEMP2)*DX(I)*DY(J) +  MPG_CYCLIC(3)*DX(I)*DY(J)*DZ(K )
                  ELSEIF(K.EQ.KMAX1) THEN
                     TEMP2 = (P_G(IJKM)+P_G(IJK))/2
                     TEMP1 = (P_G(IJK)*DZ(K) - TEMP2*DZ(K)/2)/(DZ(K)/2)
                     P_FORCE(IJK,3) = (TEMP2 - TEMP1)*DX(I)*DY(J) +  MPG_CYCLIC(3)*DX(I)*DY(J)*DZ(K )
                  ELSEIF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
                     TEMP2 = (P_G(IJK)+P_G(IJKP))/2
                     TEMP1 = (P_G(IJKM)+P_G(IJK))/2
                     P_FORCE(IJK,3) = (TEMP1 - TEMP2)*DX(I)*DY(J) +  MPG_CYCLIC(3)*DX(I)*DY(J)*DZ(K )
                  ENDIF
               ENDIF

               IF (.NOT. DES_INTERP_ON) THEN
! average fluid velocity at scalar cell center
                  UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
                  VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
                  DO M = 1, DES_MMAX
                     EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
                     IF(EP_SM.GT.ZERO) THEN
                        SOLID_DRAG(IJK,M,1) = -F_GS(IJK,M)*&
                           (DES_U_S(IJK,M)-UGC)
                        SOLID_DRAG(IJK,M,2) = -F_GS(IJK,M)*&
                           (DES_V_S(IJK,M)-VGC)
                        IF(DIMN.EQ.3) THEN
                           WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
                           SOLID_DRAG(IJK,M,3) = -F_GS(IJK,M)*&
                              (DES_W_S(IJK,M)-WGC)
                        ENDIF
                        OEPS = ONE/EP_SM
                        SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF      ! end if(pinc(ijk).gt.0)
         ENDDO         ! end do loop (ijk=ijkstart3,ijkend3)
      ENDIF            ! end if(calc_fc)
!----------------------------------------------- 


! update the contact forces (FC) on the particle to include gas pressure
! and gas-solids drag (this section is performed when drag is not interpolated -
! otherwise code further down is used for updating FC)
!-----------------------------------------------
      IF(CALC_FC .AND. .NOT.DES_INTERP_ON) THEN
         PC = 1              
         DO NP = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF (.NOT.PEA(NP,1)) CYCLE
            IJK = PIJK(NP,4)
            M = PIJK(NP,5)
            OVOL = ONE/VOL(IJK)
            FC(NP,:) = FC(NP,:) + SOLID_DRAG(IJK,M,:)*PVOL(NP)
            IF(MODEL_B) THEN    !Do not add the pressure gradient force
            ELSE !Add the pressure gradient force 
               FC(NP,:) = FC(NP,:) + (P_FORCE(IJK,:)*OVOL)*PVOL(NP)
            ENDIF
            PC = PC + 1
         ENDDO
      ENDIF
!----------------------------------------------- 


! this section is used to calculate the gas solids drag force on each particle
! using particle velocity and the fluid velocity interpolated to particle
! position
!-----------------------------------------------      
      IF(DES_INTERP_ON) THEN 

! when periodic boundaries are used some care needs to be taken when
! the velocity adjacent to a ghost cell is interpolated as the velocity
! in a ghost cell is undefined. therefore the velocity in ghost cells
! along the west, south and bottom edges is temporarily set to the 
! velocity of the fluid cell adjacent to the opposing boundary so that
! the velocity will be interpolated correctly.  the velocity in ghost
! cells along the east, north and top edges does not need to be set
! since these particular cells are never explictly indexed in the
! interpolation routines.  the velocities are reset at the end.
         IF (DES_PERIODIC_WALLS) THEN
            IF(DES_PERIODIC_WALLS_X) THEN 
               I = 1   ! ghost cell
               DO J = 1, JMAX2
                  DO K = 1, KMAX2
                     IJK = funijk(I,J,K)
                     IJK_IMAX1 = funijk(IMAX1,J,K)   ! imax1=fluid cell
                     TEMP_U_G_X(J,K) = U_G(IJK) 
                     U_G(IJK) = U_G(IJK_IMAX1)
                     TEMP_V_G_X(J,K) = V_G(IJK) 
                     V_G(IJK) = V_G(IJK_IMAX1)
                     IF(DIMN.EQ.3) THEN
                        TEMP_W_G_X(J,K) = W_G(IJK) 
                        W_G(IJK) = W_G(IJK_IMAX1)
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
            IF(DES_PERIODIC_WALLS_Y) THEN
               J = 1   ! ghost cell
               DO I = 1, IMAX2
                  DO K = 1, KMAX2
                     IJK = funijk(I,J,K)
                     IJK_JMAX1 = funijk(I,JMAX1,K)   ! jmax1=fluid cell
                     TEMP_U_G_Y(I,K) = U_G(IJK) 
                     U_G(IJK) = U_G(IJK_JMAX1)
                     TEMP_V_G_Y(I,K) = V_G(IJK) 
                     V_G(IJK) = V_G(IJK_JMAX1)
                     IF(DIMN.EQ.3) THEN 
                        TEMP_W_G_Y(I,K) = W_G(IJK) 
                        W_G(IJK) = W_G(IJK_JMAX1)
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF 
            IF(DES_PERIODIC_WALLS_Z .AND. DIMN .EQ. 3) THEN
               K = 1   ! ghost cell
               DO J = 1, JMAX2
                  DO I = 1, IMAX2
                     IJK = funijk(I,J,K)
                     IJK_KMAX1 = funijk(I,J,KMAX1)   ! kmax1=fluid cell
                     TEMP_U_G_Z(I,J) = U_G(IJK) 
                     U_G(IJK) = U_G(IJK_KMAX1)
                     TEMP_V_G_Z(I,J) = V_G(IJK) 
                     V_G(IJK) = V_G(IJK_KMAX1)
                     TEMP_W_G_Z(I,J) = W_G(IJK) 
                     W_G(IJK) = W_G(IJK_KMAX1)
                  ENDDO
               ENDDO
            ENDIF 
         ENDIF   ! end if des_periodic_walls


! sets several quantities including interp_scheme, scheme, and order and
! allocates arrays necessary for interpolation
         CALL SET_INTERPOLATION_SCHEME(2)

         drag_am = ZERO
         drag_bm = ZERO
         WTBAR = zero

         PC = 1
         DO NP = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF (.NOT.PEA(NP,1)) CYCLE
            FOCUS = .FALSE.
        
            I = PIJK(NP, 1)
            J = PIJK(NP, 2)
            K = PIJK(NP, 3)

! generally a particle may not exist in a ghost cell. however, if the
! particle is adjacent to the west, south or bottom boundary, then pcell
! may be assigned indices of a ghost cell which will be passed to
! set_interpolation_stencil
            PCELL(1) = I-1
            PCELL(2) = J-1
            PCELL(3) = (3-DIMN)*1+(DIMN-2)*(K-1)  ! =K-1 (in 3D) or =1 (in 2D)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
            CALL SET_INTERPOLATION_STENCIL(PCELL, IW, IE, JS, JN, KB,&
               KTP, INTERP_SCHEME, DIMN, ORDERNEW = ONEW) 

            IF(NP.EQ.FOCUS_PARTICLE) THEN 
               FOCUS = .TRUE.
               PRINT*, 'PCELL : I-1, J-1, K-1 = ',&
                  PCELL(1), PCELL(2), PCELL(3)
               PRINT*, 'IW, IE, JS, JN, KB, KTP = ', &
                  IW, IE, JS, JN, KB, KTP
            ENDIF

! compute the fluid velocity interpolated to the particle position
! (VEL_FP)
            CALL INTERPOLATE_QUANTS(IW, IE, JS, JN, KB, KTP, ONEW, &
               DES_POS_NEW(NP,1:DIMN), VEL_FP(NP,1:DIMN), FOCUS)
 

! calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and gas velocity interpolated to the particle position
            CALL DES_DRAG_GS(NP, VEL_FP(NP,1:DIMN), &
               DES_VEL_NEW(NP,1:DIMN))


! update the contact forces (FC) to include gas pressure and gas-solids
! drag forces on the particle  
            IF(CALC_FC) THEN
               IJK = PIJK(NP, 4)
               OVOL = ONE/VOL(IJK)
               D_FORCE(1:DIMN) = F_GP(NP)*&
                  ( VEL_FP(NP,1:DIMN)-DES_VEL_NEW(NP,1:DIMN) )
               FC(NP,:) = FC(NP,:) + D_FORCE(:)
               
            IF(MODEL_B) THEN    !Do not add the pressure gradient force
            ELSE !Add the pressure gradient force 
               FC(NP,:) = FC(NP,:) + (P_FORCE(IJK,:)*OVOL)*PVOL(NP)
            ENDIF
            ENDIF


! if callfromdes is true, then the pertinent mean fields (in this case
! ROP_S and F_GS) are not computed/updated in this call. this is done to
! speed up the simulation.  so the following section will only be called
! at the end of a given DEM simulation for a given fluid time step
! (i.e., only called once per fluid time step)
            IF(.NOT.CALLFROMDES) THEN 
               M = PIJK(NP,5)

               DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
                  DO J = 1, ONEW
                     DO I = 1, ONEW

! Shift loop index to new variables for manipulation                     
                        II = IW + I-1
                        JJ = JS + J-1
                        KK = KB + K-1

! adjust for periodic boundaries
                        IF (DES_PERIODIC_WALLS) THEN
                           IF (DES_PERIODIC_WALLS_X) THEN
                              IF(II.LT.1)     II = IMAX1+II-1
                              IF(II.GT.IMAX1) II = II-IMAX1+1
                           ENDIF
                           IF (DES_PERIODIC_WALLS_Y) THEN
                              IF(JJ.LT.1)     JJ = JMAX1+JJ-1
                              IF(JJ.GT.JMAX1) JJ = JJ-JMAX1+1
                           ENDIF
                           IF (DIMN.EQ.3 .AND. DES_PERIODIC_WALLS_Z) THEN
                              IF(KK.LT.1)     KK = KMAX1+KK-1
                              IF(KK.GT.KMAX1) KK = KK-KMAX1+1
                           ENDIF
                        ENDIF

! should this volume be for current ijk index or always particle index?                        
                        VCELL = VOL(PIJK(NP,4))
! todo: adjust this for non-uniform cell sizes
                        IF(II.EQ.1.or.II.EQ.IMAX1) VCELL = 0.5d0*VCELL
                        IF(JJ.EQ.1.or.JJ.EQ.JMAX1) VCELL = 0.5d0*VCELL
                        IF(DIMN.EQ.3)THEN
                           IF(KK.EQ.1.or.KK.EQ.KMAX1) VCELL = 0.5d0*VCELL
                        ENDIF

                        OVOL = ONE/VCELL
                  
                        drag_am(II,JJ,KK,M) = drag_am(II,JJ,KK,M) + &
                           F_GP(NP) * WEIGHTP(I,J,K)*OVOL
                        !drag_am(II,JJ,KK,M) = zero

! first remove the velocity component at this grid point from the vel_fp
                        drag_bm_tmp(1:DIMN) = VEL_FP(NP, 1:DIMN) - &
                           WEIGHTP(I,J,K)*vstencil(I,J,K, 1:DIMN)
! now find the remaning drag force
                        drag_bm_tmp(1:DIMN) = DES_VEL_NEW(NP,1:DIMN) !- DRAG_BM_TMP(1:DIMN)
                  
                        drag_bm(II,JJ,KK, 1:DIMN,M) = &
                           drag_bm(II,JJ,KK,1:DIMN,M) + &
                           F_GP(NP) * drag_bm_tmp(1:DIMN) * &
                           WEIGHTP(I,J,K)*OVOL
                        !drag_bm(II,JJ,KK, 1:DIMN,M) = zero 

                        WTBAR(II,JJ,KK,M) = WTBAR(II,JJ,KK,M) + &
                           WEIGHTP(I,J,K) *RO_Sol(NP)*OVOL*PVOL(NP)

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF          ! end if(.not.callfromdes)
            PC = PC + 1
         ENDDO             ! end loop over NP = 1, PARTICLES


         IF(.NOT.CALLFROMDES) THEN 

! adjust for periodic boundaries
            IF (DES_PERIODIC_WALLS) THEN
               IF (DES_PERIODIC_WALLS_X) THEN
                  I = 1
                  DO K = 1, KMAX1
                     DO J = 1, JMAX1
                        drag_bm(I,J,K,1, 1:DES_MMAX) = HALF*&
                           ( drag_bm(I,J,K,1, 1:DES_MMAX) + &
                             drag_bm(IMAX1,J,K,1, 1:DES_MMAX) )
                        drag_bm(IMAX1,J,K,1,1:DES_MMAX) = drag_bm(1,J,K,1, 1:DES_MMAX)
                        
                        WTBAR(I,J,K,:) = HALF*&
                           ( WTBAR(I,J,K,:) + WTBAR(IMAX1,J,K,:) )
                        WTBAR(IMAX1,J,K,:) = WTBAR(1,J,K,:)
   
                        drag_am(I,J,K,:) = HALF*&
                           ( drag_am(I,J,K,:) + drag_am(IMAX1,J,K, :) )
                        drag_am(IMAX1,J,K,:) = drag_am(1,J,K,:)
                     ENDDO
                  ENDDO
               ENDIF
               IF(DES_PERIODIC_WALLS_Y) THEN
                  J = 1
                  DO K = 1, KMAX1
                     DO I = 1, IMAX1
                        drag_bm(I,J,K,2, 1:DES_MMAX) = HALF*&
                           ( drag_bm(I,J,K,2, 1:DES_MMAX)+&
                             drag_bm(I,JMAX1,K,2, 1:DES_MMAX) )
                        drag_bm(I,JMAX1,K,2, 1:DES_MMAX) = drag_bm(I,1,K,2, 1:DES_MMAX)
   
                        WTBAR(I,J,K,:) = HALF*&
                           ( WTBAR(I,J,K,:) + WTBAR(I,JMAX1,K,:) )
                        WTBAR(I,JMAX1,K,:) = WTBAR(I,1,K,:)
                  
                        drag_am(I,J,K,:) = HALF*&
                           ( drag_am(I,J,K,:) + drag_am(I,JMAX1,K,:) )
                        drag_am(I,JMAX1,K,:) = drag_am(I,1,K,:)
                     ENDDO
                  ENDDO
               ENDIF
               IF(DES_PERIODIC_WALLS_Z .AND. DIMN .EQ. 3) THEN
                  K = 1
                  DO J = 1, JMAX1
                     DO I = 1, IMAX1
                        drag_bm(I,J,K,3, 1:DES_MMAX) = HALF*&
                           ( drag_bm(I,J,K,3, 1:DES_MMAX) + &
                             drag_bm(I,J,KMAX1,3, 1:DES_MMAX) )
                        drag_bm(I,J,KMAX1,3,1:DES_MMAX) = drag_bm(I,J,1,3,1:DES_MMAX)
   
                        WTBAR(I,J,K,:) = HALF*&
                           ( WTBAR(I,J,K,:) + WTBAR(I,J,KMAX1,:) )
                        WTBAR(I,J,KMAX1,:) = WTBAR(I,J,1,:)
   
                        drag_am(I,J,K,:) = HALF*&
                           ( drag_am(I,J,K,:) + drag_am(I,J,KMAX1,:) )
                        drag_am(I,J,KMAX1,:) = drag_am(I,J,1,:)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

            DO M = 1, DES_MMAX
               DO IJK = IJKSTART3, IJKEND3
                  IF(FLUID_AT(IJK)) THEN
                     I = I_of(IJK)
                     J = J_of(IJK)
                     K = K_of(IJK)
               
                     F_GS(IJK,M) = AVG_FACTOR*(drag_am(I,J,K,M) +&
                        drag_am(I,J-1,K,M) + drag_am(I-1,J-1,K,M) +&
                        drag_am(I-1,J,K,M))
               
                     DES_ROP_S(IJK,M) = AVG_FACTOR*(WTBAR(I,J,K,M) +&
                        WTBAR(I,J-1,K,M) + WTBAR(I-1,J-1,K,M) +&
                        WTBAR(I-1,J,K,M))
               
                     IF(DIMN.EQ.3) THEN 
                        F_GS(IJK,M) = F_GS(IJK,M) + AVG_FACTOR*&
                           (drag_am(I,J,K-1,M) + drag_am(I,J-1,K-1,M) +&
                           drag_am(I-1,J-1,K-1,M)+drag_am(I-1,J,K-1,M) )
                  
                        DES_ROP_S(IJK,M) = DES_ROP_S(IJK,M) + AVG_FACTOR*&
                           (WTBAR(I,J,K-1,M) + WTBAR(I,J-1,K-1,M) + &
                           WTBAR(I-1,J-1,K-1,M)+WTBAR(I-1,J,K-1,M) )
                     ENDIF
                  ENDIF
               ENDDO  ! end do loop over ijk
            ENDDO   ! end do loop over m=1,DES_MMAX

         ENDIF        ! end if(.not.callfromdes)


! reset the velocities which were previously adjusted for periodic boundaries
         IF (DES_PERIODIC_WALLS) THEN      
            IF(DES_PERIODIC_WALLS_X) THEN 
               I = 1
               DO J = 1, JMAX2
                  DO K = 1, KMAX2
                     IJK = funijk(I,J,K)
                     U_G(IJK) = TEMP_U_G_X(J,K)
                     V_G(IJK) = TEMP_V_G_X(J,K)
                     IF(DIMN.EQ.3) W_G(IJK) = TEMP_W_G_X(J,K)
                  ENDDO
               ENDDO
            ENDIF
            IF(DES_PERIODIC_WALLS_Y) THEN 
               J = 1
               DO K = 1, KMAX2
                  DO I = 1, IMAX2
                     IJK = funijk(I,J,K)
                     U_G(IJK) = TEMP_U_G_Y(I,K)
                     V_G(IJK) = TEMP_V_G_Y(I,K)
                     IF(DIMN.EQ.3) W_G(IJK) = TEMP_W_G_Y(I,K)
                  ENDDO
               ENDDO
            ENDIF 
            IF(DES_PERIODIC_WALLS_Z.AND.DIMN.EQ.3) THEN 
               K = 1
               DO J = 1, JMAX2
                  DO I = 1, IMAX2
                     IJK = funijk(I,J,K)
                     U_G(IJK) = TEMP_U_G_Z(I,J)
                     V_G(IJK) = TEMP_V_G_Z(I,J)
                     W_G(IJK) = TEMP_W_G_Z(I,J)
                  ENDDO
               ENDDO
            ENDIF 
         ENDIF   ! end if des_periodic_walls

      ENDIF        ! end if(des_interp_on)
!-----------------------------------------------  


      RETURN
      END SUBROUTINE DRAG_FGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      SUBROUTINE INTERPOLATE_QUANTS(IW, IE, JS, JN, KB, KTP, ONEW, &
         POSP, FVEL, FOCUS)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      USE param
      USE param1
      USE constant
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE compar
      USE discretelement
      USE interpolation
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!----------------------------------------------- 
! starting grid indices used in interpolation 
! appears that only IW, JS, and KB are necessary
      INTEGER, INTENT(IN) :: IW, IE, JS, JN, KB, KTP
! interpolation order      
      INTEGER, INTENT(IN) :: onew
! the x,y,z position of particle
      DOUBLE PRECISION, DIMENSION(DIMN) , INTENT(IN) :: POSP
! the x,y,z velocity of the gas at the particle position (interpolated
! to particle position after this call)
      DOUBLE PRECISION, DIMENSION(DIMN) , INTENT(OUT) :: FVEL
! flag to tell whether local debugging was called in drag_fgs      
      LOGICAL :: FOCUS
! constant whose value depends on dimension of system            
      DOUBLE PRECISION :: AVG_FACTOR
! indices
      INTEGER I, J, K, II, JJ, KK, IJK
      INTEGER IPJK, IJPK, IJKP, IPJPK, IPJKP, IJPKP, &
              IPJPKP
!-----------------------------------------------  

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

! avg_factor=0.25 (in 3D) or =0.50 (in 2D)  
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)

      DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
         DO J = 1, ONEW
            DO I = 1, ONEW
               II = IW + I-1
               JJ = JS + J-1
               KK = KB + K-1
               GSTENCIL(I,J,K,1) = XE(II)
               GSTENCIL(I,J,K,2) = YN(JJ)
               GSTENCIL(I,J,K,3) = ZT(KK)*(DIMN-2) + DZ(1)*(3-DIMN)
              
               IF (DES_PERIODIC_WALLS) THEN 
                  IF (DES_PERIODIC_WALLS_X) THEN
                     IF(II.LT.1)     II = IMAX1+II
                     IF(II.GT.IMAX1) II = II-IMAX1
                  ENDIF
                  IF (DES_PERIODIC_WALLS_Y) THEN
                     IF(JJ.LT.1)     JJ = JMAX1+JJ
                     IF(JJ.GT.JMAX1) JJ = JJ-JMAX1
                  ENDIF
                  IF (DES_PERIODIC_WALLS_Z .AND. DIMN .EQ. 3) THEN
                     IF(KK.LT.1)     KK = KMAX1+KK
                     IF(KK.GT.KMAX1) KK = KK-KMAX1
                  ENDIF
               ENDIF
               
               ijk = funijk(II,JJ,KK)
               
               ipjk = ip_of (ijk)    
               ijpk = jp_of (ijk)
               ijkp = kp_of (ijk)
               ijpkp = kp_of(ijpk)
               ipjkp = kp_of(ipjk)
               ipjpk = jp_of(ipjk)
               ipjpkp = kp_of(ipjpk)
               
               vstencil(i,j,k,1) = AVG_FACTOR*( u_g(ijk) + u_g(ijpk) + &
                  (u_g(ijkp) + u_g(ijpkp)) * (DIMN-2) )
               
               vstencil(i,j,k,2) = AVG_FACTOR*( v_g(ijk) + v_g(ipjk) + &
                  (v_g(ijkp) + v_g(ipjkp)) * (DIMN-2) )
               
               IF(FOCUS) THEN 
                  print*, 'II, JJ, KK = ', ii, jj, kk
                  PRINT*, 'v_g = ', v_g(ijk), v_g(ipjk), &
                     v_g(ijkp), v_g(ipjkp), vstencil(i,j,k,2)
               ENDIF
               
               IF(DIMN.eq.3) THEN 
                  vstencil(i,j,k,3) = AVG_FACTOR*(w_g(ijk) +&
                  & w_g(ijpk) + w_g(ipjk) + w_g(ipjpk) )
               ELSE 
! doesn't matter what value is put here
                  vstencil(i,j,k,3) = 0.d0
               ENDIF               
            ENDDO
         ENDDO
      ENDDO
      
      IF(DIMN.EQ.2) THEN 
         CALL interpolator( GSTENCIL(1:onew, 1:onew, 1, 1:DIMN), &
            VSTENCIL(1:onew, 1:onew, 1, 1:2), &
            POSP(1:2), FVEL(1:2), ONEW, INTERP_SCHEME, &
            WEIGHTP )
      ELSE 
         CALL interpolator( GSTENCIL(1:onew, 1:onew, 1:onew, 1:DIMN), &
            VSTENCIL(1:onew, 1:onew, 1:onew, 1:DIMN), &
            POSP(1:3), FVEL(1:3), ONEW, INTERP_SCHEME, &
            WEIGHTP )   
      ENDIF

! local debugging
      IF(FOCUS) THEN 
         PRINT*, ' '
         PRINT*, 'IN INTERP'
         PRINT*, 'IW, IE, JS, JN, KB, KTP = ', IW, IE, JS, JN, KB, KTP

         DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
            DO J = 1, ONEW
               DO I = 1, ONEW
                  PRINT*, 'V and Weight= ',&
                     VSTENCIL(I,J,K, 2), WEIGHTP(I,J,K)
               ENDDO
            ENDDO
         ENDDO

         PRINT*, 'FVEL(Y) = ', fvel(2)
         PRINT*, 'PPOSY) = ', POSP(2) !, gstencil(1:onew,1:onew,1,1)
      ENDIF

      END SUBROUTINE INTERPOLATE_QUANTS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_DRAG_GS                                            C
!  Purpose: Calculate the gas-particle drag coefficient for des        C
!  calculation such that in drag correlation the gas velocity
!  interpolated to the particle position and the particle velocity
!  are used
!      
!  Comments: No BVK drag model in this subroutine. BVK requires an
!  average particle diameter which needs to be defined for DEM case
!  diameter which needs
!
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
! Arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) ::  KK
! fluid velocity interpolated to particle position      
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: fvel
! particle velocity
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: des_vel

!-----------------------------------------------
! Local Parameters 
!-----------------------------------------------
! Parameters in the Cluster-effect model
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
!     a1 depends upon solids flux.  It has been represented by C(1) 
!     defined in the data file.
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0 
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0 
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0 
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0 

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices 
      INTEGER I, IJK, IMJK, IJMK, IJKM, IM, M
     
! Cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION USCM, VSCM, WSCM

! Cell center value of U_g , V_g, W_g
      DOUBLE PRECISION UGC, VGC, WGC
     
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION VREL 

! Gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION Mu

! Reynolds number 
      DOUBLE PRECISION Re 

! Reynolds number that includes void fraction
      DOUBLE PRECISION Re_g
    
! Ratio of settling velocity of a multiparticle system to that of a
! single particle 
      DOUBLE PRECISION V_rm 
     
! Single sphere drag coefficient 
      DOUBLE PRECISION C_d 

! Single sphere drag coefficient x Re 
      DOUBLE PRECISION C_DsxRe, C_DsxReT 
     
! Drag coefficient 
      DOUBLE PRECISION DgA  

! CEAF: 03/23/06 Gidaspow switch function variables
      DOUBLE PRECISION Ergun
      DOUBLE PRECISION WenYu
      DOUBLE PRECISION PHI_gs

! Function of EP_g in Syam O'Brien Model
      DOUBLE PRECISION A_SO, B_SO

! total solids volume fraction
      DOUBLE PRECISION phis

! solids volume fraction of phase M in fluid cell of interest
      DOUBLE PRECISION EP_SM    
! tmp variable for particle diameter      
      DOUBLE PRECISION PART_DIAM
! tmp variable for particle volume   
      DOUBLE PRECISION PART_VOL   

! Current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp

! Koch and Hill drag correlation variables
!**************************************
! Stokes Drag Force
      DOUBLE PRECISION F_STOKES
     
! zero Re function for low Reynolds number
      DOUBLE PRECISION F_0
     
! inertial function for low Reynolds number
      DOUBLE PRECISION F_1
     
! zero Re function for high Reynolds number
      DOUBLE PRECISION F_2
     
! inertial function for high Reynolds number
      DOUBLE PRECISION F_3
     
! dimensionless drag force F
      DOUBLE PRECISION F
     
! transition Reynolds numbers
      DOUBLE PRECISION Re_Trans_1, Re_Trans_2
     
! weighting factor to compute F_0 and F_2
      DOUBLE PRECISION w, D_p_av, Y_i
     
! Hill and Koch Reynolds number
      DOUBLE PRECISION Re_kh
!**************************************
     
!-----------------------------------------------    

      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER) ! Tsuji drag
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2 ! Dalla Valle (1948) 
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0 * Re**0.657D0)      ! Turton and
!     &          + 0.413D0 * Re**2.09D0 / (Re**1.09D0 + 16300.D0) ! Levenspiel (1986)


      IJK = PIJK(KK,4)
      M = PIJK(KK,5)
      EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
      PART_DIAM = 2.D0*DES_RADIUS(KK)
      PART_VOL = PVOL(KK)

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

! Reynolds numbers
      IF(Mu > ZERO)THEN
         RE = PART_DIAM * VREL*RO_G(IJK)/Mu 
! Note the presence of gas volume fraction in ROP_G
         RE_G = (PART_DIAM * VREL*ROP_G(IJK))/Mu
! Note the presence of gas volume fraction in ROP_G and additional
! factor of 1/2
         RE_kh = (0.5D0*PART_DIAM*VREL*ROP_G(IJK))/Mu
         
      ELSE
         RE = LARGE_NUMBER 
         RE_G = LARGE_NUMBER
         RE_kh = LARGE_NUMBER
      ENDIF
!     f_gp() =  single particle drag excluding vector(v_g - v_p)
!     


!---------------Begin Syamlal and O'Brien ---------------------------

      IF(TRIM(DRAG_TYPE).EQ.'SYAM_OBRIEN') then
         IF (EP_sM <= ZERO) THEN 
            F_gstmp = ZERO 
         ELSEIF (EP_G(IJK) == ZERO) THEN 
            F_gstmp = ZERO 
         ELSE 
            A_SO = EP_G(IJK)**4.14D0 
            IF (EP_G(IJK) <= 0.85D0) THEN 
               B_SO = drag_c1*EP_G(IJK)**1.28D0 
            ELSE 
               B_SO = EP_G(IJK)**drag_d1
            ENDIF 

! Calculate V_rm
            V_RM=HALF*(A_SO - 0.06D0*RE + SQRT(3.6D0-3.D0*RE*RE + &
               0.12D0*RE*(2.D0*B_SO-A_SO) + A_SO*A_SO) ) 
!------------------Begin cluster correction --------------------------
!     uncomment the following four lines and comment the above line V_RM=... 
!     V_RM=HALF*(A_SO - 0.06D0*RE + SQRT(3.6D-3.D0*RE*RE + & 
!        0.12D0*RE*(2.D0*B_SO-A_SO) + A_SO*A_SO) ) * (ONE + &
!        C(1) * exp( -a2*(Re - Re_c)**2 - a3*(EP_g(IJK)-ep_c)**2 ) &
!        * Re * (1. - EP_g(IJK))  )
!------------------End cluster correction ----------------------------
     
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
            IF(TSUJI_DRAG) THEN
               IF(EP_G(IJK) >= 0.8D0) THEN
                  F_gstmp = (Mu*PART_VOL/(PART_DIAM**2))*&
                     (150.D0*(EP_SM/EP_G(IJK)) + 1.75D0*RE)
               ELSEIF(EP_G(IJK) > 0.8D0) THEN
                  IF(RE*EP_G(IJK) > 1000.D0) THEN
                     F_gstmp = 0.75D0*0.43D0*Mu*PART_VOL*RE/&
                        (PART_DIAM**2 * EP_G(IJK)**1.7D0)
                  ELSEIF(RE*EP_G(IJK).LE.1000D0) THEN
                     F_gstmp = 0.75D0*C_DSXRET(RE*EP_G(IJK))*Mu*&
                       PART_VOL*RE/(PART_DIAM**2 *EP_G(IJK)**1.7D0)
                  ENDIF
               ENDIF 
            ELSEIF(MODEL_B) THEN 
               F_gstmp = 0.75D0*Mu*(PART_VOL)*&
                  C_DSXRE(RE/V_RM) / (V_RM*PART_DIAM*PART_DIAM) 
            ELSE
               F_gstmp = 0.75D0*Mu*(PART_VOL)*EP_G(IJK)*&
                  C_DSXRE(RE/V_RM) / (V_RM*PART_DIAM*PART_DIAM) 
            ENDIF
         ENDIF
!---------------End Syamlal and O'Brien ---------------------------


!--------------------------Begin Gidaspow --------------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') then
         IF(EP_g(IJK) <= 0.8D0) THEN
            DgA = 150D0 * (ONE - EP_g(IJK)) * Mu / &
               (EP_g(IJK) * PART_DIAM**2) + &
               1.75D0*RO_g(IJK)*VREL/PART_DIAM
         ELSE
            IF(Re_G .LE. 1000D0)THEN
               C_d = (24.D0/(Re_G+SMALL_NUMBER)) * &
                  (ONE + 0.15D0 * Re_G**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF
            DgA = 0.75D0 * C_d*VREL*ROP_g(IJK)*EP_g(IJK)**(-2.65D0) / &
               PART_DIAM
         ENDIF
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
! F_gstmp = DgA * EP_sM/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL !EP_s(IJK, M)
         ENDIF
!--------------------------End Gidaspow --------------------------


!-----------------------Begin Gidaspow_blend ---------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND') then

! Dense phase - EP_g < 0.8
         Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu / &
            ( EP_g(IJK) * PART_DIAM**2 ) + &
            + 1.75D0 * RO_g(IJK) * VREL / PART_DIAM
     
! Dilute phase - EP_g >= 0.8
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * &
               (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         WenYu = 0.75D0 * C_d*VREL*ROP_g(IJK)*EP_g(IJK)**(-2.65D0) / &
            PART_DIAM
     
! Switch function
         PHI_gs = ATAN(150.D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
! Blend the models
         DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
! F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL !EP_s(IJK, M)
         ENDIF
!-----------------------End Gidaspow_blend -----------------------


!--------------------------Begin WEN_YU --------------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'WEN_YU') then
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         DgA = 0.75D0 * C_d*VREL*ROP_g(IJK)*EP_g(IJK)**(-2.65D0) / &
            PART_DIAM
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
! F_gstmp = DgA * EP_s(IJK, M)/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
!--------------------------End WEN_YU ----------------------------


!--------------------Begin Koch & Hill (2001) --------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') then
     
         F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)*EP_g(IJK)/PART_DIAM**2
         
         phis = ONE-EP_G(IJK)   ! EP_s(IJK,M) for polydisperse systems 
         w = EXP(-10.0D0*(0.4D0-phis)/phis)
         
         IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
            F_0 = (1.0D0-w) *                                          &
            (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis     &
            *LOG(phis) + 17.14D0*phis) / (1.0D0 + 0.681D0*             &
            phis - 8.48D0*phis*phis + 8.16D0*phis**3) + w *            &
            10.0D0*phis/(1.0D0-phis)**3
            
         ELSEIF(phis >= 0.4D0) THEN
            F_0 = 10.0D0*phis/(1.0D0-phis)**3
         ENDIF
         
         IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
            F_1 = dsqrt(2.0D0/phis) / 40.0D0
         ELSE IF(phis > 0.1D0) THEN
            F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
         ENDIF
         
         IF(phis < 0.4D0) THEN
            F_2 = (1.0D0-w) *                                          &
               (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + 135.0D0/64.0D0*phis  &
               *LOG(phis) + 17.89D0*phis) / (1.0D0 + 0.681D0*          &
               phis - 11.03D0*phis*phis + 15.41D0*phis**3)+ w *        &
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
         ELSEIF(phis > 0.01D0 .AND. Re_kh <= Re_Trans_2) THEN
            F = F_0 + F_1*Re_kh*Re_kh
         ELSEIF(phis <= 0.01D0 .AND. Re_kh > Re_Trans_1 .OR.         &
                phis >  0.01D0 .AND. Re_kh > Re_Trans_2) THEN
            F = F_2 + F_3*Re_kh
         ELSE
            F = zero
         ENDIF
         
! This is a check for phis (or eps_(ijk,m)) to be within physical range
         IF(phis <= ZERO .OR. phis > ONE) F = zero
         
         DgA = F * F_STOKES

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
!-------------------End Koch and Hill (2001)----------------------


      ELSEIF((DRAG_TYPE).EQ.'DILUTE_CASE') then
         C_d =  C_DSXRET(RE)
         DgA = (0.75D0 * C_d*VREL * RO_g(IJK))/PART_DIAM
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF

!--------------------------End Dilue case -------------------------
      ELSE
         CALL START_LOG 
!IF(.not.DMP_LOG)call open_pe_log(ier)
         if(mype == pe_io) WRITE (*, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         CALL END_LOG 
         call mfix_exit(myPE)  

      ENDIF

      F_gp(KK) = (ONE - UR_F_gs) * F_gp(KK) + UR_F_gs * F_gstmp

      END SUBROUTINE DES_DRAG_GS

    
