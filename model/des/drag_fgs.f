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
!  Revision Number 4                                  Date: 31-jan-11  C
!  Authour: Pradeep G
!  Revision: 1.Modified drag_interpolation routines to change three    C
!            dimensional array (drag_am,drag_bm)into one dimensional   C
!            array (IJK)                                               C
!            2.Modified treatment for periodic boundary conditions     C
!            3.modified the loop structure from particle to IJK        C
!            4.added volume at node to include effect of volume in to  C
!              the backward interpolation                              C
!            
!  Purpose: Now the drag_fgs routine is called from calc_drag in model 
!  directory as well as by calc_forces_des. Calling arguments have     
!  also changed. Depending on the choice, once can obtain drag force   
!  based on local velocities or averaged velocities
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      
      SUBROUTINE COMPUTE_PG_GRAD_CG
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      implicit none 
      
      ! general i, j, k indices
      INTEGER I, J, K, IJK, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      INTEGER X_COUNT, Y_COUNT, Z_COUNT
      
      ! temporary variables used to calculate pressure at scalar cell edge      
      DOUBLE PRECISION TEMP1, TEMP2

!mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION MPG_CYCLIC(DIMN)
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      MPG_CYCLIC(1:DIMN) = ZERO 
      
      IF(CYCLIC_X_PD) MPG_CYCLIC(1) = DELP_X/XLENGTH
      IF(CYCLIC_Y_PD) MPG_CYCLIC(2) = DELP_Y/YLENGTH
      IF(CYCLIC_Z_PD.AND.DIMN.EQ.3) MPG_CYCLIC(3) = DELP_Z/ZLENGTH

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         P_FORCE(IJK, :) = ZERO

         IF(.NOT.FLUID_AT(IJK).OR..NOT.IS_ON_myPE_owns(I, J, K)) CYCLE 

         X_COUNT = 0
         Y_COUNT = 0
         Z_COUNT = 0
         
         IJKE = IP_OF(IJK)
         IJKW = IM_OF(IJK)
         IJKN = JP_OF(IJK)
         IJKS = JM_OF(IJK)
         IJKT = KP_OF(IJK)
         IJKB = KM_OF(IJK)

         IF(FLUID_AT(IJKE)) THEN 
            X_COUNT = X_COUNT + 1
            P_FORCE(IJK, 1) = P_FORCE(IJK, 1) + 2.d0*(P_G(IJKE) - P_G(IJK))/(DX(I) + DX(I_OF(IJKE)))
         ENDIF
            
         IF(FLUID_AT(IJKW)) THEN 
            X_COUNT = X_COUNT + 1
            P_FORCE(IJK, 1) = P_FORCE(IJK, 1) + 2.d0*(P_G(IJK) - P_G(IJKW))/(DX(I) + DX(I_OF(IJKW)))
         ENDIF

         X_COUNT = MAX(1, X_COUNT) !to prvent division from zero 
         P_FORCE(IJK, 1) = MPG_CYCLIC(1) - P_FORCE(IJK,1)/REAL(X_COUNT) 
         !P_FORCE (by convention) is stored as -dp/dx. MPG_CYCLIC is already -dp/dx. 
         !therefore, P_force is multiplied by "-" for consistency

         
         IF(FLUID_AT(IJKN)) THEN 
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) + 2.d0*(P_G(IJKN) - P_G(IJK))/(DY(J) + DY(J_OF(IJKN)))
         ENDIF
            
         IF(FLUID_AT(IJKS)) THEN 
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) + 2.d0*(P_G(IJK) - P_G(IJKS))/(DY(J) + DY(J_OF(IJKS)))
         ENDIF
         Y_COUNT = MAX(1, Y_COUNT) !to prvent division from zero 

         P_FORCE(IJK, 2) = MPG_CYCLIC(2) - P_FORCE(IJK,2)/REAL(Y_COUNT) 


         IF(DIMN.eq.3) then 
            IF(FLUID_AT(IJKT)) THEN 
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) + 2.d0*(P_G(IJKT) - P_G(IJK))/(DZ(K) + DZ(K_OF(IJKT)))
            ENDIF
            
            IF(FLUID_AT(IJKB)) THEN 
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) + 2.d0*(P_G(IJK) - P_G(IJKB))/(DZ(K) + DZ(K_OF(IJKB)))
            ENDIF
            Z_COUNT = MAX(1, Z_COUNT) !to prvent division from zero 

            P_FORCE(IJK, 3) = MPG_CYCLIC(3) - P_FORCE(IJK,3)/REAL(Z_COUNT) 
         ENDIF
         
      ENDDO
      END SUBROUTINE COMPUTE_PG_GRAD_CG

      SUBROUTINE COMPUTE_PG_GRAD
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE cutcell 
      implicit none 
      
      ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
      
      ! temporary variables used to calculate pressure at scalar cell edge      
      DOUBLE PRECISION TEMP1, TEMP2

!mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION MPG_CYCLIC(DIMN)
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      IF(CARTESIAN_GRID) THEN 
         CALL COMPUTE_PG_GRAD_CG
         RETURN
      ENDIF

      MPG_CYCLIC(1:DIMN) = ZERO 
      
      IF(CYCLIC_X_PD) MPG_CYCLIC(1) = DELP_X/XLENGTH
      IF(CYCLIC_Y_PD) MPG_CYCLIC(2) = DELP_Y/YLENGTH
      IF(CYCLIC_Z_PD.AND.DIMN.EQ.3) MPG_CYCLIC(3) = DELP_Z/ZLENGTH

      
      DO IJK = IJKSTART3, IJKEND3
         P_FORCE(IJK, :) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE 

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)

         IF(IMIN1.EQ.IMAX1) THEN
            P_FORCE(IJK,1) = MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMIN1) THEN
            TEMP2 = AVG_X(P_G(IJK), P_G(IPJK), I) 
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,1) = 2.d0*(TEMP1-TEMP2)/DX(I)  +  MPG_CYCLIC(1)
         ELSEIF(I.EQ.IMAX1) THEN
            TEMP2 = AVG_X(P_G(IMJK), P_G(IJK), I-1)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,1) = 2.d0*(TEMP2 - TEMP1)/DX(I) +  MPG_CYCLIC(1)
         ELSEIF((I.GT.IMIN1).AND.(I.LT.IMAX1)) THEN
            TEMP2 = AVG_X(P_G(IJK),  P_G(IPJK), I)
            TEMP1 = AVG_X(P_G(IMJK), P_G(IJK),  I-1)
            P_FORCE(IJK,1) = (TEMP1 - TEMP2)/DX(I) + MPG_CYCLIC(1)
         ENDIF

         IF(JMIN1.EQ.JMAX1) THEN
            P_FORCE(IJK,2) = MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMIN1) THEN
            TEMP2 = AVG_Y(P_G(IJK), P_G(IJPK), J)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,2) = 2.d0*(TEMP1 - TEMP2)/DY(J)  + MPG_CYCLIC(2)
         ELSEIF(J.EQ.JMAX1) THEN
            TEMP2 = AVG_Y(P_G(IJMK), P_G(IJK), J-1)
            TEMP1 = P_G(IJK)
            P_FORCE(IJK,2) = 2.d0*(TEMP2 - TEMP1)/DY(J) + MPG_CYCLIC(2)
         ELSEIF((J.GT.JMIN1).AND.(J.LT.JMAX1)) THEN
            TEMP2 = AVG_Y(P_G(IJK),  P_G(IJPK), J)
            TEMP1 = AVG_Y(P_G(IJMK), P_G(IJK),  J-1)
            P_FORCE(IJK,2) = (TEMP1 - TEMP2)/DY(J) +  MPG_CYCLIC(2)
         ENDIF

         IF(DIMN.EQ.3) THEN
            IF(KMIN1.EQ.KMAX1) THEN
               P_FORCE(IJK,3) = MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMIN1) THEN
               TEMP2 = AVG_Z(P_G(IJK), P_G(IJKP), K)
               TEMP1 = P_G(IJK)
               P_FORCE(IJK,3) = 2.d0*(TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF(K.EQ.KMAX1) THEN
               TEMP2 = AVG_Z(P_G(IJKM), P_G(IJK), K-1)
               TEMP1 = P_G(IJK)
               P_FORCE(IJK,3) = 2.d0*(TEMP2 - TEMP1)/DZ(K) +  MPG_CYCLIC(3)
            ELSEIF((K.GT.KMIN1).AND.(K.LT.KMAX1)) THEN
               TEMP2 = AVG_Z(P_G(IJK),  P_G(IJKP), K)
               TEMP1 = AVG_Z(P_G(IJKM), P_G(IJK),  K-1)
               P_FORCE(IJK,3) = (TEMP1 - TEMP2)/DZ(K) +  MPG_CYCLIC(3)
            ENDIF
         ENDIF
      ENDDO         ! end do loop over ijk
    END SUBROUTINE COMPUTE_PG_GRAD

    SUBROUTINE DES_CALC_PART_DRAG_FORCE_INTERP_OFF
      
      USE param
      USE param1
      USE parallel
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      use desmpi 
      USE cutcell 
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------         
! local variable used for debugging
      LOGICAL FOCUS 

      DOUBLE PRECISION DRAG_FORCE(DIMN)

! temporary variables used to calculate pressure at scalar cell edge      
      DOUBLE PRECISION TEMP1, TEMP2

! average fluid velocity in x, y, z direction at scalar cell center      
      DOUBLE PRECISION UGC, VGC, WGC

          
! indices used with periodic boundaries
      INTEGER IJK_IMAX1, IJK_JMAX1, IJK_KMAX1

! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
              IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK, &
! Pradeep added following 
              IMJMK,IMJKM,IJMKM,IMJMKM


! index of solid phase that particle NP belongs to      
      INTEGER M

! one over the solids volume fraction and one over the volume       
      DOUBLE PRECISION OEPS, OVOL 

! particle number index, used for looping      
      INTEGER NP

! index to track accounted for particles 
      INTEGER PC 

<<<<<<< drag_fgs.f
! for error messages      
      INTEGER IER

! Statistical weight of the particle. Equal to one for DEM 
   
      DOUBLE PRECISION WTP, EPS
=======
! for error messages      
      INTEGER IER
>>>>>>> 1.20.2.9

<<<<<<< drag_fgs.f
      double precision  VELG_ARR(DIMN), VELS_ARR(DIMN, MMAX)
=======
! Statistical weight of the particle. Equal to one for DEM 
   
      DOUBLE PRECISION WTP, EPS

      double precision  VELG_ARR(DIMN), VELS_ARR(DIMN, MMAX)
>>>>>>> 1.20.2.9

! see the discussion for IJK_U ..... in comments       
      INTEGER  IJK_U, IJK_V, IJK_W, ICUR, JCUR, KCUR 
      
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, velf_part(dimn), ymid

      
!
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      DO M = 1, SMAX
         !compute the F_gs field with the latest average fluid and solid 
         !velocity fields 
      CALL DRAG_GS (M, IER)
      ENDDO

      IF(CALC_FC) THEN 
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IJK_U = IMJK 
            IJK_V = IJMK
            IJK_W = IJKM
!rahul: UGC (and likewise for VGC and WGC) are computed at the
!center of the scalar cell. The center of the Ith scalar cell 
!can also be thought of as east face of the (I-1)th U- cell. It is 
!better to think in terms of U, V, W grids/cell as the interpolation
!arays (like Theta_Ue_bar) for the case of cut-cell are based 
!on respective grids, i.e., theta_Ue_bar is based on U- grid. See 
!conv_diff_u_g for an example of this usage. 
!For U uncut grid, the average at the east face of IJK_U will be
!U_AVG(at East face of IJK_U) = HALF*(U(IJK_U) + U(IP_OF(IJK_U)))
!It can be verified that the above formula and old formula of 
!U_AVG(center of IJK) = HALF*(U(IJK) + U(IM_OF(IJK))) are identical
!since IJK_U = IM_OF(IJK)
!

! Pradeep no changes to indices (imin,imax) is required for parallel implementation
            IF(PINC(IJK).GT.0) THEN

! average fluid velocity at scalar cell center
               IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
                  UGC = (Theta_Ue_bar(IJK_U)*U_G(IJK_U) + Theta_Ue(IJK_U)*U_G(IP_OF(IJK_U)))
               ELSE 
                  UGC = HALF * (U_G(IJK_U) + U_G(IP_OF(IJK_U)))
               ENDIF
                  
               IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
                  VGC = (Theta_Vn_bar(IJK_V)*V_G(IJK_V) + Theta_Vn(IJK_V)*V_G(JP_OF(IJK_V)))
               ELSE
                  VGC = HALF * (V_G(IJK_V) + V_G(JP_OF(IJK_V)))
               ENDIF

               IF(DIMN.EQ.3) THEN
                  IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                     WGC = (Theta_Wt_bar(IJK_W)*W_G(IJK_W) + Theta_Wt(IJK_W) * W_G(KP_OF(IJK_W)))
                  ELSE 
                     WGC = HALF * (W_G(IJK_W) + W_G(KP_OF(IJK_W)))
                  ENDIF
               ENDIF
<<<<<<< drag_fgs.f
<<<<<<< drag_fgs.f

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
=======
                  !original terms 
                  !UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
                  !VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
                  !WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
                  
               VELG_ARR(1) = UGC
               VELG_ARR(2) = VGC
               IF(DIMN.eq.3) VELG_ARR(3) = WGC 
                  
               VELS_ARR(1,:) = DES_U_S(IJK, :)
               VELS_ARR(2,:) = DES_V_S(IJK, :)
               IF(DIMN.eq.3) VELS_ARR(3,:) = DES_W_S(IJK, :)
               DO M = 1, MMAX
                  IF(EP_S(IJK,M).GT.ZERO) THEN
                     SOLID_DRAG(IJK,M,1) = -F_GS(IJK,M)*&
                     (DES_U_S(IJK,M)-UGC)
                     SOLID_DRAG(IJK,M,2) = -F_GS(IJK,M)*&
                     (DES_V_S(IJK,M)-VGC)
                     
                     IF(DIMN.EQ.3) THEN
                        SOLID_DRAG(IJK,M,3) = -F_GS(IJK,M)*&
                        (DES_W_S(IJK,M)-WGC)
                     ENDIF
                     OEPS = ONE/EP_S(IJK,M)
                     SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
                  ENDIF
                     
                  !rahul: temp start 
                  IF(MPPIC) THEN 
                     IF(EP_S(IJK,M).GT.ZERO) THEN
                        
                        if(MPPIC_PDRAG_IMPLICIT) THEN 
                           
                        !implicit treatment for drag term 
                           SOLID_DRAG(IJK, M, :) = F_GS(IJK,M)*VELG_ARR(:)
                           
                        ELSE
                        !explicit treatment 
                           SOLID_DRAG(IJK, M, :) = F_GS(IJK,M)*(VELG_ARR(:)-VELS_ARR(:, M))
                        endif
                     
                        EPs = MIN(EP_s(IJK,M), 1.d0-EP_STAR)
                        !EPs = EP_s(IJK,M)
                        OEPS = ONE/EPs
>>>>>>> 1.20.2.9
=======
                  !original terms 
                  !UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
                  !VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
                  !WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
                  
               VELG_ARR(1) = UGC
               VELG_ARR(2) = VGC
               IF(DIMN.eq.3) VELG_ARR(3) = WGC 
                  
               VELS_ARR(1,:) = DES_U_S(IJK, :)
               VELS_ARR(2,:) = DES_V_S(IJK, :)
               IF(DIMN.eq.3) VELS_ARR(3,:) = DES_W_S(IJK, :)
               DO M = 1, MMAX
                  IF(EP_S(IJK,M).GT.ZERO) THEN
                     SOLID_DRAG(IJK,M,1) = -F_GS(IJK,M)*&
                     (DES_U_S(IJK,M)-UGC)
                     SOLID_DRAG(IJK,M,2) = -F_GS(IJK,M)*&
                     (DES_V_S(IJK,M)-VGC)
                     
                     IF(DIMN.EQ.3) THEN
                        SOLID_DRAG(IJK,M,3) = -F_GS(IJK,M)*&
                        (DES_W_S(IJK,M)-WGC)
                     ENDIF
                     OEPS = ONE/EP_S(IJK,M)
                     SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
                  ENDIF
                     
                  !rahul: temp start 
                  IF(MPPIC) THEN 
                     IF(EP_S(IJK,M).GT.ZERO) THEN
                        
                        if(MPPIC_PDRAG_IMPLICIT) THEN 
                           
                        !implicit treatment for drag term 
                           SOLID_DRAG(IJK, M, :) = F_GS(IJK,M)*VELG_ARR(:)
                           
                        ELSE
                        !explicit treatment 
                           SOLID_DRAG(IJK, M, :) = F_GS(IJK,M)*(VELG_ARR(:)-VELS_ARR(:, M))
                        endif
                     
                        EPs = MIN(EP_s(IJK,M), 1.d0-EP_STAR)
                        !EPs = EP_s(IJK,M)
                        OEPS = ONE/EPs
>>>>>>> 1.20.2.9
                        SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
                     ENDIF
                  ENDIF
                  !rahul: temp end
               ENDDO
            ENDIF      ! end IF(PINC(IJK).GT.0)
         ENDDO         ! end do loop over ijk
      ENDIF            ! end IF(CALC_FC)
!----------------------------------------------- 


! update the contact forces (FC) on the particle to include gas pressure
! and gas-solids drag (this section is performed when drag is not interpolated -
! otherwise code further down is used for updating FC)
!-----------------------------------------------
      PC = 1              
      DO NP = 1, MAX_PIP
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(np,1)) cycle 
         pc = pc+1
         if(pea(np,4)) cycle 
         IJK = PIJK(NP,4)
         M = PIJK(NP,5)
         OVOL = ONE/VOL(IJK)
         OEPS = ONE/EP_S(IJK,M)

         IF(MPPIC) THEN 
            !EPs = MIN(EP_s(IJK,M), 1.d0-EP_STAR)
            EPs = EP_s(IJK,M)
            OEPS = ONE/EPs
            
            IF(EPs.gt.zero) then 
               OEPS = ONE/EPs
            ELSE
               OEPS = ZERO 
            ENDIF
            
            if(MPPIC_PDRAG_IMPLICIT) THEN 
               F_gp(NP) = (F_GS(IJK,M)*PVOL(NP))*OEPS
            ELSE
               F_gp(NP) = ZERO 
            endif
         ENDIF
         
         DRAG_FORCE(:) = (SOLID_DRAG(IJK,M,:)*PVOL(NP)) 

         FC(NP,:) = FC(NP,:) + DRAG_FORCE(:) 
         
         IF(MODEL_B) THEN    !Do not add the pressure gradient force
         ELSE                !Add the pressure gradient force 
            FC(NP,:) = FC(NP,:) + (P_FORCE(IJK,:))*PVOL(NP)
            
               !now P_force is in fact -dp/dx and not -(dp/dx)*VOL(IJK) as in old implementation
         ENDIF
               
      ENDDO

    end SUBROUTINE DES_CALC_PART_DRAG_FORCE_INTERP_OFF
    
    SUBROUTINE DRAG_FGS
      
      USE param
      USE param1
      USE parallel
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      use desmpi 
      USE cutcell 
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------         
! local variable used for debugging
      LOGICAL FOCUS 

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
              IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK, &
! Pradeep added following 
              IMJMK,IMJKM,IJMKM,IMJMKM

! i,j,k indices of the fluid cell the particle resides in minus 1 
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL

! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER IW, IE, JS, JN, KB, KTP

! order of interpolation set in the call to set_interpolation_scheme unless it
! is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW

! constant whose value depends on dimension of system      
      DOUBLE PRECISION AVG_FACTOR     

! index of solid phase that particle NP belongs to      
      INTEGER M

! volume of fluid cell particle resides in
      DOUBLE PRECISION VCELL, VCELL2 
! one over the solids volume fraction and one over the volume       
      DOUBLE PRECISION OEPS, OVOL 

! particle number index, used for looping      
      INTEGER NP

! index to track accounted for particles 
      INTEGER PC 

! for error messages      
      INTEGER IER

! Statistical weight of the particle. Equal to one for DEM 
   
      DOUBLE PRECISION WTP, EPS

! Pradeep temporary indices for periodic boundary adjustments
      integer korder,cur_ijk,nindx
! Pradeep introducing volume at grid nodes for backward interpolation 
      double precision  VELG_ARR(DIMN), VELS_ARR(DIMN, MMAX)

! see the discussion for IJK_U ..... in comments       
      INTEGER  IJK_U, IJK_V, IJK_W, ICUR, JCUR, KCUR 
      
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, velf_part(dimn)

      
!
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!Rahul:
!The drag force computation for interpolation off has been moved to 
!a separate subroutine to make this routine more readable. 

      IF(.NOT.DES_INTERP_ON) THEN 
         CALL DES_CALC_PART_DRAG_FORCE_INTERP_OFF
         RETURN 
      ENDIF

! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
      AVG_FACTOR = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)   

! if calc_fc is true, then the contact forces (FC) will be updated 
! to include gas-solids drag and gas pressure force on the particle
! in this section the gas pressure force on the particle is computed.
! if the user did not specify using an interpolated drag force (i.e.,
! des_interp_on=t), then the gas solid drag force on the particle is
! also computed based on cell average quantities
!-----------------------------------------------      


      
! this section is used to calculate the gas solids drag force on each particle
! using particle velocity and the fluid velocity interpolated to particle
! position
!-----------------------------------------------      
      call set_interpolation_scheme(2)
      korder = 1+(dimn-2)
      drag_am = ZERO
      drag_bm = ZERO
      wtbar = zero

!Pradeep Changing from particles loop to cell loop 
      do ijk = ijkstart3,ijkend3
         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle 
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = (3-dimn)*1+(dimn-2)*(k-1)  ! =k-1 (in 3d) or =1 (in 2d)
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
              ktp,interp_scheme,dimn,ordernew = onew) 

!Compute velocity at grid nodes and set the geometric stencil 
         avg_factor = 0.25d0*(dimn-2) + 0.5d0*(3-dimn)
         do k = 1,(3-dimn)*1+(dimn-2)*onew
            do j = 1,onew
               do i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1
                  cur_ijk = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk))
                  ipjk    = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk))
                  ijpk    = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk))
                  ipjpk   = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk))
               
                  gstencil(i,j,k,1) = xe(ii)
                  gstencil(i,j,k,2) = yn(jj)
                  gstencil(i,j,k,3) = zt(kk)*(dimn-2) + dz(1)*(3-dimn)
                  vstencil(i,j,k,1) = avg_factor*(u_g(cur_ijk)+u_g(ijpk))
                  vstencil(i,j,k,2) = avg_factor*(v_g(cur_ijk)+v_g(ipjk)) 
                  if(dimn.eq.3) then 
                  
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))

                     vstencil(i,j,k,1) = vstencil(i,j,k,1)+avg_factor*(u_g(ijkp) + u_g(ijpkp))
                     vstencil(i,j,k,2) = vstencil(i,j,k,2)+avg_factor*(v_g(ijkp) + v_g(ipjkp))
                     vstencil(i,j,k,3) = avg_factor*(w_g(cur_ijk)+&
                          w_g(ijpk)+w_g(ipjk)+w_g(ipjpk))
                  else 
                     vstencil(i,j,k,3) = 0.d0
                  endif
               enddo
            enddo
         enddo

!loop through particles in the cell  
         do nindx = 1,pinc(ijk)
            focus = .false.
            np = pic(ijk)%p(nindx)
            
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            
            if (dimn .eq. 2) then 
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
                    vstencil(1:onew,1:onew,1,1:dimn), &
                    des_pos_new(np,1:dimn),vel_fp(np,1:dimn),  &
                    onew,interp_scheme,weightp)
            else 
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
                    vstencil(1:onew,1:onew,1:onew,1:dimn), &
                    des_pos_new(np,1:dimn),vel_fp(np,1:dimn),  &
                    onew,interp_scheme,weightp)
            end if
            call des_drag_gs(np,vel_fp(np,1:dimn),des_vel_new(np,1:dimn))

! Drag force on each particle is equal to \beta *(u_g -u_s)*Vol_p/eps,
! where Vol_p is the particle volume
! f_gp obtained from des_drag_fgs subroutine below is
! equal to \beta*Vol_p/eps 
! Therefore, drag force = f_gp*(u_g - u_s)


!Add drag to the contact force
            if(calc_fc) then
               ovol = one/vol(ijk)
               d_force(:) = f_gp(np)*(vel_fp(np,:)-des_vel_new(np,:))
               if(mppic) then 
                  if(MPPIC_PDRAG_IMPLICIT) THEN 
                     
               !implicit treatment of the drag term for mppic 
                     d_force(:) = f_gp(np)*(vel_fp(np,:))
                  ELSE
                     d_force(:) = f_gp(np)*(vel_fp(np,:)-des_vel_new(np,:))
                  endif
               endif 
               fc(np,:)=fc(np,:)+d_force(:)
               if(.not.model_b) then 
                  fc(np,:) =fc(np,:)+p_force(ijk,:)*pvol(np)
               endif
            endif

! if callfromdes is true, then the pertinent mean fields (in this case
! rop_s and f_gs) are not computed/updated in this call. this is done to
! speed up the simulation.  so the following section will only be called
! at the end of a given dem simulation for a given fluid time step
! (i.e., only called once per fluid time step)
            if(.not.callfromdes) then 
               m = pijk(np,5)
               do k = 1, (3-dimn)*1+(dimn-2)*onew
                  do j = 1, onew
                     do i = 1, onew
! shift loop index to new variables for manipulation                     
                        ii = iw + i-1
                        jj = js + j-1
                        kk = kb + k-1
!Pradeep: The interpolation is done using node. so one should use consistent numbering system
!in the current version imap_c is used instead of ip_of or im_of
                        
                        icur = imap_c(ii)
                        jcur = jmap_c(jj)
                        kcur = kmap_c(kk)
                     
                        cur_ijk = funijk(icur, jcur, kcur) !imap_c(ii),jmap_c(jj),kmap_c(kk))
                     
! should this volume be for current ijk index or always particle index?                        
! Pradeep replacing the volume of cell to volume at the node
                        vcell = des_vol_node(cur_ijk)
                        
                        ovol = one/vcell
                        drag_am(cur_ijk,m) = drag_am(cur_ijk,m) + &
                             f_gp(np)*weightp(i,j,k)*ovol*wtp
! first remove the velocity component at this grid point from the vel_fp
                        drag_bm_tmp(1:dimn) = vel_fp(np, 1:dimn) - &
                             weightp(i,j,k)*vstencil(i,j,k, 1:dimn)
! now find the remaning drag force
                        drag_bm_tmp(1:dimn) = des_vel_new(np,1:dimn) !- drag_bm_tmp(1:dimn)
                        drag_bm(cur_ijk, 1:dimn,m) = &
                             drag_bm(cur_ijk,1:dimn,m) + &
                             f_gp(np) * drag_bm_tmp(1:dimn) * &
                             weightp(i,j,k)*ovol*wtp 
                        wtbar(cur_ijk,m) = wtbar(cur_ijk,m) + &
                             weightp(i,j,k) *ro_s(m)*ovol*pvol(np)*WTP
                     enddo
                  enddo
               enddo
            endif       ! if(.not.callfromdes)
         enddo          ! pinc(ijk) loop 
      end do            ! ijk loop


      if(.not.callfromdes) then 
! Pradeep at the interface drag_am,drag_bm,wtbar has to be added
! send recv will be called and the node values will be added 
! at the junction 
         call des_addnodevalues
         avg_factor = 0.125d0*(dimn-2) + 0.25d0*(3-dimn)
         do ijk = ijkstart3, ijkend3
            if(fluid_at(ijk)) then
               i = i_of(ijk)
               j = j_of(ijk)
               k = k_of(ijk)
               if (i.lt.istart2 .or. i.gt.iend2) cycle
               if (j.lt.jstart2 .or. j.gt.jend2) cycle
               if (k.lt.kstart2 .or. k.gt.kend2) cycle
               imjk = funijk(imap_c(i-1),jmap_c(j),kmap_c(k))
               ijmk = funijk(imap_c(i),jmap_c(j-1),kmap_c(k))
               imjmk = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k))
               f_gs(ijk,:) = avg_factor*(drag_am(ijk,:) +&
                    drag_am(ijmk,:) + drag_am(imjmk,:) +&
                    drag_am(imjk,:))
               rop_s(ijk,:) = avg_factor*(wtbar(ijk,:) +&
                    wtbar(ijmk,:) + wtbar(imjmk,:) +&
                    wtbar(imjk,:))
               if(dimn.eq.3) then 
                  ijkm = funijk(imap_c(i),jmap_c(j),kmap_c(k-1))
                  imjkm = funijk(imap_c(i-1),jmap_c(j),kmap_c(k-1))
                  ijmkm = funijk(imap_c(i),jmap_c(j-1),kmap_c(k-1))
                  imjmkm = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k-1))
                  f_gs(ijk,:) = f_gs(ijk,:) + avg_factor*&
                       (drag_am(ijkm,:) + drag_am(ijmkm,:) +&
                       drag_am(imjmkm,:)+drag_am(imjkm,:) )
                  rop_s(ijk,:) = rop_s(ijk,:) + avg_factor*&
                       (wtbar(ijkm,:) + wtbar(ijmkm,:) + &
                       wtbar(imjmkm,:)+wtbar(imjkm,:) )
               endif
            endif
         enddo  ! ijk loop 
      endif     ! if(.not.callfromdes)
!-----------------------------------------------  

      RETURN
   
    END SUBROUTINE DRAG_FGS


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
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
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
      INTEGER          I,  IJK, IMJK, IJMK, IJKM, IM, M
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
     
<<<<<<< drag_fgs.f
      double precision :: epg
!     Current value of F_gs (i.e., without underrelaxation)

      DOUBLE PRECISION F_gstmp
=======
      double precision :: epg
!     Current value of F_gs (i.e., without underrelaxation)

      DOUBLE PRECISION F_gstmp

      DOUBLE PRECISION:: EPS, DIAMETER
>>>>>>> 1.20.2.9

<<<<<<< drag_fgs.f
      DOUBLE PRECISION:: EPS, DIAMETER

      INCLUDE 'ep_s1.inc'
=======
      INCLUDE 'ep_s1.inc'
>>>>>>> 1.20.2.9
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!     
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER) ! Tsuji drag
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2 ! Dalla Valle (1948) 
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0 * Re**0.657D0)      ! Turton and
!     &          + 0.413D0 * Re**2.09D0 / (Re**1.09D0 + 16300.D0) ! Levenspiel (1986)


      IJK = PIJK(KK,4)
      M = PIJK(KK,5)
      
      !EPS  = EP_S(IJK, M)
      EPs = MIN(EP_s(IJK,M), 1.d0-EP_STAR)
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

<<<<<<< drag_fgs.f
!     Reynolds number
      if(Mu > ZERO)then
         RE = diameter * VREL*RO_G(IJK)/Mu 
         
!     Note the presence of gas volume fraction in ROP_G
         RE_G = (diameter * VREL*ROP_G(IJK))/Mu
         
!     Note Reynolds' number for Hill and Koch has an additional factor of 1/2 & ep_g
         RE_kh = (0.5D0*diameter*VREL*ROP_G(IJK))/Mu
=======
!     Reynolds number
      if(Mu > ZERO)then
         RE = diameter * VREL*RO_G(IJK)/Mu 
         
!     Note the presence of gas volume fraction in ROP_G
         RE_G = (diameter * VREL*ROP_G(IJK))/Mu
>>>>>>> 1.20.2.9
         
<<<<<<< drag_fgs.f
      else 
=======
!     Note Reynolds' number for Hill and Koch has an additional factor of 1/2 & ep_g
         RE_kh = (0.5D0*diameter*VREL*ROP_G(IJK))/Mu
         
      else 
>>>>>>> 1.20.2.9
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
<<<<<<< drag_fgs.f
<<<<<<< drag_fgs.f

! Calculate V_rm
            V_RM=HALF*(A_SO - 0.06D0*RE + SQRT( (3.6D-3)*RE*RE + &
               0.12D0*RE*(2.D0*B_SO-A_SO) + A_SO*A_SO) ) 
=======
            V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) 
>>>>>>> 1.20.2.9
=======
            V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) 
>>>>>>> 1.20.2.9
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
<<<<<<< drag_fgs.f
<<<<<<< drag_fgs.f
               IF(EP_G(IJK) <= 0.8D0) THEN
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
=======
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
>>>>>>> 1.20.2.9
=======
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
>>>>>>> 1.20.2.9
            ELSE
               F_gstmp = 0.75D0*Mu*(PVOL(KK))*EP_G(IJK)*C_DSXRE(RE&
               /V_RM)/(V_RM*DIAMETER*DIAMETER) 
            ENDIF
         ENDIF
!---------------End Syamlal and O'Brien ---------------------------
!     
!--------------------------Begin Gidaspow --------------------------
      ELSE IF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') then
         !EPg = MAX(EP_STAR, EP_G(IJK))
         EPg = EP_G(IJK)
         IF(EPg .LE. 0.8D0) THEN
            DgA = 150D0 * (ONE - EPg) * Mu &
            / ( EPg * diameter**2 ) &
            + 1.75D0 * RO_g(IJK) * VREL / diameter
         ELSE
            IF(Re_G .LE. 1000D0)THEN
               C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF
            DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EPg**(-2.65D0) &
            /diameter
         ENDIF
         
!     Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         
         IF(Model_B)THEN
            F_gstmp = DgA*(PVOL(KK))/EPg
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

    
