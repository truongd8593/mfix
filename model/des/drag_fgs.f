! COMMENTS: 
! - The current methods for calculating the drag force on the continuum
!   phase and on the discrete particles will result in differences.

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: COMPUTE_PG_GRAD_CG                                      C
!  Purpose: Calculate cell centered pressure force exerted on the      C
!           particles in the cell by the gas/fluid phase               C
!           (cut-cell version)                                         C
!                                                                      C
!  Notes: This pressure force needs to be calculated once in a DEM     C
!         time step (at the beggining) since the gas/fluid phase is    C
!         not updated (is static) during the DEM portion of the        C
!         simulation.                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE COMPUTE_PG_GRAD_CG

!-----------------------------------------------
! Modules
!-----------------------------------------------
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
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! counters
      INTEGER :: X_COUNT, Y_COUNT, Z_COUNT
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: MPG_CYCLIC(DIMN)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

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
         X_COUNT = MAX(1, X_COUNT) !to prevent division from zero 
! P_FORCE (by convention) is stored as -dp/dx. MPG_CYCLIC is already -dp/dx. 
! therefore, P_force is multiplied by "-" for consistency
         P_FORCE(IJK, 1) = MPG_CYCLIC(1) - P_FORCE(IJK,1)/REAL(X_COUNT) 
         
         IF(FLUID_AT(IJKN)) THEN 
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) + 2.d0*(P_G(IJKN) - P_G(IJK))/(DY(J) + DY(J_OF(IJKN)))
         ENDIF
            
         IF(FLUID_AT(IJKS)) THEN 
            Y_COUNT = Y_COUNT + 1
            P_FORCE(IJK, 2) = P_FORCE(IJK, 2) + 2.d0*(P_G(IJK) - P_G(IJKS))/(DY(J) + DY(J_OF(IJKS)))
         ENDIF
         Y_COUNT = MAX(1, Y_COUNT) !to prevent division from zero 
         P_FORCE(IJK, 2) = MPG_CYCLIC(2) - P_FORCE(IJK,2)/REAL(Y_COUNT) 

         IF(DIMN.EQ.3) THEN
            IF(FLUID_AT(IJKT)) THEN 
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) + 2.d0*(P_G(IJKT) - P_G(IJK))/(DZ(K) + DZ(K_OF(IJKT)))
            ENDIF            
            IF(FLUID_AT(IJKB)) THEN 
               Z_COUNT = Z_COUNT + 1
               P_FORCE(IJK, 3) = P_FORCE(IJK, 3) + 2.d0*(P_G(IJK) - P_G(IJKB))/(DZ(K) + DZ(K_OF(IJKB)))
            ENDIF
            Z_COUNT = MAX(1, Z_COUNT) !to prevent division from zero 
            P_FORCE(IJK, 3) = MPG_CYCLIC(3) - P_FORCE(IJK,3)/REAL(Z_COUNT) 
         ENDIF
         
      ENDDO

      RETURN
      END SUBROUTINE COMPUTE_PG_GRAD_CG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: COMPUTE_PG_GRAD                                         C
!  Purpose: Calculate cell centered pressure force exerted on the      C
!           particles in the cell by the gas/fluid phase               C
!           (cut-cell version)                                         C
!                                                                      C
!  Notes: This pressure force only needs to be calculated once during  C
!         the DEM loop (at the beginning) since the gas/fluid phase    C
!         is essentially static at that point (i.e., gas field is not  C
!         updated during DEM loop                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE COMPUTE_PG_GRAD

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
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
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! temporary variables used to calculate pressure at scalar cell edge
      DOUBLE PRECISION :: TEMP1, TEMP2
! mean pressure gradient for the case of periodic boundaries
      DOUBLE PRECISION :: MPG_CYCLIC(DIMN)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

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

      RETURN
      END SUBROUTINE COMPUTE_PG_GRAD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_NONINTERP                                      C
!  Purpose: This subroutine is called from the subroutine des_drag,    C
!           which is only called from the DISCRETE phase.              C
!           This subroutine calculates the drag force exerted on the   C
!           particles by the gas/fluid phase using cell average        C
!           quantities (i.e., non-interpolated version). The drag      C
!           coefficient (F_GS) is calculated from the subroutine       C
!           drag_gs which is called during the continuum time step     C
!           (in calc_drag) and here during the discrete time step(s).  C
!           Accordingly, the drag coefficient in each call will be     C
!           based on the most currently available values (see notes)   C
!           The subroutine then adds the gas-solids drag force and     C
!           gas pressure force to the total contact force on the       C
!           particle                                                   C
!                                                                      C
!  Notes:                                                              C
!  During the continuum time step:                                     C
!     - F_GS will be based on updated gas velocity fields but          C
!       static solids velocity fields                                  C
!     - the field variable ep_g is not updated as the particles        C
!       are essentially static during the continuum loop               C
!  During the dem time step:                                           C
!     - F_GS will be based on updated solids velocity fields but       C
!       static gas velocity fields                                     C
!     - the field variable ep_g will be updated as particles move      C
!                                                                      C
!  Comments:                                                           C
!       Since calculation of the the drag coefficient is made from     C
!       both the continuum and discrete sides it will be different     C
!       Thus the total drag force acting on each side will be          C
!       different...                                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_DRAG_NONINTERP

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
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
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! see the discussion for IJK_U ..... in comments       
      INTEGER :: IJK_U, IJK_V, IJK_W
! average fluid velocity in x, y, z direction at scalar cell center      
      DOUBLE PRECISION :: UGC, VGC, WGC
! average continuum solids velocity in x, y, z direction at scalar
! cell center
      DOUBLE PRECISION :: USC, VSC, WSC
! average fluid and solid velocity in array form
      DOUBLE PRECISION :: VELG_ARR(DIMN), &
                          VELDS_ARR(DIMN, DIM_M), &
                          VELCS_ARR(DIMN, DIM_M)
! local drag force
      DOUBLE PRECISION :: SOLID_DRAG (DIMENSION_3, DIM_M, DIMN)
      DOUBLE PRECISION :: D_FORCE(DIMN)
! index of solid phase that particle NP belongs to      
      INTEGER :: M
! continuous solids phase index
      INTEGER :: CM
! particle number index, used for looping      
      INTEGER :: NP
! solids volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM      
! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION :: OEPS, OVOL 
! for error messages      
      INTEGER :: IER
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------


! Calculate the gas-solid drag force exerted on the particle 
!----------------------------------------------------------------->>> 

! initializations              
      SOLID_DRAG(:,:,:) = ZERO

! computing F_gs (drag coefficient) with the latest average fluid 
! and solid velocity fields 
      DO M = 1, SMAX
         CALL DRAG_GS (M, IER)
      ENDDO

!$omp parallel do default(shared)                                 &
!$omp private(ijk,i,j,k,imjk,ijmk,ijkm,ijk_u,ijk_v,ijk_w,         &
!$omp         ugc,vgc,wgc,velg_arr,velds_arr,                     &
!$omp         usc,vsc,wsc,velcs_arr,                              &
!$omp         m,cm,oeps,ep_sm,solid_drag) schedule (guided,50)
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

! UGC (and likewise for VGC and WGC) are computed at the center of 
! the scalar cell. The center of the Ith scalar cell can also be 
! thought of as east face of the (I-1)th U- cell. 
! For cut-cell, it is easier to think in terms of U, V, W grids/cell
! as the interpolation arays (like Theta_Ue_bar) are based on the
! respective grids, i.e., theta_Ue_bar is based on U- grid. See 
! conv_diff_u_g for an example of this usage. 
! For U uncut grid, the average at the east face of IJK_U will be
! U_AVG(at East face of IJK_U) = HALF*(U(IJK_U) + U(IP_OF(IJK_U)))
! It can be verified that the above formula and old formula of 
! U_AVG(center of IJK) = HALF*(U(IJK) + U(IM_OF(IJK))) are identical
! since IJK_U = IM_OF(IJK)

         IF(PINC(IJK).GT.0) THEN

! average fluid velocity at scalar cell center
            IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
               UGC = (Theta_Ue_bar(IJK_U)*U_G(IJK_U) + Theta_Ue(IJK_U)*U_G(IP_OF(IJK_U)))
            ELSE 
               UGC = HALF * (U_G(IJK_U) + U_G(IP_OF(IJK_U)))
               !UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)                  
            ENDIF
               
            IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
               VGC = (Theta_Vn_bar(IJK_V)*V_G(IJK_V) + Theta_Vn(IJK_V)*V_G(JP_OF(IJK_V)))
            ELSE
               VGC = HALF * (V_G(IJK_V) + V_G(JP_OF(IJK_V)))
               !VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))                  
            ENDIF

            IF(DIMN.EQ.3) THEN
               IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                  WGC = (Theta_Wt_bar(IJK_W)*W_G(IJK_W) + Theta_Wt(IJK_W) * W_G(KP_OF(IJK_W)))
               ELSE 
                  WGC = HALF * (W_G(IJK_W) + W_G(KP_OF(IJK_W)))
                  !WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
               ENDIF
            ENDIF

            VELG_ARR(1) = UGC
            VELG_ARR(2) = VGC
            VELDS_ARR(1,:) = DES_U_S(IJK,:)
            VELDS_ARR(2,:) = DES_V_S(IJK,:)
            IF(DIMN.EQ.3) THEN
               VELG_ARR(3) = WGC
               VELDS_ARR(3,:) = DES_W_S(IJK,:)
            ENDIF
            
! average continuum solids velocity at scalar cell center            
            IF(DES_CONTINUUM_HYBRID) THEN
               DO CM = 1, SMAX
                  IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
                  ELSE 
                     USC = AVG_X_E(U_S(IMJK,CM),U_S(IJK,CM),I) 
                  ENDIF
                     
                  IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
                  ELSE
                     VSC = AVG_Y_N(V_S(IJMK,CM),V_S(IJK,CM))
                  ENDIF
      
                  IF(DIMN.EQ.3) THEN
                     IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                     ELSE 
                        WSC = AVG_Z_T(W_S(IJKM,CM),W_S(IJK,CM))
                     ENDIF
                  ENDIF
               ENDDO    ! do loop (cm=1,smax)
               VELCS_ARR(1,CM) = USC
               VELCS_ARR(2,CM) = VSC
               IF (DIMN.EQ.3) VELCS_ARR(3,CM) = WSC
            ENDIF   ! end if(des_continuum_hybrid)

            DO M = 1, DES_MMAX
!               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               EP_SM = EP_S(IJK,M)
               IF(EP_SM.GT.ZERO) THEN
                  IF (.NOT.DES_CONTINUUM_HYBRID) THEN
                     SOLID_DRAG(IJK,M,:) = -F_GS(IJK,M)*&
                        (VELDS_ARR(:,M)-VELG_ARR(:))
                  ELSE   ! des_continuum_hybrid branch
                     SOLID_DRAG(IJK,M,:) = -F_GDS(IJK,M)*&
                        (VELDS_ARR(:,M)-VELG_ARR(:))
                     DO CM = 1, SMAX
                        SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:) + &
                           -F_SDS(IJK,CM,M)*&
                           (VELDS_ARR(:,M)-VELCS_ARR(:,CM))
                     ENDDO
                  ENDIF   ! end if/else (des_continuum_hybrid)

                  OEPS = ONE/EP_SM
                  SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
               ENDIF  ! end if ep_sm>0

               IF(MPPIC) THEN
!                 EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
                  EP_SM = EP_S(IJK,M)
                  IF(EP_SM.GT.ZERO) THEN
                     IF(MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment for drag term
                        SOLID_DRAG(IJK,M, :) = F_GS(IJK,M)*VELG_ARR(:)
                     ELSE
! explicit treatment 
                        SOLID_DRAG(IJK,M, :) = F_GS(IJK,M)*&
                           (VELG_ARR(:)-VELDS_ARR(:,M))
                     ENDIF
                     OEPS = ONE/EP_SM
                     SOLID_DRAG(IJK,M,:) = SOLID_DRAG(IJK,M,:)*OEPS
                  ENDIF
               ENDIF   ! end if(mppic))

            ENDDO   ! end do loop (dm=1,des_mmax)
         ENDIF      ! end if(pinc(ijk).gt.0)
      ENDDO         ! end do loop (ijk=ijkstart3, ijkend3)
!$omp end parallel do          
!-----------------------------------------------------------------<<<


! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag 
!----------------------------------------------------------------->>>
!$omp parallel do private(np,ijk,m,ovol,oeps,ep_sm,  &
!$omp                     solid_drag,d_force) 
!$omp schedule (guided,100)    
      DO NP = 1, MAX_PIP
! skipping indices that do not represent particles and ghost particles
         if(.not.pea(np,1)) cycle 
         if(pea(np,4)) cycle 

         IJK = PIJK(NP,4)
         M = PIJK(NP,5)
         OVOL = ONE/VOL(IJK)

         IF(MPPIC) THEN 
            EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
            OEPS = ONE/EP_SM
            IF(EP_SM.GT.zero) THEN
               OEPS = ONE/EP_SM
            ELSE
               OEPS = ZERO 
            ENDIF            
            IF(MPPIC_PDRAG_IMPLICIT) THEN 
               F_gp(NP) = (F_GS(IJK,M)*PVOL(NP))*OEPS
            ELSE
               F_gp(NP) = ZERO 
            ENDIF
         ENDIF   ! end if(mppic)
         
         D_FORCE(:) = (SOLID_DRAG(IJK,M,:)*PVOL(NP)) 

         FC(NP,:) = FC(NP,:) + D_FORCE(:) 
         
         IF(.NOT.MODEL_B) THEN    
! Add the pressure gradient force.
! P_force is in fact -dp/dx and not -(dp/dx)*VOL(IJK) as in old implementation                 
            FC(NP,:) = FC(NP,:) + (P_FORCE(IJK,:))*PVOL(NP)
         ENDIF
               
      ENDDO   ! end do loop (np=1,max_pip)
!$omp end parallel do
!-----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE DES_DRAG_NONINTERP
    

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_FGS                                                C
!  Purpose: This routine is called from the DISCRETE/CONTINUUM sides.  C
!           It performs the following functions:                       C
!     - If non-interpolated then execution of the code is directed     C
!       to the subroutine des_drag_noninterp for the appropriate       C
!       calculations.                                                  C
!     - If interpolated, then it calculates the particle centered      C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. This F_GP is used to:             C
!     - Adds the gas-solids drag and gas pressure force to the total   C
!       contact force on the particle (when callfromdes=.false.)       C
!     - Determines the contribution of particle centered drag to the   C
!       center coefficient of the A matrix and the b (source) vector   c
!       in the matrix equation (A*VEL_FP=b) equation for the gas       C
!       phase x, y and z momentum balances. (when callfromdes)         C
!                                                                      C
!  Notes:                                                              C
!     - It does not make sense to calculate/update rop_s in this       C
!       routine since particles will not have moved during the         C
!       continuum step.  Moreover, an accurate value of rop_s would    C
!       not be thoroughly integrated into the continuum side since     C
!       ep_g is only updated based on rop_s from the call to           C
!       particles_in_cell. The value of rop_s from this subroutine,    C
!       which is carried into the continuum side (via epg), would      C
!       not reflect the most recent particle position (particles       C
!       positions/velocities are updated after this subroutine).       C
!                                                                      C
!  Comments:                                                           C
!       Since calculation of the the drag coefficient is made from     C
!       both the continuum and discrete sides it will be different     C
!       Thus the total drag force acting on each side will be          C
!       different...                                                   C
!                                                                      C
!  Author/Revision: Pradeep G                                          C
!  Revisions:                                                          C
!      1. modified drag_interpolation routines to change three         C
!         dimensional arrays with separate indices i,j,k               C
!         (drag_am,drag_bm) into one dimensional arrays with the       C
!         composite index ijk                                          C
!      2. modified treatment for periodic boundary conditions          C
!      3. modified the loop structure from particle to IJK             C
!      4. added volume at node to include effect of volume in to       C
!         the backward interpolation                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DRAG_FGS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
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
! local drag forces
      DOUBLE PRECISION D_FORCE(DIMN)
      DOUBLE PRECISION drag_bm_tmp(DIMN)
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: II, JJ, KK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, &
                 IMJMK, IMJKM, IJMKM, IMJMKM
      INTEGER :: ICUR, JCUR, KCUR 
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! temporary indices for periodic boundary adjustments
      INTEGER :: cur_ijk, nindx
! i,j,k indices of the fluid cell the particle resides in minus 1 
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3) :: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW
! index of solid phase that particle NP belongs to      
      INTEGER :: M
! particle number index, used for looping      
      INTEGER :: NP
! one over the volume of fluid cell
      DOUBLE PRECISION :: OVOL 
! volume of fluid cell particle resides in
      DOUBLE PRECISION :: VCELL
! constant whose value depends on dimension of system 
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
! Statistical weight of the particle. Equal to one for DEM 
      DOUBLE PRECISION :: WTP
! for error messages      
      INTEGER :: IER
!-----------------------------------------------   
! Include statement functions
!-----------------------------------------------   
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!----------------------------------------------- 

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime	      

!!$      omp_start=omp_get_wtime()
      IF(.NOT.DES_INTERP_ON) THEN 
         CALL DES_DRAG_NONINTERP
         RETURN 
      ENDIF
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp_off:',omp_end - omp_start 


! initializations 
      drag_am = ZERO
      drag_bm = ZERO
      wtbar = zero

! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.250d0*(DIMN-2) + 0.50d0*(3-DIMN)

! sets several quantities including interp_scheme, scheme, and 
! order and allocates arrays necessary for interpolation
      call set_interpolation_scheme(2)

! There is some issue associated to gstencil, vstencil which are
! allocatable variables


!!$      omp_start=omp_get_wtime()
!!$omp parallel do default(shared)                                 &
!!$omp private(ijk,i,j,k,pcell,iw,ie,js,jn,kb,ktp,onew,            &
!!$omp         avg_factor,ii,jj,kk,cur_ijk,ipjk,ijpk,ipjpk,        &
!!$omp         gstencil,vstencil,ijpkp,ipjkp,ipjpkp,ijkp,nindx,    &
!!$omp         focus,np,wtp,weightp,ovol,d_force,m,icur,jcur,kcur,vcell,   &
!!$omp         drag_bm_tmp) schedule (guided,50)           	
      DO ijk = ijkstart3,ijkend3
         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle 
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)

! generally a particle may not exist in a ghost cell. however, if the
! particle is adjacent to the west, south or bottom boundary, then pcell
! may be assigned indices of a ghost cell which will be passed to
! set_interpolation_stencil
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = (3-dimn)*1+(dimn-2)*(k-1)  ! =k-1 (in 3d) or =1 (in 2d)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
              ktp,interp_scheme,dimn,ordernew = onew) 

! Compute velocity at grid nodes and set the geometric stencil
         DO k = 1,(3-dimn)*1+(dimn-2)*onew
            DO j = 1,onew
               DO i = 1,onew
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
                  IF(DIMN.EQ.3) THEN
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))

                     vstencil(i,j,k,1) = vstencil(i,j,k,1)+avg_factor*(u_g(ijkp) + u_g(ijpkp))
                     vstencil(i,j,k,2) = vstencil(i,j,k,2)+avg_factor*(v_g(ijkp) + v_g(ipjkp))
                     vstencil(i,j,k,3) = avg_factor*(w_g(cur_ijk)+&
                          w_g(ijpk)+w_g(ipjk)+w_g(ipjpk))
                  ELSE
                     vstencil(i,j,k,3) = 0.d0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

! loop through particles in the cell
! interpolate the fluid velocity (VEL_FP) to the particle's position.
         DO nindx = 1,PINC(IJK)
            NP = PIC(ijk)%p(nindx)
 
            IF (DIMN .EQ. 2) THEN
               CALL interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
                    vstencil(1:onew,1:onew,1,1:dimn), &
                    des_pos_new(np,1:dimn),vel_fp(np,1:dimn),  &
                    onew,interp_scheme,weightp)
            ELSE
               CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
                    vstencil(1:onew,1:onew,1:onew,1:dimn), &
                    des_pos_new(np,1:dimn),vel_fp(np,1:dimn),  &
                    onew,interp_scheme,weightp)
            ENDIF

! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the interpolated gas velocity.  Note F_GP
! obtained from des_drag_gs subroutine is given as f_gp=beta*vol_p/eps
! where vol_p is the particle volume.  The drag force on each particle
! is equal to beta(u_g-u_s)*vol_p/eps. 
! Therefore, the drag force = f_gp*(u_g - u_s)
            CALL DES_DRAG_GS(NP, VEL_FP(NP,1:DIMN), &
               DES_VEL_NEW(NP,1:DIMN))
!-----------------------------------------------------------------<<<


! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag if calc_fc
!----------------------------------------------------------------->>>
            if(calc_fc) then
               d_force(:) = f_gp(np)*(vel_fp(np,:)-des_vel_new(np,:))

               if(mppic) then 
                  if(MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment of the drag term for mppic 
                     d_force(:) = f_gp(np)*(vel_fp(np,:))
                  ELSE
                     d_force(:) = f_gp(np)*(vel_fp(np,:)-des_vel_new(np,:))
                  endif
               endif 

               FC(NP,:) = FC(NP,:) + d_force(:)
               IF(.NOT.MODEL_B) THEN
! P_force is in fact -dp/dx and not -(dp/dx)*VOL(IJK) as in old implementation
                  FC(NP,:) = FC(NP,:) + p_force(ijk,:)*pvol(NP)
               ENDIF
            endif
!-----------------------------------------------------------------<<<


! Calculate the corresponding gas solids drag force that is used in 
! the gas phase momentum balances.               
!----------------------------------------------------------------->>>
            if(.not.callfromdes) then
! invoke this section at the end of given dem simulation for the given
! fluid time step and every outer iteration in the fluid time step 

               focus = .false.
               WTP = ONE
               IF(MPPIC) WTP = DES_STAT_WT(NP)

               M = pijk(np,5)
               DO k = 1, (3-dimn)*1+(dimn-2)*onew
                  DO j = 1, onew
                     DO i = 1, onew
! shift loop index to new variables for manipulation                     
                        ii = iw + i-1
                        jj = js + j-1
                        kk = kb + k-1
! The interpolation is done using node. so one should use consistent 
! numbering system. in the current version imap_c is used instead of 
! ip_of or im_of
                        icur = imap_c(ii)
                        jcur = jmap_c(jj)
                        kcur = kmap_c(kk)
                        cur_ijk = funijk(icur, jcur, kcur) 
                     
! Replacing the volume of cell to volume at the node
                        vcell = des_vol_node(cur_ijk)
                        ovol = one/vcell

! first remove the velocity component at this grid point from the vel_fp
                        drag_bm_tmp(1:dimn) = vel_fp(np, 1:dimn) - &
                             weightp(i,j,k)*vstencil(i,j,k, 1:dimn)
! now find the remaning drag force
                        drag_bm_tmp(1:dimn) = des_vel_new(np,1:dimn) !- drag_bm_tmp(1:dimn)

!!$omp critical
                        drag_am(cur_ijk,m) = drag_am(cur_ijk,m) + &
                             f_gp(np)*weightp(i,j,k)*ovol*wtp
                        drag_bm(cur_ijk, 1:dimn,m) = &
                             drag_bm(cur_ijk,1:dimn,m) + &
                             f_gp(np) * drag_bm_tmp(1:dimn) * &
                             weightp(i,j,k)*ovol*wtp 
                        wtbar(cur_ijk,m) = wtbar(cur_ijk,m) + &
                             weightp(i,j,k) *ro_s(m)*ovol*pvol(np)*WTP
!!$omp end critical
                     ENDDO
                  ENDDO
               ENDDO
             ENDIF       ! if(.not.callfromdes)
         ENDDO   ! end do (nindx = 1,pinc(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)!!$omp end parallel do
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp:',omp_end - omp_start  

      if(.not.callfromdes) then 
! At the interface drag_am, drag_bm, and wtbar have to be added
! send recv will be called and the node values will be added 
! at the junction. drag_am, drag_bm and wtbar are altered by the 
! routine when periodic boundaries are invoked. so all three
! quantities are needed at the time of this call.      
         call des_addnodevalues
!-----------------------------------------------------------------<<<


! Calculate the mean fields rop_s and f_gs. ROP_s is used to update
! ep_g and f_gs is needed for the pressure correctione quation.
!----------------------------------------------------------------->>>
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
         AVG_FACTOR = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)

!$omp parallel do default(shared)                               &
!$omp private(ijk,i,j,k,imjk,ijmk,imjmk,ijkm,imjkm,ijmkm,       & 
!$omp         imjmkm) schedule (guided,20)           
         DO ijk = ijkstart3, ijkend3
            IF(fluid_at(ijk)) THEN
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
               IF(dimn.EQ.3) THEN
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

               ENDIF
            ENDIF   ! end if (fluid_at(ijk))
         ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do 
      endif     ! if(.not.callfromdes)
!-----------------------------------------------------------------<<<
    

      RETURN
      END SUBROUTINE DRAG_FGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GS                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_fgs and des_drag                     C
!                                                                      C
!  Comments: No BVK drag model in this subroutine. BVK requires an     C
!            average particle diameter which needs to be defined for   C
!            DEM case diameter.  Similarly no drag models with the     C
!            polydisperse correction factor (i.e., _PCF suffix) are    C
!            in this subroutine.                                       C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_DRAG_GS(KK, fvel, des_vel) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
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
! Dummy arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) ::  KK
! fluid velocity interpolated to particle position      
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: fvel
! particle velocity
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: des_vel
!-----------------------------------------------
! Local parameters 
!-----------------------------------------------
! Parameters in the Cluster-effect model
! a1 depends upon solids flux.  It has been represented 
! by C(1) defined in the data file.
!     PARAMETER (a1 = 250.)   !for G_s = 98 kg/m^2.s
!     PARAMETER (a1 = 1500.)  !for G_s = 147 kg/m^2.s
      DOUBLE PRECISION, PARAMETER :: A2 = 0.005D0 
      DOUBLE PRECISION, PARAMETER :: A3 = 90.0D0 
      DOUBLE PRECISION, PARAMETER :: RE_C = 5.D0 
      DOUBLE PRECISION, PARAMETER :: EP_C = 0.92D0 
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices 
      INTEGER :: I, IJK, IMJK, IJMK, IJKM, IM
! Solides phase index
      INTEGER :: M
! Cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION :: USCM, VSCM, WSCM
! Cell center value of U_g , V_g, W_g
      DOUBLE PRECISION :: UGC, VGC, WGC
! Magnitude of gas-solids relative velocity 
      DOUBLE PRECISION :: VREL 
! Gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! Reynolds number with and without void fraction in its definition,
! respectively 
      DOUBLE PRECISION :: RE_g, RE
! Single sphere drag coefficient 
      DOUBLE PRECISION :: C_d 
! Single sphere drag coefficient x Re 
      DOUBLE PRECISION :: C_DsxRe, C_DsxReT 
! Drag coefficient 
      DOUBLE PRECISION :: DgA  
! Current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION :: F_gstmp
! total solids volume fraction
      DOUBLE PRECISION :: phis
! solids volume fraction of phase M in fluid cell of interest
      DOUBLE PRECISION :: EP_SM    
! tmp variable for particle diameter
      DOUBLE PRECISION :: PART_DIAM
! tmp variable for particle volume
      DOUBLE PRECISION :: PART_VOL

!***********************************************************
! Blended Gidaspow drag correlation variables
!***********************************************************
! Gidaspow switch function variables [ceaf 2006-03-23]
      DOUBLE PRECISION :: Ergun
      DOUBLE PRECISION :: WenYu
      DOUBLE PRECISION :: PHI_gs
!*********************************************************** 

!***********************************************************
! Syam and O'Brien drag correlation variables
!***********************************************************
! Variables which are function of EP_g
      DOUBLE PRECISION :: A, B
! Ratio of settling velocity of a multiparticle system to 
! that of a single particle
      DOUBLE PRECISION :: V_rm 
!***********************************************************

!***********************************************************
! Koch and Hill drag correlation variables
!***********************************************************
! Stokes Drag Force
      DOUBLE PRECISION :: F_STOKES
! zero Re function for low Reynolds number
      DOUBLE PRECISION :: F_0
! inertial function for low Reynolds number
      DOUBLE PRECISION :: F_1
! zero Re function for high Reynolds number
      DOUBLE PRECISION :: F_2
! inertial function for high Reynolds number
      DOUBLE PRECISION :: F_3
! dimensionless drag force F
      DOUBLE PRECISION :: F
! transition Reynolds numbers
      DOUBLE PRECISION :: Re_Trans_1, Re_Trans_2
! weighting factor to compute F_0 and F_2
      DOUBLE PRECISION :: w
! Hill and Koch Reynolds number
      DOUBLE PRECISION :: Re_kh
!***********************************************************  

!-----------------------------------------------    
! Include statement functions
!-----------------------------------------------  
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'      
!-----------------------------------------------  

! Tsuji Drag      
      C_DSXRET(RE) = 24D0*(1 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
! Dalla Valle (1948)      
      C_DSXRE(RE) = (0.63D0*SQRT(RE) + 4.8D0)**2
! Turton and Levenspiel (1986)      
!     C_DsxRe (Re) = 24.D0 * (1.D0 + 0.173D0*Re**0.657D0) + &
!        0.413D0*Re**2.09D0 / (Re**1.09D0 + 16300.D0)

      IJK = PIJK(KK,4)
      M = PIJK(KK,5)
!      EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)      
      EP_SM = EP_S(IJK,M)
      PART_DIAM = 2.D0*DES_RADIUS(KK)
      PART_VOL = PVOL(KK)      
      
      IF(DIMN == 2)THEN
         VREL = SQRT((FVEL(1)- DES_VEL(1))**2 +&
                     (FVEL(2) - DES_VEL(2))**2)
      ELSE
         VREL = SQRT((FVEL(1) - DES_VEL(1))**2 +&
                     (FVEL(2) - DES_VEL(2))**2 +&
                     (FVEL(3) - DES_VEL(3))**2)
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
! Note the presence of gas volume fraction in ROP_G and
! an additional factor of 1/2
         RE_kh = (0.5D0*PART_DIAM*VREL*ROP_G(IJK))/Mu
      ELSE
         RE = LARGE_NUMBER 
         RE_G = LARGE_NUMBER
         RE_kh = LARGE_NUMBER
      ENDIF

  

!---------------Begin Syamlal and O'Brien ---------------------------
      IF(TRIM(DRAG_TYPE).EQ.'SYAM_OBRIEN') THEN
         IF (EP_sM <= ZERO) THEN 
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

! Calculate V_rm
            V_RM=HALF*(A - 0.06D0*RE + &
               SQRT( (3.6D-3)*RE*RE + 0.12D0*RE*(2.D0*B-A) + A*A)) 
!------------------Begin cluster correction --------------------------
!     uncomment the following four lines and comment the above line V_RM=... 
!     V_RM=HALF*(A-0.06D0*RE+SQRT(3.6D-3*RE*RE+0.12D0*RE*(2.D0*B-A)+A*A)) & 
!     * ( ONE + C(1) * exp( -a2*(Re - Re_c)**2 &
!     - a3*(EP_g(IJK)-ep_c)**2 &
!     )       * Re * (1. - EP_g(IJK))                )
!------------------End cluster correction ----------------------------
    
            DgA =  0.75D0*Mu*EP_G(IJK)* C_DSXRE(RE/V_RM) / &
              (V_RM*PART_DIAM*PART_DIAM) 

            IF(TSUJI_DRAG) THEN
               IF(EP_G(IJK) <= 0.8D0) THEN
                  F_gstmp = (Mu*PART_VOL/(PART_DIAM**2))*&
                  (150D0*(EP_SM/EP_G(IJK)) + 1.75D0*RE)
               ELSEIF(EP_G(IJK) > 0.8D0) THEN
                  IF(RE*EP_G(IJK) > 1000.D0) THEN
                     F_gstmp = 0.75D0*0.43D0*Mu*PART_VOL*RE/&
                        (PART_DIAM**2 *EP_G(IJK)**1.7D0)
                  ELSEIF(RE*EP_G(IJK) <= 1000.D0) THEN
                     F_gstmp = 0.75D0*C_DSXRET(RE*EP_G(IJK))*Mu*&
                        PART_VOL*RE/(PART_DIAM**2 *EP_G(IJK)**1.7D0)
                  ENDIF
               ENDIF

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
            ELSE IF(MODEL_B) THEN 
               F_gstmp = DgA * PART_VOL/EP_G(IJK)
            ELSE
               F_gstmp = DgA * PART_VOL
            ENDIF
         ENDIF
!---------------End Syamlal and O'Brien ---------------------------


!--------------------------Begin Gidaspow --------------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'GIDASPOW') THEN
         IF(EP_g(IJK) <= 0.8D0) THEN
            DgA = 150D0 * (ONE - EP_g(IJK)) * Mu /&
               (EP_g(IJK) * PART_DIAM**2 ) + &
               1.75D0 * RO_g(IJK) * VREL / PART_DIAM
         ELSE
            IF(Re_G .LE. 1000D0)THEN
               C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
            ELSE
               C_d = 0.44D0
            ENDIF
            DgA = 0.75D0 * C_d*VREL*ROP_g(IJK)*EP_g(IJK)**(-2.65D0) / &
               PART_DIAM
         ENDIF
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL 
         ENDIF
!--------------------------End Gidaspow --------------------------


!-----------------------Begin Gidaspow_blend ---------------------
      ELSEIF(TRIM(DRAG_TYPE).EQ.'GIDASPOW_BLEND') THEN
! Dense phase - EP_g < 0.8
         Ergun = 150D0 * (ONE - EP_g(IJK)) * Mu &
         / ( EP_g(IJK) * PART_DIAM**2 ) &
         + 1.75D0 * RO_g(IJK) * VREL / PART_DIAM
     
! Dilute phase - EP_g >= 0.8
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         WenYu = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
         /PART_DIAM
     
! Switch function
         PHI_gs = ATAN(150D0 * 1.75D0 * (EP_g(IJK) - 0.8D0)) / PI + 0.5D0
! Blend the models
         DgA = (1D0 - PHI_gs) * Ergun + PHI_gs * WenYu
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
!-----------------------End Gidaspow_blend -----------------------


!--------------------------Begin WEN_YU --------------------------
      ELSE IF(TRIM(DRAG_TYPE).EQ.'WEN_YU') THEN
         IF(Re_G .LE. 1000D0)THEN
            C_d = (24.D0/(Re_G+SMALL_NUMBER)) * (ONE + 0.15D0 * Re_G**0.687D0)
         ELSE
            C_d = 0.44D0
         ENDIF
         DgA = 0.75D0 * C_d * VREL * ROP_g(IJK) * EP_g(IJK)**(-2.65D0) &
         /PART_DIAM
         
! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
         
!--------------------------End WEN_YU ----------------------------


!--------------------Begin Koch & Hill (2001) --------------------
! Added by Clay Sutton (Lehigh University) 7-14-04
! Clay's implementation was modified by Sof (01-21-2005)
! for a report explaining these changes contact s.benyahia

      ELSE IF(TRIM(DRAG_TYPE).EQ.'KOCH_HILL') THEN
    
         F_STOKES = 18D0*MU_g(IJK)*EP_g(IJK)*EP_g(IJK)/PART_DIAM**2
         
         phis = ONE-EP_G(IJK)   
         w = EXP(-10.0D0*(0.4D0-phis)/phis)
         
         IF(phis > 0.01D0 .AND. phis < 0.4D0) THEN
            F_0 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
               135.0D0/64.0D0*phis*LOG(phis) + 17.14D0*phis) / &
               (1.0D0 + 0.681D0*phis - 8.48D0*phis*phis + &
               8.16D0*phis**3) + w*10.0D0*phis/(1.0D0-phis)**3
         ELSE IF(phis >= 0.4D0) THEN
            F_0 = 10.0D0*phis/(1.0D0-phis)**3
         ENDIF
         
         IF(phis > 0.01D0 .AND. phis <= 0.1D0) THEN
            F_1 = dsqrt(2.0D0/phis) / 40.0D0
         ELSE IF(phis > 0.1D0) THEN
            F_1 = 0.11D0 + 5.1D-04 * exp(11.6D0*phis)
         ENDIF
         
         IF(phis < 0.4D0) THEN
            F_2 = (1.0D0-w) * (1.0D0 + 3.0D0*dsqrt(phis/2.0D0) + &
               135.0D0/64.0D0*phis*LOG(phis) + 17.89D0*phis) / &
               (1.0D0 + 0.681D0*phis - 11.03D0*phis*phis + &
               15.41D0*phis**3) + w*10.0D0*phis/(1.0D0-phis)**3
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
         ELSEIF(phis <= 0.01D0 .AND. Re_kh > Re_Trans_1 .OR.   &
                phis >  0.01D0 .AND. Re_kh > Re_Trans_2) THEN
            F = F_2 + F_3*Re_kh            
         ELSE
            F = zero
         ENDIF
         
! This is a check for phis (or ep_sm) to be within physical range
         IF(phis <= ZERO .OR. phis > ONE) F = zero
         
         DgA = F * F_STOKES

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
!-------------------End Koch and Hill (2001)-----------------------


!------------------------Begin Dilue case -------------------------
      ELSEIF((DRAG_TYPE).EQ.'DILUTE_CASE') THEN
         C_d =  C_DSXRET(RE)
         DgA = (0.75D0 * C_d*VREL * RO_g(IJK))/(PART_DIAM)
         IF(Model_B)THEN
            F_gstmp = DgA * PART_VOL/EP_g(IJK)
         ELSE
            F_gstmp = DgA * PART_VOL
         ENDIF
!--------------------------End Dilue case --------------------------


      ELSE   ! value for (trim(drag_type)) not a valid option
         CALL START_LOG 
         if(mype == pe_io) WRITE (*, '(A,A)') &
            'Unknown DRAG_TYPE: ', DRAG_TYPE
         WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
         CALL END_LOG 
         call mfix_exit(myPE)  
      ENDIF   ! end if/else selection of drag_type

! f_gp() =  single particle drag excluding vector(v_g - v_p)   
      F_gp(kk) = (ONE - UR_F_gs) * F_gp(KK) + UR_F_gs * F_gstmp

      END SUBROUTINE DES_DRAG_GS

