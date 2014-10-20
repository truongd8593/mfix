! Reorganizing was performed for better calculation of interpolated
! rop_s.  Contact (j.galvin) if you forsee issues/improvements/
! problems.
!  1) The calculation of rop_s (when des_interp_on) should be done
!     in the subroutine particles_in_cell to reflect the most current
!     values of particle position.  It should also be done every dem
!     time step if it is to be reflected in calculation of the drag
!     coefficient.
!     In the old method, an accurate value of rop_s would not
!     be carried into the continuum side since ep_g is only updated
!     based on rop_s from the call to particles_in_cell. Thus the
!     value of rop_s calculated in the drag subroutine would not
!     reflect the most recent particle position (as particles
!     positions/velocities are updated after the drag subroutine).
!     To avoid this issue calculation of the interpolated rop_s
!     was placed into its own subroutine.


! Comments:
! 1) The interpolation mechanism/code should probably be streamlined
!    in these subroutines but this requires more work at this point.

! 2) The drag coefficient is calculated during the continuum time step
!    and re-calculated during the subsequent discrete time steps.  As a
!    result the total drag force acting on each phase may be differing.
!    - During iterations in the continuum time step the particles are
!      treated as static entities:
!      - F_GS (or F_GP) will be based on updated gas velocity fields but
!        static solids velocity fields (or particle velocities)
!      - The field variable ep_g does not change since the particles
!        are essentially static during the continuum step
!    - During the dem time step (within the continuum time step):
!      - F_GS (or F_GP) will be based on updated solids velocity fields
!        (or particle velocities) but essentially static gas velocity
!        (except for changes when interpolation is used)
!      - The field variable ep_g will be updated as particles

! Commented by Handan Liu on August 2012, revised at June 2013.
! 1. Added two new subroutines in the end to calculate interpolation,
!        i.e. DRAG_INTERPLATION_2D and DRAG_INTERPLATION_3D,
!        to replace invoking 'interpolator' interface module,
!        which will introduce global variables; thus cause datarace with OpenMP.
! 2. Set the intermediate variables 'desposnew(:),velfp(:)' (small arrays) as private
!        to replace the global variables des_pos_new(:,np),vel_fp(:,np) (big arrays)
!        to avoid datarace; and save the calculating time.
! 3. Set the intermediate array 'weight_ft' to replace the pointer 'weightp' (global)
!        in 'interpolator', which leads errors with OpenMP.
! 4. Set the intermediate variables 'gst_tmp & vst_tmp' as private to replace global variables
!        'gstencil,vstencil' to avoid datarace and segmentation fault.

! 2014-02-24 mmeredith combined the above into one DRAG_INTERPOLATION() subroutine

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS_NONINTERP                              C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions.                             C
!     - Calculates the fluid-solids drag force exerted on the          C
!       particles by the fluid phase using cell average quantities.    C
!       The gas solids drag coefficient (F_GS) is calculated from      C
!       the subroutine drag_gs, which is called during the continuum   C
!       time step (from calc_drag) and during the discrete time(s)     C
!       (from here).                                                   C
!     - The total contact force on the particle is then updated to     C
!       include the gas-solids drag force and gas pressure force       C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_DES_DRAG_GS_NONINTERP

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
      use desmpi
      USE cutcell
      USE mfix_pic
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! see the discussion for IJK_U ..... in comments
      INTEGER :: IJK_U, IJK_V, IJK_W
! average fluid and solid velocity at scalar cell center in array form
      DOUBLE PRECISION :: VELG_ARR(3), &
                          VELDS_ARR(DES_MMAX, 3)
! local drag force
      DOUBLE PRECISION :: GS_DRAG (DIMENSION_3, DES_MMAX, 3)
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP
! solids volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM
! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION :: OEPS
! for error messages
      INTEGER :: IER
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.
      LOGICAL :: DISCRETE_FLAG
!-----------------------------------------------


! initializing
      GS_DRAG(:,:,:) = ZERO

! computing F_gs (gas-solids drag coefficient) with the latest average
! fluid and solid velocity fields.  see comments below
      DISCRETE_FLAG = .TRUE.   ! only matters if des_continuum_hybrid
      DO M = 1, DES_MMAX
         IF (RO_G0/=ZERO) THEN
! this call can not be readily replaced with a call to des_drag_gp
! due to difficulties going from particle phase and ijk loops to a
! particle loop & vice versa
            CALL DRAG_GS (M, DISCRETE_FLAG, IER)
         ENDIF
      ENDDO

!$omp parallel do default(shared)                                 &
!$omp private(ijk,i,j,k,imjk,ijmk,ijkm,ijk_u,ijk_v,ijk_w,         &
!$omp         velg_arr,velds_arr,                                 &
!$omp         m,oeps,ep_sm)                               &
!$omp schedule (guided,50)
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

! average fluid velocity at scalar cell center in array form
            IF(CUT_U_TREATMENT_AT(IJK_U)) THEN
               VELG_ARR(1) = (Theta_Ue_bar(IJK_U)*U_G(IJK_U) + &
                              Theta_Ue(IJK_U)    *U_G(IP_OF(IJK_U)))
            ELSE
               VELG_ARR(1) = HALF * (U_G(IJK_U) + U_G(IP_OF(IJK_U)))
            ENDIF

            IF(CUT_V_TREATMENT_AT(IJK_V)) THEN
               VELG_ARR(2) = (Theta_Vn_bar(IJK_V)*V_G(IJK_V) + &
                              Theta_Vn(IJK_V)    *V_G(JP_OF(IJK_V)))
            ELSE
               VELG_ARR(2) = HALF * (V_G(IJK_V) + V_G(JP_OF(IJK_V)))
            ENDIF

            VELDS_ARR(:,1) = DES_U_S(IJK,:)
            VELDS_ARR(:,2) = DES_V_S(IJK,:)

            IF(DO_K) THEN
               IF(CUT_W_TREATMENT_AT(IJK_W)) THEN
                  VELG_ARR(3) = (Theta_Wt_bar(IJK_W)*W_G(IJK_W) + &
                                 Theta_Wt(IJK_W)    * W_G(KP_OF(IJK_W)))
               ELSE
                  VELG_ARR(3) = HALF * (W_G(IJK_W) + W_G(KP_OF(IJK_W)))
               ENDIF
               VELDS_ARR(:,3) = DES_W_S(IJK,:)
            ELSE
               VELG_ARR(3) = ZERO
               VELDS_ARR(:,3) = ZERO
            ENDIF


            DO M = 1, DES_MMAX
! the call to drag coefficient should probably be here for cut-cell
! since cell center velocities are now known and these would be used
! in vrel within drag correlations). this would require some
! rearrangement/reworking?  alternatively these calculations should be
! moved into drag_gs?

               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)

               IF(EP_SM.GT.ZERO) THEN
                  IF (MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
                     GS_DRAG(IJK,M, :) = F_GS(IJK,M)*VELG_ARR(:)
                  ELSEIF (DES_CONTINUUM_HYBRID) THEN
                     GS_DRAG(IJK,M,:) = -F_GDS(IJK,M)*&
                        (VELDS_ARR(M,:)-VELG_ARR(:))
                  ELSE ! default cuase
                     GS_DRAG(IJK,M,:) = -F_GS(IJK,M)*&
                        (VELDS_ARR(M,:)-VELG_ARR(:))
                  ENDIF   ! end if/else (des_continuum_hybrid)
                  OEPS = ONE/EP_SM
                  GS_DRAG(IJK,M,:) = GS_DRAG(IJK,M,:)*OEPS
               ENDIF  ! end if ep_sm>0

            ENDDO   ! end do loop (dm=1,des_mmax)
         ENDIF      ! end if(pinc(ijk).gt.0)

      ENDDO         ! end do loop (ijk=ijkstart3, ijkend3)
!$omp end parallel do


!$omp parallel do private(np,ijk,m,oeps,ep_sm)             &
!$omp schedule (guided,100)
      DO NP = 1, MAX_PIP
! skipping indices that do not represent particles and ghost particles
         if(.not.pea(np,1)) cycle
         if(pea(np,4)) cycle

         IJK = PIJK(NP,4)
         M = PIJK(NP,5)

! Update the contact forces (FC) on the particle to include
! gas pressure and gas-solids drag
!----------------------------------------------------------------->>>
         FC(:,NP) = FC(:,NP) + GS_DRAG(IJK,M,:)*PVOL(NP)

         IF(.NOT.MODEL_B) THEN
! Add the pressure gradient force
! P_force is evaluated as -dp/dx
            FC(:,NP) = FC(:,NP) + (P_FORCE(IJK,:))*PVOL(NP)
         ENDIF
!-----------------------------------------------------------------<<<

! For mppic calculate the gas-particle drag coefficient (F_GP).
!----------------------------------------------------------------->>>
         IF(MPPIC) THEN
            EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
            IF(EP_SM.GT.ZERO) THEN
               OEPS = ONE/EP_SM
            ELSE
               OEPS = ZERO
            ENDIF
            IF(MPPIC_PDRAG_IMPLICIT) THEN
               F_gp(NP) = (F_GS(IJK,M)*PVOL(NP))*OEPS
            ELSE
               F_gp(NP) = ZERO
            ENDIF
         ENDIF   ! end if (mppic)
!-----------------------------------------------------------------<<<

      ENDDO   ! end do loop (np=1,max_pip)
!$omp end parallel do

      RETURN
      END SUBROUTINE CALC_DES_DRAG_GS_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS                                        C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions:                             C
!     - If non-interpolated then execution of the code is directed     C
!       to the subroutine calc_des_drag_gs_noninterp for the           C
!       appropriate calculations (i.e., calculation of the fluid-      C
!       solids drag force on the particle)                             C
!     - If interpolated, then it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. This F_GP is used to calculate    C
!       the fluid-solids drag on the particle.                         C
!     - The total contact force on the particle is then updated to     C
!       include the gas-solids drag force and gas pressure force       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_DES_DRAG_GS

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
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL FOCUS
! general i, j, k indices
      INTEGER :: I, J, K, IJK, cur_ijk
      INTEGER :: II, JJ, KK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3) :: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP, nindx
! constant whose value depends on dimension of system
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
! for error messages
      INTEGER :: IER

!Handan Liu added temporary variables on April 20 2012
          DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
          DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
          DOUBLE PRECISION :: velfp(3), desposnew(3)
          DOUBLE PRECISION :: D_FORCE(3)
          DOUBLE PRECISION, DIMENSION(3) :: VEL_NEW
!
!-----------------------------------------------
!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
! NON-INTERPOLATED fluid-solids drag:
! Calculate the gas solids drag force on a particle using the cell
! averaged particle velocity and the cell average fluid velocity
!----------------------------------------------------------------->>>
      IF(.NOT.DES_INTERP_ON) THEN
         CALL CALC_DES_DRAG_GS_NONINTERP
         RETURN
      ENDIF
!-----------------------------------------------------------------<<<


! INTERPOLATED fluid-solids drag (the rest of this routine):
! Calculate the gas solids drag coefficient using the particle
! velocity and the fluid velocity interpolated to particle
! position.
!----------------------------------------------------------------->>>
! initializing
! RG: I do not understand the need for this large array. It is not
! used anywhere else except the same subroutine it is calculated.
! A local single rank array (like in the old implementation) was just fine
      vel_fp = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = merge(0.50d0, 0.25d0, NO_K)

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
!!$omp         focus,np,weightp,m)                                 &
!!$omp schedule (guided,50)

!Handan Liu modified the following do-loop on Jan 20 2013.
!$omp parallel do default(shared)                                       &
!$omp private(ijk, i, j, k, pcell, iw, ie, js, jn, kb, ktp,             &
!$omp         onew, ii, jj, kk,cur_ijk, ipjk, ijpk, ipjpk,              &
!$omp         gst_tmp, vst_tmp, velfp, desposnew, ijpkp, ipjkp, &
!$omp         ipjpkp,ijkp,nindx,np,weight_ft,d_force)
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
         pcell(3) = merge(1, k-1, NO_K)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
              ktp,interp_scheme,dimn,ordernew = onew)

! Compute velocity at grid nodes and set the geometric stencil
         DO k = 1,merge(1, ONEW, NO_K)
            DO j = 1,onew
               DO i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1
                  cur_ijk = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk))
                  ipjk    = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk))
                  ijpk    = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk))
                  ipjpk   = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk))

                  gst_tmp(i,j,k,1) = xe(ii)
                  gst_tmp(i,j,k,2) = yn(jj)
                  gst_tmp(i,j,k,3) = merge(DZ(1), zt(kk), NO_K)
                  vst_tmp(i,j,k,1) = avg_factor*(u_g(cur_ijk)+u_g(ijpk))
                  vst_tmp(i,j,k,2) = avg_factor*(v_g(cur_ijk)+v_g(ipjk))

                  if(DO_K) then
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))
                                         vst_tmp(i,j,k,1) = vst_tmp(i,j,k,1) + avg_factor*(u_g(ijkp) + u_g(ijpkp))
                     vst_tmp(i,j,k,2) = vst_tmp(i,j,k,2) + avg_factor*(v_g(ijkp) + v_g(ipjkp))
                     vst_tmp(i,j,k,3) = avg_factor*(w_g(cur_ijk)+&
                          w_g(ijpk)+w_g(ipjk)+w_g(ipjpk))
                  else
                     vst_tmp(i,j,k,3) = 0.d0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!===================================================================>> Handan Liu
! loop through particles in the cell
! interpolate the fluid velocity (VEL_FP) to the particle's position.
         DO nindx = 1,PINC(IJK)
            NP = PIC(ijk)%p(nindx)
! skipping indices that do not represent particles and ghost particles
            if(.not.pea(np,1)) cycle
            if(pea(np,4)) cycle

            desposnew(:) = des_pos_new(:,np)
            call DRAG_INTERPOLATION(gst_tmp,vst_tmp,desposnew,velfp,weight_ft)
            vel_fp(1:3,np) = velfp(1:3)




!**********************************************************************!
!                      TEST CASE MODIFICATIONS                         !
!``````````````````````````````````````````````````````````````````````!
!          Massless particles for the circle advection test            !
!----------------------------------------------------------------------!
                   DES_VEL_NEW(:,NP) = VEL_FP(:,NP)
!......................................................................!



! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the interpolated gas velocity.  Note F_GP
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            VEL_NEW(:) = DES_VEL_NEW(:,NP)
            CALL DES_DRAG_GP(NP, velfp(1:3), VEL_NEW)

! Calculate the gas-solids drag force on the particle
            IF(MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment of the drag term for mppic
!------------------------------------------------------------------<<<< Handan Liu
               D_FORCE(1:3) = F_GP(NP)*(VEL_FP(1:3,NP))
            ELSE
! default case
               D_FORCE(1:3) = F_GP(NP)*(VEL_FP(1:3,NP)-VEL_NEW)
            ENDIF

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
            FC(:3,NP) = FC(:3,NP) + D_FORCE(:3)

            IF(.NOT.MODEL_B) THEN
! P_force is evaluated as -dp/dx
               FC(:3,NP) = FC(:3,NP) + p_force(ijk,1:3)*pvol(NP)
            ENDIF
!------------------------------------------------------------------>>>> Handan Liu
         ENDDO       ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!$omp end parallel do
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp:',omp_end - omp_start

!-----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE CALC_DES_DRAG_GS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GS                                             C
!  Purpose: This subroutine is only called from the CONTINUUM side.    C
!     It performs the following functions:                             C
!     - If non-interpolated, then execution of the code is directed    C
!       to the subroutine drag_gs for the appropriate calculations     C
!       (i.e., calculation of the fluid-solids drag coefficient)       C
!     - If interpolated then, it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. It then determines the            C
!       the contributions of fluid-particle drag to the center         C
!       coefficient of the A matrix and the b (source) vector in the   C
!       matrix equation (A*VEL_FP=b) equation for the fluid phase      C
!       x, y and z momentum balances using F_GP.                       C
!                                                                      C
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

      SUBROUTINE DES_DRAG_GS

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
      USE desmpi
      USE cutcell
      USE mfix_pic
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL :: FOCUS
! general i, j, k indices
      INTEGER :: I, J, K, IJK, cur_ijk
      INTEGER :: II, JJ, KK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, &
                 IMJMK, IMJKM, IJMKM, IMJMKM
      INTEGER :: ICUR, JCUR, KCUR
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3) :: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to set_interpolation_stencil
      INTEGER :: ONEW
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP, nindx
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
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.
      LOGICAL :: DISCRETE_FLAG
!Handan Liu added temporary variables on April 20 2012
      DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
      DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
      DOUBLE PRECISION :: velfp(3), desposnew(3)
      DOUBLE PRECISION, DIMENSION(3) :: VEL_NEW

!-----------------------------------------------

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime

!!$      omp_start=omp_get_wtime()

! NON-INTERPOLATED fluid-solids DRAG:
! Calculate the fluid solids drag coefficient (F_GS) using the cell
! averaged particle velocity and the cell average fluid velocity
!----------------------------------------------------------------->>>
      IF (.NOT.DES_INTERP_ON) THEN
         DISCRETE_FLAG = .TRUE.   ! only matters if des_continuum_hybrid
         DO M = 1, DES_MMAX
            IF (RO_G0/=ZERO) THEN
! this call can not be readily replaced with a call to des_drag_gp
! due to difficulties going from particle phase and ijk loops to a
! particle loop & vice versa
               CALL DRAG_GS (M, DISCRETE_FLAG, IER)
            ENDIF
         ENDDO
         RETURN
      ENDIF
!-----------------------------------------------------------------<<<
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp_off:',omp_end - omp_start

! INTERPOLATED fluid-solids drag (the rest of this routine):
! Calculate the fluid solids drag coefficient using the particle
! velocity and the fluid velocity interpolated to particle
! position.
!----------------------------------------------------------------->>>
! initializations
      drag_am = ZERO
      drag_bm = ZERO
      vel_fp = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = merge(0.50d0, 0.25d0, NO_K)

! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      call set_interpolation_scheme(2)
! There is some issue associated to gstencil, vstencil which are
! allocatable variables

!!$      omp_start=omp_get_wtime()
!!$omp parallel do default(shared)                                 &
!!$omp private(ijk,i,j,k,iw,ie,js,jn,kb,ktp,                       &
!!$omp         ii,jj,kk,cur_ijk,icur,jcur,kcur,                    &
!!$omp         ijpkp,ipjkp,ipjpkp,ijkp,ipjk,ijpk,ipjpk,            &
!!$omp         avg_factor,pcell,vcell,onew,gstencil,vstencil,      &
!!$omp         focus,np,wtp,weightp,ovol,m,nindx)                  &
!!$omp schedule (guided,50)

!Handan Liu modified the following do-loop on Jan 15 2013,
!       again modified on May 9 2013
!!$omp parallel default(shared)                                          &
!!$omp private(ijk,i,j,k,pcell,iw,ie,js,jn,kb,ktp,onew,                          &
!!$omp         ii,jj,kk,cur_ijk,ipjk,ijpk,ipjpk,                                         &
!!$omp         gst_tmp,vst_tmp,velfp,desposnew,ijpkp,ipjkp,                      &
!!$omp         ipjpkp,ijkp,nindx,focus,np,wtp,m,weight_ft,                       &
!!$omp             icur,jcur,kcur,vcell,ovol)
!!$omp do reduction(+:drag_am) reduction(+:drag_bm)
      DO IJK = IJKSTART3,IJKEND3
         IF(.NOT.FLUID_AT(IJK) .OR. PINC(IJK)==0) cycle
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)

! generally a particle may not exist in a ghost cell. however, if the
! particle is adjacent to the west, south or bottom boundary, then pcell
! may be assigned indices of a ghost cell which will be passed to
! set_interpolation_stencil
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = merge(1, k-1, NO_K)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
              ktp,interp_scheme,dimn,ordernew = onew)

! Compute velocity at grid nodes and set the geometric stencil
         DO k = 1, merge(1, ONEW, NO_K)
            DO j = 1,onew
               DO i = 1,onew
                  ii = iw + i-1
                  jj = js + j-1
                  kk = kb + k-1
                  cur_ijk = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk))
                  ipjk    = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk))
                  ijpk    = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk))
                  ipjpk   = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk))
                  GST_TMP(I,J,K,1) = XE(II)
                  GST_TMP(I,J,K,2) = YN(JJ)
                  GST_TMP(I,J,K,3) = merge(DZ(1), ZT(KK), NO_K)
                  VST_TMP(I,J,K,1) = AVG_FACTOR*(U_G(CUR_IJK)+U_G(IJPK))
                  VST_TMP(I,J,K,2) = AVG_FACTOR*(V_G(CUR_IJK)+V_G(IPJK))

                  IF(DO_K) THEN
                     IJPKP   = FUNIJK(IMAP_C(II),JMAP_C(JJ+1),KMAP_C(KK+1))
                     IPJKP   = FUNIJK(IMAP_C(II+1),JMAP_C(JJ),KMAP_C(KK+1))
                     IPJPKP  = FUNIJK(IMAP_C(II+1),JMAP_C(JJ+1),KMAP_C(KK+1))
                     IJKP    = FUNIJK(IMAP_C(II),JMAP_C(JJ),KMAP_C(KK+1))
                     VST_TMP(I,J,K,1) = VST_TMP(I,J,K,1) + &
                     AVG_FACTOR*(U_G(IJKP) + U_G(IJPKP))

                     VST_TMP(I,J,K,2) = VST_TMP(I,J,K,2) + &
                     AVG_FACTOR*(V_G(IJKP) + V_G(IPJKP))

                     VST_TMP(I,J,K,3) = AVG_FACTOR*(W_G(CUR_IJK)+&
                          W_G(IJPK)+W_G(IPJK)+W_G(IPJPK))
                  ELSE
                     VST_TMP(I,J,K,3) = 0.D0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!===================================================================>> Handan Liu
! loop through particles in the cell
! interpolate the fluid velocity (VEL_FP) to the particle's position.
         DO nindx = 1,PINC(IJK)
            NP = PIC(ijk)%p(nindx)
! skipping indices that do not represent particles and ghost particles
            if(.not.pea(np,1)) cycle
            if(pea(np,4)) cycle
            desposnew(:) = des_pos_new(:,np)
            call DRAG_INTERPOLATION(gst_tmp,vst_tmp,desposnew,velfp,weight_ft)
            vel_fp(1:3,np) = velfp(1:3)
!===================================================================>> Handan Liu
!
! Calculate the particle centered drag coefficient (F_GP) using the
! particle velocity and the interpolated gas velocity.  Note F_GP
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            VEL_NEW(:) = DES_VEL_NEW(:,NP)
            CALL DES_DRAG_GP(NP, velfp(1:3), &
               VEL_NEW)
!-----------------------------------------------------------------<<<
! Calculate the corresponding gas solids drag force that is used in
! the gas phase momentum balances.
!----------------------------------------------------------------->>>
            focus = .false.
            M = pijk(np,5)
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)

            DO k = 1, merge(1, ONEW, NO_K)
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

!!$omp critical
                     drag_am(cur_ijk,m) = drag_am(cur_ijk,m) + &
                     !f_gp(np)*weightp(i,j,k)*ovol*wtp          !Handan Liu revised on Jan 15 2013
                     f_gp(np)*weight_ft(i,j,k)*ovol*wtp

                     drag_bm(cur_ijk, 1:3,m) = &
                     drag_bm(cur_ijk,1:3,m) + &
                     f_gp(np) * vel_new(1:3) * &
                     !weightp(i,j,k)*ovol*wtp           !Handan Liu revised on Jan 15 2013
                     weight_ft(i,j,k)*ovol*wtp
!!$omp end critical
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp:',omp_end - omp_start

! At the interface drag_am and drag_bm have to be added
! send recv will be called and the node values will be added
! at the junction. drag_am are drag_bm are altered by the
! routine when periodic boundaries are invoked. so both
! quantities are needed at the time of this call.
      call des_addnodevalues
!-----------------------------------------------------------------<<<
! Calculate/update the cell centered drag coefficient F_GS for use
! in the pressure correction equation
!----------------------------------------------------------------->>>
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
      AVG_FACTOR = merge(0.25d0, 0.125D0, NO_K)

!$omp parallel do default(shared)                               &
!$omp private(ijk,i,j,k,imjk,ijmk,imjmk,ijkm,imjkm,ijmkm,       &
!$omp         imjmkm)                                           &
!$omp schedule (guided,20)
      DO ijk = ijkstart3, ijkend3
         IF(FLUID_AT(IJK)) THEN

            i = i_of(ijk)
            j = j_of(ijk)
            k = k_of(ijk)
            if (i.lt.istart2 .or. i.gt.iend2) cycle
            if (j.lt.jstart2 .or. j.gt.jend2) cycle
            if (k.lt.kstart2 .or. k.gt.kend2) cycle
            imjk = funijk(imap_c(i-1),jmap_c(j),kmap_c(k))
            ijmk = funijk(imap_c(i),jmap_c(j-1),kmap_c(k))
            imjmk = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k))

            IF (.NOT.DES_CONTINUUM_HYBRID) THEN
               f_gs(ijk,:) = avg_factor*&
                  (drag_am(ijk,:)   + drag_am(ijmk,:) +&
                   drag_am(imjmk,:) + drag_am(imjk,:))
            ELSE
               f_gds(ijk,:) = avg_factor*&
                  (drag_am(ijk,:)   + drag_am(ijmk,:) +&
                   drag_am(imjmk,:) + drag_am(imjk,:))
            ENDIF
            IF(DO_K) THEN
               ijkm = funijk(imap_c(i),jmap_c(j),kmap_c(k-1))
               imjkm = funijk(imap_c(i-1),jmap_c(j),kmap_c(k-1))
               ijmkm = funijk(imap_c(i),jmap_c(j-1),kmap_c(k-1))
               imjmkm = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k-1))
               IF(.NOT.DES_CONTINUUM_HYBRID) THEN
                  f_gs(ijk,:) = f_gs(ijk,:) + avg_factor*&
                       (drag_am(ijkm,:) + drag_am(ijmkm,:) +&
                        drag_am(imjmkm,:)+drag_am(imjkm,:) )
               ELSE
                  f_gds(ijk,:) = f_gds(ijk,:) + avg_factor*&
                       (drag_am(ijkm,:) + drag_am(ijmkm,:) +&
                        drag_am(imjmkm,:)+drag_am(imjkm,:) )
               ENDIF
            ENDIF   ! end if
         ELSE   ! else branch of if (fluid_at(ijk))
            IF (DES_CONTINUUM_HYBRID) THEN
               F_GDS(IJK,:) = ZERO
            ELSE
               F_GS(IJK,:) = ZERO
            ENDIF
         ENDIF   ! end if/else (fluid_at(ijk))

      ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do
!-----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE DES_DRAG_GS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GP                                             C
!  Purpose: Calculate the gas-particle drag coefficient using          C
!           the gas velocity interpolated to the particle position     C
!           and the particle velocity.                                 C
!           Invoked from des_drag_gs and calc_des_drag_gs              C
!                                                                      C
!  Comments: The BVK drag model and all drag models with the           C
!            polydisperse correction factor (i.e., suffix _PCF)        C
!            require an average particle diameter. This has been       C
!            loosely defined for discrete particles based on their     C
!            solids phase                                              C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_DRAG_GP(LL, FLUID_VEL, PARTICLE_VEL)

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
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) :: LL
! fluid velocity interpolated to particle position
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: FLUID_VEL
! particle velocity
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: PARTICLE_VEL
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices, associated with current particle
      INTEGER :: IJK
! solids phase index, associated with current particle
      INTEGER :: M
! magnitude of gas-solids relative velocity
      DOUBLE PRECISION :: VREL
! gas laminar viscosity redefined here to set viscosity at pressure
! boundaries
      DOUBLE PRECISION :: Mu
! drag coefficient
      DOUBLE PRECISION :: DgA
! current value of F_gs (i.e., without underrelaxation)
      DOUBLE PRECISION F_gstmp
! indices of solids phases (continuous, discrete)
      INTEGER :: CM, DM, L
! temporary shift of total number of solids phases to account for both
! discrete and continuous solids phases used for the hybrid mdoel
      INTEGER :: MAXM
! tmp local variable for the particle diameter of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: DP_loc(2*DIM_M)
! tmp local variable for the solids volume fraction of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: EPs_loc(2*DIM_M)
! tmp local variable for the particle density of solids
! phase M (continuous or discrete)
      DOUBLE PRECISION :: ROs_loc(2*DIM_M)
! correction factors for implementing polydisperse drag model
! proposed by van der Hoef et al. (2005)
      DOUBLE PRECISION :: F_cor, tmp_sum, tmp_fac
! average particle diameter in polydisperse systems
      DOUBLE PRECISION :: DPA
! diameter ratio in polydisperse systems
      DOUBLE PRECISION :: Y_i
! total solids volume fraction
      DOUBLE PRECISION :: phis
! aliases for void fraction, gas density, gas bulk density,
! solids volume fraction, particle diameter, particle density
      DOUBLE PRECISION :: EPG, ROg, ROPg, EP_SM, DPM, ROs
!-----------------------------------------------

! values based on current particle
      IJK = PIJK(LL,4)
! solids phase index of current particle
      M = PIJK(LL,5)

! Assign local variables DP_loc, EPs_loc, and MAXM.  These
! represent arrays for the particle diameter, solids volume
! fraction, and number of particle types (i.e., phases).
      IF (.NOT.DES_CONTINUUM_HYBRID) THEN
         MAXM = DES_MMAX
         DO DM = 1,MAXM
            DP_loc(DM) = DES_D_p0(DM)
            EPs_loc(DM) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
            ROs_loc(DM) = DES_RO_S(DM)
         ENDDO
      ELSE   ! des_continuum_hybrid branch
! For the hybrid model the diameters and solids volume fractions of
! of both discrete and continuous are stored in this single quantity.
! Any loops of solids phases will include all solids phases (discrete
! and continuum)
         MAXM = SMAX + DES_MMAX
! populate DP, EPS starting with discrete phases
         DO DM = 1,DES_MMAX
            DP_loc(DM) = DES_D_p0(DM)
            EPs_loc(DM) = DES_ROP_S(IJK,DM)/DES_RO_S(DM)
            ROs_loc(DM) = DES_RO_S(DM)
         ENDDO
         DO CM = 1,SMAX
            L = DES_MMAX + CM
            DP_loc(L) = D_P(IJK,CM)
            EPs_loc(L) = EP_S(IJK,CM)
            ROs_loc(L) = RO_S(IJK,CM)
         ENDDO
      ENDIF   ! end if/else (.not.des_continuum_hybrid)


! magnitude of gas-particle relative velocity
      IF(NO_K)THEN
         VREL = SQRT((FLUID_VEL(1) - PARTICLE_VEL(1))**2 +&
                     (FLUID_VEL(2) - PARTICLE_VEL(2))**2)
      ELSE
         VREL = SQRT((FLUID_VEL(1) - PARTICLE_VEL(1))**2 +&
                     (FLUID_VEL(2) - PARTICLE_VEL(2))**2 +&
                     (FLUID_VEL(3) - PARTICLE_VEL(3))**2)
      ENDIF

! Laminar viscosity at a pressure boundary is given the value of the
! fluid cell next to it. This applies just to the calculation of the
! drag, in other routines the value of viscosity at a pressure boundary
! always has a zero value.
! This will never happen since this subroutine is currently only called
! for fluid_at cells (does not include flow boundaries)
! This points to an inconsitency in calculation of drag between
! continuum and discrete models that is probably not addressed in the
! solution of the gas phase momentum balances
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

! calculate the total solids volume fraction
      phis = ZERO
      DO L = 1, MAXM
! this is slightly /= one-ep_g due to round-off
         phis = phis + EPs_loc(L)
      ENDDO

! calculate the average paricle diameter and particle ratio
      DPA = ZERO
      tmp_sum = ZERO
      tmp_fac = ZERO
      DO L = 1, MAXM
         IF (phis .GT. ZERO) THEN
            tmp_fac = EPs_loc(L)/phis
            tmp_sum = tmp_sum + tmp_fac/DP_loc(L)
          ELSE
            tmp_sum = tmp_sum + ONE/DP_loc(L) ! not important, but will avoid NaN's in empty cells
          ENDIF
      ENDDO
      DPA = ONE / tmp_sum
      Y_i = DP_loc(M) * tmp_sum

! assign variables for short dummy arguments
      EPg = EP_G(IJK)
      ROg = RO_G(IJK)
      ROPg = ROP_G(IJK)
      EP_SM = EPs_loc(M)
      DPM = DP_loc(M)
      ROs = ROs_loc(M)

! determine the drag coefficient
      IF (EP_SM <= ZERO) THEN
! this won't happen in DEM case since routine is performed over
! particles not cells as in continuum case
         DgA = ZERO
      ELSEIF (EPg == ZERO) THEN
! this case will already be caught in most drag subroutines whenever
! RE==0 (for correlations in which RE includes EPg). however, this will
! prevent potential divisions by zero in some models by setting it now.
         DgA = ZERO
      ELSE
! determine the drag coefficient
         SELECT CASE(DRAG_TYPE_ENUM)
         CASE (SYAM_OBRIEN)
            CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,DPM)
         CASE (GIDASPOW)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)
         CASE (GIDASPOW_PCF)
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (GIDASPOW_BLEND)
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPM)
         CASE (GIDASPOW_BLEND_PCF)
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,DPA)
         CASE (WEN_YU)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPM)
         CASE (WEN_YU_PCF)
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,DPA)
         CASE (KOCH_HILL)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPM,phis)
         CASE (KOCH_HILL_PCF)
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (BVK)
            CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,DPM,DPA,phis)
         CASE (USER_DRAG)
            CALL DRAG_USR(IJK, M, DgA, EPg, Mu, ROg, VREL, DPM, ROs)
         CASE DEFAULT
            CALL START_LOG
            IF(DMP_LOG) WRITE (*, '(A,A)') &
               'Unknown DRAG_TYPE: ', DRAG_TYPE
            WRITE (UNIT_LOG, '(A,A)') 'Unknown DRAG_TYPE: ', DRAG_TYPE
            CALL END_LOG
            CALL mfix_exit(myPE)
         END SELECT   ! end selection of drag_type
      ENDIF   ! end if/elseif/else (ep_sm <= zero, ep_g==0)


! Modify drag coefficient to account for possible corrections and
! for differences between Model B and Model A
      IF(DRAG_TYPE_ENUM == GIDASPOW_PCF .OR. &
         DRAG_TYPE_ENUM == GIDASPOW_BLEND_PCF .OR. &
         DRAG_TYPE_ENUM == WEN_YU_PCF .OR. &
         DRAG_TYPE_ENUM == KOCH_HILL_PCF .OR. &
         DRAG_TYPE_ENUM == BVK) THEN
! see erratum by Beetstra et al. (2007) : the correction factor differs
! for model A versus model B.
! application of the correction factor for model A is found from
! the correction factor for model B and neglects the Y_i**3 term
         IF(Model_B) THEN
            IF (M == 1) THEN
               F_cor = (EPg*Y_i + phis*Y_i**2)
            ELSE
               F_cor = (EPg*Y_i + phis*Y_i**2 + &
                  0.064d0*EPg*Y_i**3)
            ENDIF
         ELSE
            F_cor = Y_i
         ENDIF
         DgA = ONE/(Y_i*Y_i) * DgA * F_cor
      ENDIF

! Calculate the drag coefficient (Model B coeff = Model A coeff/EP_g)
      IF(MODEL_B) THEN
         F_gstmp = DgA * PVOL(LL)/EP_G(IJK)
      ELSE
         F_gstmp = DgA * PVOL(LL)
      ENDIF

! Determine drag force coefficient accounting for any under relaxation
! f_gp() =  single particle drag excluding vector(v_g - v_p)
      F_gp(LL) = (ONE - UR_F_gs) * F_gp(LL) + UR_F_gs * F_gstmp

      RETURN
      END SUBROUTINE DES_DRAG_GP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Subroutine: DRAG_INTERPOLATION                                       C
!  Purpose: DES - Calculate the fluid velocity interpolated at the      C
!           particle's location and weights. Replace 'interpolator'     C
!                       interface for OpenMP implementation.            C
!                                                                       C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_INTERPOLATION(GSTEN,VSTEN,DESPOS,VELFP,WEIGHTFACTOR)

      use geometry, only: NO_K

        IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: GSTEN
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: VSTEN
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: DESPOS
        DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: VELFP
        DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(OUT) :: WEIGHTFACTOR
        INTEGER :: II, JJ, KK

        DOUBLE PRECISION, DIMENSION(2) :: XXVAL, YYVAL, ZZVAL
        DOUBLE PRECISION :: DXX, DYY, DZZ
        DOUBLE PRECISION, DIMENSION(3) :: ZETAA

        DXX = GSTEN(2,1,1,1) - GSTEN(1,1,1,1)
        DYY = GSTEN(1,2,1,2) - GSTEN(1,1,1,2)

        ZETAA(1:2) = DESPOS(1:2) - GSTEN(1,1,1,1:2)

        ZETAA(1) = ZETAA(1)/DXX
        ZETAA(2) = ZETAA(2)/DYY

        XXVAL(1)=1-ZETAA(1)
        YYVAL(1)=1-ZETAA(2)
        XXVAL(2)=ZETAA(1)
        YYVAL(2)=ZETAA(2)

        VELFP(:) = 0.D0

        IF(NO_K) THEN
           DO JJ=1,2
              DO II=1,2
                 WEIGHTFACTOR(II,JJ,1) = XXVAL(II)*YYVAL(JJ)
                 VELFP(1:2) = VELFP(1:2) + VSTEN(II,JJ,1,1:2)*WEIGHTFACTOR(II,JJ,1)
              ENDDO
           ENDDO
        ELSE
           DZZ = GSTEN(1,1,2,3) - GSTEN(1,1,1,3)
           ZETAA(3) = DESPOS(3) - GSTEN(1,1,1,3)
           ZETAA(3) = ZETAA(3)/DZZ
           ZZVAL(1)=1-ZETAA(3)
           ZZVAL(2)=ZETAA(3)
           DO KK=1,2
              DO JJ=1,2
                 DO II=1,2
                    WEIGHTFACTOR(II,JJ,KK) = XXVAL(II)*YYVAL(JJ)*ZZVAL(KK)
                    VELFP(1:3) = VELFP(1:3) + VSTEN(II,JJ,KK,1:3)*WEIGHTFACTOR(II,JJ,KK)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

      END SUBROUTINE DRAG_INTERPOLATION
