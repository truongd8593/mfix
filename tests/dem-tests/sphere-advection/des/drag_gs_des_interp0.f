!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS                                        C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions:                             C
!     - If interpolated, then it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. This F_GP is used to calculate    C
!       the fluid-solids drag on the particle.                         C
!     - The total contact force on the particle is then updated to     C
!       include the gas-solids drag force and gas pressure force       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_GS_DES_INTERP0

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
         ENDDO       ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!$omp end parallel do


      RETURN
      END SUBROUTINE DRAG_GS_DES_INTERP0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DRAG_GS                                             C
!  Purpose: This subroutine is only called from the CONTINUUM side.    C
!     It performs the following functions:                             C
!     - If interpolated then, it calculates the fluid-particle         C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. It then determines the            C
!       the contributions of fluid-particle drag to the center         C
!       coefficient of the A matrix and the b (source) vector in the   C
!       matrix equation (A*VEL_FP=b) equation for the fluid phase      C
!       x, y and z momentum balances using F_GP.                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_GS_GAS_INTERP0

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


!!$omp parallel default(shared)                                        &
!!$omp private(ijk,i,j,k,pcell,iw,ie,js,jn,kb,ktp,onew,                &
!!$omp         ii,jj,kk,cur_ijk,ipjk,ijpk,ipjpk,                       &
!!$omp         gst_tmp,vst_tmp,velfp,desposnew,ijpkp,ipjkp,            &
!!$omp         ipjpkp,ijkp,nindx,focus,np,wtp,m,weight_ft,             &
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
                     drag_am(cur_ijk) = drag_am(cur_ijk) + &
                        f_gp(np)*weight_ft(i,j,k)*ovol*wtp

                     drag_bm(cur_ijk,1:3) = &
                        drag_bm(cur_ijk,1:3) + &
                        f_gp(np) * vel_new(1:3) * &
                        weight_ft(i,j,k)*ovol*wtp
!!$omp end critical
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel


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
                  (drag_am(ijk)   + drag_am(ijmk) +&
                   drag_am(imjmk) + drag_am(imjk))
            ELSE
               f_gds(ijk,:) = avg_factor*&
                  (drag_am(ijk)   + drag_am(ijmk) +&
                   drag_am(imjmk) + drag_am(imjk))
            ENDIF
            IF(DO_K) THEN
               ijkm = funijk(imap_c(i),jmap_c(j),kmap_c(k-1))
               imjkm = funijk(imap_c(i-1),jmap_c(j),kmap_c(k-1))
               ijmkm = funijk(imap_c(i),jmap_c(j-1),kmap_c(k-1))
               imjmkm = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k-1))
               IF(.NOT.DES_CONTINUUM_HYBRID) THEN
                  f_gs(ijk,:) = f_gs(ijk,:) + avg_factor*&
                       (drag_am(ijkm) + drag_am(ijmkm) +&
                        drag_am(imjmkm)+drag_am(imjkm) )
               ELSE
                  f_gds(ijk,:) = f_gds(ijk,:) + avg_factor*&
                       (drag_am(ijkm) + drag_am(ijmkm) +&
                        drag_am(imjmkm)+drag_am(imjkm) )
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

      RETURN
      END SUBROUTINE DRAG_GS_GAS_INTERP0
