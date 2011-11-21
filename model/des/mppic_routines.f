      SUBROUTINE MPPIC_COMPUTE_PS_GRAD
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
      USE constant 
      USE cutcell 
      USE interpolation
      USE mfix_pic
      implicit none 
      
      ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IDIM, M
      
      ! temporary variables used to calculate pressure at scalar cell edge      
      DOUBLE PRECISION TEMP1, TEMP2, avg_factor 

      integer :: korder, iw,ie,js,jn,kb,ktp, onew, pcell(3), cur_ijk, NP, nindx
      
      integer :: ii,jj,kk, ipjpk, ijpkp, ipjkp, ipjpkp, I1, J1, K1
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      
      if(MPPIC_SOLID_STRESS_SNIDER) then 

         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN
               P_S(IJK,1) = PSFAC_FRIC_PIC * ((1.d0 - EP_G(IJK))**FRIC_EXP_PIC)/ MAX(EP_G(IJK) - EP_STAR, FRIC_NON_SING_FAC*EP_G(IJK))
               !write(102,'(2(2x,i4),2(2x,g17.8))') I_OF(IJK), J_OF(IJK), EP_S(IJK,1), P_STAR(IJK) 
            ELSE
               !So that ghost cells have higher pressure 
               P_S(IJK,1) = PSFAC_FRIC_PIC * ((1.d0 - 0.1d0)**FRIC_EXP_PIC)/ MAX(0.1d0 - EP_STAR, FRIC_NON_SING_FAC*0.1d0)

            ENDIF
         ENDDO


      ELSE 
         DO IJK = IJKSTART3, IJKEND3
            PS_FORCE_PIC(IJK,:) = ZERO 
            
            IF(FLUID_AT(IJK)) THEN
               if(EP_G(IJK).lt.ep_star) then 
                  P_S(IJK,1) = one*(one-ep_g(ijk)) 
               else 
               
                  P_S(IJK,1) = ZERO
               endif
               
            ELSE
               P_S(IJK,1) = 1.!\*0.d0
            ENDIF
         ENDDO
         
      ENDIF
      
      CALL SEND_RECV(P_S,1)

      !Since EP_G is already shared across the processors, the pressure gradient calculation 
      !can be made a function call so that the extra communication of P_S can be avoided. 

      !DO k = kstart2, kend1 
      !do j = jstart2, jend1
      !do i = istart2, iend1 
      DO IJK = IJKSTART3, IJKEND3
         PS_FORCE_PIC(IJK,:) = ZERO 
         IF(.NOT.FLUID_AT(IJK)) CYCLE 

               !I = I_OF(IJK)
               !J = J_OF(IJK)
               !K = K_OF(IJK)
         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)
         PS_FORCE_PIC(IJK,1) = 2.d0*(P_S(IPJK,1) - P_S(IJK,1))/(DX(I)+DX(I_of(ipjk)))
         
         PS_FORCE_PIC(IJK,2) = 2.d0*(P_S(IJPK,1) - P_S(IJK,1))/(DY(j)+Dy(j_of(ijpk)))
         if(dimn.eq.3) PS_FORCE_PIC(IJK,3) = 2.d0*(P_S(IJKP,1) - P_S(IJK,1))/(Dz(k)+Dz(k_of(ijkp)))
         
      ENDDO
      !Rahul:
      !the above will not compute pressure gradients normal to the east, south and bottom faces 
      !which are very important
      I1 = 1 
      DO K1 = KSTART3, KEND3
         DO J1 = JSTART3, JEND3
!//
            IF(I1.NE.ISTART2)   EXIT
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
            IJK = FUNIJK(I1,J1,K1) 
            IPJK = IP_OF(IJK)
            PS_FORCE_PIC(IJK,1) = 2.d0*(P_S(IPJK,1) - P_S(IJK,1))/(DX(I)+DX(I_of(ipjk)))
                     
         ENDDO
      ENDDO
      J1 = 1 
      DO K1 = KSTART3, KEND3
         DO I1 = ISTART3, IEND3
!//
            IF(J1.NE.JSTART2)   EXIT
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
            IJK = FUNIJK(I1,J1,K1) 
            IJPK = JP_OF(IJK)
            PS_FORCE_PIC(IJK,2) = 2.d0*(P_S(IJPK,1) - P_S(IJK,1))/(DY(j)+Dy(j_of(ijpk)))

         END DO 
      END DO 

      IF(DIMN.eq.3) then 
         K1 = 1 
         DO J1 = JSTART3, JEND3
            DO I1 = ISTART3, IEND3
               IF(K1.NE.KSTART2)   EXIT
               IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
               IJK = FUNIJK(I1,J1,K1) 
               IJKP = KP_OF(IJK)
               PS_FORCE_PIC(IJK,3) = 2.d0*(P_S(IJKP,1) - P_S(IJK,1))/(Dz(k)+Dz(k_of(ijkp)))
            END DO 
         END DO 
      ENDIF

      DO IDIM = 1, DIMN
         CALL SEND_RECV(PS_FORCE_PIC(:,IDIM),1)
      ENDDO

      CALL SET_INTERPOLATION_SCHEME(2)

      KORDER = 1+(DIMN-2)

      do ijk = ijkstart3,ijkend3
         
         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle 
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = (3-dimn)*1+(dimn-2)*(k-1) ! =k-1 (in 3d) or =1 (in 2d)
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
                  if(dimn.eq.3) then 
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                  endif
                  gstencil(i,j,k,1) = xe(ii)
                  gstencil(i,j,k,2) = yn(jj)
                  gstencil(i,j,k,3) = zt(kk)*(dimn-2) + dz(1)*(3-dimn)
                  
                  
                  DO M = 1, MMAX
                     VEL_SOL_STENCIL(i,j,k,1,M) = avg_factor*(u_s(cur_ijk,M)+u_s(ijpk,M))
                     VEL_SOL_STENCIL(i,j,k,2,M) = avg_factor*(v_s(cur_ijk,M)+v_s(ipjk,M)) 
                  ENDDO
                  psgradstencil(i,j,k,1) = avg_factor*(PS_FORCE_PIC(cur_ijk,1)+PS_FORCE_PIC(ijpk,1))
                  psgradstencil(i,j,k,2) = avg_factor*(PS_FORCE_PIC(cur_ijk,2)+PS_FORCE_PIC(ipjk,2)) 

                  vstencil(i,j,k,1) = avg_factor*(u_g(cur_ijk)+u_g(ijpk))
                  vstencil(i,j,k,2) = avg_factor*(v_g(cur_ijk)+v_g(ipjk)) 
                  if(dimn.eq.3) then 
                     psgradstencil(i,j,k,1) = psgradstencil(i,j,k,1)+avg_factor*(PS_FORCE_PIC(ijkp,1) + PS_FORCE_PIC(ijpkp,1))
                     psgradstencil(i,j,k,2) = psgradstencil(i,j,k,2)+avg_factor*(PS_FORCE_PIC(ijkp,2) + PS_FORCE_PIC(ipjkp,2))
                     psgradstencil(i,j,k,3) = avg_factor*(PS_FORCE_PIC(cur_ijk,3)+&
                     PS_FORCE_PIC(ijpk,3)+PS_FORCE_PIC(ipjk,3)+PS_FORCE_PIC(ipjpk,3))
                     
                     vstencil(i,j,k,1) = vstencil(i,j,k,1)+avg_factor*(u_g(ijkp) + u_g(ijpkp))
                     vstencil(i,j,k,2) = vstencil(i,j,k,2)+avg_factor*(v_g(ijkp) + v_g(ipjkp))
                     vstencil(i,j,k,3) = avg_factor*(w_g(cur_ijk)+&
                     & w_g(ijpk)+w_g(ipjk)+w_g(ipjpk))

                     DO M = 1, MMAX     
                        VEL_SOL_STENCIL(i,j,k,1, M) = VEL_SOL_STENCIL(i,j,k,1,M)+avg_factor*(u_s(ijkp,M) + u_s(ijpkp,M))
                        VEL_SOL_STENCIL(i,j,k,2, M) = VEL_SOL_STENCIL(i,j,k,2,M)+avg_factor*(v_s(ijkp,M) + v_s(ipjkp,M))
                        VEL_SOL_STENCIL(i,j,k,3, M) = avg_factor*(w_s(cur_ijk,M)+&
                        w_s(ijpk,M)+w_s(ipjk,M)+w_s(ipjpk,M))
                     ENDDO
                  else 
                     psgradstencil(i,j,k,3) = 0.d0
                     VEL_SOL_STENCIL(i,j,k,3, 1:MMAX) = 0.d0
                     vstencil(i,j,k,3) = 0.d0
                                          
                  endif
                  
               enddo
            enddo
         enddo
         
                  !loop through particles in the cell  
         do nindx = 1,pinc(ijk)
            np = pic(ijk)%p(nindx)
            m = pijk(np,5)
                        
            if(dimn.eq.2) then 
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
               psgradstencil(1:onew,1:onew,1,1:dimn), &
               des_pos_new(np,1:dimn),PS_GRAD(np,1:dimn),  &
               onew,interp_scheme,weightp)
            else 
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
               psgradstencil(1:onew,1:onew,1:onew,1:dimn), &
               des_pos_new(np,1:dimn),PS_GRAD(np,1:dimn),  &
               onew,interp_scheme,weightp)
            endif
            
            do idim = 1, dimn
               AVGSOLVEL_P(NP,IDIM) = ARRAY_DOT_PRODUCT(VEL_SOL_STENCIL(:,:,:,IDIM,M),WEIGHTP(:,:,:))
               VEL_FP(NP,IDIM) = ARRAY_DOT_PRODUCT(VSTENCIL(:,:,:,IDIM),WEIGHTP(:,:,:))
            ENDDO
            
         END DO
      END DO
    

      END SUBROUTINE MPPIC_COMPUTE_PS_GRAD


      SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS 
      
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


 ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
              IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK, &
! Pradeep added following 
              IMJMK,IMJKM,IJMKM,IMJMKM
      INTEGER I1, I2, J1, J2, K1, K2, IDIM, IJK2
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
      DOUBLE PRECISION AVG_FACTOR, TEMP1, AVG_FACTOR_FACE

! index of solid phase that particle NP belongs to      
      INTEGER M

! volume of fluid cell particle resides in
      DOUBLE PRECISION VCELL, VCELL2 
! one over the solids volume fraction and one over the volume       
      DOUBLE PRECISION OEPS, OVOL, MASS_SOL1, MASS_SOL2, MASS_SOL3, MASS_SOL4

! particle number index, used for looping      
      INTEGER NP

! index to track accounted for particles 
      INTEGER PC 

! for error messages      
      INTEGER IER

! Statistical weight of the particle. Equal to one for DEM 
   
      DOUBLE PRECISION WTP, ZCOR, JUNK_VAL(3)

! Pradeep temporary indices for periodic boundary adjustments
      integer korder,cur_ijk,nindx
! Pradeep introducing volume at grid nodes for backward interpolation 
      double precision VELG_ARR(DIMN), VELS_ARR(DIMN, MMAX), np_m

! see the discussion for IJK_U ..... in comments       
      INTEGER  IJK_U, IJK_V, IJK_W, ICUR, JCUR, KCUR

      character*100 :: filename
      integer  FLUID_IND      
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      
      integer :: epg_min_loc(1), count_nodes_outside, count_nodes_inside, count_nodes_inside_max, COUNT_TEMP
      double precision :: epg_min2
      
      double precision :: RESID_WTBAR(MMAX), RESID_BM(DIMN, MMAX), NORM_FACTOR
!
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'


     
      MASS_SOL1 = zero
      MASS_SOL2 = zero 
      MASS_SOL3 = zero 
      MASS_SOL4 = zero 
      call set_interpolation_scheme(2)
      korder = 1+(dimn-2)
      drag_bm = ZERO
      drag_am = zero 
      wtbar = zero
      if(dimn.eq.2) count_nodes_inside_max = 4
      if(dimn.eq.3) count_nodes_inside_max = 8
      IJKLOOP: DO IJK = IJKSTART3,IJKEND3
         ROP_S(IJK,:) = zero 
         DES_U_S(IJK, :) = ZERO
         DES_V_S(IJK, :) = ZERO
         if(dimn.eq.3) DES_W_S(IJK, :) = ZERO
         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle 
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = (3-dimn)*1+(dimn-2)*(k-1) ! =k-1 (in 3d) or =1 (in 2d)
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
         ktp,interp_scheme,dimn,ordernew = onew) 

!Compute velocity at grid nodes and set the geometric stencil 
         avg_factor = 0.25d0*(dimn-2) + 0.5d0*(3-dimn)
         count_nodes_outside = 0 
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
                  if(dimn.eq.3) then 
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                  endif
                  gstencil(i,j,k,1) = xe(ii)
                  gstencil(i,j,k,2) = yn(jj)
                  gstencil(i,j,k,3) = zt(kk)*(dimn-2) + dz(1)*(3-dimn)
                  vstencil(i,j,k,:) = zero 
                  if(cartesian_grid) then 
                     if(scalar_node_atwall(cur_ijk)) count_nodes_outside = &
                     & count_nodes_outside + 1 
                  endif
               enddo
            enddo
         enddo

         count_nodes_inside = count_nodes_inside_max - count_nodes_outside

         !loop through particles in the cell  
         do nindx = 1,pinc(ijk)
            focus = .false.
            np = pic(ijk)%p(nindx)
            m = pijk(np,5)
            
            NP_M = ROP_SO(IJK,m)*VOL(IJK)/(ro_s(m)*pvol(np))

            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            
            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP
            if (dimn .eq. 2) then 
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
                    vstencil(1:onew,1:onew,1,1:dimn), &
                    des_pos_new(np,1:dimn),JUNK_VAL(1:dimn),  &
                    onew,interp_scheme,weightp)
            else 
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
                    vstencil(1:onew,1:onew,1:onew,1:dimn), &
                    des_pos_new(np,1:dimn),JUNK_VAL(1:dimn),  &
                    onew,interp_scheme,weightp)
            end if
            
            do k = 1, (3-dimn)*1+(dimn-2)*onew
               do j = 1, onew
                  do i = 1, onew
                     ! shift loop index to new variables for manipulation                     
                     ii = iw + i-1
                     jj = js + j-1
                     kk = kb + k-1

                     icur = imap_c(ii)
                     jcur = jmap_c(jj)
                     kcur = kmap_c(kk)
                     
                     cur_ijk = funijk(icur, jcur, kcur) !imap_c(ii),jmap_c(jj),kmap_c(kk))
                     
! should this volume be for current ijk index or always particle index?                        
                     vcell = des_vol_node(cur_ijk)
                     ovol = one/vcell
                     
                     temp1 = weightp(i,j,k)*ro_s(m)*pvol(np)*wtp!*ovol
                                         
                     wtbar(cur_ijk,m) = wtbar(cur_ijk,m) + temp1 
                     
                     drag_bm(cur_ijk, 1:dimn,m) = &
                     drag_bm(cur_ijk,1:dimn,m) + temp1*des_vel_new(np, 1:dimn)
                     
                     MASS_SOL3 = MASS_SOL3 +  TEMP1

                  enddo
               enddo
            enddo
            
         enddo          ! pinc(ijk) loop 
         
         
         if(count_nodes_inside.lt.count_nodes_inside_max) then 
            i = i_of(ijk)
            j = j_of(ijk)
            k = k_of(ijk)
            
            I1 = I-1
            I2 = I
            J1 = J-1
            J2 = J
            
            IF(NO_K) THEN 
               K1 = K
               K2 = K
            ELSE
               K1 = K-1
               K2 = K
            ENDIF
         !Convention used to number node numbers is described below 
         
         ! i=1, j=2           i=2, j=2
         !   _____________________
         !   |                   |
         !   |  I = 2, J = 2     |
         !   |___________________|
         ! i=1, j=1           i=2, j=1
         !first calculate the residual wtbar and drag_bm that was computed
         !on nodes that do not belong to the domain
            RESID_WTBAR(1:MMAX) = ZERO 
            RESID_BM(1:DIMN, 1:MMAX) = ZERO
            DO KK = K1, K2
               DO JJ = J1, J2
                  DO II = I1, I2
                     IJK2 = funijk(II, JJ, KK) 
                     IF(SCALAR_NODE_ATWALL(IJK2)) THEN 
                        RESID_WTBAR(1:MMAX) = RESID_WTBAR(1:MMAX) + WTBAR(IJK2,1:MMAX)
                        WTBAR(IJK2,1:MMAX) = ZERO 
                        DO IDIM = 1, DIMN
                           RESID_BM(IDIM, 1:MMAX) = RESID_BM(IDIM, 1:MMAX) +  & 
                           & DRAG_BM(IJK2,IDIM, 1:MMAX)
                           DRAG_BM(IJK2,IDIM, 1:MMAX) = ZERO 
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
                  
                     
           ! Now add this residual equally to the remaining nodes
           ! Pradeep: at the interface drag_am,drag_bm,wtbar has to be added
           ! send recv will be called and the node values will be added 
           ! at the junction 
            
            NORM_FACTOR = one/real(count_nodes_inside)          
            COUNT_TEMP = 0
            DO KK = K1, K2
               DO JJ = J1, J2
                  DO II = I1, I2
                     IJK2 = funijk(II, JJ, KK) 
                     IF(.NOT.SCALAR_NODE_ATWALL(IJK2)) THEN 
                        COUNT_TEMP = COUNT_TEMP + 1
                        WTBAR(IJK2,1:MMAX) = WTBAR(IJK2,1:MMAX) + RESID_WTBAR(1:MMAX)*NORM_FACTOR
                        DO IDIM = 1, DIMN
                           DRAG_BM(IJK2,IDIM, 1:MMAX) = DRAG_BM(IJK2,IDIM, 1:MMAX) + RESID_BM(IDIM, 1:MMAX)*NORM_FACTOR
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            
            !WRITE(*,*) 'hello: NODES INSIDE AND OUTSIDE', count_nodes_inside, count_nodes_inside_max, one/real(count_temp), norm_factor
         ENDIF
         

      END DO IJKLOOP            ! IJK LOOP
      
      MASS_SOL4 = SUM(WTBAR(:,1))
      CALL DES_ADDNODEVALUES
      

      DO IJK = IJKSTART3, IJKEND3
         DO M = 1, MMAX
            IF(WTBAR(IJK,M).GT.ZERO) then
               DRAG_BM(IJK, :, M) = DRAG_BM(IJK, :, M)/WTBAR(IJK, M)
            endif
            !now correct wtbar for nodes that fall on the boundary 
            WTBAR(IJK, M) = WTBAR(IJK, M)*DES_VOL_NODE_RATIO(IJK)
         ENDDO
      ENDDO

      AVG_FACTOR =  0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)
      AVG_FACTOR_FACE = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)

      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IMJK = FUNIJK(I-1,J,K)
            IJMK = FUNIJK(I,J-1,K)
            IMJMK = FUNIJK(I-1,J-1,K) 

            ROP_S(IJK,:) = AVG_FACTOR*(WTBAR(IJK,:)+&
            WTBAR(IJMK,:) + WTBAR(IMJMK,:) +&
            WTBAR(IMJK,:))

            
            DES_U_S(IJK,:) = avg_factor*(DRAG_BM(IJK,1,:) + &
            DRAG_BM(IJMK,1,:) + DRAG_BM(IMJMK,1,:) + &
            DRAG_BM(IMJK,1,:))

            DES_V_S(IJK,:) = avg_factor*(DRAG_BM(IJK,2,:) + &
            DRAG_BM(IJMK,2,:) + DRAG_BM(IMJMK,2,:) + &
            DRAG_BM(IMJK,2,:))

            U_S(IJK,:) = AVG_FACTOR_FACE*(DRAG_BM(IJK,1,:) + DRAG_BM(IJMK, 1,:))
            V_S(IJK,:) = AVG_FACTOR_FACE*(DRAG_BM(IJK,2,:) + DRAG_BM(IMJK, 2,:))
            if(dimn.eq.3) then 
               IJKM = FUNIJK(I,J,K-1)
               IJMKM = FUNIJK(I,J-1,K-1)
               IMJKM = FUNIJK(I-1,J,K-1)
               IMJMKM = FUNIJK(I-1,J-1,K-1)


               ROP_S(IJK,:) = ROP_S(IJK,:) + AVG_FACTOR*&
               (WTBAR(IJKM,:) + WTBAR(IJMKM,:) + &
               WTBAR(IMJMKM,:) + WTBAR(IMJKM,:))

               DES_U_S(ijk,:) = DES_U_S(ijk,:) + avg_factor*&
               (DRAG_BM(ijkm,1,:) + DRAG_BM(ijmkm,1,:) &
               + DRAG_BM(imjmkm,1,:)+DRAG_BM(imjkm,1,:))
               
               
               DES_V_S(ijk,:) = DES_V_S(ijk,:) + avg_factor*&
               (DRAG_BM(ijkm,2,:) + DRAG_BM(ijmkm,2,:) &
               + DRAG_BM(imjmkm,2,:)+DRAG_BM(imjkm,2,:) )
               
               DES_W_S(ijk,:) =  avg_factor*&
               (DRAG_BM(ijk,3,:) + DRAG_BM(ijmk,3,:) &
               + DRAG_BM(imjmk,3,:)+DRAG_BM(imjk,3,:) )
               
               DES_W_S(ijk,:) = DES_W_S(ijk,:) + avg_factor*&
               (DRAG_BM(ijkm,3,:) + DRAG_BM(ijmkm,3,:) &
               + DRAG_BM(imjmkm,3,:)+DRAG_BM(imjkm,3,:) )
               
               
               U_S(IJK, :) = U_S(IJK,:) + AVG_FACTOR_FACE*(DRAG_BM(IJKM,1,:) + DRAG_BM(IJMKM, 1,:))
               V_S(IJK, :) = V_S(IJK,:) + AVG_FACTOR_FACE*(DRAG_BM(IJKM,2,:) + DRAG_BM(IMJKM, 2,:))
               W_S(IJK, :) = AVG_FACTOR_FACE*(DRAG_BM(IJK,3,:) + DRAG_BM(IMJK, 3,:) + DRAG_BM(IJMK,3,:) + DRAG_BM(IMJMK,3,:))

               
            ENDIF
!!$            DO M = 1, MMAX
!!$               IF(ROP_S(IJK,M).GT.ZERO) then 
!!$                  DES_U_S(IJK,M) = DES_U_S(IJK,M)/ROP_S(IJK,M)
!!$                  DES_V_S(IJK,M) = DES_V_S(IJK,M)/ROP_S(IJK,M)
!!$                  if(DIMN.eq.3) DES_W_S(IJK,M) = DES_W_S(IJK,M)/ROP_S(IJK,M)
!!$
!!$               ENDIF
!!$            ENDDO

            ROP_S(IJK,:) = ROP_S(IJK,:)/VOL(IJK)

            EP_G(IJK) = ONE   
            
            !MASS_SOL3 = MASS_SOL3 + ROP_SO(IJK, 1)*VOL(IJK)
            DO M = 1, MMAX
               
               IF(ROP_S(IJK,M) > ZERO) THEN
                  
                  MASS_SOL2 = MASS_SOL2 + ROP_S(IJK, M)*VOL(IJK)
                  
                  
                  IF(.not.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
                  
                  ROP_SO(IJK,M)  = ROP_S(IJK,M) 
                  
                  IF(EP_G(IJK).LT.ZERO .AND. DES_CONTINUUM_COUPLED) THEN 
                                ! this does not matter if pure granular flow simulation (i.e. no fluid)
                     IF(DMP_LOG)  WRITE(UNIT_LOG, 2000) EP_G(IJK), IJK, &
                     & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                     & EP_S(IJK,M), &
                     & PINC(IJK), & 
                     & CUT_CELL_AT(IJK) 
                     
                     if(mype.eq.pe_IO) WRITE(*, 2000) EP_G(IJK), IJK, &
                     & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                     & EP_S(IJK,M), &
                     & PINC(IJK), & 
                     & CUT_CELL_AT(IJK) 
                     
                     IF(cartesian_grid) then 
                        CALL WRITE_DES_DATA
                        CALL WRITE_VTK_FILE
                        !if(dmp_log) write(unit_log, *) 'will not terminate the simulation here in normal compute fields'
                        !if(mype.eq.pe_io) write(*, *) 'will not terminate the simulation here in normal compute fields'
                        if(dmp_log) write(unit_log, *) 'Terminal error, stopping the simulation here'
                        if(mype.eq.pe_io) write(*, *) 'Terminal error, stopping the simulation here'                         
                        
                        call mfix_exit(myPE)
                     ELSE
                        
                        if(dmp_log) write(unit_log, *) 'Terminal error, stopping the simulation here'
                        if(mype.eq.pe_io) write(*, *) 'Terminal error, stopping the simulation here'                         
                        call mfix_exit(myPE)
                     end IF
                  ENDIF
                  ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
               ENDIF
            ENDDO
         ENDIF
                  
      ENDDO                     ! ijk loop 

      EPG_MIN2 = MINVAL(EP_G(:))
      epg_min_loc = MINLOC(EP_G(:))
      IJK = epg_min_loc(1)
      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      WRITE(*,1014) epg_min2, & 
      & I_OF(IJK), j_of(ijk), k_of(ijk), &
      & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
      & PINC(IJK), & 
      & cut_cell_at(ijk),fluid_at(ijk)

! 1014 FORMAT(1x,70('*') , / &
 1014 FORMAT(1x,70('*') , / &
      &      10x,'EPGMIN NORMAL = ', 2x,g17.8, / & 
      &      10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2x,i5),/, &
      &      10x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8),/ & 
      &      10x,'No of paricles in cell = ',I10, & 
      &      10x,'CUT CELL, FLUID AT IJK ?    ', 2(2x, L2)) !,/& 
!      &      1X,70('*')/)
      
 2000 format(/5X, 'Message from normal compute mean fields', & 
      & /,5X,'Warning, EP_G = ', g17.8, 2x, 'LT Zero at IJK',I10, & 
      & /,5X,'I,J,K = ', I10, 2X,I10, 2x, I10, & 
      & /,5X,'EP_S = ', ES15.9, & 
      & /,5X,'No of paricles in cell = ',I10, & 
      & /,5X,'Cut cell ? ', L2,/)

      CALL SET_WALL_BC(IER) 
      !the above routine will apply noslip or free slip BC as per the mfix  convention. 
      !currently, this implies NSW or FSW wall BC's will be re-applied to gas-phase 
      !field as well. This can be changed later on to be more specific to MPPIC case    
      WRITE(*,'(10x,A,4(2x,g17.8))') 'NORM: SOLIDS MASS 1 AND 2 =  ', MASS_SOL1, MASS_SOL2, MASS_SOL3, MASS_SOL4

      IF(.not.cartesian_grid) return 

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK)
         U_so(IJK, :) = U_s(IJK, :)
         V_so(IJK, :) = V_s(IJK, :)
         W_so(IJK, :) = W_s(IJK, :)
         U_S(IJK, :) = ZERO
         V_S(IJK, :) = ZERO
         W_S(IJK, :) = ZERO

         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	 
         
         IF(WALL_U_AT(IJK)) THEN
            U_S(IJK, :) = ZERO 
            !currently only No slip BC is being set on this mean 
            !solid's velocity field. Later this wall part can be 
            !treated separately and U_S set only for scalar cells 
            !where FLUID_AT(IJK) is true. 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IPJK = IP_OF(IJK) 
            IF(FLUID_AT(IPJK)) THEN 
               DO M = 1, MMAX

                  U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
               ENDDO
            ELSE
               U_S(IJK,:) = DES_U_S(IJK, :)
            ENDIF
         ENDIF

         
         
         IF(WALL_V_AT(IJK)) THEN
            V_S(IJK, :) = ZERO 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IJPK = JP_OF(IJK) 
            IF(FLUID_AT(IJPK)) THEN 
               DO M = 1, MMAX
                  V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
               ENDDO
            ELSE
               V_S(IJK,:) = DES_V_S(IJK, :)
            ENDIF
         ENDIF

         
         IF(DIMN.EQ.3) THEN 
            IF(WALL_W_AT(IJK)) THEN
               W_S(IJK, :) = ZERO 
            ELSE
               if(.not.FLUID_AT(IJK)) cycle 
               
               IJKP = KP_OF(IJK) 
               IF(FLUID_AT(IJKP)) THEN 
                  DO M = 1, MMAX
                     W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
                  ENDDO
               ELSE
                  W_S(IJK,:) = DES_W_S(IJK, :)
               ENDIF
            ENDIF
         ENDIF
      ENDDO


   
!!$      DO IJK = IJKSTART3, IJKEND3
!!$         I = I_OF(IJK)
!!$         J = J_OF(IJK)
!!$         K = K_OF(IJK)
!!$         IF(.not.IS_ON_myPE_wobnd(I,J,K)) CYCLE 
!!$         IPJK = IP_OF(IJK) 
!!$         IJPK = JP_OF(IJK) 
!!$         IJKP = KP_OF(IJK) 
!!$         
!!$         DO M = 1, MMAX
!!$            U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
!!$            V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
!!$            if(dimn.eq.3) W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
!!$         ENDDO
!!$      ENDDO

      !RETURN 
      
      WRITE(filename,'(A,"_",I5.5,".dat")') 'EPS_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', status='unknown')
      write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "EPS" ', ' "EPSO" ',' "MEANUS" ', ' "MEANVS" ', ' "FLUID" ', ' "VOL" ' , ' "MEANU_S" ', ' "MEANV_S" '
      write(1000,*)'ZONE F=POINT, I=', (IEND2-ISTART2)+1,  ', J=', JEND2-JSTART2+1, ', K=', KEND2-KSTART2 + 1
      
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               IJK  = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK)) THEN 
                  FLUID_IND = 1
               ELSE 
                  FLUID_IND = 0
               END IF
               
               if(dimn.eq.2) zcor = zt(k)
               if(dimn.eq.3) zcor = zt(k-1) + dz(k)
               M=1
               write(1000,'(7(2x,g17.8),(2x,i4),3( 2x, g17.8))') XE(I-1)+DX(I), YN(J-1)+DY(J),ZCOR, ROP_S(IJK,m)/ro_s(m), ROP_SO(IJK,m)/ro_s(m), des_u_s(ijk,m), des_v_s(ijk,m), fluid_ind, vol(ijk), U_S(IJK,1), V_S(IJK,1)
                  
            enddo
         enddo
      enddo
      close(1000, status='keep')
      !stop
      END SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS


      
      SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS2
      
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


 ! general i, j, k indices
      INTEGER I, J, K, IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
              IPJPK, IPJKP, IJPKP, IPJPKP, II, JJ, KK, &
! Pradeep added following 
              IMJMK,IMJKM,IJMKM,IMJMKM
      INTEGER I1, I2, J1, J2, K1, K2, IDIM, IJK2
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
      DOUBLE PRECISION AVG_FACTOR, TEMP1, AVG_FACTOR_FACE

! index of solid phase that particle NP belongs to      
      INTEGER M

! volume of fluid cell particle resides in
      DOUBLE PRECISION VCELL, VCELL2 
! one over the solids volume fraction and one over the volume       
      DOUBLE PRECISION OEPS, OVOL, MASS_SOL1, MASS_SOL2, MASS_SOL3, MASS_SOL4

! particle number index, used for looping      
      INTEGER NP

! index to track accounted for particles 
      INTEGER PC 

! for error messages      
      INTEGER IER

! Statistical weight of the particle. Equal to one for DEM 
   
      DOUBLE PRECISION WTP, ZCOR, JUNK_VAL(3)

! Pradeep temporary indices for periodic boundary adjustments
      integer korder,cur_ijk,nindx
! Pradeep introducing volume at grid nodes for backward interpolation 
      double precision VELG_ARR(DIMN), VELS_ARR(DIMN, MMAX), np_m

! see the discussion for IJK_U ..... in comments       
      INTEGER  IJK_U, IJK_V, IJK_W, ICUR, JCUR, KCUR

      character*100 :: filename
      integer  FLUID_IND      
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      
      integer :: epg_min_loc(1), count_nodes_outside, count_nodes_inside, count_nodes_inside_max, COUNT_TEMP
      DOUBLE PRECISION :: EPG_MIN2, VOL_SURR
      
      double precision :: RESID_WTBAR(MMAX), RESID_BM(DIMN, MMAX), NORM_FACTOR
!
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'


     
      MASS_SOL1 = zero
      MASS_SOL2 = zero 
      MASS_SOL3 = zero 
      MASS_SOL4 = zero 
      call set_interpolation_scheme(2)
      korder = 1+(dimn-2)
      drag_bm = ZERO
      drag_am = zero 
      wtbar = zero
      if(dimn.eq.2) count_nodes_inside_max = 4
      if(dimn.eq.3) count_nodes_inside_max = 8
      IJKLOOP: DO IJK = IJKSTART3,IJKEND3
         
         ROP_S(IJK,:) = zero 
         DES_U_S(IJK, :) = ZERO
         DES_V_S(IJK, :) = ZERO
         if(dimn.eq.3) DES_W_S(IJK, :) = ZERO
         if(.not.fluid_at(ijk) .or. pinc(ijk).eq.0) cycle 
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         pcell(1) = i-1
         pcell(2) = j-1
         pcell(3) = (3-dimn)*1+(dimn-2)*(k-1) ! =k-1 (in 3d) or =1 (in 2d)
         call set_interpolation_stencil(pcell,iw,ie,js,jn,kb,&
         ktp,interp_scheme,dimn,ordernew = onew) 

!Compute velocity at grid nodes and set the geometric stencil 
         avg_factor = 0.25d0*(dimn-2) + 0.5d0*(3-dimn)
         count_nodes_outside = 0 
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
                  if(dimn.eq.3) then 
                     ijkp    = funijk(imap_c(ii),jmap_c(jj),kmap_c(kk+1))
                     ijpkp   = funijk(imap_c(ii),jmap_c(jj+1),kmap_c(kk+1))
                     ipjkp   = funijk(imap_c(ii+1),jmap_c(jj),kmap_c(kk+1))
                     ipjpkp  = funijk(imap_c(ii+1),jmap_c(jj+1),kmap_c(kk+1))
                  endif
                  gstencil(i,j,k,1) = xe(ii)
                  gstencil(i,j,k,2) = yn(jj)
                  gstencil(i,j,k,3) = zt(kk)*(dimn-2) + dz(1)*(3-dimn)
                  vstencil(i,j,k,:) = zero 
                  if(cartesian_grid) then 
                     if(scalar_node_atwall(cur_ijk)) count_nodes_outside = &
                     & count_nodes_outside + 1 
                  endif
               enddo
            enddo
         enddo

         count_nodes_inside = count_nodes_inside_max - count_nodes_outside

         !loop through particles in the cell  
         do nindx = 1,pinc(ijk)
            focus = .false.
            np = pic(ijk)%p(nindx)
            m = pijk(np,5)
            
            NP_M = ROP_SO(IJK,m)*VOL(IJK)/(ro_s(m)*pvol(np))

            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            
            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP
            if (dimn .eq. 2) then 
               call interpolator(gstencil(1:onew,1:onew,1,1:dimn), &
                    vstencil(1:onew,1:onew,1,1:dimn), &
                    des_pos_new(np,1:dimn),JUNK_VAL(1:dimn),  &
                    onew,interp_scheme,weightp)
            else 
               call interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn), &
                    vstencil(1:onew,1:onew,1:onew,1:dimn), &
                    des_pos_new(np,1:dimn),JUNK_VAL(1:dimn),  &
                    onew,interp_scheme,weightp)
            end if
            
            do k = 1, (3-dimn)*1+(dimn-2)*onew
               do j = 1, onew
                  do i = 1, onew
                     ! shift loop index to new variables for manipulation                     
                     ii = iw + i-1
                     jj = js + j-1
                     kk = kb + k-1

                     icur = imap_c(ii)
                     jcur = jmap_c(jj)
                     kcur = kmap_c(kk)
                     
                     cur_ijk = funijk(icur, jcur, kcur) !imap_c(ii),jmap_c(jj),kmap_c(kk))
                     
! should this volume be for current ijk index or always particle index?                        
                     vcell = des_vol_node(cur_ijk)
                     ovol = one/vcell
                     
                     temp1 = weightp(i,j,k)*ro_s(m)*pvol(np)*wtp!*ovol
                                         
                     wtbar(cur_ijk,m) = wtbar(cur_ijk,m) + temp1 
                     
                     drag_bm(cur_ijk, 1:dimn,m) = &
                     drag_bm(cur_ijk,1:dimn,m) + temp1*des_vel_new(np, 1:dimn)
                     
                     MASS_SOL3 = MASS_SOL3 +  TEMP1

                  enddo
               enddo
            enddo
            
         enddo          ! pinc(ijk) loop 
         
         
         if(count_nodes_inside.lt.count_nodes_inside_max) then 
            i = i_of(ijk)
            j = j_of(ijk)
            k = k_of(ijk)
            
            I1 = I-1
            I2 = I
            J1 = J-1
            J2 = J
            
            IF(NO_K) THEN 
               K1 = K
               K2 = K
            ELSE
               K1 = K-1
               K2 = K
            ENDIF
         !Convention used to number node numbers is described below 
         
         ! i=1, j=2           i=2, j=2
         !   _____________________
         !   |                   |
         !   |  I = 2, J = 2     |
         !   |___________________|
         ! i=1, j=1           i=2, j=1
         !first calculate the residual wtbar and drag_bm that was computed
         !on nodes that do not belong to the domain
            RESID_WTBAR(1:MMAX) = ZERO 
            RESID_BM(1:DIMN, 1:MMAX) = ZERO
            DO KK = K1, K2
               DO JJ = J1, J2
                  DO II = I1, I2
                     IJK2 = funijk(II, JJ, KK) 
                     IF(SCALAR_NODE_ATWALL(IJK2)) THEN 
                        RESID_WTBAR(1:MMAX) = RESID_WTBAR(1:MMAX) + WTBAR(IJK2,1:MMAX)
                        WTBAR(IJK2,1:MMAX) = ZERO 
                        DO IDIM = 1, DIMN
                           RESID_BM(IDIM, 1:MMAX) = RESID_BM(IDIM, 1:MMAX) +  & 
                           & DRAG_BM(IJK2,IDIM, 1:MMAX)
                           DRAG_BM(IJK2,IDIM, 1:MMAX) = ZERO 
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
                  
                     
           ! Now add this residual equally to the remaining nodes
           ! Pradeep: at the interface drag_am,drag_bm,wtbar has to be added
           ! send recv will be called and the node values will be added 
           ! at the junction 
            
            NORM_FACTOR = one/real(count_nodes_inside)          
            COUNT_TEMP = 0
            DO KK = K1, K2
               DO JJ = J1, J2
                  DO II = I1, I2
                     IJK2 = funijk(II, JJ, KK) 
                     IF(.NOT.SCALAR_NODE_ATWALL(IJK2)) THEN 
                        COUNT_TEMP = COUNT_TEMP + 1
                        WTBAR(IJK2,1:MMAX) = WTBAR(IJK2,1:MMAX) + RESID_WTBAR(1:MMAX)*NORM_FACTOR
                        DO IDIM = 1, DIMN
                           DRAG_BM(IJK2,IDIM, 1:MMAX) = DRAG_BM(IJK2,IDIM, 1:MMAX) + RESID_BM(IDIM, 1:MMAX)*NORM_FACTOR
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            
            !WRITE(*,*) 'hello: NODES INSIDE AND OUTSIDE', count_nodes_inside, count_nodes_inside_max, one/real(count_temp), norm_factor
         ENDIF
         

      END DO IJKLOOP            ! IJK LOOP
      
      MASS_SOL4 = SUM(WTBAR(:,1))

      CALL DES_ADDNODEVALUES
      
      
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IJK = funijk(I,J,K)

               IF(SCALAR_NODE_ATWALL(IJK)) CYCLE

               !Now going from node to scalar center. Same convention
               !as sketched earlier 
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               IF(NO_K) THEN 
                  K1 = K
                  K2 = K
               ELSE
                  K1 = K
                  K2 = K+1
               ENDIF
               
               VOL_SURR = ZERO 

               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK) 
                        IF(FLUID_AT(IJK2)) VOL_SURR = VOL_SURR+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO
                 
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK) 
                        IF(FLUID_AT(IJK2)) THEN 
                           DO M = 1, MMAX
                              ROP_S(IJK2, M) = ROP_S(IJK2, M) + WTBAR(IJK,M)*VOL(IJK2)/VOL_SURR
                              DES_U_S(IJK2, M) = DES_U_S(IJK2, M) + DRAG_BM(IJK, 1, M)*VOL(IJK2)/VOL_SURR
                              DES_V_S(IJK2, M) = DES_V_S(IJK2, M) + DRAG_BM(IJK, 2, M)*VOL(IJK2)/VOL_SURR
                              IF(DIMN.eq.3) DES_W_S(IJK2, M) = DES_W_S(IJK2, M) + DRAG_BM(IJK, 3, M)*VOL(IJK2)/VOL_SURR
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      

      DO IJK = IJKSTART3, IJKEND3
         IF(.not.FLUID_AT(IJK)) cycle
         DO M = 1, MMAX
            IF(ROP_S(IJK, M).GT.ZERO) THEN 
               DES_U_S(IJK, M) = DES_U_S(IJK,M)/ROP_S(IJK, M)
               DES_V_S(IJK, M) = DES_V_S(IJK,M)/ROP_S(IJK, M)
               IF(DIMN.eq.3) DES_W_S(IJK, M) = DES_W_S(IJK,M)/ROP_S(IJK, M)
               
               ROP_S(IJK, M) = ROP_S(IJK, M)/VOL(IJK)
            ENDIF
         ENDDO
                  
         EP_G(IJK) = ONE   
            
            !MASS_SOL3 = MASS_SOL3 + ROP_SO(IJK, 1)*VOL(IJK)
         DO M = 1, MMAX
            
            IF(ROP_S(IJK,M) > ZERO) THEN
               
               MASS_SOL2 = MASS_SOL2 + ROP_S(IJK, M)*VOL(IJK)
               
                  
               IF(.not.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
               
               ROP_SO(IJK,M)  = ROP_S(IJK,M) 
               
               IF(EP_G(IJK).LT.ZERO .AND. DES_CONTINUUM_COUPLED) THEN 
                                ! this does not matter if pure granular flow simulation (i.e. no fluid)
                  IF(DMP_LOG)  WRITE(UNIT_LOG, 2000) EP_G(IJK), IJK, &
                  & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                  & EP_S(IJK,M), &
                  & PINC(IJK), & 
                  & CUT_CELL_AT(IJK) 
                     
                  if(mype.eq.pe_IO) WRITE(*, 2000) EP_G(IJK), IJK, &
                  & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                  & EP_S(IJK,M), &
                  & PINC(IJK), & 
                  & CUT_CELL_AT(IJK) 
                  
                  IF(cartesian_grid) then 
                     CALL WRITE_DES_DATA
                     CALL WRITE_VTK_FILE
                        !if(dmp_log) write(unit_log, *) 'will not terminate the simulation here in normal compute fields'
                        !if(mype.eq.pe_io) write(*, *) 'will not terminate the simulation here in normal compute fields'
                     if(dmp_log) write(unit_log, *) 'Terminal error, stopping the simulation here'
                     if(mype.eq.pe_io) write(*, *) 'Terminal error, stopping the simulation here'                         
                     
                     call mfix_exit(myPE)
                  ELSE
                        
                     if(dmp_log) write(unit_log, *) 'Terminal error, stopping the simulation here'
                     if(mype.eq.pe_io) write(*, *) 'Terminal error, stopping the simulation here'                         
                     call mfix_exit(myPE)
                  end IF
               ENDIF
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
            ENDIF
         ENDDO
         
      ENDDO                     ! ijk loop 
         
      EPG_MIN2 = MINVAL(EP_G(:))
      epg_min_loc = MINLOC(EP_G(:))
      IJK = epg_min_loc(1)
      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      WRITE(*,1014) epg_min2, & 
      & I_OF(IJK), j_of(ijk), k_of(ijk), &
      & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
      & PINC(IJK), & 
      & cut_cell_at(ijk),fluid_at(ijk)

! 1014 FORMAT(1x,70('*') , / &
 1014 FORMAT( /, &
      &      10x,'EPGMIN NORMAL = ', 2x,g17.8, / & 
      &      10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2x,i5),/, &
      &      10x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8),/ & 
      &      10x,'No of paricles in cell = ',I10, & 
      &      10x,'CUT CELL, FLUID AT IJK ?    ', 2(2x, L2)) !,/& 
!      &      1X,70('*')/)
      
 2000 format(/5X, 'Message from normal compute mean fields', & 
      & /,5X,'Warning, EP_G = ', g17.8, 2x, 'LT Zero at IJK',I10, & 
      & /,5X,'I,J,K = ', I10, 2X,I10, 2x, I10, & 
      & /,5X,'EP_S = ', ES15.9, & 
      & /,5X,'No of paricles in cell = ',I10, & 
      & /,5X,'Cut cell ? ', L2,/)

      CALL SET_WALL_BC(IER) 
      !the above routine will apply noslip or free slip BC as per the mfix  convention. 
      !currently, this implies NSW or FSW wall BC's will be re-applied to gas-phase 
      !field as well. This can be changed later on to be more specific to MPPIC case    
      WRITE(*,'(10x,A,4(2x,g17.8))') 'NORM: SOLIDS MASS 1 AND 2 =  ', MASS_SOL1, MASS_SOL2, MASS_SOL3, MASS_SOL4

      IF(.not.cartesian_grid) return 

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK)
         U_so(IJK, :) = U_s(IJK, :)
         V_so(IJK, :) = V_s(IJK, :)
         W_so(IJK, :) = W_s(IJK, :)
         U_S(IJK, :) = ZERO
         V_S(IJK, :) = ZERO
         W_S(IJK, :) = ZERO

         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	 
         
         IF(WALL_U_AT(IJK)) THEN
            U_S(IJK, :) = ZERO 
            !currently only No slip BC is being set on this mean 
            !solid's velocity field. Later this wall part can be 
            !treated separately and U_S set only for scalar cells 
            !where FLUID_AT(IJK) is true. 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IPJK = IP_OF(IJK) 
            IF(FLUID_AT(IPJK)) THEN 
               DO M = 1, MMAX

                  U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
               ENDDO
            ELSE
               U_S(IJK,:) = DES_U_S(IJK, :)
            ENDIF
         ENDIF

         
         
         IF(WALL_V_AT(IJK)) THEN
            V_S(IJK, :) = ZERO 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IJPK = JP_OF(IJK) 
            IF(FLUID_AT(IJPK)) THEN 
               DO M = 1, MMAX
                  V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
               ENDDO
            ELSE
               V_S(IJK,:) = DES_V_S(IJK, :)
            ENDIF
         ENDIF

         
         IF(DIMN.EQ.3) THEN 
            IF(WALL_W_AT(IJK)) THEN
               W_S(IJK, :) = ZERO 
            ELSE
               if(.not.FLUID_AT(IJK)) cycle 
               
               IJKP = KP_OF(IJK) 
               IF(FLUID_AT(IJKP)) THEN 
                  DO M = 1, MMAX
                     W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
                  ENDDO
               ELSE
                  W_S(IJK,:) = DES_W_S(IJK, :)
               ENDIF
            ENDIF
         ENDIF
      ENDDO


   
!!$      DO IJK = IJKSTART3, IJKEND3
!!$         I = I_OF(IJK)
!!$         J = J_OF(IJK)
!!$         K = K_OF(IJK)
!!$         IF(.not.IS_ON_myPE_wobnd(I,J,K)) CYCLE 
!!$         IPJK = IP_OF(IJK) 
!!$         IJPK = JP_OF(IJK) 
!!$         IJKP = KP_OF(IJK) 
!!$         
!!$         DO M = 1, MMAX
!!$            U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
!!$            V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
!!$            if(dimn.eq.3) W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
!!$         ENDDO
!!$      ENDDO

      !RETURN 
      
      WRITE(filename,'(A,"_",I5.5,".dat")') 'EPS_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', status='unknown')
      write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "EPS" ', ' "EPSO" ',' "MEANUS" ', ' "MEANVS" ', ' "FLUID" ', ' "VOL" ' , ' "MEANU_S" ', ' "MEANV_S" '
      write(1000,*)'ZONE F=POINT, I=', (IEND2-ISTART2)+1,  ', J=', JEND2-JSTART2+1, ', K=', KEND2-KSTART2 + 1
      
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               IJK  = FUNIJK(I,J,K)
               IF(FLUID_AT(IJK)) THEN 
                  FLUID_IND = 1
               ELSE 
                  FLUID_IND = 0
               END IF
               
               if(dimn.eq.2) zcor = zt(k)
               if(dimn.eq.3) zcor = zt(k-1) + dz(k)
               M=1
               write(1000,'(7(2x,g17.8),(2x,i4),3( 2x, g17.8))') XE(I-1)+DX(I), YN(J-1)+DY(J),ZCOR, ROP_S(IJK,m)/ro_s(m), ROP_SO(IJK,m)/ro_s(m), des_u_s(ijk,m), des_v_s(ijk,m), fluid_ind, vol(ijk), U_S(IJK,1), V_S(IJK,1)
                  
            enddo
         enddo
      enddo
      close(1000, status='keep')
      !stop
      END SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS2



      SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS_CG
      
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE vtk
      
      USE fldvar


      USE physprop
      USE bc
      USE discretelement
      USE drag
      USE interpolation
      use desmpi 

      IMPLICIT NONE

      INTEGER :: I, J, K,IJK,NINDX,NP,M,NODE,NC, IER, IPJK, IJPK, IJKP

      DOUBLE PRECISION :: WTP

      DOUBLE PRECISION :: MASS_SOL1,MASS_SOL2, ZCOR 

      DOUBLE PRECISION :: XND,YND,ZND,XP,YP,ZP,SUM_O_DIST,DIST_TO_NODE

      DOUBLE PRECISION :: SUM_OF_WEIGHT,SUM_CONTRIBUTIONS, SUM_CONTRIB_VEL(DIMN)

      DOUBLE PRECISION, DIMENSION(15) :: O_DIST_TO_NODE
      DOUBLE PRECISION :: WEIGHT_ON_NODE, TEMP1, OVOL
      DOUBLE PRECISION, DIMENSION(DIMENSION_3+DIMENSION_MAX_CUT_CELL,DIM_M) :: CONTRIBUTION_OF_NODE
      
      DOUBLE PRECISION, DIMENSION(DIMENSION_3+DIMENSION_MAX_CUT_CELL,DIMN, DIM_M) :: CONTRIBUTION_OF_NODE_SOLVEL
      INTEGER :: CUTCELL_IND
      LOGICAL ::FOCUS
      
      character*100 :: filename
      integer  FLUID_IND      
      double precision :: epg_min2
      integer :: epg_min_loc(1)
!
!-----------------------------------------------   

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'


      CONTRIBUTION_OF_NODE = ZERO 
      CONTRIBUTION_OF_NODE_SOLVEL = ZERO 
      MASS_SOL1 = ZERO 
      DO IJK = IJKSTART3,IJKEND3
         
         IF(.NOT.FLUID_AT(IJK) .OR. PINC(IJK).EQ.0) CYCLE 
         
         !print*,'==========================================================='
         !print*,'IJK,  I,J=',IJK,I_OF(IJK),J_OF(IJK)


!         print*,'NUMBER_OF_NODES=',NUMBER_OF_NODES(IJK)
!         DO NODE = 1,NUMBER_OF_NODES(IJK)
!            print*,'CNCT=',NODE,CONNECTIVITY(IJK,NODE),SCALAR_NODE_XYZ(CONNECTIVITY(IJK,NODE),1),SCALAR_NODE_XYZ(CONNECTIVITY(IJK,NODE),2)
!         ENDDO
!         print*,''

         !print*,'NUMBER OF PARTICLES IN CELL =',PINC(IJK)

        !LOOP THROUGH PARTICLES IN THE CELL  
         DO NINDX = 1,PINC(IJK)
            FOCUS = .FALSE.
            NP = PIC(IJK)%P(NINDX)
            M = PIJK(NP,5)
            
            
            XP = des_pos_new(np,1)
            YP = des_pos_new(np,2)
            
            IF(NO_K) THEN
               ZP = des_pos_new(np,3)
            ELSE
               ZP = ZERO
            ENDIF

            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            
            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP

            SUM_O_DIST = ZERO

            !print*,'NUMBER_OF_NODES=',NUMBER_OF_NODES(IJK)
            
            DO NODE = 1,NUMBER_OF_NODES(IJK) ! Distribute the contribution of each particle on the cell nodes
                                                                ! Based on inverse distance from the particle and each node
                                                                ! First loop to prepare distribution by computing particle-node distance

               NC = CONNECTIVITY(IJK,NODE)
               
               XND = SCALAR_NODE_XYZ(NC,1)
               YND = SCALAR_NODE_XYZ(NC,2)

               IF(NO_K) THEN
                  ZND = SCALAR_NODE_XYZ(NC,3)
               ELSE
                  ZND = ZERO
               ENDIF
               
               DIST_TO_NODE = DSQRT((XND-XP)**2 + (YND-YP)**2 + (ZND-ZP)**2)

!               IF(DIST_TO_NODE < HALF*D_P0(M)) THEN
               IF(DIST_TO_NODE < 1.0D-12) THEN ! Hard-wired tolerance, maybe need to make it user-defined
                  O_DIST_TO_NODE(NODE) = UNDEFINED
               ELSE
                  O_DIST_TO_NODE(NODE) = ONE / DIST_TO_NODE
               ENDIF    

               SUM_O_DIST = SUM_O_DIST + O_DIST_TO_NODE(NODE)

!            print*,'NODE,DIST=',NODE,DIST_TO_NODE

            ENDDO

            SUM_OF_WEIGHT = ZERO
            
            DO NODE = 1,NUMBER_OF_NODES(IJK)                  ! Second loop to distribute particle contribution onto each node

               NC = CONNECTIVITY(IJK,NODE)
               
               WEIGHT_ON_NODE = O_DIST_TO_NODE(NODE)/SUM_O_DIST
               
               SUM_OF_WEIGHT = SUM_OF_WEIGHT + WEIGHT_ON_NODE

               !print*,'NODE,WEIGHT=',NODE,WEIGHT_ON_NODE

               !WTBAR(NC,M) = WTBAR(NC,M) + WEIGHT_OF_NODE(NC)*RO_S(M)*OVOL*PVOL(NP)*WTP
               OVOL = MAX(ZERO, Ovol_around_node(NC)) 

               TEMP1 =  WEIGHT_ON_NODE * RO_S(M)*PVOL(NP)*WTP*OVOL 
               CONTRIBUTION_OF_NODE(NC,M) = CONTRIBUTION_OF_NODE(NC,M) + TEMP1 
               
               CONTRIBUTION_OF_NODE_SOLVEL(NC, 1:DIMN, M) = CONTRIBUTION_OF_NODE_SOLVEL(NC, 1:DIMN, M) + TEMP1*DES_VEL_NEW(NP, 1:DIMN)
            ENDDO

            !print*,'SUM OF WEIGHTH (should be equal to ONE) =',SUM_OF_WEIGHT      
            !read(*,*)

            
         ENDDO                  ! PINC(IJK) LOOP 
      END DO                    ! IJK LOOP



!!$      DO IJK = IJKSTART3, IJKEND3 + NUMBER_OF_NEW_POINTS ! Loop over all nodes 
!!$         IF(Ovol_around_node(IJK)>ZERO) THEN
!!$            CONTRIBUTION_OF_NODE(IJK,M) = CONTRIBUTION_OF_NODE(IJK,M) * Ovol_around_node(IJK) ! Divide by volume surrounding a given node
!!$         ENDIF
!!$      ENDDO

      MASS_SOL2 = ZERO

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK)
         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	 
         IF (.NOT.FLUID_AT(IJK)) CYCLE 
         ROP_SO(IJK,:) = ROP_S(IJK,:)
         ROP_S(IJK,:) = ZERO 
!         IF(.NOT.FLUID_AT(IJK) .OR. PINC(IJK).EQ.0) CYCLE 
         IF(.NOT.FLUID_AT(IJK)) CYCLE 
         EP_G(IJK) = ONE   

         DO M = 1, MMAX
            SUM_CONTRIBUTIONS = ZERO
            SUM_CONTRIB_VEL(1:DIMN) = ZERO 
            DO NODE = 1,NUMBER_OF_NODES(IJK)
               NC = CONNECTIVITY(IJK,NODE)
               SUM_CONTRIBUTIONS = SUM_CONTRIBUTIONS + CONTRIBUTION_OF_NODE(NC,M)
               SUM_CONTRIB_VEL(1:DIMN) = SUM_CONTRIB_VEL(1:DIMN) +  CONTRIBUTION_OF_NODE_SOLVEL(NC, 1:DIMN, M)
            ENDDO
            ROP_S(IJK,M) = SUM_CONTRIBUTIONS / NUMBER_OF_NODES(IJK) ! ROP_S in a cellis the average of the ROP_S at the nodes
            
            !below, it is not <Us>, but EP_S*RO_S*<U_s>
            DES_U_S(IJK, M) = SUM_CONTRIB_VEL(1)/NUMBER_OF_NODES(IJK)
            DES_V_S(IJK, M) = SUM_CONTRIB_VEL(2)/NUMBER_OF_NODES(IJK)
            IF(DIMN.eq.3) DES_W_S(IJK, M) = SUM_CONTRIB_VEL(3)/NUMBER_OF_NODES(IJK)
            
         ENDDO

!Rahul:
!At this point, a routine will need to be called for MPI runs that adds the 
!contibutions along the processor boundaries. This will be similar to the 
!routine des_addnodevalues, except that the summation needs to be done at the 
!scalar cell centers, unlike the scalar cell nodes in des_addnodevalues

         !convert EP_S*RO_S*<U_s> to <U_s> by dividing by ROP_S
         DO M = 1, MMAX
            IF(ROP_S(IJK,M).GT.ZERO) then 
               DES_U_S(IJK,M) = DES_U_S(IJK,M)/ROP_S(IJK,M)
               DES_V_S(IJK,M) = DES_V_S(IJK,M)/ROP_S(IJK,M)
               if(DIMN.eq.3) DES_W_S(IJK,M) = DES_W_S(IJK,M)/ROP_S(IJK,M)
            ELSE
               DES_U_S(IJK,M) = ZERO 
               DES_V_S(IJK,M) = ZERO 
               if(DIMN.eq.3) DES_W_S(IJK,M) = ZERO 
               
            ENDIF
         ENDDO
         
         
         DO M = 1, MMAX
            
            MASS_SOL2 = MASS_SOL2 + ROP_S(IJK,M) * VOL(IJK)
            
            IF(.not.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
                  
                  
            IF(EP_G(IJK).LT.ZERO .AND. DES_CONTINUUM_COUPLED) THEN 
               
               CALL WRITE_DES_DATA
               CALL WRITE_VTK_FILE
            ! this does not matter if pure granular flow simulation (i.e. no fluid)
               IF(DMP_LOG)  WRITE(UNIT_LOG, 2000) EP_G(IJK), IJK, &
               & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
               & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
               & EP_S(IJK,M), &
               & PINC(IJK), & 
               & CUT_CELL_AT(IJK) , &
               & VOL(IJK), DX(I_OF(IJK))*DY(J_OF(IJK))*DZ(K_OF(IJK)), (VOL(IJK)/DX(I_OF(IJK))/(DY(J_OF(IJK))*DZ(K_OF(IJK))))*100.
                     
               if(mype.eq.pe_IO) WRITE(*, 2000) EP_G(IJK), IJK, &
               & I_OF(IJK), J_OF(IJK), K_OF(IJK), &
               & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
               & EP_S(IJK,M), &
               & PINC(IJK), & 
               & CUT_CELL_AT(IJK), &
               & VOL(IJK), DX(I_OF(IJK))*DY(J_OF(IJK))*DZ(K_OF(IJK)), (VOL(IJK)/(DX(I_OF(IJK))*DY(J_OF(IJK))*DZ(K_OF(IJK))))*100.
               
                        
               if(dmp_log) write(unit_log, *) 'Terminal error, stopping the simulation here'
               if(mype.eq.pe_io) write(*, *) 'Terminal error, stopping the simulation here'                         
               call mfix_exit(myPE)
            
            ENDIF
            
            ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

            !print*,'IJK,EPS(IJK,M),EPG=',IJK,EP_S(IJK,M),ONE-EP_S(IJK,M)

         ENDDO
      ENDDO

 2000 format(/10X, 'Message from CG compute mean fields case', & 
      & /,10X,'Warning, EP_G    = ', g17.8, 2x, 'LT Zero at IJK',I10, & 
      & /,10X,'I,J,K            = ', I10, 2X,I10, 2x, I10, & 
      & /,10x,'XMID, YMID, ZMID = ', 3(2x,g17.8), & 
      & /,10X,'EP_S             = ', ES15.9, & 
      & /,10X,'No of paricles in cell = ', I10, & 
      & /,10X,'Cut cell ? ', L2 , & 
      & /,10X,'Cell, Uncut cell volume, and % ratio = ', 3(2x,g17.8), /)


      EPG_MIN2 = MINVAL(EP_G(:))
      epg_min_loc = MINLOC(EP_G(:))
      IJK = epg_min_loc(1)
      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      WRITE(*,1014) epg_min2, & 
      & I_OF(IJK), j_of(ijk), k_of(ijk), &
      & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
      & cut_cell_at(ijk),fluid_at(ijk)

 1014 FORMAT(1x,70('*') , / &
      &      10x,'EPGMIN NORMAL = ', 2x,g17.8, / & 
      &      10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2x,i5),/, &
      &      10x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8),/ & 
      &      10x,'CUT CELL, FLUID AT IJK ?    ', 2(2x, L2),/& 
      &      1X,70('*')/)      
      call send_recv(ep_g,2)
      call send_recv(rop_g,2)
      call send_recv(des_u_s,2)
      call send_recv(des_v_s,2) 
      if(dimn.eq.3) call send_recv(des_w_s,2) 
      call send_recv(rop_s,2)

      
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK)
         U_so(IJK, :) = U_s(IJK, :)
         V_so(IJK, :) = V_s(IJK, :)
         W_so(IJK, :) = W_s(IJK, :)
         
         IF (.NOT.IS_ON_myPE_plus1layer(I,J,K)) CYCLE	 
         
         
         IF(WALL_U_AT(IJK)) THEN
            U_S(IJK, :) = ZERO 
            !currently only No slip BC is being set on this mean 
            !solid's velocity field. Later this wall part can be 
            !treated separately and U_S set only for scalar cells 
            !where FLUID_AT(IJK) is true. 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IPJK = IP_OF(IJK) 
            IF(FLUID_AT(IPJK)) THEN 
               DO M = 1, MMAX

                  U_S(IJK,M) = 0.5d0*(DES_U_S(IJK, M) + DES_U_S(IPJK,M))
               ENDDO
            ELSE
               U_S(IJK,:) = DES_U_S(IJK, :)
            ENDIF
         ENDIF

         
         
         IF(WALL_V_AT(IJK)) THEN
            V_S(IJK, :) = ZERO 
         ELSE
            if(.not.FLUID_AT(IJK)) cycle 

            IJPK = JP_OF(IJK) 
            IF(FLUID_AT(IJPK)) THEN 
               DO M = 1, MMAX
                  V_S(IJK,M) = 0.5d0*(DES_V_S(IJK, M) + DES_V_S(IJPK,M))
               ENDDO
            ELSE
               V_S(IJK,:) = DES_V_S(IJK, :)
            ENDIF
         ENDIF

         
         IF(DIMN.EQ.3) THEN 
            IF(WALL_W_AT(IJK)) THEN
               W_S(IJK, :) = ZERO 
            ELSE
               if(.not.FLUID_AT(IJK)) cycle 
               
               IJKP = KP_OF(IJK) 
               IF(FLUID_AT(IJKP)) THEN 
                  DO M = 1, MMAX
                     W_S(IJK,M) = 0.5d0*(DES_W_S(IJK, M) + DES_W_S(IJKP,M))
                  ENDDO
               ELSE
                  W_S(IJK,:) = DES_W_S(IJK, :)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      !print*,'Sum of all PINC:'
      !print*,PINC(9)+PINC(14)+PINC(19)+PINC(8)+PINC(13)+PINC(18)+PINC(7)+PINC(12)
      !print*,'=============================================================='
      !print*,'PINC in each cell'
      !print*,PINC(9),PINC(14),PINC(19)
      !print*,PINC(8),PINC(13),PINC(18)
      !print*,PINC(7),PINC(12)
      !print*,'=============================================================='
      !print*,'EP_G with new method'
      !print*,ONE-EP_S(9,1),ONE-EP_S(14,1),ONE-EP_S(19,1)
      !print*,ONE-EP_S(8,1),ONE-EP_S(13,1),ONE-EP_S(18,1)
      !print*,ONE-EP_S(7,1),ONE-EP_S(12,1)
      !print*,'=============================================================='
      !print*,'EP_G with old method'
      !print*,EP_G(9),EP_G(14),EP_G(19)
      !print*,EP_G(8),EP_G(13),EP_G(18)
      !print*,EP_G(7),EP_G(12)
      !print*,'=============================================================='

    
      WRITE(*,'(10x,A30,4(2x,g17.8))') 'CG: SOLIDS MASS 1 AND 2 =  ', MASS_SOL1, MASS_SOL2
      
      !RETURN 
      
      CALL SET_WALL_BC(IER) 
      !the above routine will apply noslip or free slip BC as per the mfix  convention. 
      !currently, this implies NSW or FSW wall BC's will be re-applied to gas-phase 
      !field as well. This can be changed later on to be more specific to MPPIC case       
      
      WRITE(filename,'(A,"_",I5.5,".dat")') 'CG_EPS_',myPE
      OPEN(1000, file = TRIM(filename), form ='formatted', status='unknown')
      write(1000,*)'VARIABLES= ',' "X" ',' "Y" ',' "Z" ',' "EPS" ', ' "EPSO" ',' "MEANUS" ', ' "MEANVS" ', ' "FLUID" ', ' "CUTCELL" ' , ' "VOL" ' , ' "MEANU_S" ', ' "MEANV_S" ', ' "MEANU_SO" ' , ' "MEANV_SO" '
      write(1000,*)'ZONE F=POINT, I=', (IEND2-ISTART2)+1,  ', J=', JEND2-JSTART2+1, ', K=', KEND2-KSTART2 + 1
      
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               IJK  = FUNIJK(I,J,K)
               FLUID_IND = 0 
               CUTCELL_IND = 0 
               IF(FLUID_AT(IJK)) FLUID_IND = 1
               IF(CUT_CELL_AT(IJK)) CUTCELL_IND = 1
               
               if(dimn.eq.2) zcor = zt(k)
               if(dimn.eq.3) zcor = zt(k-1) + dz(k)
               M=1
               write(1000,'(7(2x,g17.8),2((2x,i4)),3( 2x, g17.8))') XE(I-1)+DX(I), YN(J-1)+DY(J),ZCOR, ROP_S(IJK,m)/ro_s(m), ROP_SO(IJK,m)/ro_s(m), des_u_s(ijk,m), des_v_s(ijk,m), fluid_ind, cutcell_ind, vol(ijk), U_S(IJK,1), V_S(IJK,1), U_so(IJK,1), V_so(IJK, 1) 
                  
            enddo
         enddo
      enddo
      close(1000, status='keep')
      !stop
      END SUBROUTINE MPPIC_COMPUTE_MEAN_FIELDS_CG


    
      SUBROUTINE MPPIC_BC_U_S  
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
      USE output
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IPJK, IP, IJK_WALL
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
   	       IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
               IJK_WALL = FUNIJK(I1,J1,K1) 
               IJK = FUNIJK(I1,J1,K1+1) 
               U_S(IJK_WALL, :) = -U_S(IJK,:)
            END DO 
         END DO 
	 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
               
   	       IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
               IJK_WALL = FUNIJK(I1,J1,K1) 
               IJK = FUNIJK(I1,J1,K1-1) 
               U_S(IJK_WALL, :) = -U_S(IJK,:)
            END DO 
         END DO 
      ENDIF 

      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
            
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
            IJK_WALL = FUNIJK(I1,J1,K1) 
            IJK = FUNIJK(I1,J1+1,K1) 
            U_S(IJK_WALL, :) = -U_S(IJK,:)
         END DO 
      END DO 
      J1 = JMAX2 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	    
            IJK_WALL = FUNIJK(I1,J1,K1) 
            IJK = FUNIJK(I1,J1-1,K1) 
            U_S(IJK_WALL, :) = -U_S(IJK,:)
         END DO 
      END DO 
       
      RETURN  
      END SUBROUTINE MPPIC_BC_U_S  
      
      SUBROUTINE MPPIC_BC_V_S  
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
      USE output
      USE compar  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, JM, I1, I2, J1, J2, K1, K2, IJK,& 
                       IM, KM, IJKS, IJMK, IJPK, IJK_WALL 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE		 
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1+1)
               V_S(IJK_WALL, :) = -V_S(IJK,:)
            END DO 
         END DO 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE		 
               IJK_WALL = FUNIJK(I1,J1,K1)
               IJK = FUNIJK(I1,J1,K1-1)
               V_S(IJK_WALL, :) = -V_S(IJK,:)
            END DO 
         END DO 
      ENDIF 
!
      I1 = 1
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE		 
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1+1,J1,K1)
            V_S(IJK_WALL, :) = -V_S(IJK,:)
            
         END DO 
      END DO 
      I1 = IMAX2 
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
             
            IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE		 
            IJK_WALL = FUNIJK(I1,J1,K1)
            IJK = FUNIJK(I1-1,J1,K1)
            V_S(IJK_WALL, :) = -V_S(IJK,:)
            
         END DO 
      END DO
      END SUBROUTINE MPPIC_BC_V_S  

      

      SUBROUTINE MPPIC_BC_W_S  
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
      USE output
      USE compar  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, KM, I1, I2, J1, J2, K1, K2, IJK,& 
                       IM, JM, IJKB, IJKM, IJKP, IJK_WALL 
! 
      INCLUDE 'function.inc'
!
!  Set the default boundary conditions
!
      J1 = 1 
      DO K1 = kmin3,kmax3 
         DO I1 = imin3,imax3 
   	    IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	 	    	    
            
            IJK = FUNIJK(I1,J1+1,K1) 
            IJK_WALL = FUNIJK(I1,J1,K1) 
            W_S(IJK_WALL,:) = -W_S(IJK,:)
         END DO 
      END DO 
      J1 = JMAX2 
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3 
            
   	    IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	 	    	    
            
            IJK = FUNIJK(I1,J1-1,K1) 
            IJK_WALL = FUNIJK(I1,J1,K1) 
            W_S(IJK_WALL,:) = -W_S(IJK,:)
         END DO 
      END DO 
      I1 = 1 
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3 
            
   	    IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	 	    	    
            
            IJK = FUNIJK(I1+1,J1,K1) 
            IJK_WALL = FUNIJK(I1,J1,K1) 
            W_S(IJK_WALL,:) = -W_S(IJK,:)
         END DO 
      END DO 
      I1 = IMAX2 
      DO K1 = kmin3,kmax3 
         DO J1 = jmin3,jmax3 
            
   	    IF (.NOT.IS_ON_myPE_owns(I1,J1,K1)) CYCLE	 	    	    
            
            IJK = FUNIJK(I1-1,J1,K1) 
            IJK_WALL = FUNIJK(I1,J1,K1) 
            W_S(IJK_WALL,:) = -W_S(IJK,:)
         END DO 
      END DO 
      
      
      END SUBROUTINE MPPIC_BC_W_S  

      
      subroutine MPPIC_ADD_FRIC_FORCE(NP)
      USE param
      USE param1
      USE parallel
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE cutcell 

      IMPLICIT NONE 
      integer, intent(in) :: np
      INTEGER :: COUNT, COUNT_BC, IJK_WALL, IJK, idir, idim, I, J, K 
      character*80 :: wall_type 
      double precision :: vel_norm(dimn), vel_tang(dimn), normal(dimn), tangent(dimn)
      double precision :: vel_tang_mod, WALL_COOR(DIMN), DIST, WALLCOR_MIN(DIMN), WALLCOR_MAX(DIMN) 

      double precision ::  NORM_CF(3), XPOS, YPOS, ZPOS, max_dist, dist_fun, ramp_fun, FCN 
      logical :: doit 
      INCLUDE 'function.inc'
      IJK = PIJK(NP, 4) 

      FC(NP, :) = FC(NP,:) + PMASS(NP) * GRAV(:)

      max_dist = 8.d0*des_radius(NP)

      COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC 
      DO COUNT = 1, COUNT_BC
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IJK_WALL  = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%IJK_SCAL 
         if(cut_cell_at(ijk)) then 
            if(ijk_wall.ne.ijk) cycle
         else
            doit = .false. 
            doit = i.eq.imin1.or.i.eq.imax1
            doit = doit.or.j.eq.jmin1.or.j.eq.jmax1
            if(dimn.eq.3) doit = doit.or.k.eq.kmin1.or.k.eq.kmax1
            if(.not.doit) cycle 
         endif
         
         WALL_TYPE = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%DES_BC_TYPE

         SELECT CASE (TRIM(WALL_TYPE)) 
         
         CASE('NORMAL_WALL')
            WALLCOR_MIN(1) = XE(I-1)
            WALLCOR_MAX(1) = XE(I)
            WALLCOR_MIN(2) = YN(J-1)
            WALLCOR_MAX(2) = YN(J)
            
            IF(DIMN.EQ.3) THEN 
               WALLCOR_MIN(3) = ZT(K-1)
               WALLCOR_MAX(3) = ZT(K)
            END IF
            NORMAL(1:DIMN) = DES_CELLWISE_BCDATA(IJK)%BDRY_LIST(COUNT)%NORMAL(1:DIMN)
            !Find the direction of the normal 
            IDIR = 0
            DO IDIM = 1, DIMN
               IDIR = IDIR + ABS(NORMAL(IDIM))*IDIM 
            end DO

            WALL_COOR(1:DIMN)  = WALLCOR_MIN(1:DIMN)
            
            IF(NORMAL(IDIR).GT.0) WALL_COOR(IDIR) = WALLCOR_MAX(IDIR)
!let's say if the wall is the east wall of this scalar cell, then the wall cordinate 
!will be XE(I) which corresponds to WALLCOR_MAX 
               
            DIST =  NORMAL(IDIR)*(DES_POS_NEW(NP, IDIR) - WALL_COOR(IDIR))
            
            !according to the above convention for normal for 'normal_walls', 
!distance will be negative for particles inside the domain and 
!positive for particles outside the domain. This is becuase
!the wall normal points away from the physical domain, and any 
!point inside the domain will have negative distance. 
            DIST = -DIST
            NORMAL(:) = -NORMAL(:) 

         CASE('CUT_FACE')
            XPOS = DES_POS_NEW(NP,1) 
            YPOS = DES_POS_NEW(NP,2)
            ZPOS = ZERO 
            IF (DIMN .EQ. 3) THEN
               ZPOS = DES_POS_NEW(NP,3)
            ENDIF
            
            
            CALL GET_DEL_H_DES(IJK_wall,'SCALAR',XPOS , YPOS, ZPOS,& 
            & DIST, NORM_CF(1), NORM_CF(2), NORM_CF(3), .true.)

            NORMAL(1:DIMN) = NORM_CF(1:DIMN) 
            
         CASE DEFAULT
            !this might happen for small cells. do nothing in this case. 

         END SELECT
         FCN = (DOT_PRODUCT(FC(np,1:dimn), normal(1:dimn)))

         vel_norm(:) = (DOT_PRODUCT(des_vel_new(np,1:dimn), normal(1:dimn)))*normal(1:dimn)
         vel_tang(:) = des_vel_new(np,:) - vel_norm(:)
         
         vel_tang_mod = dot_product(vel_tang(1:dimn), vel_tang(1:dimn))
         vel_tang_mod = sqrt(vel_tang_mod)
         if(vel_tang_mod.gt.zero) then
            tangent(:) = vel_tang(:)/vel_tang_mod
         else
            tangent(:) = zero 
         endif
         
         !currently only treating those walls for friction that are native 
         !to this cell 
         dist_fun = min(dist/max_dist, 1.d0)
         dist_fun = dist_fun - 1.d0 
         ramp_fun = (1.d0 - exp(dist_fun))/(1.d0-exp(-1.d0))
         
         FC(NP, :) = FC(NP, :) - MEW_W*FCN*TANGENT(:)*ramp_fun
         !write(*,'(A,9(2x,g17.8))') 'vel, norm, tangent', des_vel_new(NP,:), normal(:), tangent(:)
         if(ramp_fun.lt.zero) write(*,'(A,3(2x,g17.8))') 'dist/maxdist, dist, ramp_fun ', dist_fun+1.d0, dist_fun, ramp_fun
         if(normal(2).eq.1.d0) then 
            write(*,'(A,9(2x,g17.8))') 'vel, norm, tangent', des_vel_new(NP,:), normal(:), tangent(:)
            write(*,'(A,3(2x,g17.8))') 'dist/maxdist, dist, ramp_fun ', dist_fun+1.d0, dist_fun, ramp_fun
            write(*,'(A,9(2x,g17.8))') 'FC, FCN, FCT', FC(np,:),FCN,  MEW_W*FCN*TANGENT(:)*ramp_fun
         !read(*,*)
      endif
      enddo
      
      FC(NP, :) = FC(NP,:) - PMASS(NP) * GRAV(:)

      end subroutine MPPIC_ADD_FRIC_FORCE
    
      
      
      SUBROUTINE MPPIC_APPLY_PS_GRAD

      USE param
      USE param1
      USE parallel
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mppic_wallbc
      USE physprop
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM 
      INTEGER I, J, K, IJK, IJK_C, IJK_OLD, IJK2, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, DP_BAR, COEFF_EN, MEANVEL(DIMN), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles  
      INTEGER PC , epg_min_loc(1)

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC 
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case 
      
      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, MEANUS(DIMN, MMAX), RELVEL(DIMN)
      DOUBLE PRECISION MEANUS_e(DIMN, MMAX), MEANUS_w(DIMN, MMAX),MEANUS_n(DIMN, MMAX),MEANUS_s(DIMN, MMAX),MEANUS_t(DIMN, MMAX), MEANUS_b(DIMN, MMAX)
      DOUBLE PRECISION :: DPS_DXE, DPS_DXW, DPS_DYN, DPS_DYS, DPS_DZT, DPS_DZB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, epg_min2, velf_part(dimn)
      INTEGER :: TOT_CASE, case1_count, case2_count, case3_count, case4_count 
      
!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      PC = 1
      FOCUS_PARTICLE = -1
      !DTPIC_CFL = LARGE_NUMBER 
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      
      EPG_MIN2 = MINVAL(EP_G(:))
      epg_min_loc = MINLOC(EP_G)
      IJK = epg_min_loc(1)
      !WRITE(*,*) 'IN APPLY_MPPIC_GRAD_PS '
      !WRITE(*,*) 'EPG MIN  = ', epg_min2, PINC(IJK)
      !WRITE(*,*) 'LOCATION = ',IJK, I_OF(IJK), j_of(ijk), k_of(ijk)

      case1_count = 0 ; case2_count = 0 ; case3_count = 0 ; case4_count = 0 

      J = JMIN2
      K = KMIN1
      !DO I = IMIN1, IMAX1
      !   IJK = funijk(i,j,k) 
      !   write(*,'(L2, 2x, A,2x,10(2x,g17.8))') FLUID_AT(IJK), 'GC VELS = ', U_G(IJK), V_G(IJK), P_G(IJK), RO_G(IJK)
      !enddo
!      STOP 


      DO L = 1, MAX_PIP
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 

         DES_LOC_DEBUG = .FALSE.
         ! If a particle is classified as new, then forces are ignored. 
! Classification from new to existing is performed in routine
! des_check_new_particle.f


         VEL_ORIG(:) = DES_VEL_NEW(L,:)
         IF(L.EQ.FOCUS_PARTICLE) THEN 
                  
            WRITE(*,'(A20,2x,3(2x,i4))') 'CELL ID = ', PIJK(L,1:3)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'EPS = ', 1.d0 - EP_g(PIJK(L,4))
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', DES_VEL_NEW(L,:)
            
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'FC = ', FC(L,:)
         ENDIF

!comment for Pradeep: Pradeep, Im not using this interpolaitonf or mean us
!right now but I want to keep it here as it was written with cut-cell in mind
!so it might be useful when i start looking at cut-cell again.
               
         M = PIJK(L,5)
         IJK = PIJK(L,4)

         IJK_OLD = IJK
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKE = IP_OF(IJK)
         IJKW = IM_OF(IJK)
         IJKN = JP_OF(IJK)
         IJKS = JM_OF(IJK)
         IJKT = KP_OF(IJK)
         IJKB = KM_OF(IJK)
         
         COEFF_EN = MPPIC_COEFF_EN

         XI_EAST = ZERO
         XI_WEST = ZERO
         XI_NORTH = ZERO
         XI_SOUTH = ZERO
         XI_TOP = ZERO
         XI_BOTTOM = ZERO
         
         !if(mod(L,500).eq.0) write(*,*) 'coeff_en = ', coeff_en, mppic_coeff_en
         IF(FLUID_AT(IJKE)) THEN 
            !DPS_DXE =  2.d0*(P_S(IJKE,1) - P_S(IJK,1))/(DX(I) + DX(I_OF(IJKE)))
            MEANUS_E(1,:) = (DES_U_S(IJKE,:)*DX(I) + DES_U_S(IJK,:)*DX(I_OF(IJKE)))/(DX(I) + DX(I_OF(IJKE)))
            XI_EAST = (DES_POS_NEW(L,1) - XE(I_OF(IJKW)))/DX(I)
         ELSE
            !DPS_DXE = 2.d0*(P_S(IJKE,1) - P_S(IJK,1))/(DX(I) + DX(I_OF(IJKE)))
            !DPS_DXE = zero!
            MEANUS_E(1,:) = zero !DES_U_S(IJK,:) 
            XI_EAST = (DES_POS_NEW(L,1) - XE(I_OF(IJKW)))/DX(I)
         ENDIF
         
         IF(FLUID_AT(IJKW)) THEN 
            !DPS_DXW =  2.d0*(P_S(IJK,1) - P_S(IJKW,1))/(DX(I) + DX(I_OF(IJKW)))
            MEANUS_W(1,:) = (DES_U_S(IJKW,:)*DX(I) + DES_U_S(IJK,:)*DX(I_OF(IJKW)))/(DX(I) + DX(I_OF(IJKW)))
            XI_WEST = (XE(I) - DES_POS_NEW(L,1))/DX(I)
         ELSE
            !DPS_DXW =  2.d0*(P_S(IJK,1) - P_S(IJKW,1))/(DX(I) + DX(I_OF(IJKW)))
            !DPS_DXW = zero!
            MEANUS_W(1,:) = zero !DES_U_S(IJK,:) 
            XI_WEST = (XE(I) - DES_POS_NEW(L,1))/DX(I)
         ENDIF
         DPS_DXE = PS_FORCE_PIC(IJK,1)
         DPS_DXW = PS_FORCE_PIC(IJKW,1)
         
         
         PS_FORCE(1) = XI_EAST*DPS_DXE + XI_WEST*DPS_DXW
         VELF_PART(1) = XI_EAST*U_G(IJK) + XI_WEST*U_G(IJKW)
         MEANUS(1,:) = XI_EAST*MEANUS_E(1,:) + XI_WEST*MEANUS_W(1,:)
         !MEANUS(1,:) = XI_EAST*U_S(IJK,:) + XI_WEST*U_S(IJKW,:)
         IF(FLUID_AT(IJKN)) THEN 
            !DPS_DYN =  2.d0*(P_S(IJKN,1) - P_S(IJK,1))/(DY(J) + DY(J_OF(IJKN)))
            XI_NORTH = (DES_POS_NEW(L,2) - YN(J_OF(IJKS)))/DY(J) 
            MEANUS_N(2,:) = (DES_V_S(IJKN,:)*DY(J) + DES_V_S(IJK,:)*DY(J_OF(IJKN)))/(DY(J) + DY(J_OF(IJKN)))
         ELSE
            !DPS_DYN = 2.d0*(P_S(IJKN,1) - P_S(IJK,1))/(DY(J) + DY(J_OF(IJKN)))
            !DPS_DYN = zero!
            XI_NORTH = (DES_POS_NEW(L,2) - YN(J_OF(IJKS)))/DY(J) 
            MEANUS_N(2,:) = zero !DES_V_S(IJK,:)  
         ENDIF
         
         IF(FLUID_AT(IJKS)) THEN 
            !DPS_DYS =  2.D0*(P_S(IJK,1) - P_S(IJKS,1))/(DY(J) + DY(J_OF(IJKS)))
            XI_SOUTH = (YN(J_OF(IJK)) - DES_POS_NEW(L,2))/DY(J) 
            MEANUS_S(2,:) = (DES_V_S(IJKS,:)*DY(J) + DES_V_S(IJK,:)*DY(J_OF(IJKS)))/(DY(J) + DY(J_OF(IJKS)))
         ELSE
            
            !DPS_DYS =  2.D0*(P_S(IJK,1) - P_S(IJKS,1))/(DY(J) + DY(J_OF(IJKS)))
            !DPS_DYS = zero! 2.D0*(P_S(IJK,1) - P_S(IJKS,1))/(DY(J) + DY(J_OF(IJKS)))
            XI_SOUTH = (YN(J_OF(IJK)) - DES_POS_NEW(L,2))/DY(J) 
            MEANUS_S(2,:) = zero !DES_V_S(IJK,:)  
         ENDIF

         DPS_DYN = PS_FORCE_PIC(IJK,2)
         DPS_DYS = PS_FORCE_PIC(IJKS,2)
         
         PS_FORCE(2) = XI_NORTH*DPS_DYN + XI_SOUTH*DPS_DYS
               
         VELF_PART(2) = XI_NORTH*V_G(IJK) + XI_SOUTH*V_G(IJKS)
         
         MEANUS(2,:) = XI_NORTH*MEANUS_N(2,:) + XI_SOUTH*MEANUS_S(2,:)
         !MEANUS(2,:) = XI_NORTH*V_S(IJK,:) + XI_SOUTH*V_S(IJKS,:)
         
         IF(DIMN.eq.3) then 
            IF(FLUID_AT(IJKT)) THEN 
               !DPS_DZT =  2.d0*(P_S(IJKT,1) - P_S(IJK,1))/(DZ(K) + DZ(K_OF(IJKT)))
               XI_TOP = (DES_POS_NEW(L,3) - ZT(K_OF(IJKB)))/DZ(K)
               MEANUS_T(3,:) = (DES_W_S(IJKT,:)*DZ(K) + DES_W_S(IJK,:)*DZ(K_OF(IJKT)))/(DZ(K) + DZ(K_OF(IJKT)))
            ELSE
               !DPS_DZT = zero!
               !DPS_DZT =  2.d0*(P_S(IJKT,1) - P_S(IJK,1))/(DZ(K) + DZ(K_OF(IJKT)))
               XI_TOP = (DES_POS_NEW(L,3) - ZT(K_OF(IJKB)))/DZ(K)
               MEANUS_T(3,:) = zero !DES_W_S(IJKT,:)
            ENDIF
            
            IF(FLUID_AT(IJKB)) THEN 
               !DPS_DZB =  2.d0*(P_S(IJK,1) - P_S(IJKB,1))/(DZ(K) + DZ(K_OF(IJKB)))
               XI_BOTTOM = (ZT(K_OF(IJK)) - DES_POS_NEW(L,3))/DZ(K)
               MEANUS_B(3,:) = (DES_W_S(IJKB,:)*DZ(K) + DES_W_S(IJK,:)*DZ(K_OF(IJKB)))/(DZ(K) + DZ(K_OF(IJKB)))
            ELSE
               !DPS_DZB =  2.d0*(P_S(IJK,1) - P_S(IJKB,1))/(DZ(K) + DZ(K_OF(IJKB)))
               !DPS_DZB = ZERO
               XI_BOTTOM = (ZT(K_OF(IJK)) - DES_POS_NEW(L,3))/DZ(K)
               MEANUS_B(3,:) = zero !DES_W_S(IJK,:)
            ENDIF
            
            DPS_DZT = PS_FORCE_PIC(IJK,3)
            DPS_DZB = PS_FORCE_PIC(IJKB,3)
            PS_FORCE(3) = XI_TOP*DPS_DZT + XI_BOTTOM*DPS_DZB
            
            VELF_PART(3) = XI_TOP*W_G(IJK) + XI_BOTTOM*W_G(IJKB)

            MEANUS(3,:) = XI_TOP*MEANUS_T(3,:) + XI_BOTTOM*MEANUS_B(3,:)
            !MEANUS(3,:) = XI_TOP*W_S(IJK,:) + XI_BOTTOM*W_S(IJKB,:)
            
         ENDIF
         
         PS_FORCE(:) = PS_GRAD(L, :) 
         !IF(ABS(PS_FORCE(2)).GT.ZERO)  WRITE(*,*) 'PS_FORCE = ', PS_FORCE
         dp_bar = zero 
         DELUP(:) =-( DTSOLID*PS_FORCE(:))

         DELUP(:) = DELUP(:)/ROP_S(IJK_OLD,M)
         
         MEANVEL(1) = DES_U_S(IJK_OLD,M)
         MEANVEL(2) = DES_V_S(IJK_OLD,M)
         IF(DIMN.EQ.3) MEANVEL(3) = DES_W_S(IJK_OLD,M)
         
          !if(MOD(L,500).eq.0) write(*,*) 'mean sol vel =', meanus(:,M)
         MEANUS(:,M) =  AVGSOLVEL_P (L,:)
         
         !MEANUS(:,M) = MEANVEL(:)
         !RELVEL(:) = DES_VEL_NEW(L,:) - MEANUS(:,M)
         !DO IDIM = 1, DIMN
        !    IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle 

         !   IF(RELVEL(IDIM)*DELUP(IDIM).GT.ZERO) THEN
               !do nothing
        !    ELSE
        !       DES_VEL_NEW(L,IDIM) = MEANUS(IDIM,M) - 0.4d0*RELVEL(IDIM)
               
               !IF(DES_VEL_NEW(L,IDIM)*DELUP(IDIM).LT.ZERO) DES_VEL_NEW(L,IDIM) = -0.5d0*DES_VEL_NEW(L,IDIM)
        !    ENDIF
        ! ENDDO
        ! CYCLE 

         RELVEL(:) = DES_VEL_NEW(L,:) - MEANUS(:,M)
         
         DO IDIM = 1, DIMN

            IF(ABS(PS_FORCE(IDIM)).eq.zero) cycle 
            IF(DES_VEL_NEW(L,IDIM).GE.ZERO) then 
               signvel = 1.d0
            else
               signvel = -1.d0
            endif
            IF(DES_VEL_NEW(L,IDIM)*MEANUS(IDIM,M).GT.ZERO) THEN 
               
               IF(DES_VEL_NEW(L,IDIM)*DELUP(IDIM).GT.ZERO) THEN 
                  
                  IF(ABS(MEANUS(IDIM,M)) .GT. ABS(DES_VEL_NEW(L,IDIM))) THEN 
                       IJK_C = IJK
                     IF(IDIM.eq.1) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                     ELSEIF(IDIM.eq.2) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                     ELSEIF(IDIM.eq.3) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                     ENDIF
                     if(fluid_at(IJK_C)) then 
                        DES_VEL_NEW(L,IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                        !DES_VEL_NEW(L,IDIM) = (1.D0+COEFF_EN)*DES_VEL_NEW(L,IDIM) 
                     endif
                      case4_count = case4_count + 1
                  ENDIF
               ELSE 
                  IF(ABS(DES_VEL_NEW(L,IDIM)).GT.ABS(MEANUS(IDIM,M))) then 
                     !DES_VEL_NEW(L,IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                     IJK_C = IJK
                     IF(IDIM.eq.1) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = IM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = IP_OF(IJK)
                     ELSEIF(IDIM.eq.2) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = JM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = JP_OF(IJK)
                     ELSEIF(IDIM.eq.3) then 
                        if(DES_VEL_NEW(L,IDIM).LT.ZERO) IJK_C = KM_OF(IJK)
                        if(DES_VEL_NEW(L,IDIM).GT.ZERO) IJK_C = KP_OF(IJK)
                     ENDIF

                     if((IDIM.EQ.2.AND.DES_VEL_NEW(L,IDIM).LT.ZERO).or.(.not.fluid_at(IJK_C))) then 
                     !if(.not.fluid_at(IJK_C)) then 
                        DES_VEL_NEW(L,IDIM) = -COEFF_EN*DES_VEL_NEW(L,IDIM)
                        !DES_VEL_NEW(L,IDIM) = coeff_en*des_vel_new(L, IDIM)

                     else
                        DES_VEL_NEW(L,IDIM) = (MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM))
                        !DES_VEL_NEW(L,IDIM) = coeff_en*des_vel_new(L, IDIM)
                     endif
                     
                     case1_count = case1_count + 1
                  ELSE
                     !do nothing
                                !DES_VEL_NEW(L,IDIM) = DES_VEL_NEW(L,IDIM)
                     case1_count = case1_count + 1
                     
                  ENDIF
               ENDIF
            ELSE
               IF(MEANUS(IDIM,M)*DELUP(IDIM).GT.ZERO) THEN 
                  DES_VEL_NEW(L,IDIM) = MEANUS(IDIM,M) - COEFF_EN*RELVEL(IDIM)
                  !DES_VEL_NEW(L,IDIM) = -COEFF_EN*DES_VEL_NEW(L,IDIM)
                  
                  case2_count = case2_count + 1
               ELSE 
                  case3_count = case3_count + 1
                  !DO NOTHING 
               ENDIF
            ENDIF
         ENDDO

         !
         if(L.eq.FOCUS_PARTICLE) THEN 
         !iF((IJK.eq.epg_min_loc(1).or.IJK_OLD.eq.epg_min_loc(1)).and.epg_min2.lt.0.38) then
                                !if(j.ne.2) cycle 
            WRITE(*,'(A20,2x,i6, 4(2x,g17.8))') 'L,XIE, XIW, XIN, XIS', L, XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH
            
            WRITE(*,'(A20,2x,3(2x,i5))') 'ORIGINAL I, J, K =', I_OF(IJK_OLD),J_OF(IJK_OLD),K_OF(IJK_OLD) 
            WRITE(*,'(A20,2x,3(2x,i5))') 'PIJK I, J, K =', I_OF(IJK),J_OF(IJK),K_OF(IJK) 
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', VEL_ORIG(:)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL NEW = ', DES_VEL_NEW(L,:)
            !WRITE(*,'(A20,2x,3(2x,g17.8))') 'MEANUS N and S = ', MEANUS_N(2,:), MEANUS_S(2,:)
            
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'MEANVEL = ', MEANVEL(:)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'MEANUS = ', MEANUS(:,1)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'FVEL = ', VELF_PART(:)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_POS_NEW = ', DES_POS_NEW(L,:)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'GRAD PS = ', PS_FORCE(:)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DPS_DYN, S = ', DPS_DYN, DPS_DYS

            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DELUP =  ', DELUP(:)

            !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU_INT = ', UPRIMETAU_INT(:)
            !WRITE(*,'(A20,2x,3(2x,g17.8))') 'UPRIMETAU = ', UPRIMETAU(:)
            !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) WRITE(*,'(A20,2x,3(2x,g17.8))') 'U*DT, MFP =', UPRIMEMOD*DTSOLID, MEAN_FREE_PATH
            read(*,*)
         ENDIF
      end DO
      
      TOT_CASE = case1_count + case2_count + case3_count + case4_count 
      !IF(TOT_CASE.GT.0) THEN
      !WRITE(*,'(A, 4(2x,i10))') 'CASE COUNT NUMBERS  = ', case1_count ,case2_count ,case3_count ,case4_count 
      !WRITE(*,'(A, 4(2x,g12.7))') 'CASE COUNT %AGE = ', real(case1_count)*100./real(tot_case),real(case2_count)*100./real(tot_case), real(case3_count)*100./real(tot_case), real(case4_count)*100./real(tot_case)
      !ENDIF
      RETURN 


      end SUBROUTINE MPPIC_APPLY_PS_GRAD


      subroutine mppic_avg_eps

        !this needs to refined futrther. This was taken from somewhere in 
        !particles in cell and dumped here to keep the code orderly: Rah
!!$            MASS_SOL2 = zero 
!!$      DO IJK = ijkstart3, ijkend3
!!$         I = I_OF(IJK)
!!$         J = J_OF(IJK)
!!$         K = K_OF(IJK)
!!$         
!!$         IF(.NOT.FLUID_AT(IJK).OR..NOT.IS_ON_myPE_owns(I, J, K)) CYCLE 
!!$         !for cut-cell, it is important to check both. FLUID_AT(IJK) 
!!$         !alone is not enough as it is shared between procs and 
!!$         !fluid_at(ijkend3) might be true when in fact it does 
!!$         !not belong to that proc 
!!$         !IF(.NOT.CUT_CELL_AT(IJK)) CYCLE 
!!$         
!!$         IF(DIMN.EQ.2) THEN 
!!$            KPLUS1 = K
!!$            KMINUS1 = K
!!$         ELSE
!!$            KPLUS1 =  MIN(K+1, KEND3) !in serial kend3 = kend2. therefore, the min and max ops
!!$            KMINUS1 = MAX(K-1, KSTART3)
!!$         ENDIF
!!$                  
!!$         JPLUS1 =  MIN(J+1, JEND3)
!!$         JMINUS1 = MAX(J-1, JSTART3)
!!$         
!!$         IPLUS1 =  MIN(I+1, IEND3)
!!$         IMINUS1 = MAX(I-1, ISTART3)
!!$         !COUNT_AVG = 1          !1 for ijk itself 
!!$         !VOL_AVG = VOL(IJK) 
!!$         !ROP_S(IJK,:) = ROP_SO(IJK, :)*VOL(IJK)
!!$         VOL_AVG = ZERO 
!!$         ROP_S(IJK,:) = ZERO 
!!$         COUNT_AVG = 0
!!$         DO KK = KMINUS1, KPLUS1
!!$            DO JJ = JMINUS1, JPLUS1
!!$               DO II = IMINUS1, IPLUS1 
!!$                  IJK2 = FUNIJK(II,JJ,KK)
!!$
!!$                  !IF(FLUID_AT(IJK2).AND.(IJK2.NE.IJK)) THEN 
!!$                  IF(FLUID_AT(IJK2)) THEN 
!!$                     !if(rop_so(ijk2,1).gt.zero) then 
!!$                        VOL_AVG = VOL_AVG + VOL(IJK2) 
!!$                        COUNT_AVG = COUNT_AVG + 1
!!$                        ROP_S(IJK, :) = ROP_S(IJK, :) + ROP_SO(IJK2, :)*VOL(IJK2)/(ro_s(:)*pvol(1))
!!$                     !endif
!!$                  ENDIF
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$         ROP_S(IJK,:) = (ROP_S(IJK,:)*(ro_s(:)*pvol(1)))/VOL_AVG !REAL(COUNT_AVG)
!!$         !ROP_S(IJK,:) = ROP_S(IJK,:)/REAL(COUNT_AVG)
!!$         SUM_EPS = ZERO
!!$         DO M = 1, MMAX 
!!$            MASS_SOL2(M) = MASS_SOL2(M) + ROP_S(IJK,M)*VOL(IJK)
!!$            SUM_EPS = SUM_EPS + EP_S(IJK,M) 
!!$         ENDDO
!!$
!!$         !EPG_OLD = EP_G(IJK) 
!!$         EP_G(IJK) = 1.d0
!!$         EP_G(IJK) = EP_G(IJK) - SUM_EPS 
!!$
!!$         
!!$         !WRITE(*,*) 'COUNT_AVG = ', COUNT_AVG, EPG_OLD, EP_G(IJK), SUM_EPS
!!$         IF(EP_G(IJK).LT.ZERO .AND. DES_CONTINUUM_COUPLED) THEN 
!!$               
!!$            if(dmp_log) write(unit_log,1012)  IJK, CUT_CELL_AT(IJK), I_OF(IJK), I_OF(IJK), I_OF(IJK), EP_G(IJK), ROP_S(IJK,:)/RO_S(:) 
!!$            IF(PRINT_DES_SCREEN) write(*,1012)  IJK, CUT_CELL_AT(IJK), I_OF(IJK), I_OF(IJK), I_OF(IJK), EP_G(IJK), ROP_S(IJK,:)/RO_S(:) 
!!$            READ(*,*) 
!!$         ENDIF 
!!$         ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
!!$         
!!$      ENDDO
    end subroutine mppic_avg_eps
    
