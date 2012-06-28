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
!           Note that P_force is evaluated as -dp/dx                   C
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
!  Subroutine: CALC_DES_ROP_S                                          C
!  Purpose: Calculate the bulk density of particles belonging to the   C
!           same 'mth' solids phase associated with the ijk fluid cell C
!           based on interpolation of particle position (only invoked  C
!           when des_interp_on)                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_DES_ROP_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices      
      USE compar
      USE sendrecv      
      USE discretelement
      USE interpolation
      use desmpi 
      USE cutcell 
      USE mfix_pic
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
! volume of fluid cell, & one over the volume of fluid cell
      DOUBLE PRECISION :: VCELL, OVOL
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
!-----------------------------------------------

! initializations      
      wtbar = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.250d0*(DIMN-2) + 0.50d0*(3-DIMN)

! sets several quantities including interp_scheme, scheme, and 
! order and allocates arrays necessary for interpolation      
      call set_interpolation_scheme(2)
      
! There is some issue associated to gstencil, vstencil which are 
! allocable variables

!!$      omp_start=omp_get_wtime()
!!$omp parallel do default(shared)                                 &
!!$omp private(ijk,i,j,k,pcell,iw,ie,js,jn,kb,ktp,                 &
!!$omp         avg_factor,ii,jj,kk,cur_ijk,ipjk,ijpk,ipjpk,        &
!!$omp         ijpkp,ipjkp,ipjpkp,ijkp,icur,jcur,kcur              &
!!$omp         gstencil,vstencil,nindx,onew,                       &
!!$omp         focus,np,wtp,weightp,ovol,m,vcell)                  &
!!$omp schedule (guided,50)           	
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
                  IF(dimn.EQ.3) THEN                  
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
         DO nindx = 1,pinc(ijk)
            np = pic(ijk)%p(nindx)            
  
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

            focus = .false.
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)

! Calculate wtbar so des_rop_s, and in turn, ep_g can be updated 
!----------------------------------------------------------------->>>
            M = pijk(np,5)
            DO k = 1, (3-dimn)*1+(dimn-2)*onew
               DO j = 1, onew
                  DO i = 1, onew
! shift loop index to new variables for manipulation                     
                     ii = iw + i-1
                     jj = js + j-1
                     kk = kb + k-1
! The interpolation is done using node. so one should use consistent
! numbering system in the current version imap_c is used instead of 
! ip_of or im_of                       
                     icur = imap_c(ii)
                     jcur = jmap_c(jj)
                     kcur = kmap_c(kk)                     
                     cur_ijk = funijk(icur, jcur, kcur) !imap_c(ii),jmap_c(jj),kmap_c(kk))

! Replacing the volume of cell to volume at the node
                     vcell = des_vol_node(cur_ijk)
                     ovol = one/vcell

!!$omp critical                 
                     wtbar(cur_ijk,m) = wtbar(cur_ijk,m) + &
                          weightp(i,j,k) *DES_ro_s(m)*ovol*pvol(np)*WTP
!!$omp end critical                             
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel do
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'drag_interp:',omp_end - omp_start  

! At the interface wtbar has to be added
! send recv will be called and the node values will be added 
! at the junction. wtbar is altered by the routine when
! periodic boundaries are invoked
      call des_addnodevalues2
!-----------------------------------------------------------------<<<


! Calculate/update the cell centered bulk density 
!----------------------------------------------------------------->>>  
! avg_factor=0.125 (in 3D) or =0.25 (in 2D)
      AVG_FACTOR = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)   

!$omp parallel do default(shared)                               &
!$omp private(ijk,i,j,k,imjk,ijmk,imjmk,ijkm,imjkm,ijmkm,       & 
!$omp         imjmkm)                                           &
!$omp schedule (guided,20)           
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
            DES_rop_s(ijk,:) = avg_factor*(wtbar(ijk,:) +&
                 wtbar(ijmk,:) + wtbar(imjmk,:) +&
                 wtbar(imjk,:))
            IF(dimn.EQ.3) THEN
               ijkm = funijk(imap_c(i),jmap_c(j),kmap_c(k-1))
               imjkm = funijk(imap_c(i-1),jmap_c(j),kmap_c(k-1))
               ijmkm = funijk(imap_c(i),jmap_c(j-1),kmap_c(k-1))
               imjmkm = funijk(imap_c(i-1),jmap_c(j-1),kmap_c(k-1))
               DES_rop_s(ijk,:) = DES_rop_s(ijk,:) + avg_factor*&
                    (wtbar(ijkm,:) + wtbar(ijmkm,:) + &
                    wtbar(imjmkm,:)+wtbar(imjkm,:) )
            ENDIF
! to more closely mimic current implementation
! (should be able to remove this...)            
            IF (.NOT. DES_CONTINUUM_HYBRID) THEN
               DO M = 1,DES_MMAX
                  ROP_S(IJK,M) = DES_ROP_S(IJK,M)
               ENDDO
            ENDIF                       
         ENDIF   ! end if (fluid_at(ijk))
         
      ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do 
!-----------------------------------------------------------------<<<

      RETURN   
      END SUBROUTINE CALC_DES_ROP_S



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS_NONINTERP                              C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions.                             C
!     - Calculates the drag force exerted on the particles by the      C
!       gas/fluid phase using cell average quantities.  The drag       C
!       coefficient (F_GS) is calculated using the subroutine          C
!       drag_gs, which is called during the continuum time step        C
!       (from calc_drag) and during the discrete time(s) (from here).  C
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
                          VELDS_ARR(DIMN, DES_MMAX), &
                          VELCS_ARR(DIMN, MMAX)
! local drag force
      DOUBLE PRECISION :: GS_DRAG (DIMENSION_3, DES_MMAX, DIMN)
! index of solid phase that particle NP belongs to      
      INTEGER :: M
! particle number index, used for looping      
      INTEGER :: NP, nindx
! solids volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM      
! one over solids volume fraction in fluid cell and one over the
! volume of fluid cell
      DOUBLE PRECISION :: OEPS, OVOL 
! for error messages      
      INTEGER :: IER
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.  
      LOGICAL :: DISCRETE_FLAG
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------


! Calculate the gas-solid drag force exerted on the particle
!----------------------------------------------------------------->>>
! initializations              
      GS_DRAG(:,:,:) = ZERO
      gd_force(:,:) = ZERO

! computing F_gs (drag coefficient) with the latest average fluid 
! and solid velocity fields.  see comments below
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
!$omp         ugc,vgc,wgc,velg_arr,velds_arr,                     &
!$omp         np,nindx,m,oeps,ep_sm,gs_drag)                      &
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

            DO M = 1, DES_MMAX
! the call to drag coefficient should probably be here for cut-cell
! since cell center velocities are now known and these would be used 
! in vrel within drag correlations). this would require some 
! rearrangement/reworking?  alternatively these calculations should be
! moved into drag_gs?

               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)

               IF(EP_SM.GT.ZERO) THEN
                  OEPS = ONE/EP_SM

                  IF (MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
                     GS_DRAG(IJK,M, :) = F_GS(IJK,M)*VELG_ARR(:)
                  ELSEIF (DES_CONTINUUM_HYBRID) THEN
                     GS_DRAG(IJK,M,:) = -F_GDS(IJK,M)*&
                        (VELDS_ARR(:,M)-VELG_ARR(:))
                  ELSE ! default cuase
                     GS_DRAG(IJK,M,:) = -F_GS(IJK,M)*&
                        (VELDS_ARR(:,M)-VELG_ARR(:))
                  ENDIF   ! end if/else (des_continuum_hybrid)

                  GS_DRAG(IJK,M,:) = GS_DRAG(IJK,M,:)*OEPS

               ENDIF  ! end if ep_sm>0               

            ENDDO   ! end do loop (dm=1,des_mmax)

! loop through particles in the cell and determine the drag force on
! each particle 
            DO nindx = 1,PINC(IJK)
               NP = PIC(ijk)%p(nindx)
! skipping indices that do not represent particles and ghost particles
               if(.not.pea(np,1)) cycle 
               if(pea(np,4)) cycle

               M = PIJK(NP,5)
               GD_FORCE(NP,:) = GS_DRAG(IJK,M,:)*PVOL(NP)

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag
               FC(NP,:) = FC(NP,:) + GD_FORCE(NP,:) 
               IF(.NOT.MODEL_B) THEN    
! Add the pressure gradient force
! P_force is evaluated as -dp/dx
                  FC(NP,:) = FC(NP,:) + (P_FORCE(IJK,:))*PVOL(NP)
               ENDIF               
            ENDDO   ! end do loop (nindex=1,pinc(ijk))

         ENDIF      ! end if(pinc(ijk).gt.0)

      ENDDO         ! end do loop (ijk=ijkstart3, ijkend3)
!$omp end parallel do          


!----------------------------------------------------------------->>>
! For mppic calculate the gas-particle drag coefficient (F_GP).
      IF(MPPIC) THEN
!$omp parallel do private(np,ijk,m,oeps,ep_sm)  &
!$omp schedule (guided,100)    
         DO NP = 1, MAX_PIP
! skipping indices that do not represent particles and ghost particles
            if(.not.pea(np,1)) cycle 
            if(pea(np,4)) cycle 

            IJK = PIJK(NP,4)
            M = PIJK(NP,5)

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
         ENDDO   ! end do loop (np=1,max_pip)
!$omp end parallel do
      ENDIF   ! end if(mppic)
!-----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE CALC_DES_DRAG_GS_NONINTERP
    


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DES_DRAG_GS                                        C
!  Purpose: This subroutine is only called from the DISCRETE side.     C
!     It performs the following functions:                             C
!     - If non-interpolated then execution of the code is directed     C
!       to the subroutine calc_des_drag_gs_noninterp for the           C
!       appropriate calculations                                       C
!     - If interpolated, then it calculates the particle centered      C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity. This F_GP is used to calculate    C
!       the gas-solids drag on the particle.                           C
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
!-----------------------------------------------   
! Include statement functions
!-----------------------------------------------   
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!----------------------------------------------- 

! non-interpolated drag:
! Calculate the gas solids drag force on a particle using the cell 
! averaged particle velocity and the cell average fluid velocity
!----------------------------------------------------------------->>>
      IF(.NOT.DES_INTERP_ON) THEN
         CALL CALC_DES_DRAG_GS_NONINTERP
         RETURN 
      ENDIF
!-----------------------------------------------------------------<<<


! interpolated drag (the rest of this routine):
! Calculate the gas solids drag coefficient using the particle 
! velocity and the fluid velocity interpolated to particle
! position. 
!----------------------------------------------------------------->>>
! initializations
      gd_force(:,:) = ZERO

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
!!$omp         focus,np,weightp,m)                                 &
!!$omp schedule (guided,50)           	
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
! skipping indices that do not represent particles and ghost particles
            if(.not.pea(np,1)) cycle
            if(pea(np,4)) cycle

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
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            CALL DES_DRAG_GP(NP, VEL_FP(NP,1:DIMN), &
               DES_VEL_NEW(NP,1:DIMN))

! Calculate the gas-solids drag force on the particle
            IF(MPPIC .AND. MPPIC_PDRAG_IMPLICIT) THEN
! implicit treatment of the drag term for mppic 
               GD_FORCE(NP,:) = F_GP(NP)*(VEL_FP(NP,:))
            ELSE
! default case
               GD_FORCE(NP,:) = F_GP(NP)*(VEL_FP(NP,:)-DES_VEL_NEW(NP,:))
            ENDIF

! Update the contact forces (FC) on the particle to include gas
! pressure and gas-solids drag 
            FC(NP,:) = FC(NP,:) + GD_FORCE(NP,:)
            IF(.NOT.MODEL_B) THEN
! P_force is evaluated as -dp/dx 
               FC(NP,:) = FC(NP,:) + p_force(ijk,:)*pvol(NP)
            ENDIF               

         ENDDO       ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel do
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
!     - If interpolated then, it calculates the particle centered      C
!       drag coefficient (F_GP) based on the particle velocity and     C
!       interpolated fluid velocity.                                   C
!     - It determines the contribution of particle centered drag to    C
!       the center coefficient of the A matrix and the b (source)      C
!       vector in the matrix equation (A*VEL_FP=b) equation for the    C
!       gas phase x, y and z momentum balances using F_GP              C
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
      use desmpi 
      USE cutcell 
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable used for debugging
      LOGICAL :: FOCUS 
! local drag forces
      DOUBLE PRECISION :: drag_bm_tmp(DIMN)
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
!-----------------------------------------------   
! Include statement functions
!-----------------------------------------------   
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime	      

!!$      omp_start=omp_get_wtime()

! NON-INTERPOLATED DRAG:
! Calculate the gas solids drag coefficient (F_GS) using the cell 
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


! INTERPOLATED DRAG (the rest of this routine):
! Calculate the gas solids drag coefficient using the particle 
! velocity and the fluid velocity interpolated to particle
! position. 
!----------------------------------------------------------------->>>
! initializations 
      drag_am = ZERO
      drag_bm = ZERO

! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.250d0*(DIMN-2) + 0.50d0*(3-DIMN)

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
!!$omp         focus,np,wtp,weightp,ovol,m,nindx,                  &
!!$omp         drag_bm_tmp)                                        &
!!$omp schedule (guided,50)
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
! skipping indices that do not represent particles and ghost particles
            if(.not.pea(np,1)) cycle
            if(pea(np,4)) cycle

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
! obtained from des_drag_gp subroutine is given as:
!    f_gp=beta*vol_p/eps where vol_p is the particle volume.
! The drag force on each particle is equal to:
!    beta(u_g-u_s)*vol_p/eps.
! Therefore, the drag force = f_gp*(u_g - u_s)
            CALL DES_DRAG_GP(NP, VEL_FP(NP,1:DIMN), &
               DES_VEL_NEW(NP,1:DIMN))
!-----------------------------------------------------------------<<<


! Calculate the corresponding gas solids drag force that is used in 
! the gas phase momentum balances.               
!----------------------------------------------------------------->>>
            focus = .false.
            M = pijk(np,5)
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)

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
!!$omp end critical
                  ENDDO
               ENDDO
            ENDDO

         ENDDO   ! end do (nindx = 1,pinc(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
!!$omp end parallel do
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
      AVG_FACTOR = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)

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
            IF(DIMN.EQ.3) THEN
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
            ENDIF   ! end if (dimn.eq.3)
         ELSE   ! else branch of if (fluid_at(ijk))
            IF (DES_CONTINUUM_HYBRID) THEN
               F_GDS(IJK,M) = ZERO
            ELSE
               F_GS(IJK,M) = ZERO
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

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle number id.
      INTEGER , INTENT(IN) :: LL
! fluid velocity interpolated to particle position      
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: FLUID_VEL
! particle velocity
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: PARTICLE_VEL
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
! solids volume fraction, particle diameter      
      DOUBLE PRECISION :: EPG, ROg, ROPg, EP_SM, DPM
!-----------------------------------------------    
! Include statement functions
!-----------------------------------------------  
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
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
         ENDDO
         DO CM = 1,SMAX
            L = DES_MMAX + CM
            DP_loc(L) = D_P(IJK,CM)
            EPs_loc(L) = EP_S(IJK,CM)
         ENDDO
      ENDIF   ! end if/else (.not.des_continuum_hybrid)


! magnitude of gas-particle relative velocity
      IF(DIMN == 2)THEN
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
      Y_i = DP_loc(M)/DPA

! assign variables for short dummy arguments
      EPg = EP_G(IJK)
      ROg = RO_G(IJK)
      ROPg = ROP_G(IJK)
      EP_SM = EPs_loc(M)
      DPM = DP_loc(M)

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
         SELECT CASE(TRIM(DRAG_TYPE))
         CASE ('SYAM_OBRIEN')
            CALL DRAG_SYAM_OBRIEN(DgA,EPG,Mu,ROg,VREL,&
                DPM)
         CASE ('GIDASPOW') 
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,&
                 DPM)
         CASE ('GIDASPOW_PCF')
            CALL DRAG_GIDASPOW(DgA,EPg,Mu,ROg,ROPg,VREL,&
                 DPA)
         CASE ('GIDASPOW_BLEND')
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,&
                 DPM)
         CASE ('GIDASPOW_BLEND_PCF')
            CALL DRAG_GIDASPOW_BLEND(DgA,EPg,Mu,ROg,ROPg,VREL,&
                 DPA)
         CASE ('WEN_YU')
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,&
                 DPM)
         CASE ('WEN_YU_PCF')
            CALL DRAG_WEN_YU(DgA,EPg,Mu,ROPg,VREL,&
                 DPA)
         CASE ('KOCH_HILL')
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,&
                 DPM,DPM,phis)
         CASE ('KOCH_HILL_PCF')
            CALL DRAG_KOCH_HILL(DgA,EPg,Mu,ROPg,VREL,&
                 DPM,DPA,phis)
         CASE ('BVK')
            CALL DRAG_BVK(DgA,EPg,Mu,ROPg,VREL,&
                 DPM,DPA,phis)
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
      IF(TRIM(DRAG_TYPE) == 'GIDASPOW_PCF' .OR. &
         TRIM(DRAG_TYPE) == 'GIDASPOW_BLEND_PCF' .OR. &
         TRIM(DRAG_TYPE) == 'WEN_YU_PCF' .OR. &
         TRIM(DRAG_TYPE) == 'KOCH_HILL_PCF' .OR. &
         TRIM(DRAG_TYPE) == 'BVK') THEN
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

