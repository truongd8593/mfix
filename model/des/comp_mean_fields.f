!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS

      use discretelement, only: DES_INTERP_MEAN_FIELDS

      IF(DES_INTERP_MEAN_FIELDS) THEN
         CALL COMP_MEAN_FIELDS_INTERP
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF


      END SUBROUTINE COMP_MEAN_FIELDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle index
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase index
      INTEGER M, CM
! ijk indices
      INTEGER I, J, K, IJK
! Variable to distribute particle volume
      DOUBLE PRECISION ::  WTP
! 1 over volume of fluid cell
      DOUBLE PRECISION :: OVOL
! total volume of mth phase solids in cell and 1 over that value
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,DES_MMAX), OSOLVOL
! solids volume fraction of mth solids phase
      DOUBLE PRECISION EP_SM
! total solids volume fraction of continuum solids phases
      DOUBLE PRECISION SUM_EPS
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

! Add particle values (volume, velocity) to ongoing summations
!!$      omp_start1=omp_get_wtime()
!!$omp single private(l,wtp,i,j,k,ijk,m)
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) EXIT
! skipping indices with no particles (non-existent particles)
         IF(.NOT.PEA(L,1)) CYCLE
! incrementing particle account when particle exists
         PC = PC + 1
! skipping ghost particles
         IF(PEA(L,4)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         M = PIJK(L,5)

         WTP = ONE
         IF(MPPIC) WTP = DES_STAT_WT(L)
! adding particle volume to ongoing summation of solids volume
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + PVOL(L)*WTP
! adding particle velocity to ongoing summation of solids velocity
         DES_U_S(IJK,M) = DES_U_S(IJK,M) + PVOL(L)*DES_VEL_NEW(1,L)*WTP
         DES_V_S(IJK,M) = DES_V_S(IJK,M) + PVOL(L)*DES_VEL_NEW(2,L)*WTP
         IF(DO_K) DES_W_S(IJK,M) = DES_W_S(IJK,M) + &
            PVOL(L)*DES_VEL_NEW(3,L)*WTP
      ENDDO      ! end loop over L = 1,particles


! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,i,j,k,cm,m,sum_eps,ep_sm,                    &
!$omp         osolvol,ovol)
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF (.NOT.DES_CONTINUUM_HYBRID) THEN
            EP_G(IJK) = ONE
         ELSE
! summing together total continuum solids volume
            SUM_EPS = ZERO
            DO CM = 1,SMAX
               SUM_EPS = SUM_EPS + EP_S(IJK,CM)
            ENDDO
            EP_G(IJK) = ONE - SUM_EPS
         ENDIF  ! end if/else (.not.des_continuum_hybrid)


! calculating the cell average solids velocity for each solids phase
          DO M = 1, DES_MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DO_K) DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            IF(VOL(IJK).GT.0) THEN
               OVOL = ONE/(VOL(IJK))
               DES_ROP_S(IJK,M) = DES_RO_S(M)*SOLVOLINC(IJK,M)*OVOL
            ENDIF

! calculating void fraction in fluid cell based on value of bulk density
! calculated above
            IF(DES_ROP_S(IJK,M) >= ZERO) THEN
! calculating solids volume fraction based on bulk density
               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               !IF(.NOT.DES_ONEWAY_COUPLED)
               EP_G(IJK) = EP_G(IJK) - EP_SM
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

! ep_g does not matter if granular flow simulation (i.e. no fluid)
               IF(EP_G(IJK)<ZERO .AND.DES_CONTINUUM_COUPLED) THEN
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1000)
                     WRITE(UNIT_LOG,1004) IJK, I_OF(IJK), J_OF(IJK), &
                        EP_SM, PINC(IJK)
                     WRITE(UNIT_LOG,1001)
                  ENDIF

                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDDO   ! end loop over M=1,DES_MMAX

      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do


 1000 FORMAT(3X,'---------- FROM COMPUTE_MEAN_FIELDS_ZERO_ORDER ',&
             '---------->')
 1001 FORMAT(3X,'<--------- END COMPUTE_MEAN_FIELDS_ZERO_ORDER ',&
             '----------')

 1004 FORMAT(5X,'WARNING: EP_G < 0 at IJK=', I10,' I=', I10, &
         ' J=', I10,/5X,'EP_S=', ES15.9, ' & PINC (number of ',&
         'particles in cell)= ',I10)


      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_INTERP

!-----------------------------------------------
! Modules
!-----------------------------------------------
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
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK, &
                 IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, &
                 II, JJ, KK
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: IDIM, IJK2
      INTEGER :: ICUR, JCUR, KCUR, CUR_IJK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to
! set_interpolation_stencil
      INTEGER :: ONEW
! constant whose value depends on dimension of system
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
! index of solid phase that particle NP belongs to
      INTEGER :: M
! particle number index, used for looping
      INTEGER :: NP, NINDX
! index to track accounted for particles
      INTEGER :: PC
! Statistical weight of the particle. Equal to one for DEM
      DOUBLE PRECISION :: WTP
! one over the solids volume fraction
      DOUBLE PRECISION :: OEPS

      DOUBLE PRECISION :: VOL_SURR
      DOUBLE PRECISION :: MASS_SOL1, MASS_SOL2
! sum of mass_sol1 and mass_sol2 across all processors
      DOUBLE PRECISION :: MASS_SOL1_ALL, MASS_SOL2_ALL

      DOUBLE PRECISION :: TEMP1, TEMP2

! for error messages
      INTEGER :: IER

      DOUBLE PRECISION :: JUNK_VAL(3)

      INTEGER :: COUNT_NODES_OUTSIDE, COUNT_NODES_INSIDE, &
                 COUNT_NODES_INSIDE_MAX, COUNT_TEMP
      double precision :: RESID_ROPS(DES_MMAX), &
                          RESID_VEL(3, DES_MMAX)
      double precision :: NORM_FACTOR
!Handan Liu added on Jan 17 2013
          DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
          DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
          DOUBLE PRECISION :: desposnew(3)
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

! initializing
      MASS_SOL1 = ZERO
      MASS_SOL2 = ZERO
      MASS_SOL1_ALL = ZERO
      MASS_SOL2_ALL = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

! cartesian_grid related quantities
      COUNT_NODES_INSIDE_MAX = merge(4, 8, NO_K)


! Initialize entire arrays to zero
      DES_VEL_NODE = ZERO
      DES_ROPS_NODE = ZERO
      DES_ROP_S = zero
      DES_U_S = ZERO
      DES_V_S = ZERO
      IF(DO_K) DES_W_S = ZERO


! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      CALL SET_INTERPOLATION_SCHEME(2)

!Handan Liu added on Jan 17 2013; again on June 2013
!$omp   parallel default(shared)                                        &
!$omp   private(IJK,I,J,K,PCELL,COUNT_NODES_INSIDE,II,JJ,KK,IW, &
!$omp           IE,JS,JN,KB,KTP,ONEW,CUR_IJK,IPJK,IJPK,IPJPK,IJKP,      &
!$omp           IJPKP,IPJKP,IPJPKP,gst_tmp,vst_tmp,nindx,np,wtp,m,      &
!$omp           JUNK_VAL,desposnew,weight_ft,icur,jcur,kcur,            &
!$omp           I1, I2, J1, J2, K1, K2, IDIM,IJK2,NORM_FACTOR,          &
!$omp           RESID_ROPS,RESID_VEL,COUNT_NODES_OUTSIDE, TEMP1)
!$omp do reduction(+:MASS_SOL1) reduction(+:DES_ROPS_NODE,DES_VEL_NODE)
      !IJKLOOP: DO IJK = IJKSTART3,IJKEND3      ! Removed by Handan Liu
      DO IJK = IJKSTART3,IJKEND3

! Cycle this cell if not in the fluid domain or if it contains no
! particle/parcel
         IF(.NOT.FLUID_AT(IJK) .OR. PINC(IJK).EQ.0) CYCLE

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(.NOT.IS_ON_myPE_owns(I,J,K)) THEN
! PINC array count reflects only the acutal particles and does not
! include the ghost particles. Therefore, PINC can only be non-zero
! for scalar cells that belong to this processor. So this is just
! a sanity check that will ensure an error is flagged if the logic
! is broken elsewhere in the code.
            IF(DMP_LOG) THEN
               WRITE(UNIT_LOG,*) &
                  'Critical Error in compute_mean_fields_interp:', &
                  'This cell does not belong to this proc ', myPE
               WRITE(UNIT_LOG,'(A,10(2x,i5))') &
                  'PINC in I, J, K, IJP, NP = ', &
                  I_OF(IJK), J_OF(IJK), K_OF(IJK), PINC(IJK), &
                  PINC(IM_OF(IJK)), PINC(IP_OF(IJK))
            ENDIF
            WRITE(*,*) 'Critical Error in compute_mean_fields_interp:', &
               'This cell does not belong to this proc ', myPE, &
               'Exiting the code, check the log file for details'
            CALL MFIX_EXIT(myPE)
         ENDIF

         PCELL(1) = I-1
         PCELL(2) = J-1
         PCELL(3) = merge(1, K-1, NO_K)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         CALL SET_INTERPOLATION_STENCIL(PCELL,IW,IE,JS,JN,KB,&
            KTP,INTERP_SCHEME,DIMN,ORDERNEW = ONEW)

         COUNT_NODES_OUTSIDE = 0
! Computing/setting the geometric stencil
         DO K = 1,merge(1, ONEW, NO_K)
            DO J = 1,ONEW
               DO I = 1,ONEW
                  II = IW + I-1
                  JJ = JS + J-1
                  KK = KB + K-1
                  CUR_IJK = funijk(IMAP_C(II),JMAP_C(JJ),KMAP_C(KK))
                  IPJK    = funijk(IMAP_C(II+1),JMAP_C(JJ),KMAP_C(KK))
                  IJPK    = funijk(IMAP_C(II),JMAP_C(JJ+1),KMAP_C(KK))
                  IPJPK   = funijk(IMAP_C(II+1),JMAP_C(JJ+1),KMAP_C(KK))
                  IF(DO_K) THEN
                     IJKP    = funijk(IMAP_C(II),JMAP_C(JJ),KMAP_C(KK+1))
                     IJPKP   = funijk(IMAP_C(II),JMAP_C(JJ+1),KMAP_C(KK+1))
                     IPJKP   = funijk(IMAP_C(II+1),JMAP_C(JJ),KMAP_C(KK+1))
                     IPJPKP  = funijk(IMAP_C(II+1),JMAP_C(JJ+1),KMAP_C(KK+1))
                  ENDIF

                  GST_TMP(I,J,K,1) = XE(II)
                  GST_TMP(I,J,K,2) = YN(JJ)
                  GST_TMP(I,J,K,3) = merge(DZ(1), ZT(KK), NO_K)
                  VST_TMP(I,J,K,:) = ZERO
!===================================================================>>> Handan Liu

                  IF(CARTESIAN_GRID) THEN
                     IF(SCALAR_NODE_ATWALL(CUR_IJK)) COUNT_NODES_OUTSIDE = &
                     & COUNT_NODES_OUTSIDE + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO


! Calculate des_rops_node so des_rop_s, and in turn, ep_g can be updated
!----------------------------------------------------------------->>>

! looping through particles in the cell
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)

! conducting some sanity checks
            IF(PEA(NP, 4)) THEN
               IF(DMP_LOG) THEN
                  WRITE(UNIT_LOG,*) 'Encountered a ghost ',&
                     'particle in compute_mean_fields_interp'
                  WRITE(UNIT_LOG,'(A,10(2x,i5))') &
                     'PINC in I, J, K, IJP, NP = ', I_OF(IJK), &
                     J_OF(IJK), K_OF(IJK), PINC(IJK), PINC(IM_OF(IJK)),&
                     PINC(IP_OF(IJK))
               ENDIF
               WRITE(*,*) 'Encountered a ghost particle ',&
                  'in compute_mean_fields_interp'
               WRITE(*,'(A,10(2x,i5))') &
                  'PINC in I, J, K, IJP, NP = ', I_OF(IJK), J_OF(IJK), &
                  K_OF(IJK), PINC(IJK), PINC(IM_OF(IJK)), &
                  PINC(IP_OF(IJK))
               CALL MFIX_EXIT(myPE)
            ENDIF

            desposnew(:) = des_pos_new(:,np)
            call DRAG_INTERPOLATION(gst_tmp,vst_tmp,desposnew,JUNK_VAL,weight_ft)
!===================================================================>>> Handan Liu

            M = PIJK(NP,5)
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP

            TEMP2 = DES_RO_S(M)*PVOL(NP)*WTP

            DO K = 1, merge(1, ONEW, NO_K)
               DO J = 1, ONEW
                  DO I = 1, ONEW
! shift loop index to new variables for manipulation
                     II = IW + I-1
                     JJ = JS + J-1
                     KK = KB + K-1

! The interpolation is done using node. so one should use consistent
! numbering system in the current version imap_c is used instead of
! ip_of or im_of
                     ICUR = IMAP_C(II)
                     JCUR = JMAP_C(JJ)
                     KCUR = KMAP_C(KK)
                     CUR_IJK = funijk(ICUR, JCUR, KCUR)

                     !TEMP1 = WEIGHTP(I,J,K)*DES_RO_S(M)*PVOL(NP)*WTP
! Changed TEMP1 as an array TEMP1(NP) to ensure different TEMP1
! for each particle in an ijk cell <June 18 2013>
                     TEMP1 = WEIGHT_FT(I,J,K)*TEMP2
                     DES_ROPS_NODE(CUR_IJK,M) = DES_ROPS_NODE(CUR_IJK,M) + TEMP1
                     DES_VEL_NODE(CUR_IJK,:,M) = &
                        DES_VEL_NODE(CUR_IJK,:,M) + TEMP1*DES_VEL_NEW(:,NP)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx=1,pinc(ijk))
!-----------------------------------------------------------------<<<


! Only for cutcell cases may count_nodes_inside become less than its
! original set value. In such an event, the contribution of scalar nodes
! that do not reside in the domain is added to a residual array. This
! array is then redistribited equally to the nodes that are in the fluid
! domain. These steps are done to conserve mass.
!----------------------------------------------------------------->>>
         IF (CARTESIAN_GRID) THEN

! only for cartesian_grid will count_nodes_outside be modified from zero
            COUNT_NODES_INSIDE = &
               COUNT_NODES_INSIDE_MAX - COUNT_NODES_OUTSIDE

            IF(COUNT_NODES_INSIDE.LT.COUNT_NODES_INSIDE_MAX) THEN

! initializing
               RESID_ROPS(1:DES_MMAX) = ZERO
               RESID_VEL(:, 1:DES_MMAX) = ZERO

! Convention used to number node numbers
! i=1, j=2           i=2, j=2
!   _____________________
!   |                   |
!   |  I = 2, J = 2     |
!   |___________________|
! i=1, j=1           i=2, j=1
! setting indices based on convention
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               I1 = I-1
               I2 = I
               J1 = J-1
               J2 = J
               K1 = merge(K, K-1, NO_K)
               K2 = K
! first calculate the residual des_rops_node and des_vel_node that was
! computed on nodes that do not belong to the domain

               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK)

                        IF(SCALAR_NODE_ATWALL(IJK2)) THEN
                           RESID_ROPS(1:DES_MMAX) = &
                              RESID_ROPS(1:DES_MMAX) +&
                              DES_ROPS_NODE(IJK2,1:DES_MMAX)
                           DES_ROPS_NODE(IJK2,1:DES_MMAX) = ZERO
                           DO IDIM = 1, merge(2,3,NO_K)
                              RESID_VEL(IDIM, 1:DES_MMAX) = &
                                 RESID_VEL(IDIM, 1:DES_MMAX) + &
                                 DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX)
                              DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) = ZERO
                           ENDDO
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO

! now add this residual equally to the remaining nodes
               NORM_FACTOR = ONE/REAL(COUNT_NODES_INSIDE)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK)

                        IF(.NOT.SCALAR_NODE_ATWALL(IJK2)) THEN
                           DES_ROPS_NODE(IJK2,1:DES_MMAX) = &
                              DES_ROPS_NODE(IJK2,1:DES_MMAX) + &
                              RESID_ROPS(1:DES_MMAX)*NORM_FACTOR
                           DO IDIM = 1, merge(2,3,NO_K)
                              DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) = &
                                 DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) + &
                                 RESID_VEL(IDIM, 1:DES_MMAX)*NORM_FACTOR
                           ENDDO
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF   ! end if (cartesian_grid)
!-----------------------------------------------------------------<<<
      !ENDDO IJKLOOP   ! end do ijkloop (ijk=ijkstart3,ijkend3)
      ENDDO
!$omp end parallel

! At the interface des_rops_node has to be added since particles
! across the processors will contribute to the same scalar node.
! sendrecv will be called and the node values will be added
! at the junction. des_rops_node is altered by the routine when
! periodic boundaries are invoked
      CALL DES_ADDNODEVALUES_MEAN_FIELDS


! Now go from node to scalar center. Same convention as sketched
! earlier
!----------------------------------------------------------------->>>
! Explanation by RG: 08/17/2012
! the approach used here is to make it general enough for cutcells to be
! included as well. The new changes do not alter earlier calculations
! but make the technique general as to include cartesian grid (cut-cell)
! simulations.
! Previously, the volume of the node (by array des_vol_node) was used to
! first scale the nodal values. Subsequently, these nodal values were
! equally added to compute the cell centered values for the scalar cell.

! Consider an internal node next to an edge node (a node adjacent to a
! boundary). In 2D, the volume of an edge node will be half that of an
! internal node. And an edge node will contribute double compared to
! an internal node to the value of the scalar cell they share. These
! calculations were previously accomplished via the variable volume of
! node.  Now this is accomplished by the ratio vol(ijk2)/vol_sur, where
! vol(ijk2) is the volume of the scalar cell in consideration and
! vol_sur is the sum of all the scalar cell volumes that have this node
! as the common node.

! looping over all fluid cells
!Handan Liu added here on Feb. 28 2013
!$omp   parallel do default(shared)             &
!$omp   private(K,J,I,IJK,I1,I2,J1,J2,K1,K2,    &
!$omp                   II,JJ,KK,IJK2,M,VOL_SURR)   collapse (3)
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = funijk(I,J,K)
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               K1 = K
               K2 = merge(K, K+1, NO_K)

               VOL_SURR = ZERO

! looping over stencil points (node values)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells
                        IJK2 = funijk(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK))
                        IF(FLUID_AT(IJK2)) VOL_SURR = VOL_SURR+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO

! looping over stencil points (NODE VALUES)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells

                        IJK2 = funijk(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK))
                        IF(FLUID_AT(IJK2).and.(IS_ON_myPE_wobnd(II, JJ, KK))) THEN
! Since the data in the ghost cells is spurious anyway and overwritten during
! subsequent send receives, do not compute any value here as this will
! mess up the total mass value that is computed below to ensure mass conservation
! between Lagrangian and continuum representations
                           DO M = 1, DES_MMAX
                              DES_ROP_S(IJK2, M) = DES_ROP_S(IJK2, M) + &
                                 DES_ROPS_NODE(IJK,M)*VOL(IJK2)/VOL_SURR
                              DES_U_S(IJK2, M) = DES_U_S(IJK2, M) + &
                                 DES_VEL_NODE(IJK, 1, M)*VOL(IJK2)/VOL_SURR
                              DES_V_S(IJK2, M) = DES_V_S(IJK2, M) + &
                                 DES_VEL_NODE(IJK, 2, M)*VOL(IJK2)/VOL_SURR
                              IF(DO_K) DES_W_S(IJK2, M) = DES_W_S(IJK2, M) + &
                                 DES_VEL_NODE(IJK, 3, M)*VOL(IJK2)/VOL_SURR
                           ENDDO
                        ENDIF
                     ENDDO  ! end do (ii=i1,i2)
                  ENDDO  ! end do (jj=j1,j2)
               ENDDO  ! end do (kk=k1,k2)

            ENDDO   ! end do (i=istart2,iend1)
         ENDDO   ! end do (j=jstart2,jend1)
      ENDDO   ! end do (k=kstart2,kend1)

!-----------------------------------------------------------------<<<

!Handan Liu added here on Feb. 28 2013
!$omp   parallel do default(shared)             &
!$omp   private(K,J,I,IJK,M) reduction(+:MASS_SOL2)
      DO IJK = IJKSTART3, IJKEND3
         IF(.not.FLUID_AT(IJK)) cycle
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         DO M = 1, DES_MMAX
            IF(DES_ROP_S(IJK, M).GT.ZERO) THEN
               DES_U_S(IJK, M) = DES_U_S(IJK,M)/DES_ROP_S(IJK, M)
               DES_V_S(IJK, M) = DES_V_S(IJK,M)/DES_ROP_S(IJK, M)
               IF(DO_K) DES_W_S(IJK, M) = DES_W_S(IJK,M)/DES_ROP_S(IJK, M)

! Finally divide by scalar cell volume to obtain \eps * \rho_s
               DES_ROP_S(IJK, M) = DES_ROP_S(IJK, M)/VOL(IJK)

! Note that for cut-cell, it is important to check both fluid_at and
! is_on_mype_wobnd.  Fluid_at is not enough as it is shared between procs
! and fluid_at(ijkend3) might be true when in fact it does not belong to
! that proc
               IF(IS_ON_myPE_wobnd(I, J, K)) MASS_SOL2 = MASS_SOL2 + &
               DES_ROP_S(IJK, M)*VOL(IJK)

               !WRITE(*,*) 'ROP_S = ', DES_ROP_S(IJK, M), DES_ROP_S(IJK, M)/DES_RO_S(M)

            ENDIF
         ENDDO   ! end loop over M=1,DES_MMAX
      ENDDO                     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do

! RG: Aug, 20, 2012.
! -------------------------------------------------------------------
! Decoupling the computation of ep_g and rop_g from the computation
! of des_rop_s. In the current implementation, des_rop_s (or rop_s)
! and epg and also rop_g are first calculated and then sent
! received. However, ep_g (and also rop_g) are simply functions of
! des_rop_s. Therefore, once des_rop_s is sent received, ep_g and
! rop_g can be correctly computed even for the ghost cells; thus, not
! requiring additional send_recv calls for ep_g and rop_g.

! DES_ROP_S is not exchanged across the interface.
! ROP_S and EP_G are done at the end of the DEM time step in
! des_time_march.f.
! -------------------------------------------------------------------

      IF (.NOT.MPPIC) THEN
! According to the current DEM implementation, the mean fields are
! communicated (i.e., sent received) following the inner DEM
! time steps. Since the number of these inner time steps could be large
! in DEM and the cell centered DES mean fields (such as des_rop_s,
! des_u_s, etc.) in the ghost cells (across processor boundaries)
! are not used during these substeps, the sent receive of these mean
! fields at the end of DEM substeps does not result in any error.

! It should be noted above that the ep_g values in the ghost cells
! are not correct, and if any future development requires the correct
! value of ep_g in ghost cells during DEM substeps, a send_recv
! of des_rop_s should be done before calling the following routine.
! This will impose a substantial penalty on computational time since
! communication overheads will go up.
         CALL COMP_EPG_ROP_G_FROM_ROP_S


      ELSE

! share the des_rop_s array across the processors
         CALL SEND_RECV(DES_ROP_S,2)
! compute the arrays epg and rop_g. This negates the need for further
! communication of these arrays.
         CALL COMP_EPG_ROP_G_FROM_ROP_S

! Now calculate Eulerian mean velocity fields like U_S, V_S, and W_S.
         CALL SEND_RECV(DES_U_S,2)
         CALL SEND_RECV(DES_V_S,2)
         IF(DO_K) CALL SEND_RECV(DES_W_S,2)

! The Eulerian velocity field is used to set up the stencil to interpolate
! mean solid velocity at the parcel's location. DES_U_S could have also been
! used, but that also would have require the communication at this stage.
! The final interpolated value does not change if the stencil is formed by
! first obtaining face centered Eulerian velocities (U_S, etc.)
! and then computing the node velocities from them or directly computing
! the node velocities from cell centered average velocity field (DES_U_S,
! etc.). We are using the first approach as it is more natural to set
! BC's on solid velocity field in the face centered represenation (U_S,
! etc.)

         IF(.NOT.CARTESIAN_GRID) THEN
            CALL MPPIC_COMP_EULERIAN_VELS_NON_CG
         ELSE
            CALL MPPIC_COMP_EULERIAN_VELS_CG
         ENDIF
      ENDIF   ! end if (.not.mppic)

! turn on the below statements to check if the mass is conserved
! between discrete and continuum representations. Should be turned to
! false for any production runs.
      IF(DES_REPORT_MASS_INTERP) THEN
         CALL GLOBAL_SUM(MASS_SOL1, MASS_SOL1_ALL)
         CALL GLOBAL_SUM(MASS_SOL2, MASS_SOL2_ALL)
         if(myPE.eq.pe_IO) THEN
            WRITE(*,'(/,5x,A,4(2x,g17.8),/)') &
                 'SOLIDS MASS DISCRETE AND CONTINUUM =  ', &
                 MASS_SOL1_ALL, MASS_SOL2_ALL
         ENDIF
      ENDIF
      END SUBROUTINE COMP_MEAN_FIELDS_INTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE COMP_EPG_ROP_G_FROM_ROP_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE constant
      USE geometry
      USE funits
      USE indices
      USE compar
      USE physprop
      USE fldvar
      USE discretelement
      USE cutcell
      USE mfix_pic
      use desmpi_wrapper
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: IJK
! solids phase indices
      INTEGER :: CM, M
! solids volume fraction of mth solids phase
      DOUBLE PRECISION EP_SM
! total solids volume fraction of continuum solids phases
      DOUBLE PRECISION SUM_EPS
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,cm,m,sum_eps,ep_sm)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF (.NOT.DES_CONTINUUM_HYBRID) THEN
            EP_G(IJK) = ONE
         ELSE
! summing together total continuum solids volume
            SUM_EPS = ZERO
            DO CM = 1,SMAX
               SUM_EPS = SUM_EPS + EP_S(IJK,CM)
            ENDDO
            EP_G(IJK) = ONE - SUM_EPS
         ENDIF  ! end if/else (.not.des_continuum_hybrid)

         DO M = 1, DES_MMAX
! calculating void fraction in fluid cell based on value of bulk density
            IF(DES_ROP_S(IJK, M).GT.ZERO) THEN
               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               !IF(.NOT.DES_ONEWAY_COUPLED)
               EP_G(IJK) = EP_G(IJK) - EP_SM
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

! ep_g does not matter if granular flow simulation (i.e. no fluid)
              IF(EP_G(IJK)<ZERO .AND.DES_CONTINUUM_COUPLED) THEN
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1004) EP_G(IJK), IJK, I_OF(IJK), &
                        J_OF(IJK), K_OF(IJK), EP_SM, PINC(IJK), &
                        CUT_CELL_AT(IJK)
                  ENDIF
                  WRITE(*,1004) EP_G(IJK), IJK, &
                     I_OF(IJK), J_OF(IJK), K_OF(IJK), EP_SM, &
                     PINC(IJK), CUT_CELL_AT(IJK)

                  IF(CARTESIAN_GRID) THEN
                     CALL WRITE_DES_DATA
!                     CALL WRITE_VTU_FILE
                  ENDIF

!                 CALL MFIX_EXIT(myPE)
                  call des_mpi_stop
               ENDIF
            ENDIF
         ENDDO                  ! end loop over M=1,DES_MMAX
      ENDDO                     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do
 1004 FORMAT(/5X, 'Message from comp_epg_rop_g_from_rop_s', &
      /,5X,'Warning, EP_G = ', g17.8, 2x, 'LT Zero at IJK', I20, &
      /,5X,'I,J,K = ', I10, 2X, I10, 2x, I10, &
      /,5X,'EP_S=', ES15.9, &
      /,5X,'PINC (number of particles in cell)= ',I10, &
      /,5X,'Cut cell ? ', L2,/)


      END SUBROUTINE COMP_EPG_ROP_G_FROM_ROP_S

