!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS

      use discretelement, only: DES_INTERP_MEAN_FIELDS

      IF(DES_INTERP_MEAN_FIELDS) THEN
         CALL COMP_MEAN_FIELDS_INTERP1
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
      USE functions
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
         ' J=', I10,/5X,'EP_S=', ES16.9, ' & PINC (number of ',&
         'particles in cell)= ',I10)


      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER
