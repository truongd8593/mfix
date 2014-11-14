!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_gs_Z                                                 C
!  Purpose: Calculate the average drag coefficient at i, j, k+1/2 and  C
!           multiply with W-momentum cell volume.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE VF_GS_Z(VXF_GS, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE drag
      USE discretelement
      USE fun_avg
      USE functions

      use run, only: SOLIDS_MODEL

      use run, only: DEM_SOLIDS, PIC_SOLIDS
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Volume x Drag
      DOUBLE PRECISION, INTENT(INOUT) :: VxF_gs(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: K, IJK, IJKT
! Index of continuum solids phases
      INTEGER :: M
! Index of discrete solids 'phases'
      INTEGER :: DM, MTOT
!-----------------------------------------------

      MTOT = merge(MMAX, MMAX+DES_MMAX, DES_CONTINUUM_HYBRID)

      M_LP: DO M = 1, MTOT
         IF(SOLIDS_MODEL(M) == 'DEM' .AND. DES_ONEWAY_COUPLED) THEN
            VxF_GS(:,M) = ZERO
            CYCLE M_LP
         ENDIF

         DO IJK = ijkstart3, ijkend3
            IF(IP_AT_T(IJK)) THEN
               VXF_GS(IJK,M) = ZERO
            ELSE
               K = K_OF(IJK)
               IJKT = TOP_OF(IJK)
               VXF_GS(IJK,M) = VOL_W(IJK) * &
                  AVG_Z(F_GS(IJK,M),F_GS(IJKT,M),K)
            ENDIF
         ENDDO      ! end do loop (ijk=ijkstart3,ijkend3)
      ENDDO M_LP


      IF (DES_CONTINUUM_HYBRID) THEN
         DO DM = 1, DES_MMAX
            DO IJK = ijkstart3, ijkend3
               IF (IP_AT_T(IJK)) THEN
                  VXF_GDS(IJK,DM) = ZERO
               ELSE
                  K = K_OF(IJK)
                  IJKT = TOP_OF(IJK)
                  VXF_GDS(IJK,DM) = VOL_W(IJK) * &
                     AVG_Z(F_GDS(IJK,DM),F_GDS(IJKT,DM),K)
               ENDIF
            ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)
         ENDDO   ! end do loop (dm=1,des_mmax)
      ENDIF


      RETURN
      END SUBROUTINE VF_GS_Z


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: VF_SS_Z                                                 C
!  Purpose: Calculate the average Solid-Solid drag coefficient at      C
!           i, j, k+1/2 and multiply with W-momentum cell volume.      C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEB-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE VF_SS_Z(VXF_SS, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE drag
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
! Volume x Drag
      DOUBLE PRECISION, INTENT(INOUT) :: VxF_SS(DIMENSION_3, DIMENSION_LM)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: K, IJK, IJKT
! Index of continuum solids phases
      INTEGER :: L, M, LM
! Index of discrete solids phases
      INTEGER :: DM
!-----------------------------------------------

! initialize every call
      VXF_SS(:,:) = ZERO

      DO M = 1, MMAX
         DO L = 1, MMAX
            LM = FUNLM(L,M)
            IF (L .NE. M) THEN
!!$omp  parallel do private(K,IJK,IJKT)
               DO IJK = ijkstart3, ijkend3
                  IF (.NOT.IP_AT_T(IJK)) THEN
                     K = K_OF(IJK)
                     IJKT = TOP_OF(IJK)
                     VXF_SS(IJK,LM) = AVG_Z(F_SS(IJK,LM),F_SS(IJKT,LM),K)*VOL_W(IJK)
                  ELSE     !Impermeable wall
                     VXF_SS(IJK,LM) = ZERO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO   ! end do loop (l=1,mmax)
      ENDDO   ! end do loop (m=1,mmax)

      IF (DES_CONTINUUM_HYBRID) THEN
! initialize every call
         VXF_SDS(:,:,:) = ZERO
         DO M = 1, MMAX
            DO DM = 1, DES_MMAX
!!$omp  parallel do private(K,IJK,IJKT)
               DO IJK = ijkstart3, ijkend3
                  IF (.NOT.IP_AT_T(IJK)) THEN
                     K = K_OF(IJK)
                     IJKT = TOP_OF(IJK)
                     VXF_SDS(IJK,M,DM) = AVG_Z(F_SDS(IJK,M,DM),F_SDS(IJKT,M,DM),K)*VOL_W(IJK)
                  ELSE
                     VXF_SDS(IJK,M,DM) = ZERO
                  ENDIF
               ENDDO      ! end do loop (ijk=ijkstart3,ijkend3)
            ENDDO   ! end do loop (dm=1,des_mmax)
         ENDDO   ! end do loop (m=1,mmax)
      ENDIF   ! end if (discrete_element and des_continuum_hybrid)

      RETURN
      END SUBROUTINE VF_SS_Z


