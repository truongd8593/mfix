!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GAS_DRAG                                                C
!  Purpose: Account for the equal and opposite drag force on the gas   C
!           phase due to particles by introducing the drag as a        C
!           source term.  Face centered.                               C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GAS_DRAG(A_M, B_M, VXF_GS, VELDIR, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Volume times drag coefficient at cell face
      DOUBLE PRECISION, INTENT(IN) :: VXF_GS(DIMENSION_3, DIMENSION_M)
! flag/index used to specify which gas phase momentum equation
! (x=1,y=2,z=3) is currently being solved
      INTEGER, INTENT(IN) :: VELDIR
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
! Global indices
      INTEGER :: IMJK, IJMK, IJKM, &
                 IMJMK, IMJKM, IJMKM, IMJMKM
! Solids phase index
      INTEGER :: M
! Face center values of u_sm (i+1/2), v_sm (j+1/2) and w_sm (k+1/2)
      DOUBLE PRECISION :: USFCM, VSFCM, WSFCM
! Volume of x, y, or z cell on staggered grid
      DOUBLE PRECISION :: VCELL
! temporary variables for matrix A_M and vector B_M that are
! used for local calculations
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
! (=0.25 in 3D and =0.5 in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!-----------------------------------------------

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
!!$      omp_start=omp_get_wtime()

      IF(DES_ONEWAY_COUPLED) THEN !do nothing
         RETURN
      ENDIF

! (=0.25 in 3D and =0.5 in 2D)
      AVG_FACTOR = merge(0.5d0, 0.25D0, NO_K)

      IF(VELDIR .EQ.1) THEN

!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,usfcm) schedule (guided,50)
         DO IJK = IJKSTART3, IJKEND3
            DO M = 1, DES_MMAX
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))

                  IF(DES_INTERP_ON) THEN
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(IJK,M) + &
                                            DRAG_AM(IJMK,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(IJK,1,M) + &
                                            DRAG_BM(IJMK,1,M))
                     IF(DO_K) THEN
                        IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                        IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                        IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                        IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))
                        tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(IJKM,M) + &
                                                    DRAG_AM(IJMKM,M))
                        tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(IJKM,1,M) +&
                                                    DRAG_BM(IJMKM,1,M))
                     ENDIF
! this volume comes from integration of the gas phase momentum equation
! over the corresponding control volume (in this case the U CV)
                     VCELL = VOL_U(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                  ELSE   ! non-interpolated
                     USFCM = AVG_X(DES_U_S(IJK,M),DES_U_S(EAST_OF(IJK),M),I)
                     IF (DES_CONTINUUM_HYBRID) THEN
                        tmp_A =  - VXF_GDS(IJK,M)
                        tmp_B =  - VXF_GDS(IJK,M)*USFCM
                     ELSE
                        tmp_A =  - VXF_GS(IJK,M)
                        tmp_B =  - VXF_GS(IJK,M)*USFCM
                     ENDIF
                  ENDIF   ! end if/else(des_interp_on)

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B

               ENDIF   ! end if (fluid_at(ijk))
            ENDDO   ! end do loop (m=1,des_mmax)
         ENDDO   !   end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do

      ELSEIF(VELDIR.EQ.2) THEN

!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,vsfcm) schedule (guided,50)
         DO IJK = IJKSTART3, IJKEND3
            DO M = 1, DES_MMAX
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))

                  IF(DES_INTERP_ON) THEN
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(IJK,M) + &
                                            DRAG_AM(IMJK,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(IJK,2,M) + &
                                            DRAG_BM(IMJK,2,M))
                     IF(DO_K) THEN
                        IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                        IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                        IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                        IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))
                        tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(IJKM,M) +&
                                                    DRAG_AM(IMJKM,M))
                        tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(IJKM,2,M) +&
                                                    DRAG_BM(IMJKM,2,M))
                     ENDIF
                     VCELL = VOL_V(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                  ELSE
                    VSFCM = AVG_Y(DES_V_S(IJK,M),DES_V_S(NORTH_OF(IJK),M),J)
                    IF (DES_CONTINUUM_HYBRID) THEN
                        tmp_A =  - VXF_GDS(IJK,M)
                        tmp_B =  - VXF_GDS(IJK,M)*VSFCM
                    ELSE
                       tmp_A =  - VXF_GS(IJK,M)
                       tmp_B =  - VXF_GS(IJK,M)*VSFCM
                    ENDIF
                  ENDIF   ! end if/else (des_interp_on)

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B

               ENDIF   ! end if (fluid_at(ijk))
            ENDDO   ! end do loop (m=1,des_mmax)
         ENDDO   !   end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do

      ELSEIF(VELDIR.EQ.3) THEN

!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,wsfcm) schedule (guided,50)
         DO IJK = IJKSTART3, IJKEND3
            DO M = 1, DES_MMAX
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))
                  IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                  IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                  IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                  IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))

                  IF(DES_INTERP_ON) THEN
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(IJK,M) + &
                                            DRAG_AM(IMJK,M) + &
                                            DRAG_AM(IJMK,M) +&
                                            DRAG_AM(IMJMK,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(IJK,3,M) + &
                                            DRAG_BM(IMJK,3,M) + &
                                            DRAG_BM(IJMK,3,M) + &
                                            DRAG_BM(IMJMK,3,M))
                     VCELL = VOL_W(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                  ELSE
                     WSFCM = AVG_Z(DES_W_S(IJK,M),DES_W_S(TOP_OF(IJK),M),K)
                     IF (DES_CONTINUUM_HYBRID) THEN
                        tmp_A =  - VXF_GDS(IJK,M)
                        tmp_B =  - VXF_GDS(IJK,M)*WSFCM
                     ELSE
                        tmp_A = - VXF_GS(IJK,M)
                        tmp_B = - VXF_GS(IJK,M)*WSFCM
                     ENDIF
                  ENDIF

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B

               ENDIF   ! end if (fluid_at(ijk))
            ENDDO   ! end do loop (m=1,des_mmax)
         ENDDO   !   end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel do

      ENDIF   ! end if/else (veldir=1,2,3)
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'gas_drag:',omp_end - omp_start

      RETURN
      END SUBROUTINE GAS_DRAG
