!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)            C
!>
!!  Purpose: DES - Accounting for the equal and opposite drag force   
!!           on gas due to particles by introducing the drag           
!!           as a source term. Face centered                            
!<
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)
!-----------------------------------------------
!     M o d u l e s 
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
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     
!     Error index
      INTEGER          IER
     
!     Indices
      INTEGER          IJK, I, J, K, NC
     
!     Phase index
      INTEGER          M, UV, VV, WV

!     Averaging Factor 
      DOUBLE PRECISION :: AVG_FACTOR, COUPLING_FAC

      DOUBLE PRECISION USFCM, VSFCM, WSFCM, VCELL
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
      DOUBLE PRECISION B_M(DIMENSION_3, 0:DIMENSION_M)
      DOUBLE PRECISION VXF_GS(DIMENSION_3, DIMENSION_M)  
      DOUBLE PRECISION tmp_A, tmp_B
! Pradeep global indices 
      integer IMJK,IJMK,IJKM,IMJMK,IMJKM,IJMKM,IMJMKM 
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime	      
!       by Tingwen      
!!$      omp_start=omp_get_wtime()
      
      IF(DES_ONEWAY_COUPLED) THEN !do nothing 
         RETURN 
      ENDIF

      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)

      IF(UV.EQ.1) THEN
!!         DO M = 1, MMAX
!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,usfcm) schedule (guided,50)       
            DO IJK = IJKSTART3, IJKEND3
         DO M = 1, MMAX    
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)
                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE
!pradeep changing drag_am from three to one dimensional
                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))

                  IF(DES_INTERP_ON) then 
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(IJK,M) + DRAG_AM(IJMK,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(IJK,1,M) + DRAG_BM(IJMK,1,M))
                         
                     IF(DIMN.eq.3) THEN 
                        IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                        IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                        IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                        IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))
                        tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(IJKM,M) + DRAG_AM(IJMKM,M))
                        tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(IJKM,1,M) + DRAG_BM(IJMKM,1,M))
                     ENDIF
                     
                     VCELL = VOL(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                     
                  ELSE 
                     USFCM = AVG_X(DES_U_S(IJK,M),DES_U_S(EAST_OF(IJK),M),I)
                     tmp_A =  - VXF_GS(IJK,M)
                     tmp_B =  - VXF_GS(IJK,M)*USFCM
                  ENDIF

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B

               ENDIF
            ENDDO
         ENDDO
!$omp end parallel do 
	
      ELSEIF(VV.EQ.1) THEN
!!         DO M = 1, MMAX
!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,vsfcm) schedule (guided,50)   
            DO IJK = IJKSTART3, IJKEND3
         DO M = 1, MMAX    
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)
                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE
!pradeep changing drag_am from three to one dimensional
                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))
                  
                  IF(DES_INTERP_ON) THEN                      
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(IJK,M) + DRAG_AM(IMJK,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(IJK,2,M) + DRAG_BM(IMJK,2,M))
                 
                     IF(DIMN.eq.3) THEN
                        IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                        IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                        IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                        IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))
                        tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(IJKM,M) + DRAG_AM(IMJKM,M))
                        tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(IJKM,2,M) + DRAG_BM(IMJKM,2,M))
                     ENDIF

                     VCELL = VOL(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                 
                  ELSE 
                    VSFCM = AVG_Y(DES_V_S(IJK,M),DES_V_S(NORTH_OF(IJK),M),J)
                    tmp_A =  - VXF_GS(IJK,M)
                    tmp_B =  - VXF_GS(IJK,M)*VSFCM
                  ENDIF

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B
                 
               ENDIF
             ENDDO
          ENDDO
!$omp end parallel do 
          
      ELSE IF(WV.EQ.1) THEN
!!         DO M = 1, MMAX
!$omp parallel do default(shared)                               &
!$omp private(ijk,m,i,j,k,imjk,ijmk,imjmk,tmp_a,tmp_b,          &
!$omp         ijkm,imjkm,ijmkm,imjmkm,vcell,wsfcm) schedule (guided,50)   
            DO IJK = IJKSTART3, IJKEND3
         DO M = 1, MMAX    
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)
                  IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
                  IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
                  IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE
!pradeep changing drag_am from three to one dimensional
                  IMJK = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K))
                  IJMK = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K))
                  IMJMK = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K))
                  IJKM = FUNIJK(IMAP_C(I),JMAP_C(J),KMAP_C(K-1))
                  IMJKM = FUNIJK(IMAP_C(I-1),JMAP_C(J),KMAP_C(K-1))
                  IJMKM = FUNIJK(IMAP_C(I),JMAP_C(J-1),KMAP_C(K-1))
                  IMJMKM = FUNIJK(IMAP_C(I-1),JMAP_C(J-1),KMAP_C(K-1))
                  
                  IF(DES_INTERP_ON) THEN
                     tmp_A =  - 0.25d0*(DRAG_AM(IJK,M) + &
                               DRAG_AM(IMJK,M)+DRAG_AM(IJMK,M)+DRAG_AM(IMJMK,M))
                     tmp_B =  - 0.25d0*(DRAG_BM(IJK,3,M) + &
                              DRAG_BM(IMJK,3,M)+ DRAG_BM(IJMK,3,M)+ DRAG_BM(IMJMK,3,M))
                     
                     VCELL = VOL(IJK)
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL

                  ELSE 
                     WSFCM = AVG_Z(DES_W_S(IJK,M),DES_W_S(TOP_OF(IJK),M),K)
                     tmp_A = - VXF_GS(IJK,M)
                     tmp_B = - VXF_GS(IJK,M)*WSFCM
                     
                  ENDIF

                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B
                  
               ENDIF
            ENDDO
         ENDDO
!$omp end parallel do 
	
      ENDIF
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'gas_drag:',omp_end - omp_start
!!       get worse!!!!!
      RETURN
      END SUBROUTINE GAS_DRAG
