!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_FLAGS                                              C
!  Purpose: This module assigns a flag to a cell to identify its type. C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: In cylindrical geometry, set the b.c at X=0 to free-slip   C
!  Author: M. Syamlal                                 Date: 09-APR-92  C
!  Revision Number: 2                                                  C
!  Purpose: Initialize FLAG_E, FLAG_N, FLAG_T for IS specifications.   C
!           Change definition of Flags.                                C
!  Author: M. Syamlal                                 Date: 21-OCT-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!  Revision Number: 3                                                  C
!  Purpose: Define FLAG using info from ICBC_FLAG                      C
!  Author: M. Syamlal                                 Date: 21-APR-93  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, BC_DEFINED, BC_TYPE,     C
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        IS_K_b, IS_K_t, IS_J_s, IS_J_n, IS_I_w,       C
!                        BC_I_e, NO_I, NO_J, NO_K                      C
!  Variables modified: I, J, K, IJK, FLAG, FLAG_E, FLAG_N, FLAG_T,     C
!                                                                      C
!  Local variables: L, FLAGX                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_FLAGS 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE bc
      USE is
      USE indices
      USE physprop
      USE funits 
      USE compar        !//d
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
! 
!                      Indices 
      INTEGER          I, J, K, IJK 
! 
!                      Local DO loop index for b.c. specification 
      INTEGER          L 
! 
!                      Temporary storage for FLAG value 
      INTEGER          FLAGX 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Cell flag definitions
!  FLAG  ICBC_FLAG BC_TYPE        Cell type
!  ----- --------- -------        ---------
!   1       .        -            Cell containing gas or solids or both
!  10       p      P_INFLOW       Specified pressure inflow cell
!  11       P      P_OUTFLOW      Specified pressure outflow cell
!  20       I      MASS_INFLOW    Specified mass flux inflow cell
!  21       O      MASS_OUTFLOW   Specified mass flux outflow cell
!  31       o      OUTFLOW        outflow cell
! 100       W      NO_SLIP_WALL   Internal/external wall with no-slip b.c.
! 101       S      FREE_SLIP_WALL Internal/external wall with free-slip
! 102       s      PAR_SLIP_WALL  Internal/external wall with partial-slip b.c.
! 106       c      CYCLIC         Cyclic b.c.
! 107       C      CYCLIC_PD      Cyclic b.c. with pressure drop
! Flag values greater than 100 are considered to be wall cells
! (see function.inc).
!
!!$omp  parallel do private( IJK) &
!!$omp&  schedule(static)

!//AIKEPARDBG dump the ICBC_FLAG in matrix form to verify with serial version
!      write(*,"('(PE ',I2,'): IJKSTART3,IJKEND3 = ',2(I5,2X)))") & !//AIKEPARDBG
!              myPE,IJKSTART3,IJKEND3                    !//AIKEPARDBG
!      DO K = KMIN3, KMAX3                               !//AIKEPARDBG
!         write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
!	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
!         DO J = JMIN3, JMAX3                            !//AIKEPARDBG
!           write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	 
!           DO I = IMIN3, IMAX3                          !//AIKEPARDBG
!             IJK = FUNIJK_GL(I,J,K)                     !//AIKEPARDBG
!             write(UNIT_LOG,"(2X,A3,$)") ICBC_FLAG(IJK) !//AIKEPARDBG
!           END DO                                       !//AIKEPARDBG
!           write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
!         END DO                                         !//AIKEPARDBG
!      END DO                                            !//AIKEPARDBG
!      write(*,"('(PE ',I2,'): reached debug stop in set_flags')") myPE    !//AIKEPARDBG
!      call mfix_exit(myPE)   !//AIKEPARDBGSTOP
      
!// 200 1002 change the DO loop limits from i,ijkmax2 --> ijkstart3, ijkend3
      DO IJK = IJKSTART3, IJKEND3
!
         SELECT CASE (ICBC_FLAG(IJK)(1:1))  
         CASE ('.')  
            FLAG(IJK) = 1 
         CASE ('p')  
            FLAG(IJK) = 10 
         CASE ('P')  
            FLAG(IJK) = 11 
         CASE ('I')  
            FLAG(IJK) = 20 
         CASE ('O')  
            FLAG(IJK) = 21 
         CASE ('o')  
            FLAG(IJK) = 31 
         CASE ('W')  
            FLAG(IJK) = 100 
         CASE ('S')  
            FLAG(IJK) = 101 
         CASE ('s')  
            FLAG(IJK) = 102 
         CASE ('c')  
            FLAG(IJK) = 106 
         CASE ('C')  
            FLAG(IJK) = 107 
         CASE DEFAULT 

!Access to only one thread at a time
!!$omp       critical
            WRITE (UNIT_LOG, 1000) IJK, ICBC_FLAG(IJK) 
            call mfix_exit(myPE)
!!$omp       end critical
         END SELECT 
!
!  Initialize flags for internal surfaces.  If the flag is greater than
!  or equal to 2000, there is no internal surface.  UNDEFINED_I should
!  be a large +ve value.
!
         FLAG_E(IJK) = UNDEFINED_I 
         FLAG_N(IJK) = UNDEFINED_I 
         FLAG_T(IJK) = UNDEFINED_I 
      END DO 



      DO L = 1, DIMENSION_IS 
!
!       Make sure an IS has been specified
!
         IF (IS_DEFINED(L)) THEN 
            IF (IS_TYPE(L)=='IMPERMEABLE' .OR. IS_TYPE(L)(3:13)=='IMPERMEABLE'&
               ) THEN 
               FLAGX = 0 
            ELSE IF (IS_TYPE(L)=='SEMIPERMEABLE' .OR. IS_TYPE(L)(3:15)==&
                  'SEMIPERMEABLE') THEN 
               FLAGX = 1000 + L 
            ELSE 
               WRITE (UNIT_LOG, 1100) L 
               call mfix_exit(myPE)  
            ENDIF 
            IF (IS_X_W(L)==IS_X_E(L) .AND. DO_I) THEN 
               IS_PLANE(L) = 'E' 
               I = IS_I_W(L) 
               DO K = IS_K_B(L), IS_K_T(L) 
                  DO J = IS_J_S(L), IS_J_N(L) 
!//? Check the following filter
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!// 220 1004 Need to use local FUNIJK, size(FLAG_E) = DIMENSION_3L 
                     IJK = FUNIJK(I,J,K) 
                     FLAG_E(IJK) = FLAGX 
                  END DO 
               END DO 
            ELSE IF (IS_TYPE(L)(1:1) == 'X') THEN 
               IS_PLANE(L) = 'E' 
               DO I = IS_I_W(L), IS_I_E(L) 
                  DO K = IS_K_B(L), IS_K_T(L) 
                     DO J = IS_J_S(L), IS_J_N(L) 
!//? Check the following filter
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!// 220 1004 Need to use local FUNIJK, size(FLAG_E) = DIMENSION_3L 
                        IJK = FUNIJK(I,J,K) 
                        FLAG_E(IJK) = FLAGX 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
!
            IF (IS_Y_S(L)==IS_Y_N(L) .AND. DO_J) THEN 
               IS_PLANE(L) = 'N' 
               J = IS_J_S(L) 
               DO K = IS_K_B(L), IS_K_T(L) 
                  DO I = IS_I_W(L), IS_I_E(L) 
!//? Check the following filter		  
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		  
!// 220 1004 Need to use local FUNIJK, size(FLAG_E) = DIMENSION_3L 		
                     IJK = FUNIJK(I,J,K) 
                     FLAG_N(IJK) = FLAGX 
                  END DO 
               END DO 
            ELSE IF (IS_TYPE(L)(1:1) == 'Y') THEN 
               IS_PLANE(L) = 'N' 
               DO J = IS_J_S(L), IS_J_N(L) 
                  DO K = IS_K_B(L), IS_K_T(L) 
                     DO I = IS_I_W(L), IS_I_E(L) 
!//? Check the following filter		     
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
!// 220 1004 Need to use local FUNIJK, size(FLAG_E) = DIMENSION_3L 
                        IJK = FUNIJK(I,J,K) 
                        FLAG_N(IJK) = FLAGX 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
!
            IF (IS_Z_B(L)==IS_Z_T(L) .AND. DO_K) THEN 
               IS_PLANE(L) = 'T' 
               K = IS_K_B(L) 
               DO J = IS_J_S(L), IS_J_N(L) 
                  DO I = IS_I_W(L), IS_I_E(L) 
!//? Check the following filter		  
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		  
!// 220 1004 Need to use local FUNIJK		       
                     IJK = FUNIJK(I,J,K) 
                     FLAG_T(IJK) = FLAGX 
                  END DO 
               END DO 
            ELSE IF (IS_TYPE(L)(1:1) == 'Z') THEN 
               IS_PLANE(L) = 'T' 
               DO K = IS_K_B(L), IS_K_T(L) 
                  DO J = IS_J_S(L), IS_J_N(L) 
                     DO I = IS_I_W(L), IS_I_E(L) 
!//? Check the following filter		     
!// 360 1025 Check if current i,j,k resides on this PE		     
   		       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
!// 220 1004 Need to use local FUNIJK		     
                        IJK = FUNIJK(I,J,K) 
                        FLAG_T(IJK) = FLAGX 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
!
         ENDIF 
      END DO 
      RETURN  
 1000 FORMAT(/1X,70('*')//' From: SET_FLAGS',/' Message: ICBC_FLAG(',I3,') = ',&
         A3,' is illegal',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: SET_FLAGS',/' Message: Unknown IS_TYPE(',I3,&
         ')',/1X,70('*')/) 
      END SUBROUTINE SET_FLAGS 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_FLAGS1                                             C
!  Purpose: Assign IP flag to the faces of wall cells                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified: FLAG_E, FLAG_N, FLAG_T                          C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_FLAGS1 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar 
      USE geometry
      USE bc 
      USE is 
      USE indices
      USE physprop 
      USE funits 
      USE compar        !//d
      USE mpi_utility   !//SP
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      Indices 
      INTEGER          IJK,  IMJK, IJMK, IJKM, IPJK, IJPK, IJKP 
      INTEGER          I, J, K 
! 
!//SP
      INTEGER, DIMENSION(:), allocatable :: FLAG_TEMP
      integer :: flag_size
!-----------------------------------------------
      INCLUDE 'function.inc'
!//SP allocate storage for temporary flag arrays

      flag_size = ijkmax3
      if (myPE.eq.root) then
          flag_size = ijkmax3
      endif

      allocate( flag_temp(flag_size) )

!
!// 350 1025 change do loop limits: 1,ijkmax2->ijkstart3,ijkend3
      DO IJK = ijkstart3,ijkend3
         IMJK = IM_OF(IJK) 
         IJMK = JM_OF(IJK) 
         IJKM = KM_OF(IJK) 
         IPJK = IP_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKP = KP_OF(IJK) 
         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 
	 IF(.NOT.IS_ON_myPE_OWNS(I,J,K)) CYCLE
         IF (WALL_AT(IJK)) THEN 
            FLAG_E(IJK) = 0 
            FLAG_N(IJK) = 0 
            FLAG_T(IJK) = 0 
            FLAG_E(IMJK) = 0 
            FLAG_N(IJMK) = 0 
            FLAG_T(IJKM) = 0 
!
            IF (CYCLIC_AT(IJK)) THEN             ! make the upper (E, N, T) bdry 
!
!
               IF (I == IMAX2) THEN              ! permeable
                  IF (J/=1 .AND. J/=JMAX2) THEN
                     IF (NO_K) THEN
                        FLAG_E(IMJK) = 2000
                     ELSE IF (K/=1 .AND. K/=KMAX2) THEN
                        FLAG_E(IMJK) = 2000
                     ENDIF
                  ENDIF
               ENDIF
!
               IF (J == JMAX2) THEN
                  IF (I/=1 .AND. I/=IMAX2) THEN
                     IF (NO_K) THEN
                        FLAG_N(IJMK) = 2000
                     ELSE IF (K/=1 .AND. K/=KMAX2) THEN
                        FLAG_N(IJMK) = 2000
                     ENDIF
                  ENDIF
               ENDIF
!
               IF (K == KMAX2) THEN
                  IF (J/=1 .AND. J/=JMAX2) THEN
                     IF (I/=1 .AND. I/=IMAX2) FLAG_T(IJKM) = 2000
                  ENDIF
               ENDIF
!
            ENDIF
!//SP
!              IF (I == IMAX2.OR.I == IMAX3) THEN              ! permeable
!                 IF ((J/=1.OR.J/=0.) .AND. (J/=JMAX2.OR.J/=JMAX3)) THEN
!                    IF (NO_K) THEN
!                       FLAG_E(IMJK) = 2000
!                    ELSE IF ((K/=1.OR.K/=0) .AND. (K/=KMAX2.OR.K/=KMAX3)) THEN
!                       FLAG_E(IMJK) = 2000
!                    ENDIF
!                 ENDIF
!              ENDIF
!
!              IF (J == JMAX2.OR.J == JMAX3) THEN
!                 IF ((I/=1.OR.I/=0) .AND. (I/=IMAX2.OR.I/=IMAX3)) THEN
!                    IF (NO_K) THEN
!                       FLAG_N(IJMK) = 2000
!                    ELSE IF ((K/=1.OR.K/=0) .AND. (K/=KMAX2.OR.K/=KMAX3)) THEN
!                       FLAG_N(IJMK) = 2000
!                    ENDIF
!                 ENDIF
!              ENDIF
!
!              IF (K == KMAX2.OR.K == KMAX3) THEN
!                 IF ((J/=1.OR.J/=0.) .AND. (J/=JMAX2.OR.J/=JMAX3)) THEN
!                    IF ((I/=1.OR.I/=0) .AND. (I/=IMAX2.OR.I/=IMAX3)) FLAG_T(IJKM) = 2000
!                 ENDIF
!              ENDIF
!
!           ENDIF

!
         ELSE IF (FLUID_AT(IJK)) THEN 
            IF ( .NOT.WALL_AT(IMJK) .AND. FLAG_E(IMJK)==UNDEFINED_I) FLAG_E(&
               IMJK) = 2000 + FLAG(IMJK) 
            IF ( .NOT.WALL_AT(IJMK) .AND. FLAG_N(IJMK)==UNDEFINED_I) FLAG_N(&
               IJMK) = 2000 + FLAG(IJMK) 
            IF ( .NOT.WALL_AT(IJKM) .AND. FLAG_T(IJKM)==UNDEFINED_I) FLAG_T(&
               IJKM) = 2000 + FLAG(IJKM) 
            IF ( .NOT.WALL_AT(IPJK) .AND. FLAG_E(IJK)==UNDEFINED_I) FLAG_E(IJK)&
                = 2000 + FLAG(IPJK) 
            IF ( .NOT.WALL_AT(IJPK) .AND. FLAG_N(IJK)==UNDEFINED_I) FLAG_N(IJK)&
                = 2000 + FLAG(IJPK) 
            IF ( .NOT.WALL_AT(IJKP) .AND. FLAG_T(IJK)==UNDEFINED) FLAG_T(IJK)&
                = 2000 + FLAG(IJKP) 
         ENDIF 
      END DO 

!//SP Fill the ghost layers using gather and scatter
      call gather( flag_e, flag_temp )
      call scatter( flag_e, flag_temp )

      call gather( flag_n, flag_temp )
      call scatter( flag_n, flag_temp )

      call gather( flag_t, flag_temp )
      call scatter( flag_t, flag_temp )

!//SP deallocate storage of temporary flag arrays

    deallocate( flag_temp )

!//AIKEPARDBG dump the FLAG_X in matrix form to verify with serial version
!      DO K = Kstart3, Kend3                               !//AIKEPARDBG
!         write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
!	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
!         DO J = Jstart3, Jend3                            !//AIKEPARDBG
!           write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	
!           DO I = Istart3, Iend3                          !//AIKEPARDBG
!             IJK = FUNIJK(I,J,K)                     !//AIKEPARDBG
!             write(UNIT_LOG,"(2X,A3,$)") FLAG_T(IJK) !//AIKEPARDBG
!           END DO                                       !//AIKEPARDBG
!           write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
!         END DO                                         !//AIKEPARDBG
!      END DO    
      
      RETURN  
      END SUBROUTINE SET_FLAGS1 
