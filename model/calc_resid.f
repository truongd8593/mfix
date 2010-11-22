!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_C(Var, A_m, B_m, M, RESID,                  C
!                            MAX_RESID, IJK_RESID, IER)                C
!  Purpose: Calculate residuals for continuity equations               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_C(VAR,A_M,B_M,M,NUM,DEN,RESID,MAX_RESID,IJK_RESID,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar
      USE mpi_utility 
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER

!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!!!HPF$ align Var(:) with TT(:)
!!!HPF$ align A_m(:,*,*) with TT(:)
!!!HPF$ align B_m(:,*) with TT(:)

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID
 
!                      Indices
      INTEGER          IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT

      INTEGER          I, J, K      

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Error message
      CHARACTER*80     LINE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!// New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc

      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 

!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO
!
!$omp  parallel do private( IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!$omp&  NUM1, DEN1) &
!$omp&  REDUCTION(+:NUM,DEN,NCELLS)  
      DO IJK = ijkstart3, ijkend3

      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE      
         IF (FLUID_AT(IJK)) THEN 
!
            IJKW = WEST_OF(IJK) 
            IJKS = SOUTH_OF(IJK) 
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK) 
!
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,E,M)*VAR(IJKE)+&
               A_M(IJK,W,M)*VAR(IJKW)+A_M(IJK,N,M)*VAR(IJKN)+A_M(IJK,S,M)*VAR(&
               IJKS)) 
            IF (DO_K) THEN 
               IJKB = BOTTOM_OF(IJK) 
               IJKT = TOP_OF(IJK) 
!
               NUM1 = NUM1 - (A_M(IJK,T,M)*VAR(IJKT)+A_M(IJK,B,M)*VAR(IJKB)) 
            ENDIF 
!

            NUM1 = ABS(NUM1) 
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK)) 
!
!efd            IF (NUM1 > MAX_RESID) THEN 
!efd               MAX_RESID = NUM1 
!efd               IJK_RESID = IJK 
!efd            ENDIF 

            RESID_IJK(IJK) = NUM1
!     
            NCELLS = NCELLS + 1 
!     
            NUM = NUM + NUM1 
            DEN = DEN + DEN1 
!
         ENDIF 
      END DO 
!efd

      if(.not.debug_resid) return

      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE      
         IF (RESID_IJK(IJK) > MAX_RESID) then
               IJK_RESID = IJK
               MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

!//  Determine max residual
      do nproc=0,NumPEs-1
	if(nproc.eq.myPE) then
	MAX_RESID_L(nproc) = MAX_RESID
	IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
	else
	MAX_RESID_L(nproc) = 0.0
	IJK_RESID_L(nproc) = 0
	endif
      enddo

!//  Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)

!//  Collect all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

!//  Call to determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2

      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo
      
      IF (DEN > ZERO) THEN 
         RESID = NUM/DEN 
         MAX_RESID = NCELLS*MAX_RESID/DEN 
      ELSE IF (NUM == ZERO) THEN 
         RESID = ZERO 
         MAX_RESID = ZERO 
         IJK_RESID = 0 
      ELSE 
         RESID = UNDEFINED 
         MAX_RESID = UNDEFINED 
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.' 
!         CALL WRITE_ERROR ('CALC_RESID_C', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CALC_RESID_C 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_S(Var, A_m, B_m, M, RESID,                  C
!                            MAX_RESID, IJK_RESID, TOL, IER)           C
!  Purpose: Calculate residuals for sclar equations                    C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_S(VAR,A_M,B_M,M,NUM,DEN,RESID,MAX_RESID,IJK_RESID,TOL, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar      
      USE mpi_utility
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER

!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID
 
!                      Ignore residual calculation for scalar values below this
      DOUBLE PRECISION TOL

!                      Indices
      INTEGER          IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP

      INTEGER          I, J, K

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Error message
      CHARACTER*80     LINE
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!// New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc

      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 

!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!
!$omp    parallel do &
!$omp&   private(   IJK,     &
!$omp&   IMJK,IJMK,IPJK,IJPK, IJKM,IJKP, &
!$omp&   NUM1, DEN1) &
!$omp&   REDUCTION(+:NUM, DEN,NCELLS)  
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK) .AND. ABS(VAR(IJK)) > TOL) THEN 
!
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
!
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*VAR(IJK)+A_M(IJK,E,M)*VAR(IPJK)+&
               A_M(IJK,W,M)*VAR(IMJK)+A_M(IJK,N,M)*VAR(IJPK)+A_M(IJK,S,M)*VAR(&
               IJMK)) 
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IJKP = KP_OF(IJK) 
               NUM1 = NUM1 - (A_M(IJK,T,M)*VAR(IJKP)+A_M(IJK,B,M)*VAR(IJKM)) 
            ENDIF 
!
            NUM1 = ABS(NUM1) 
            DEN1 = ABS(A_M(IJK,0,M)*VAR(IJK)) 
!
!efd            IF (NUM1 > MAX_RESID) THEN 
!efd               MAX_RESID = NUM1 
!efd               IJK_RESID = IJK 
!efd            ENDIF 

            RESID_IJK(IJK) = NUM1
!
            NCELLS = NCELLS + 1 
!
            NUM = NUM + NUM1 
            DEN = DEN + DEN1 
!
         ENDIF 
      END DO 

      if(.not.debug_resid) return

!// Collect all the information among all the procesors
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)
!efd
      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
         IF (RESID_IJK(IJK) > MAX_RESID) THEN               
               IJK_RESID = IJK
               MAX_RESID = RESID_IJK( IJK_RESID )
         ENDIF
      ENDDO

!// 
      do nproc=0,NumPEs-1
	if(nproc.eq.myPE) then
	MAX_RESID_L(nproc) = MAX_RESID
	IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
	else
	MAX_RESID_L(nproc) = 0.0
	IJK_RESID_L(nproc) = 0
	endif
      enddo

!//  Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)


!//  Collect all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

!//  Determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2

      do nproc=0,NumPEs-1

        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then

        IJK_RESID = IJK_RESID_GL(nproc)

        endif

      enddo

      IF (DEN > ZERO) THEN 
         RESID = NUM/DEN 
         MAX_RESID = NCELLS*MAX_RESID/DEN 
      ELSE IF (NUM == ZERO) THEN 
         RESID = ZERO 
         MAX_RESID = ZERO 
         IJK_RESID = 0 
      ELSE 
         RESID = UNDEFINED 
         MAX_RESID = UNDEFINED 
!        WRITE(LINE,*)'Message: All center coefficients are zero.'
!        CALL WRITE_ERROR('CALC_RESID_S', LINE, 1)
      ENDIF 
!
      RETURN  
      END SUBROUTINE CALC_RESID_S 
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_pp(B_m, NORM, RESID, MAX_RESID, IJK_RESID,  C
!                             IER)                                     C
!  Purpose: Calculate residuals for pressure correction equation       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_PP(B_M, NORM, NUM, DEN, RESID,MAX_RESID, IJK_RESID, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar   
      USE mpi_utility   
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!
!                      Normalization factor
      DOUBLE PRECISION NORM

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID
 
!                      Indices
      INTEGER          IJK

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Error message
      CHARACTER*80     LINE

      INTEGER          I, J, K      
!-----------------------------------------------
      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!//  New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc
      
      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 
!
      DEN1 = ONE 
!
!FIX AEOLUS 032608 Missing initialization 
      IJK_RESID = 1

      DO IJK = ijkstart3, ijkend3 
        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
      
         IF (FLUID_AT(IJK)) THEN 
!
            NUM1 = ABS(B_M(IJK,0)) 
!
            IF (NUM1 > MAX_RESID) THEN 
               MAX_RESID = NUM1 
               IJK_RESID = IJK 
            ENDIF 
!
            NCELLS = NCELLS + 1 
!
            NUM = NUM + NUM1 
            DEN = DEN + DEN1 
!
         ENDIF 
      END DO 


      if(.not.debug_resid) then

         call global_all_sum(NUM)
         call global_all_sum(DEN)

         IF (DEN*NORM > ZERO) THEN 
            RESID = NUM/(DEN*NORM) 
         ELSE IF (NUM == ZERO) THEN 
            RESID = ZERO 
         ELSE 
            RESID = LARGE_NUMBER 
         ENDIF 

      else

!//   Determine the global sum
         call global_all_sum(NUM)
         call global_all_sum(DEN)
         call global_all_sum(NCELLS)

!//   
         do nproc=0,NumPEs-1
            if(nproc.eq.myPE) then
               MAX_RESID_L(nproc) = MAX_RESID
               IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
            else
               MAX_RESID_L(nproc) = 0.0
               IJK_RESID_L(nproc) = 0
            endif
         enddo

!//   Determine the maximum among all the procesors
         call global_all_max(MAX_RESID)


!//   Collect all the information among all the procesors
         call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
         call global_all_sum(IJK_RESID_L, IJK_RESID_GL)

!//   Determine the global IJK location w.r.t. serial version
         IJK_RESID = IJKMAX2

         do nproc=0,NumPEs-1

            if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then

               IJK_RESID = IJK_RESID_GL(nproc)

            endif

         enddo            
         IF (DEN*NORM > ZERO) THEN 
            RESID = NUM/(DEN*NORM) 
            MAX_RESID = NCELLS*MAX_RESID/(DEN*NORM) 
         ELSE IF (NUM == ZERO) THEN 
            RESID = ZERO 
            MAX_RESID = ZERO 
            IJK_RESID = 0 
         ELSE 
            RESID = LARGE_NUMBER 
            MAX_RESID = LARGE_NUMBER 
!     WRITE (LINE, *) 'Warning: All center coefficients are zero.' 
!     CALL WRITE_ERROR ('CALC_RESID_pp', LINE, 1) 
         ENDIF 

      endif                     ! debug_resid_g
!
      RETURN  
      END SUBROUTINE CALC_RESID_PP 

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_mb(INIT, ErrorPercent)                      C
!  Purpose: Calculate overall mass balance error                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-DEC-02   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_MB(INIT, ErrorPercent) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE parallel 
      USE geometry
      USE indices
      USE run
      USE bc
      USE constant
      USE physprop
      USE compar   
      USE mpi_utility   
      USE residual
      USE rxns
      USE mflux
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Flag to check whether this is an initialization call
!                      0 -> initialize old accumulation; 1 -> calc residual
      INTEGER          init
      
!                      Error index
      INTEGER          IER
      
!                      Phase index
      INTEGER          M

!                      Total mass balance error as a % of inflow
      DOUBLE PRECISION ErrorPercent(0:MMAX)

!                      Locally define dt, so that this routine works when dt is not defined
      DOUBLE PRECISION dt_local

 
!                      Indices
      INTEGER          L, IJK

      DOUBLE PRECISION flux_in, flux_out, fin, fout, err, accum_new, denom

!     functions      
      DOUBLE PRECISION Accumulation
    
      INCLUDE 'function.inc'

      if(dt == UNDEFINED)then
        dt_local = ONE
      else
        dt_local = dt
      endif

      if(init == 0) then
!       Initilaize this routine
	!Accumulation
        if(dt == UNDEFINED)then
          Accum_resid_g = ZERO
        else
          Accum_resid_g = Accumulation(ROP_g)
        endif
	DO M=1, MMAX
          if(dt == UNDEFINED)then
            Accum_resid_s(M) = ZERO
	  else
            Accum_resid_s(M) = Accumulation(ROP_s(1,M))
	  endif
	END DO
	return
	
      else
        if(dt == UNDEFINED)then
	  Accum_new = - Accumulation(SUM_R_g) * dt_local
	else
	  Accum_new = Accumulation(ROP_g) - Accumulation(SUM_R_g) * dt_local
	endif 
	
	flux_out = zero
        flux_in = zero
        DO L = 1, DIMENSION_BC
          IF (BC_DEFINED(L)) THEN
!            call Calc_mass_flux(BC_I_W(L), BC_I_E(L), BC_J_S(L), & 
!            BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), U_g, V_g, W_g, &
!            ROP_g, fin, fout, IER) 
            IF(.NOT.Added_Mass) THEN
	      call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), & 
              BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), Flux_gE, Flux_gN, Flux_gT, &
              fin, fout, IER)  
            ELSE
	      call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), & 
              BC_J_N(L), BC_K_B(L), BC_K_T(L), BC_PLANE(L), Flux_gSE, Flux_gSN, Flux_gST, &
              fin, fout, IER)  
            ENDIF
	    flux_out = flux_out + fout  * dt_local
            flux_in = flux_in + fin * dt_local
          ENDIF 
        END DO
        
	Err = (accum_new - Accum_resid_g) - (flux_in - flux_out)
	denom = max(abs(accum_new), abs(Accum_resid_g), abs(flux_in), abs(flux_out))
        IF (denom /= ZERO) THEN
	   ErrorPercent(0) = err*100./denom
        ELSE
           ErrorPercent(0) = err*100./SMALL_NUMBER
        END IF
	
	DO M =1, MMAX
          if(dt == UNDEFINED)then
	    Accum_new =  - Accumulation(SUM_R_s(1,M)) * dt_local
	  else
	    Accum_new = Accumulation(ROP_s(1,M)) - Accumulation(SUM_R_s(1,M)) * dt_local
	  endif
	
	  flux_out = zero
          flux_in = zero
          DO L = 1, DIMENSION_BC
            IF (BC_DEFINED(L)) THEN
!              call Calc_mass_flux(BC_I_W(L), BC_I_E(L), BC_J_S(L), BC_J_N(L), &
!              BC_K_B(L), BC_K_T(L), BC_PLANE(L), U_s(1,M), V_s(1,M), W_s(1,M), &
!              ROP_s(1,M), fin, fout, IER)  
              IF(.NOT.Added_Mass .OR. M /= M_AM) THEN
                call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), BC_J_N(L), &
                BC_K_B(L), BC_K_T(L), BC_PLANE(L), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), &
                fin, fout, IER)  
              ELSE 
                call Calc_mass_fluxHR(BC_I_W(L), BC_I_E(L), BC_J_S(L), BC_J_N(L), &
                BC_K_B(L), BC_K_T(L), BC_PLANE(L), Flux_sSE, Flux_sSN, Flux_sST, &
                fin, fout, IER)  
              ENDIF
	      flux_out = flux_out + fout  * dt_local
              flux_in = flux_in + fin * dt_local
            ENDIF 
          END DO
        
	  Err = (accum_new - Accum_resid_s(M)) - (flux_in - flux_out)
	  denom = max(abs(accum_new), abs(Accum_resid_s(M)), abs(flux_in), abs(flux_out))
	  if(denom /= ZERO) THEN
	    ErrorPercent(M) = err*100./denom
	  else
	    ErrorPercent(M) = err*100./SMALL_NUMBER
	  endif
	END DO
	
      endif

      RETURN  
      END SUBROUTINE CALC_RESID_MB 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_U(U_m, V_m, W_m, A_m, B_m, M, RESID,        C
!                            MAX_RESID, IJK_RESID, IER)                C
!  Purpose: Calculate residuals for u-momentum equations               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_U(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, RESID,MAX_RESID, &
         IJK_RESID, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar 
      USE mpi_utility 
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER

!
!                      U
      DOUBLE PRECISION U_m(DIMENSION_3)
!
!                      V
      DOUBLE PRECISION V_m(DIMENSION_3)
!
!                      W
      DOUBLE PRECISION W_m(DIMENSION_3)

!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID
 
!
!                      Velocity magnitude
      DOUBLE PRECISION VEL

!                      Indices
      INTEGER          IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Error message
      CHARACTER*80     LINE

      INTEGER          I, J, K
      
!-----------------------------------------------
!     Local variables
!-----------------------------------------------

      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!// New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc


      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 
!efd
!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!
!$omp  parallel do private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp&   NUM1, DEN1,VEL) &
!$omp&  REDUCTION(+:NUM, DEN,NCELLS )  
      DO IJK = ijkstart3, ijkend3
        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
      
         IF (.NOT.IP_AT_E(IJK)) THEN 
!
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
!
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*U_M(IJK)+A_M(IJK,E,M)*U_M(IPJK)+&
               A_M(IJK,W,M)*U_M(IMJK)+A_M(IJK,N,M)*U_M(IJPK)+A_M(IJK,S,M)*U_M(&
               IJMK)) 
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IJKP = KP_OF(IJK) 
!
               NUM1 = NUM1 - (A_M(IJK,T,M)*U_M(IJKP)+A_M(IJK,B,M)*U_M(IJKM)) 
            ENDIF 
!
!         Ignore momentum residual in stagnant regions.  Need an alternative
!         criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2) 
            IF (VEL > SMALL_NUMBER) THEN 
               NUM1 = ABS(NUM1) 
               DEN1 = ABS(A_M(IJK,0,M)*VEL) 
!
!efd               IF (NUM1 > MAX_RESID) THEN 
!efd                  MAX_RESID = NUM1 
!efd                  IJK_RESID = IJK 
!efd               ENDIF 

               RESID_IJK(IJK) = NUM1
!
               NCELLS = NCELLS + 1 
!
               NUM = NUM + NUM1 
               DEN = DEN + DEN1 
            ENDIF 
!
         ENDIF 
      END DO 

      if(.not.debug_resid) return

!//
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)
!efd
      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )

      DO IJK = ijkstart3, ijkend3
        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
      
          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

!//
      do nproc=0,NumPEs-1
	if(nproc.eq.myPE) then
	MAX_RESID_L(nproc) = MAX_RESID
	IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
	else
	MAX_RESID_L(nproc) = 0.0
	IJK_RESID_L(nproc) = 0
	endif
      enddo

!//  Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)


!//  Collect all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)
      
!//  Determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2

      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

      IF (DEN > ZERO) THEN 
         RESID = NUM/DEN 
         MAX_RESID = NCELLS*MAX_RESID/DEN 
      ELSE IF (NUM == ZERO) THEN 
         RESID = ZERO 
         MAX_RESID = ZERO 
         IJK_RESID = 0 
      ELSE 
         RESID = LARGE_NUMBER 
         MAX_RESID = LARGE_NUMBER 
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.' 
!         CALL WRITE_ERROR ('CALC_RESID_U', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CALC_RESID_U 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_V(U_m, V_m, W_m, A_m, B_m, M, RESID,        C
!                            MAX_RESID, IJK_RESID, IER)                C
!  Purpose: Calculate residuals for v-momentum equations               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_V(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, RESID,MAX_RESID, &
         IJK_RESID, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar   
      USE mpi_utility  
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER

!
!                      U
      DOUBLE PRECISION U_m(DIMENSION_3)
!
!                      V
      DOUBLE PRECISION V_m(DIMENSION_3)
!
!                      W
      DOUBLE PRECISION W_m(DIMENSION_3)

!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID

!
!                      Velocity magnitude
      DOUBLE PRECISION VEL
 
!                      Indices
      INTEGER          IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Error message
      CHARACTER*80     LINE
      
      INTEGER          I, J, K
      
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
!efd
      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!// New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc



      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 
!efd
!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!
!$omp  parallel do private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp&  VEL,  NUM1, DEN1) &
!$omp&  REDUCTION(+:NUM, DEN, NCELLS)  
      DO IJK = ijkstart3, ijkend3 
        IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE

         IF (.NOT.IP_AT_N(IJK)) THEN 
!
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
!
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*V_M(IJK)+A_M(IJK,E,M)*V_M(IPJK)+&
               A_M(IJK,W,M)*V_M(IMJK)+A_M(IJK,N,M)*V_M(IJPK)+A_M(IJK,S,M)*V_M(&
               IJMK)) 
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IJKP = KP_OF(IJK) 
!
               NUM1 = NUM1 - (A_M(IJK,T,M)*V_M(IJKP)+A_M(IJK,B,M)*V_M(IJKM)) 
            ENDIF 
!
!         Ignore momentum residual in stagnant regions.  Need an alternative
!         criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2) 
            IF (VEL > SMALL_NUMBER) THEN 
               NUM1 = ABS(NUM1) 
               DEN1 = ABS(A_M(IJK,0,M)*VEL) 
!
!efd               IF (NUM1 > MAX_RESID) THEN 
!efd                  MAX_RESID = NUM1 
!efd                  IJK_RESID = IJK 
!efd               ENDIF 
               RESID_IJK(IJK) = NUM1

!
               NCELLS = NCELLS + 1 
!
               NUM = NUM + NUM1 
               DEN = DEN + DEN1 
            ENDIF 
!
         ENDIF 
      END DO 


      if(.not.debug_resid) return

!//
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)
      
      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE      
          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

!//
      do nproc=0,NumPEs-1
	if(nproc.eq.myPE) then
	MAX_RESID_L(nproc) = MAX_RESID
	IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
	else
	MAX_RESID_L(nproc) = 0.0
	IJK_RESID_L(nproc) = 0
	endif
      enddo

!//  Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)


!//  Collect all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)
      
!//  Determine the global IJK location w.r.t. serial version

      IJK_RESID = IJKMAX2
      
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo
      
      IF (DEN > ZERO) THEN 
         RESID = NUM/DEN 
         MAX_RESID = NCELLS*MAX_RESID/DEN 
      ELSE IF (NUM == ZERO) THEN 
         RESID = ZERO 
         MAX_RESID = ZERO 
         IJK_RESID = 0 
      ELSE 
         RESID = LARGE_NUMBER 
         MAX_RESID = LARGE_NUMBER 
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.' 
!         CALL WRITE_ERROR ('CALC_RESID_V', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CALC_RESID_V 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_RESID_W(U_m, V_m, W_m, A_m, B_m, M, RESID,        C
!                            MAX_RESID, IJK_RESID, IER)                C
!  Purpose: Calculate residuals for w-momentum equations               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_RESID_W(U_M, V_M, W_M, A_M, B_M, M, NUM, DEN, RESID,MAX_RESID, &
         IJK_RESID, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE parallel 
      USE geometry
      USE indices
      USE compar 
      USE mpi_utility 
      USE run 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER

!
!                      U
      DOUBLE PRECISION U_m(DIMENSION_3)
!
!                      V
      DOUBLE PRECISION V_m(DIMENSION_3)
!
!                      W
      DOUBLE PRECISION W_m(DIMENSION_3)

!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!                      Phase index
      INTEGER          M

!                      Average value of Residual
      DOUBLE PRECISION RESID

!                      Maximum value of Residual
      DOUBLE PRECISION MAX_RESID

!                      IJK of Maximum value of Residual
      INTEGER          IJK_RESID

!
!                      Velocity magnitude
      DOUBLE PRECISION VEL
 
!                      Indices
      INTEGER          IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP

!                      Numerators and denominators
      DOUBLE PRECISION NUM, NUM1, DEN, DEN1

!                      Number of fluid cells
      INTEGER          NCELLS 

!                      Error message
      CHARACTER*80     LINE

      INTEGER          I, J, K      
!-----------------------------------------------
!     Local variables
!-----------------------------------------------
      double precision, dimension (ijksize3_all(myPE)) :: RESID_IJK
!// New local variables for DMP version
      DOUBLE PRECISION     MAX_RESID_GL(0:numPEs-1), MAX_RESID_L(0:numPEs-1)
      INTEGER              IJK_RESID_GL(0:numPEs-1), IJK_RESID_L(0:numPEs-1)
      INTEGER              nproc


      INCLUDE 'function.inc'
!
      NUM = ZERO 
      DEN = ZERO 
      MAX_RESID = -ONE 
      NCELLS = 0 

!$omp parallel do private( IJK )
      DO IJK = ijkstart3, ijkend3
          RESID_IJK(IJK) = ZERO
      ENDDO

!
!$omp  parallel do private( IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
!$omp&  VEL,  NUM1, DEN1) &
!$omp&  REDUCTION(+:NUM, DEN,NCELLS )  
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE
      
         IF (.NOT.IP_AT_T(IJK)) THEN 
!
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
!
            NUM1 = B_M(IJK,M) - (A_M(IJK,0,M)*W_M(IJK)+A_M(IJK,E,M)*W_M(IPJK)+&
               A_M(IJK,W,M)*W_M(IMJK)+A_M(IJK,N,M)*W_M(IJPK)+A_M(IJK,S,M)*W_M(&
               IJMK)) 
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IJKP = KP_OF(IJK) 
!
               NUM1 = NUM1 - (A_M(IJK,T,M)*W_M(IJKP)+A_M(IJK,B,M)*W_M(IJKM)) 
            ENDIF 
!
!         Ignore momentum residual in stagnant regions.  Need an alternative
!         criteria for residual scaling for such cases.
            VEL = SQRT(U_M(IJK)**2+V_M(IJK)**2+W_M(IJK)**2) 
            IF (VEL > SMALL_NUMBER) THEN 
               NUM1 = ABS(NUM1) 
               DEN1 = ABS(A_M(IJK,0,M)*VEL) 
!
!efd               IF (NUM1 > MAX_RESID) THEN 
!efd                  MAX_RESID = NUM1 
!efd                  IJK_RESID = IJK 
!efd               ENDIF 
               RESID_IJK(IJK) = NUM1

!
               NCELLS = NCELLS + 1 
!
               NUM = NUM + NUM1 
               DEN = DEN + DEN1 
            ENDIF 
!
         ENDIF 
      END DO 


      if(.not.debug_resid) return

!//
      call global_all_sum(NUM)
      call global_all_sum(DEN)
      call global_all_sum(NCELLS)
!efd

      IJK_RESID = 1
      MAX_RESID = RESID_IJK( IJK_RESID )
      DO IJK = ijkstart3, ijkend3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK),J_OF(IJK), K_OF(IJK))) CYCLE      
      
          IF (RESID_IJK( IJK ) > MAX_RESID) THEN
              IJK_RESID = IJK
              MAX_RESID = RESID_IJK( IJK_RESID )
          ENDIF
      ENDDO

!//
      do nproc=0,NumPEs-1
	if(nproc.eq.myPE) then
	MAX_RESID_L(nproc) = MAX_RESID
	IJK_RESID_L(nproc) = FUNIJK_GL(I_OF(IJK_RESID), J_OF(IJK_RESID), K_OF(IJK_RESID))
	else
	MAX_RESID_L(nproc) = 0.0
	IJK_RESID_L(nproc) = 0
	endif
      enddo

!//  Determine the maximum among all the procesors
      call global_all_max(MAX_RESID)


!//  Collect all the information among all the procesors
      call global_all_sum(MAX_RESID_L, MAX_RESID_GL)
      call global_all_sum(IJK_RESID_L, IJK_RESID_GL)
      
!//  Determine the global IJK location w.r.t. serial version
      IJK_RESID = IJKMAX2
      
      do nproc=0,NumPEs-1
        if(MAX_RESID_GL(nproc).eq.MAX_RESID.and.IJK_RESID_GL(nproc).lt.IJK_RESID) then
          IJK_RESID = IJK_RESID_GL(nproc)
        endif
      enddo

      IF (DEN > ZERO) THEN 
         RESID = NUM/DEN 
         MAX_RESID = NCELLS*MAX_RESID/DEN 
      ELSE IF (NUM == ZERO) THEN 
         RESID = ZERO 
         MAX_RESID = ZERO 
         IJK_RESID = 0 
      ELSE 
         RESID = LARGE_NUMBER 
         MAX_RESID = LARGE_NUMBER 
!         WRITE (LINE, *) 'Warning: All center coefficients are zero.' 
!         CALL WRITE_ERROR ('CALC_RESID_W', LINE, 1) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE CALC_RESID_W 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization, I,J,K,MAX_RESID_XX,IJK_RESID_XX,nproc
!// 050 Replace global array size with  subdomain array size, RESID_IJK
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 360 Check if i,j,k resides on current processor
!// 400 Added mpi_utility module and other global reduction (sum,max) calls
!//     Some other changes for Determine the global IJK location w.r.t. serial version
