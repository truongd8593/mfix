!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_Mass_balance(init)                               C
!  Purpose: Check global species and elemental balances                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-Nov-02  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: ABORT, SUM                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_Mass_balance (init) 
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE fldvar
      USE rxns
      USE geometry
      USE run
      USE bc
      USE constant
      USE physprop
      USE indices
      USE funits 
      USE compar 
      USE mpi_utility  
      USE output
      USE check
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
!                      Flag to check whether this is an initialization call
!                      0 -> initialization call; 1 -> integration call
      INTEGER          init
!
!                      Indices
      INTEGER          I, J, K, IJK, IER
!
!                      Solids phase
      INTEGER          M
!
!                      Species index
      INTEGER          N
!
!                      Do-loop counter
      INTEGER          L
! 
!                      Starting and ending indices (Make sure these represent a plane and NOT a volume)
      INTEGER          I1, I2, J1, J2, K1, K2 
! 
!                       
      DOUBLE PRECISION Accumulation_old, Accumulation_delta, flux, flux_in_tot, flux_out_tot, error_percent 
      DOUBLE PRECISION fin, fout
      
!     functions      
      DOUBLE PRECISION Accumulation, Accumulation_sp
!-----------------------------------------------
      INCLUDE 'function.inc'

      if(report_mass_balance_dt == UNDEFINED) return

      if(init == 0) then
!       allocate arrays
!       Initilaize this routine
        start_time = time
        report_time = time + report_mass_balance_dt
	
	!initialize flux and reaction rate arrays
        DO L = 1, DIMENSION_BC
	  flux_out_g(L) = ZERO
	  flux_in_g(L) = ZERO
          DO N = 1, NMAX(0)
	    flux_out_X_g(L, N) = ZERO
	    flux_in_X_g(L, N) = ZERO
          END DO
	  
	  DO M = 1, MMAX
	    flux_out_s(L, M) = ZERO
	    flux_in_s(L, M) = ZERO
            DO N = 1, NMAX(M)
	      flux_out_X_s(L, M, N) = ZERO
	      flux_in_X_s(L, M, N) = ZERO
            END DO
          END DO

        END DO

	Integral_SUM_R_g = ZERO
        DO N = 1, NMAX(0)
          Integral_R_g(N) = ZERO
        END DO
	
	DO M = 1, MMAX
	  Integral_SUM_R_s(M) = ZERO
          DO N = 1, NMAX(M)
            Integral_R_s(M, N) = ZERO
          END DO
        END DO

	!Accumulation
        Accumulation_g = Accumulation(ROP_g)
        DO N = 1, NMAX(0)
	  Accumulation_X_g(N) = Accumulation_sp(ROP_g, X_g(1, N)) 
        END DO 
	
	DO M = 1, MMAX
          Accumulation_s(M) = Accumulation(ROP_s(1,M))
          DO N = 1, NMAX(M)
	    Accumulation_X_s(M, N) = Accumulation_sp(ROP_s(1, M), X_s(1, M, N)) 
          END DO 
        END DO 
 
	return
	
	
      else
!       Flow in and out of the boundaries
!       Use dt_prev for these integrations because adjust_dt may have changed the dt
        DO L = 1, DIMENSION_BC
          IF (BC_DEFINED(L)) THEN
            I1 = BC_I_W(L) 
            I2 = BC_I_E(L) 
            J1 = BC_J_S(L) 
            J2 = BC_J_N(L) 
            K1 = BC_K_B(L) 
            K2 = BC_K_T(L)
	     
            call Calc_mass_flux(I1, I2, J1, J2, K1, K2, BC_PLANE(L), U_g, V_g, W_g, ROP_g, fin, fout, IER) 
	    flux_out_g(L) = flux_out_g(L) + fout  * dt_prev
            flux_in_g(L) = flux_in_g(L) + fin * dt_prev
	    
	    DO M = 1, MMAX
              call Calc_mass_flux(I1, I2, J1, J2, K1, K2, BC_PLANE(L), U_s(1,m), V_s(1,m), W_s(1,m), ROP_s(1,m), fin, fout, IER) 
	      flux_out_s(L, M) = flux_out_s(L, M) + fout  * dt_prev
              flux_in_s(L, M) = flux_in_s(L, M) + fin * dt_prev
	    END DO
	    
          ENDIF 
        END DO
	
	Integral_SUM_R_g = Integral_SUM_R_g + Accumulation(SUM_R_g) * dt_prev
	IF(SPECIES_EQ(0))THEN 
          DO N = 1, NMAX(0)
            DO L = 1, DIMENSION_BC
              IF (BC_DEFINED(L)) THEN
                I1 = BC_I_W(L) 
                I2 = BC_I_E(L) 
                J1 = BC_J_S(L) 
                J2 = BC_J_N(L) 
                K1 = BC_K_B(L) 
                K2 = BC_K_T(L) 
                call Calc_mass_flux_sp(I1, I2, J1, J2, K1, K2, BC_PLANE(L), U_g, V_g, W_g, ROP_g, X_g(1, N), fin, fout, IER) 
	        flux_out_X_g(L, N) = flux_out_X_g(L, N) + fout  * dt_prev
                flux_in_X_g(L, N) = flux_in_X_g(L, N) + fin * dt_prev
              ENDIF 
	    
	      Integral_R_g(N) = Integral_R_g(N) + (Accumulation(R_gp(1,N)) - &
	                                Accumulation_sp(ROX_gc(1,N), X_g(1,N)) )* dt_prev 
            END DO
	
          END DO 
	ENDIF

        DO M = 1, MMAX	
	  IF(SPECIES_EQ(M))THEN 
	    Integral_SUM_R_s(M) = Integral_SUM_R_s(M) + Accumulation(SUM_R_s(1,M)) * dt_prev 
            DO N = 1, NMAX(M)
              DO L = 1, DIMENSION_BC
                IF (BC_DEFINED(L)) THEN
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  call Calc_mass_flux_sp(I1, I2, J1, J2, K1, K2, BC_PLANE(L), &
                  U_s(1,M), V_s(1,M), W_s(1,M), ROP_s(1,M), X_s(1, M, N), fin, fout, &
                  IER) 
	          flux_out_X_s(L, M, N) = flux_out_X_s(L, M, N) + fout  * dt_prev
                  flux_in_X_s(L, M, N) = flux_in_X_s(L, M, N) + fin * dt_prev
                ENDIF 
	    
	        Integral_R_s(M, N) = Integral_R_s(M, N) + (Accumulation(R_sp(1,M,N)) - &
	                                            Accumulation_sp(ROX_sc(1,M,N), X_s(1,M,N)) )* dt_prev 
              END DO
	
            END DO 
	  ENDIF
        END DO 
      
      endif 
      
      
      if (time >= report_time)then
	
        CALL START_LOG 
        WRITE(UNIT_LOG, '(/A,G12.5,A,G12.5)') 'Mass balance for interval ', start_time, ' to ', time

	Accumulation_old = Accumulation_g
        Accumulation_g = Accumulation(ROP_g)
	Accumulation_delta = Accumulation_g - Accumulation_old - Integral_SUM_R_g
        WRITE(UNIT_LOG, '(A)') 'Total Fluid Accumulation (g)'
        WRITE(UNIT_LOG, '(4(A,G12.5))') '  Old = ', Accumulation_old, ', New = ', &
	  Accumulation_g, ', Production = ', Integral_SUM_R_g, ', net accu(New - Old - Production) = ', Accumulation_delta
	  
        WRITE(UNIT_LOG, '(A)') 'Integral of boundary flux (g)'
	Write(Unit_log, '(A, T8, A, T21, A, T34, A)')'  BC#', 'in', 'out', '(in - out)'
	flux = zero
	flux_in_tot = zero
	flux_out_tot = zero
        DO L = 1, DIMENSION_BC
	  if(flux_out_g(L) /= ZERO .OR. flux_in_g(L) /= ZERO) &
	    Write(Unit_log, '(2X, I5, 1X, 3(G12.5, 1x))')L, flux_in_g(L), flux_out_g(L), (flux_in_g(L)-flux_out_g(L))
	  flux = flux + flux_in_g(L) - flux_out_g(L)
	  flux_in_tot = flux_in_tot + flux_in_g(L)
	  flux_out_tot = flux_out_tot + flux_out_g(L)
	END DO
	if((flux - Accumulation_delta) /= zero) then
	  error_percent = undefined
	  if(flux_in_tot /= zero) error_percent = (flux - Accumulation_delta)*100./flux_in_tot
	else
	  error_percent = zero
	endif
	Write(Unit_log, '(2X, A, 1X, 3(G12.5, 1x))')'Total', flux_in_tot, flux_out_tot, flux
	Write(Unit_log, '(A, G12.5, A, G12.5)')'Error (net influx - net accu) = ', &
	  (flux - Accumulation_delta), ' %Error = ', error_percent
	  
	DO M = 1, MMAX
	  Accumulation_old = Accumulation_s(M)
          Accumulation_s(M) = Accumulation(ROP_s(1, M))
	  Accumulation_delta = Accumulation_s(M) - Accumulation_old - Integral_SUM_R_s(M)
          WRITE(UNIT_LOG, '(/A, I1, A)') 'Total Solids-', M, ' Accumulation (g)'
          WRITE(UNIT_LOG, '(4(A,G12.5))') '  Old = ', Accumulation_old, ', New = ', &
	    Accumulation_s(M), ', Production = ', Integral_SUM_R_s(M), ', net accu(New - Old - Production) = ', Accumulation_delta
	  
          WRITE(UNIT_LOG, '(A)') 'Integral of boundary flux (g)'
	  Write(Unit_log, '(A, T8, A, T21, A, T34, A)')'  BC#', 'in', 'out', '(in - out)'
	  flux = zero
	  flux_in_tot = zero
	  flux_out_tot = zero
          DO L = 1, DIMENSION_BC
	    if(flux_out_s(L,M) /= ZERO .OR. flux_in_s(L,M) /= ZERO) &
	      Write(Unit_log, '(2X, I5, 1X, 3(G12.5, 1x))')L, flux_in_s(L,M), flux_out_s(L,M), (flux_in_s(L,M)-flux_out_s(L,M))
	    flux = flux + flux_in_s(L,M) - flux_out_s(L,M)
	    flux_in_tot = flux_in_tot + flux_in_s(L,M)
	    flux_out_tot = flux_out_tot + flux_out_s(L,M)
	  END DO
	  if((flux - Accumulation_delta) /= zero) then
	    if(flux_in_tot /= zero) then
	      error_percent = (flux - Accumulation_delta)*100./flux_in_tot
	    else
	      error_percent = (flux - Accumulation_delta)*100./Accumulation_old
	    endif
	  else
	    error_percent = zero
	  endif
	  Write(Unit_log, '(2X, A, 1X, 3(G12.5, 1x))')'Total', flux_in_tot, flux_out_tot, flux
	  Write(Unit_log, '(A, G12.5, A, G12.5)')'Error (net influx - net accu) = ', &
	    (flux - Accumulation_delta), ' %Error = ', error_percent
	END DO


        IF(SPECIES_EQ(0) )THEN
          DO N = 1, NMAX(0)
	
            WRITE(UNIT_LOG, '(/A,I2)') 'Gas species - ', N
	  
	    Accumulation_old = Accumulation_X_g(N)
	    Accumulation_X_g(N) = Accumulation_sp(ROP_g, X_g(1, N)) 
	    Accumulation_delta = Accumulation_X_g(N) - Accumulation_old - Integral_R_g(N)
            WRITE(UNIT_LOG, '(A)') 'Species Accumulation (g)'
            WRITE(UNIT_LOG, '(4(A,G12.5))') '  Old = ', Accumulation_old, ', New = ', &
	      Accumulation_X_g(N), ', Production = ', Integral_R_g(N), ', net accu(New - Old - Production) = ', Accumulation_delta
	  
            WRITE(UNIT_LOG, '(A)') 'Integral of boundary flux (g)'
	    Write(Unit_log, '(A, T8, A, T21, A, T34, A)')'  BC#', 'in', 'out', '(in - out)'
	    flux = zero
	    flux_in_tot = zero
	    flux_out_tot = zero
            DO L = 1, DIMENSION_BC
	      if(flux_out_X_g(L, N) /= ZERO .OR. flux_in_X_g(L, N) /= ZERO) &
	        Write(Unit_log, '(2X, I5, 1X, 3(G12.5, 1x))')L, flux_in_X_g(L, N), flux_out_X_g(L, N), &
	       (flux_in_X_g(L, N)-flux_out_X_g(L, N))
	      flux = flux + flux_in_X_g(L, N) - flux_out_X_g(L, N)
	      flux_in_tot = flux_in_tot + flux_in_X_g(L, N)
	      flux_out_tot = flux_out_tot + flux_out_X_g(L, N)
  	    END DO
	    if((flux - Accumulation_delta) /= zero) then
	      error_percent = undefined
	      if(flux_in_tot /= zero) error_percent = (flux - Accumulation_delta)*100./flux_in_tot
	    else
	      error_percent = zero
	    endif
	    Write(Unit_log, '(2X, A, 1X, 3(G12.5, 1x))')'Total', flux_in_tot, flux_out_tot, flux
	    Write(Unit_log, '(A, G12.5, A, I1, A, G12.5)')'Error (net influx - net accu) = ', &
	     (flux - Accumulation_delta), ' %Error(',N,') = ', error_percent
	  
          END DO
	ENDIF
	 
	DO M = 1, MMAX
          IF(SPECIES_EQ(M) )THEN
            DO N = 1, NMAX(M)
	
              WRITE(UNIT_LOG, '(/A,I1, A, I2)') 'Solids-', M, ' species - ', N
	  
	      Accumulation_old = Accumulation_X_s(M,N)
  	      Accumulation_X_s(M,N) = Accumulation_sp(ROP_s(1,M), X_s(1, M, N)) 
	      Accumulation_delta = Accumulation_X_s(M,N) - Accumulation_old - Integral_R_s(M,N)
              WRITE(UNIT_LOG, '(A)') 'Species Accumulation (g)'
              WRITE(UNIT_LOG, '(4(A,G12.5))') '  Old = ', Accumulation_old, ', New = ', &
	        Accumulation_X_s(M,N), ', Production = ', Integral_R_s(M,N), ', net accu(New - Old - Production) = ', Accumulation_delta
	  
              WRITE(UNIT_LOG, '(A)') 'Integral of boundary flux (g)'
	      Write(Unit_log, '(A, T8, A, T21, A, T34, A)')'  BC#', 'in', 'out', '(in - out)'
	      flux = zero
	      flux_in_tot = zero
	      flux_out_tot = zero
              DO L = 1, DIMENSION_BC
	        if(flux_out_X_s(L, M, N) /= ZERO .OR. flux_in_X_s(L, M, N) /= ZERO) &
	          Write(Unit_log, '(2X, I5, 1X, 3(G12.5, 1x))')L, flux_in_X_s(L, M, N), flux_out_X_s(L, M, N), &
	           (flux_in_X_s(L, M, N)-flux_out_X_s(L, M, N))
	        flux = flux + flux_in_X_s(L, M, N) - flux_out_X_s(L, M, N)
	        flux_in_tot = flux_in_tot + flux_in_X_s(L, M, N)
	        flux_out_tot = flux_out_tot + flux_out_X_s(L, M, N)
	      END DO
	      if((flux - Accumulation_delta) /= zero) then
	        error_percent = undefined
	        if(flux_in_tot /= zero) error_percent = (flux - Accumulation_delta)*100./flux_in_tot
	      else
	        error_percent = zero
	      endif
	      Write(Unit_log, '(2X, A, 1X, 3(G12.5, 1x))')'Total', flux_in_tot, flux_out_tot, flux
	      Write(Unit_log, '(A, G12.5, A, I1, A, G12.5)')'Error (net influx - net accu) = ', &
	        (flux - Accumulation_delta), ' %Error(',N,') = ', error_percent
	  
            END DO
	  ENDIF 
        END DO 
        WRITE(UNIT_LOG, '(/)')
	CALL END_LOG
	
        start_time = time
        report_time = time + report_mass_balance_dt
	
	!initialize flux and reaction rate arrays
        DO L = 1, DIMENSION_BC
	  flux_out_g(L) = ZERO
	  flux_in_g(L) = ZERO
          DO N = 1, NMAX(0)
	    flux_out_X_g(L, N) = ZERO
	    flux_in_X_g(L, N) = ZERO
          END DO
	  
	  DO M = 1, MMAX
	    flux_out_s(L, M) = ZERO
	    flux_in_s(L, M) = ZERO
            DO N = 1, NMAX(M)
	      flux_out_X_s(L, M, N) = ZERO
	      flux_in_X_s(L, M, N) = ZERO
            END DO
          END DO

        END DO

	Integral_SUM_R_g = ZERO
        DO N = 1, NMAX(0)
          Integral_R_g(N) = ZERO
        END DO
	
	DO M = 1, MMAX
	  Integral_SUM_R_s(M) = ZERO
          DO N = 1, NMAX(M)
            Integral_R_s(M, N) = ZERO
          END DO
        END DO
	
	
      endif
      

      RETURN  
      END SUBROUTINE CHECK_Mass_balance 




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Calc_mass_flux                                         C
!  Purpose: Calculate the mass flux (g/s) across a plane               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-Dec-02   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE  Calc_mass_flux(I1, I2, J1, J2, K1, K2, Plane, U, V, W, ROP, flux_in, flux_out, IER) 
      USE param 
      USE param1 
      IMPLICIT NONE
      
! 
!                      Starting and ending indices (Make sure these represent a plane and NOT a volume)
      INTEGER ::           I1, I2, J1, J2, K1, K2 

!                      For each (scalar) cell the plane (W, E, S, N, B, T) to be used for flux calc
      Character :: 	       Plane
 
!                      Components of Velocity
      DOUBLE PRECISION ::  U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
 
!                      Macroscopic density
      DOUBLE PRECISION ::  ROP(DIMENSION_3)

!                      Species Mass fraction
      DOUBLE PRECISION :: Xn(DIMENSION_3)

!                      flux across the plane (composed of many cell faces)
!                      flux_in: flux from the scalar cell into the rest of the domain
!                      flux_out: flux into the scalar cell from the rest of the domain
      DOUBLE PRECISION ::  flux_in, flux_out

      INTEGER ::           IER
      
      Xn = one
      call Calc_mass_flux_sp(I1, I2, J1, J2, K1, K2, PLANE, U, V, W, ROP, Xn, flux_in, flux_out, IER) 
      return
      end SUBROUTINE  Calc_mass_flux	 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Calc_mass_flux_sp                                      C
!  Purpose: Calculate the species mass flux (g/s) across a plane       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-Dec-02   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE  Calc_mass_flux_sp(I1, I2, J1, J2, K1, K2, Plane, U, V, W, ROP, Xn, flux_in, flux_out, IER) 
      USE param 
      USE param1 
      USE geometry
      USE physprop
      USE indices
      USE compar 
      USE mpi_utility  
      IMPLICIT NONE
      
! 
!                      Starting and ending indices (Make sure these represent a plane and NOT a volume)
      INTEGER ::           I1, I2, J1, J2, K1, K2 

!                      For each (scalar) cell the plane (W, E, S, N, B, T) to be used for flux calc
      Character :: 	       Plane
 
!                      Components of Velocity
      DOUBLE PRECISION ::  U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
 
!                      Macroscopic density
      DOUBLE PRECISION ::  ROP(DIMENSION_3)

!                      Species Mass fraction
      DOUBLE PRECISION :: Xn(DIMENSION_3)

!                      flux across the plane (composed of many cell faces)
!                      flux_in: flux from the scalar cell into the rest of the domain
!                      flux_out: flux into the scalar cell from the rest of the domain
      DOUBLE PRECISION ::  flux_in, flux_out, flux_in_global, flux_out_global

      INTEGER ::           IER
!
!                      Indices
      INTEGER          I, J, K, IJK
      
      INCLUDE 'function.inc'

      flux_in = zero
      flux_out = zero
      IER = 0
      
     
            DO K = K1, K2 
            DO J = J1, J2 
            DO I = I1, I2 
              IF(.NOT.IS_ON_myPE_owns(I, J, K)) cycle

              SELECT CASE (PLANE)  
              CASE ('W')  
                IJK = IM_OF(FUNIJK(I,J,K))
		IF(U(IJK) > ZERO)THEN
		  flux_out = flux_out + AYZ(IJK) * U(IJK) * ROP(IJK)  * Xn(IJK)
		ELSE
		  flux_in = flux_in - AYZ(IJK) * U(IJK) * ROP(IP_OF(IJK))  * Xn(IP_OF(IJK))
		ENDIF

              CASE ('E')  
                IJK = FUNIJK(I,J,K)
		IF(U(IJK) > ZERO)THEN
		  flux_in = flux_in + AYZ(IJK) * U(IJK) * ROP(IJK)  * Xn(IJK) 
		ELSE
		  flux_out = flux_out - AYZ(IJK) * U(IJK) * ROP(IP_OF(IJK))  * Xn(IP_OF(IJK)) 
		ENDIF
		
              CASE ('S')  
                IJK = JM_OF(FUNIJK(I,J,K))
		IF(V(IJK) > ZERO)THEN
		  flux_out = flux_out + AXZ(IJK) * V(IJK) * ROP(IJK)  * Xn(IJK) 
		ELSE
		  flux_in = flux_in - AXZ(IJK) * V(IJK) * ROP(JP_OF(IJK))  * Xn(JP_OF(IJK))
		ENDIF

              CASE ('N')  
                IJK = FUNIJK(I,J,K)
		IF(V(IJK) > ZERO)THEN
		  flux_in = flux_in + AXZ(IJK) * V(IJK) * ROP(IJK)  * Xn(IJK) 
		ELSE
		  flux_out = flux_out - AXZ(IJK) * V(IJK) * ROP(JP_OF(IJK))  * Xn(JP_OF(IJK)) 
		ENDIF

              CASE ('B')  
                IJK = KM_OF(FUNIJK(I,J,K))
		IF(W(IJK) > ZERO)THEN
		  flux_out = flux_out + AXY(IJK) * W(IJK) * ROP(IJK)  * Xn(IJK)
		ELSE
		  flux_in = flux_in - AXY(IJK) * W(IJK) * ROP(KP_OF(IJK))  * Xn(KP_OF(IJK)) 
		ENDIF

              CASE ('T')  
                IJK = FUNIJK(I,J,K)
		IF(W(IJK) > ZERO)THEN
		  flux_in = flux_in + AXY(IJK) * W(IJK) * ROP(IJK)  * Xn(IJK) 
		ELSE
		  flux_out = flux_out - AXY(IJK) * W(IJK) * ROP(KP_OF(IJK))  *  Xn(KP_OF(IJK)) 
		ENDIF
		
	      CASE DEFAULT 
!                IER = 1
!                CALL START_LOG 
!                WRITE (UNIT_LOG, '(A, A1)' ) 'From: Calc_mass_flux, Unknown Plane: ', Plane
!                CALL MFIX_EXIT(myPE) 
		
              END SELECT 
	    ENDDO
	    ENDDO
	    ENDDO

            call global_all_sum(flux_in)
            call global_all_sum(flux_out)
	 
      return
      end SUBROUTINE  Calc_mass_flux_sp	 
      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Accumulation(ro)                                       C
!  Purpose: Intergrate density accumulation over the entire domain     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-DEC-02  C
!  Local variables:  None                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION Accumulation(ro) 
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
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      density distributiom
      DOUBLE PRECISION ::  RO(DIMENSION_3)
!                      Species Mass fraction
      DOUBLE PRECISION :: Xn(DIMENSION_3)
      
!     functions      
      DOUBLE PRECISION Accumulation_sp
! 
!-----------------------------------------------

      Xn = ONE
      Accumulation = Accumulation_sp(ro, xn)
      
      RETURN  
      END FUNCTION Accumulation 
    
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Accumulation_sp(ro, Xn)                                C
!  Purpose: Intergrate density accumulation over the entire domain     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-DEC-02  C
!  Local variables:  None                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION Accumulation_sp(ro, Xn) 
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
      USE physprop
      USE geometry
      USE indices
      USE compar 
      USE mpi_utility 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Indices 
      INTEGER          I, J, K, IJK 
      
!                      density distributiom
      DOUBLE PRECISION ::  RO(DIMENSION_3)
!                      Species Mass fraction
      DOUBLE PRECISION :: Xn(DIMENSION_3)
      
      DOUBLE PRECISION SUM
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
      
      SUM = ZERO
!
!!$omp$   parallel do private(IJK) &
!!$omp&   reduction(+:SUM)
      DO IJK = ijkstart3, ijkend3 
	IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) cycle
          IF (FLUID_AT(IJK)) SUM = SUM + RO(IJK) * Xn(IJK) * VOL(IJK) 
      END DO 

      call global_all_sum(sum)
      Accumulation_sp = sum
      
      RETURN  
      END FUNCTION Accumulation_sp 
    
