!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_GMRES(Vname, Var, A_m, B_m,                        C
!                         CMETHOD, TOL, ITMAX, MAX_IT, IER )
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_GMRES(VNAME, VAR, A_M, B_M, cmethod,TOL,ITMAX,MAX_IT,IER)
      
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE PARAM
      USE PARAM1
      USE MATRIX
      USE GEOMETRY
      USE INDICES
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER          IER
!
!                      maximum number of iterations
      INTEGER          ITMAX
!
!                      maximum number of outer iterations
      INTEGER          MAX_IT
!
!
!
!
!                      convergence tolerance
      DOUBLE PRECISION TOL
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_M(DIMENSION_3)
!
!                      Variable name
      CHARACTER*(*)    VNAME
!
!                      Variable
      DOUBLE PRECISION VAR(DIMENSION_3)

!                    sweep direction
      CHARACTER*4 :: CMETHOD
!-------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL LEQ_MATVEC, LEQ_MSOLVE


      IF (B_M(IJKMAX2/2) == ZERO) THEN  !This is usually not true so expensive 
                                        !DNRM2 calculation is avoided 
        IF (DNRM2(IJKMAX2,B_M,1) == ZERO) THEN 
           VAR(:IJKMAX2) = ZERO 
           RETURN  
        ENDIF 
      ENDIF 
      
      
      CALL LEQ_GMRES0( VNAME, VAR, A_M, B_M,                        &
                       CMETHOD, TOL, ITMAX, MAX_IT, LEQ_MATVEC, LEQ_MSOLVE, IER )

      RETURN
      END SUBROUTINE LEQ_GMRES




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LEQ_GMRES0(Vname, Var, A_m, B_m,                       C
!                    CMETHOD, TOL, ITMAX, MAX_IT, MATVEC, MSOLVE, IER ) 
!  Purpose: Compute residual of linear system                          C
!                                                                      C
!  Author: Ed D'Azevedo                               Date: 21-JAN-99  C
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
      SUBROUTINE LEQ_GMRES0(VNAME, VAR, A_M, B_M,  &
                            CMETHOD, TOL, ITMAX, MAX_IT, &
                            MATVEC, MSOLVE, IER )
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE PARAM
      USE PARAM1
      USE PARALLEL
      USE MATRIX
      USE GEOMETRY
      USE INDICES
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error indicator
      INTEGER          IER
!
!                      maximum number of iterations
      INTEGER          ITMAX
!
!                      maximum number of outer iterations
      INTEGER          MAX_IT
!
!
!                      convergence tolerance
      DOUBLE PRECISION TOL
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_M(DIMENSION_3)
!
!                      Variable name
      CHARACTER*(*)    VNAME
!
!                      Variable
      DOUBLE PRECISION VAR(DIMENSION_3)

!                    sweep direction
      CHARACTER*4 :: CMETHOD

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

	DOUBLE PRECISION V(IJKMAX2,ITMAX+1)
	DOUBLE PRECISION H(ITMAX+1,ITMAX)
	DOUBLE PRECISION CS(ITMAX)
	DOUBLE PRECISION SN(ITMAX)
	DOUBLE PRECISION Y(ITMAX)

	DOUBLE PRECISION R(IJKMAX2)
	DOUBLE PRECISION TEMP(IJKMAX2)
	DOUBLE PRECISION E1(IJKMAX2)
	DOUBLE PRECISION SS(IJKMAX2)
        DOUBLE PRECISION WW(IJKMAX2)

        DOUBLE PRECISION VY(IJKMAX2)
        

	DOUBLE PRECISION BNRM2, ERROR, ERROR0, NORM_R0, NORM_R, NORM_W
        DOUBLE PRECISION YII,YJJ


        INTEGER II,JJ,I,J,K, INC 
        INTEGER RESTRT, M, LDA, NDIM, ITER, MDIM
        DOUBLE PRECISION DTEMP, DSUM, INV_NORM_R, NORM_S
        DOUBLE PRECISION MINH,MAXH,CONDH,NORM_Y

	INTEGER IDEBUG
	PARAMETER(IDEBUG=0)

	DOUBLE PRECISION TOLTRIG
	PARAMETER(TOLTRIG=1.0D-4)







	NDIM = IJKMAX2
	RESTRT = ITMAX

	ITER = 0
	IER = 0
        IF (USE_DOLOOP) THEN
           BNRM2 = ZERO

!$omp      parallel do private(ii) reduction(+:bnrm2)
           DO II=1,NDIM
             BNRM2 = BNRM2 + B_M(II)*B_M(II)
           ENDDO
           BNRM2 = SQRT( BNRM2 )
        ELSE
	    BNRM2 = SQRT( DOT_PRODUCT( B_M(1:NDIM), B_M(1:NDIM) ) )
        ENDIF
        IF (BNRM2 .EQ. ZERO) THEN
	   BNRM2 = ONE
	ENDIF

!	------------------------
!	r = M \ (b - A*x)
!	error = norm(r) / bnrm2
!	------------------------

	CALL MATVEC( VNAME, VAR, A_M, R )

        IF (USE_DOLOOP) THEN

!$omp       parallel do private(ii)
            DO II=1,NDIM
              TEMP(II) = B_M(II) - R(II)
            ENDDO

        ELSE
	    TEMP(1:NDIM) = B_M(1:NDIM) - R(1:NDIM)
        ENDIF

        
	CALL MSOLVE( VNAME, TEMP, A_M,  R, CMETHOD )

        IF (USE_DOLOOP) THEN
            NORM_R = ZERO

!$omp       parallel do private(ii) reduction(+:norm_r)
            DO II=1,NDIM
              NORM_R = NORM_R + R(II)*R(II)
            ENDDO
            NORM_R = SQRT( NORM_R )
        ELSE
	    NORM_R = SQRT( DOT_PRODUCT( R(1:NDIM), R(1:NDIM) ) )
        ENDIF
	ERROR = NORM_R/BNRM2
        ERROR0 = ERROR
	NORM_R0 = NORM_R

        IF (IDEBUG.GE.1) THEN
           PRINT*,'LEQ_GMRES: ', VNAME,' INITIAL ERROR ', ERROR0
           PRINT*,'LEQ_GMRES: ', VNAME,' INITIAL RESIDUAL ', NORM_R0
	ENDIF


	M = RESTRT

        IF (USE_DOLOOP) THEN

!$omp      parallel do private(ii,jj)
           DO JJ=1,M+1
              DO II=1,NDIM
                 V(II,JJ) = ZERO
              ENDDO
           ENDDO

!$omp      parallel do private(ii,jj)
           DO JJ=1,M
             DO II=1,M+1
                H(II,JJ) = ZERO
             ENDDO
           ENDDO

           DO JJ=1,M
              CS(JJ) = ZERO
              SN(JJ) = ZERO
           ENDDO

!$omp      parallel do private(ii)
           DO II=1,NDIM
              E1(II) = ZERO
           ENDDO
           E1(1) = ONE
        

        ELSE
	    V(:,:) = ZERO
	    H(:,:) = ZERO
	    CS(:) = ZERO
	    SN(:) = ZERO
	    E1(:) = ZERO
	    E1(1) = ONE
        ENDIF


!	----------------
!	begin iteration
!	----------------
	DO ITER=1,MAX_IT

!		---------------
!		r = M \ (b-A*x)
!		---------------

		CALL MATVEC( VNAME, VAR, A_M, R )
                IF (USE_DOLOOP) THEN

!$omp              parallel do private(ii)
                   DO II=1,NDIM
                     TEMP(II) = B_M(II) - R(II)
                   ENDDO
                ELSE
		   TEMP(1:NDIM) = B_M(1:NDIM) - R(1:NDIM)
                ENDIF
	  	CALL MSOLVE( VNAME, TEMP, A_M, R, CMETHOD )

                IF (USE_DOLOOP) THEN
		    NORM_R = ZERO
!$omp               parallel do private(ii) reduction(+:norm_r)
                    DO II=1,NDIM
                        NORM_R = NORM_R + R(II)*R(II)
                    ENDDO
		    NORM_R = SQRT(NORM_R)


!		    ----------------------------------------
!                   On most machines, multipy is much faster
!                   than divides
!		    ----------------------------------------
                    INV_NORM_R = ONE / NORM_R
!$omp               parallel do private(ii)
		    DO II=1,NDIM
			V(II,1) = R(II)*INV_NORM_R
                    ENDDO

!$omp               parallel do private(ii)
		    DO II=1,NDIM
                        SS(II) = NORM_R * E1(II)
                    ENDDO

                ELSE
		    NORM_R = SQRT( DOT_PRODUCT( R(1:NDIM), R(1:NDIM) ) )
		    V(:,1) = R(:)/NORM_R
		    SS(:) = NORM_R * E1(:)
                ENDIF


!		----------------------------
!		construct orthonormal basis
!		using Gram-Schmidt
!		----------------------------
		DO I=1,M
!		   ------------------
!		   w = M \ (A*V(:,i))
!		   ------------------
		   CALL MATVEC( VNAME, V(:,I), A_M, TEMP )
		   CALL MSOLVE( VNAME, TEMP, A_M, WW, CMETHOD )

                   IF (USE_DOLOOP) THEN
                        DO K=1,I
                            DTEMP = ZERO

!$omp                       parallel do private(ii) reduction(+:dtemp)
                            DO II=1,NDIM
                                DTEMP = DTEMP + WW(II)*V(II,K)
                            ENDDO
                            H(K,I) = DTEMP

!$omp                       parallel do private(ii)
                            DO II=1,NDIM
                              WW(II) = WW(II) - H(K,I)*V(II,K)
                            ENDDO
                        ENDDO
                   ELSE
		       DO K=1,I
			    H(K,I) = DOT_PRODUCT( WW(:), V(:,K) )
			    WW(:) = WW(:) - H(K,I)*V(:,K)
		       ENDDO
                   ENDIF

                   IF (USE_DOLOOP) THEN
                       NORM_W = ZERO
!$omp                  parallel do private(ii) reduction(+:norm_w)
                       DO II=1,NDIM
                           NORM_W = NORM_W  + WW(II)*WW(II)
                       ENDDO
                       NORM_W = SQRT( NORM_W  )
                   ELSE
		       NORM_W = SQRT( DOT_PRODUCT( WW(:), WW(:) ) )
                   ENDIF
		   H(I+1,I) = NORM_W


                   IF (USE_DOLOOP) THEN
                       DTEMP = ONE / H(I+1,I)

!$omp                  parallel do private(ii)
                       DO II=1,NDIM
                          V(II,I+1) = WW(II)*DTEMP
                       ENDDO
     
                   ELSE
		       V(:,I+1) = WW(:) / H(I+1,I)
                   ENDIF

!		   -----------------------
!		   apply Givens rotation
!		   -----------------------
		   DO K=1,I-1
			DTEMP    =  CS(K)*H(K,I) + SN(K)*H(K+1,I)
			H(K+1,I) = -SN(K)*H(K,I) + CS(K)*H(K+1,I)
			H(K,I)   = DTEMP
		   ENDDO

!		   --------------------------
!		   form i-th rotation matrix
!		   approximate residual norm
!		   --------------------------

		   CALL ROTMAT( H(I,I), H(I+1,I), CS(I), SN(I) )

		   DTEMP = CS(I)*SS(I)
		   SS(I+1) = -SN(I)*SS(I)
		   SS(I) = DTEMP

		   H(I,I) = CS(I)*H(I,I) + SN(I)*H(I+1,I)
		   H(I+1,I) = ZERO

		   ERROR = ABS( SS(I+1) ) / BNRM2
		   IF (ERROR .LE. TOL*ERROR0) THEN
!			-----------------------------
!			update approximation and exit
!			-----------------------------
!			-----------------------
!			triangular solve with
!			y = H(1:i,1:i) \ s(1:i)
!			-----------------------
		        MDIM = I

			DO II=1,MDIM
			  Y(II) = SS(II)
			ENDDO

			DO II=MDIM,1,-1
			   YII = Y(II)/H(II,II)
			   Y(II) = YII
			   DO JJ=1,II-1
			       Y(JJ) = Y(JJ) - H(JJ,II)*YII
			   ENDDO
			ENDDO

                	IF (IDEBUG.GE.1) THEN
!               	-------------
!               	double check
!               	-------------
			MAXH = ABS(H(1,1))
			MINH = ABS(H(1,1))
			DO II=1,MDIM
			  MAXH = MAX( MAXH, ABS(H(II,II)) )
			  MINH = MIN( MINH, ABS(H(II,II)) )
			ENDDO
			CONDH = MAXH/MINH

                   	TEMP(1:MDIM) = SS(1:MDIM) - &
                        	MATMUL( H(1:MDIM,1:MDIM ), Y(1:MDIM) )   
                   	DTEMP = DOT_PRODUCT( TEMP(1:MDIM), TEMP(1:MDIM) ) 
                   	DTEMP = SQRT( DTEMP )
                   	NORM_S = SQRT(DOT_PRODUCT( SS(1:MDIM),SS(1:MDIM) ))
                   	NORM_Y = SQRT(DOT_PRODUCT( Y(1:MDIM),Y(1:MDIM) ))
                   	IF (DTEMP .GT. CONDH*NORM_S) THEN
				PRINT *,'DTEMP, NORM_S ', DTEMP, NORM_S
				PRINT *,'CONDH, NORM_Y ', CONDH, NORM_Y
                        	STOP '** STOP IN LEQ_GMRES ** '
			ENDIF
                	ENDIF


     
			IF (USE_DOLOOP) THEN
!$omp               parallel do private(ii)
			    DO II=1,NDIM
			      VY(II) = ZERO
			    ENDDO

			   DO JJ=1,MDIM
			      YJJ = Y(JJ)

!$omp                         parallel do private(ii)
			      DO II=1,NDIM
				VY(II) = VY(II) +  V(II,JJ)*YJJ
			      ENDDO
			   ENDDO

!$omp              parallel do private(ii)
			   DO II=1,NDIM
			       VAR(II) = VAR(II) + VY(II)
			   ENDDO
 
                        ELSE
		           VAR(1:NDIM)=VAR(1:NDIM)+MATMUL(V(:,1:I),Y(1:I))
		        ENDIF


			EXIT
		   ENDIF
		ENDDO


		IF ( ERROR .LE. TOL*ERROR0) THEN
			EXIT
		ENDIF

!		---------------------
!		update approximations
!		---------------------

!		------------------------
!		y = H(1:m,1:m) \ s(1:m)
!		x = x + V(:,1:m)*y
!		r = M \ (b-A*x)
!		------------------------


		MDIM = M
                DO II=1,MDIM
		  Y(II) = SS(II)
                ENDDO

                DO II=MDIM,1,-1
		   YII = Y(II)/H(II,II)
                   Y(II) = YII
                   DO JJ=1,II-1
		       Y(JJ) = Y(JJ) - H(JJ,II)*YII
                   ENDDO
		ENDDO
		

                        IF (IDEBUG.GE.1) THEN
!                       -------------
!                       double check
!                       -------------
                        MAXH = ABS(H(1,1))
                        MINH = ABS(H(1,1))
                        DO II=1,MDIM
                          MAXH = MAX( MAXH, ABS(H(II,II)) )
                          MINH = MIN( MINH, ABS(H(II,II)) )
                        ENDDO
                        CONDH = MAXH/MINH

                        TEMP(1:MDIM) = SS(1:MDIM) - &
                                MATMUL( H(1:MDIM,1:MDIM ), Y(1:MDIM) )
                        DTEMP = DOT_PRODUCT( TEMP(1:MDIM), TEMP(1:MDIM) )
                        DTEMP = SQRT( DTEMP )
                        NORM_S = SQRT(DOT_PRODUCT( SS(1:MDIM),SS(1:MDIM) ))
                        NORM_Y = SQRT(DOT_PRODUCT( Y(1:MDIM),Y(1:MDIM) ))
                        IF (DTEMP .GT. CONDH*NORM_S) THEN
                                PRINT *,'DTEMP, NORM_S ', DTEMP, NORM_S
                                PRINT *,'CONDH, NORM_Y ', CONDH, NORM_Y
                                STOP '** STOP IN LEQ_GMRES ** '
                        ENDIF
                        ENDIF




                IF (USE_DOLOOP) THEN

                   MDIM = M

!$omp               parallel do private(ii)
                    DO II=1,NDIM
                      VY(II) = ZERO
		    ENDDO

                   DO JJ=1,MDIM
                      YJJ = Y(JJ)

!$omp                 parallel do private(ii)
		      DO II=1,NDIM
                        VY(II) = VY(II) +  V(II,JJ)*YJJ
                      ENDDO
                   ENDDO
    
!$omp              parallel do private(ii)
                   DO II=1,NDIM
		       VAR(II) = VAR(II) + VY(II)
                   ENDDO	
                ELSE
		 VAR(1:NDIM) = VAR(1:NDIM) + MATMUL( V(1:NDIM,1:M), Y(1:M) )
                ENDIF

		 
                CALL MATVEC( VNAME, VAR, A_M, R )

                IF (USE_DOLOOP) THEN

!$omp               parallel do private(ii)
                    DO II=1,NDIM
                      TEMP(II) = B_M(II) - R(II)
                    ENDDO
                ELSE
                    TEMP(1:NDIM) = B_M(1:NDIM) - R(1:NDIM)
                ENDIF    
                CALL MSOLVE( VNAME, TEMP, A_M, R, CMETHOD )

                IF (USE_DOLOOP) THEN
                    NORM_R = ZERO
!$omp               parallel do private(ii) reduction(+:norm_r)
                    DO II=1,NDIM
                      NORM_R = NORM_R + R(II)*R(II)
                    ENDDO
                    NORM_R = SQRT(NORM_R)
                ELSE
		    NORM_R = SQRT( DOT_PRODUCT( R(1:NDIM), R(1:NDIM) ) )
                ENDIF
		SS(I+1) = NORM_R

		ERROR = SS(I+1) / BNRM2

                IF (IDEBUG.GE.1) THEN
                   PRINT*,'LEQ_GMRES: I = ', I, ' ERROR = ', ERROR
		   PRINT*,'LEQ_GMRES: I = ', I,' NORM_R ', NORM_R
		ENDIF

		IF (ERROR .LE. TOL*ERROR0)  THEN
			EXIT
		ENDIF

	   ENDDO

	IF (ERROR .GE. TOL*ERROR0) THEN
		IER = 1
	ENDIF

	IF (IDEBUG.GE.1) THEN
	    CALL MATVEC( VNAME, VAR, A_M, R )
	    R(1:NDIM) = R(1:NDIM) - B_M(1:NDIM)
	    NORM_R=SQRT(DOT_PRODUCT(R(1:NDIM),R(1:NDIM)))
	    PRINT*,'LEQ_GMRES: ITER ', ITER,' I ',I
	    PRINT*,'LEQ_GMRES: ', VNAME, ' FINAL ERROR ',ERROR
	    PRINT*,'LEQ_GMRES: ', VNAME, ' FINAL RESIDUAL ',NORM_R
	    PRINT*,'LEQ_GMRES: ', VNAME, ' ERR RATIO ',ERROR/ERROR0
	    PRINT*,'LEQ_GMRES: ', VNAME, ' RESID RATIO ',NORM_R/NORM_R0

	ENDIF




        RETURN

	END SUBROUTINE LEQ_GMRES0

!       ==============================================
        SUBROUTINE ROTMAT(  A, B, C, S )
	IMPLICIT NONE
	DOUBLE PRECISION A,B,C,S

        DOUBLE PRECISION ONE,ZERO
        PARAMETER(ONE=1.0D0,ZERO=0.0D0)
        DOUBLE PRECISION TEMP

	IF (B.EQ.ZERO) THEN
            C = ONE
            S = ZERO
        ELSE IF (ABS(B) .GT. ABS(A)) THEN

           TEMP = A / B
           S = ONE / SQRT( ONE + TEMP*TEMP )
           C = TEMP * S

        ELSE

           TEMP = B / A
           C = ONE / SQRT( ONE + TEMP*TEMP )
           S = TEMP * C
        ENDIF


        RETURN
        END

  


