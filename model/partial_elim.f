!TO DO:
! 1. add capability to turn off certain momentum equations  -- done
! 2. Account for solids-solids transfer terms for scalars.  Will need to pass the
!    exchange coefficients to this routine (unlike momentum eqs)
! 3. repeatetive code may be eliminated by defining and using functions 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_S(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for scalar quantities               C
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
      SUBROUTINE PARTIAL_ELIM_S(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE compar  
      USE drag
      USE fldvar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKW, IJKS, IJKB, IJKE, IJKN, IJKT 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
! ! !	******GERA MODIFICATIONS************
      DOUBLE PRECISION SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
      INTEGER           L, M,LP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a(0:DIMENSION_M), BB(0:DIMENSION_M), F(0:DIMENSION_M,0:DIMENSION_M),&
                       Saxf(0:DIMENSION_M) 
!                     error message 
      CHARACTER*80     LINE 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'

!$omp  parallel do private( IJKW, IJKS, IJKB, IJKE, IJKN, IJKT,  &
!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, den) &
!$omp&  schedule(static)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN 
               IJKW = WEST_OF(IJK) 
               IJKS = SOUTH_OF(IJK)
               IJKE = EAST_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               DO M=0, MMAX
	         A(M)=A_M(IJK,0,M)
		 BB(M)=B_M(IJK,M)

		 if (m .ne. 0) then
		   F(M,0)=-VXF(IJK,M)
		   F(0,M)=F(M,0)
                 else
		   F(0,0)=0.   
	         end if
		 
		 DO L =1, MMAX
		   IF ((L .NE. M) .AND. (M .NE. 0)) THEN
		     F(M,L)=0.0  !insert solids-solids exchange coefficients here
		   ELSE
		     F(M,L)=0.
		   END IF
		   F(L,M)=0.0
		 END DO
		  
		 IF (M == 0 ) THEN
		   SAXF(M) = -(A_M(IJK,E,M)*VAR_G(IJKE)+A_M(IJK,W,M)*VAR_G(IJKW&
                          )+A_M(IJK,N,M)*VAR_G(IJKN)+A_M(IJK,S,M)*VAR_G(IJKS)) 
		 ELSE
		   SAXF(M) = -(A_M(IJK,E,M)*VAR_S(IJKE,M)+A_M(IJK,W,M)*VAR_S(&
                          IJKW,M)+A_M(IJK,N,M)*VAR_S(IJKN,M)+A_M(IJK,S,M)*VAR_S(&
                          IJKS,M))
		 ENDIF 

                 IF (DO_K) THEN 
                   IJKB = BOTTOM_OF(IJK) 
                   IJKT = TOP_OF(IJK)
		   IF ( M ==0) THEN
                     SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_G(IJKT)+A_M(IJK,B,M)*&
                               VAR_G(IJKB)) 
		   ELSE
                     SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_S(IJKT,M)+A_M(IJK,B,M)*&
                               VAR_S(IJKB,M)) 
		   ENDIF
                 ENDIF 
			
	       END DO
		  
	       Do M=0,MMAX
		 SUM_A=0.
		 SUM_B=0.0
		 do L =0,MMAX
		   SUM_A_LPRIME=0.0
		   SUM_B_LPRIME=0.0
		      
		   DO LP=0,MMAX
		     IF ( LP .NE. M) THEN
			SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
			IF (LP == 0) THEN
			  SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
			ELSE
			  SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
			 END IF
		     END IF
		   END DO
		   
                   DEN = A(L) + SUM_A_LPRIME + F(L,M)   
		   IF ( DEN .NE. ZERO) THEN
			 SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
			 SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
			            )/DEN
		   ENDIF
	         END DO
	         A_M(IJK,0,M)= SUM_A+A(M)
		 B_M(IJK,M)  = SUM_B+BB(M)
	       END DO
            ENDIF 
         END DO 

      RETURN  
      END SUBROUTINE PARTIAL_ELIM_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_U(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for X vector quantities
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
      SUBROUTINE PARTIAL_ELIM_U(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
      USE compar 
      USE drag
      USE fldvar       
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,LM,I,IJKE 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
! !	******GERA MODIFICATIONS************
      DOUBLE PRECISION SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
      INTEGER           L, M,LP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a(0:DIMENSION_M), BB(0:DIMENSION_M), F(0:DIMENSION_M,0:DIMENSION_M),&
                       Saxf(0:DIMENSION_M) 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
!      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      

!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, den) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3
               IF (FLOW_AT_E(IJK)) THEN 
                 IMJK = IM_OF(IJK) 
                 IJMK = JM_OF(IJK) 
                 IPJK = IP_OF(IJK) 
                 IJPK = JP_OF(IJK) 
		 F = ZERO
		 DO M=0, MMAX
                   IF (MOMENTUM_X_EQ(M)) THEN
		     A(M) = A_M(IJK,0,M)
		     BB(M)= B_M(IJK,M)

		     if (m .ne. 0) then
		       IF (MOMENTUM_X_EQ(0)) F(M,0)=-VXF(IJK,M)
		       F(0,M) = F(M,0)
		     end if
		     DO L =1, MMAX
		       IF ((L .NE. M) .AND. (M .NE. 0) .AND. MOMENTUM_X_EQ(L)) THEN
		         LM = FUNLM(L,M)	
		         IF (.NOT.IP_AT_E(IJK)) THEN 
                	   I = I_OF(IJK) 
               		   IJKE = EAST_OF(IJK) 

               		   F(M,L) = -AVG_X(F_SS(IJK,LM),F_SS(IJKE,LM),I)*VOL_U(IJK) 
            	         ENDIF 
		       END IF
		       F(L,M)=F(M,L)
		     END DO

		     IF (M == 0 ) THEN
			 SAXF(M) = -(A_M(IJK,E,M)*VAR_G(IPJK)+A_M(IJK,W,M)*VAR_G(IMJK&
                                )+A_M(IJK,N,M)*VAR_G(IJPK)+A_M(IJK,S,M)*VAR_G(IJMK)) 
		     ELSE
			 SAXF(M) = -(A_M(IJK,E,M)*VAR_S(IPJK,M)+A_M(IJK,W,M)*VAR_S(&
                          IMJK,M)+A_M(IJK,N,M)*VAR_S(IJPK,M)+A_M(IJK,S,M)*VAR_S(&
                          IJMK,M))
		     ENDIF

                     IF (DO_K) THEN 
                       IJKM = KM_OF(IJK) 
                       IJKP = KP_OF(IJK) 
		       IF (M ==0) THEN
                         SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_G(IJKP)+A_M(IJK,B,M)*&
                                  VAR_G(IJKM))
		       ELSE
                         SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_S(IJKP,M)+A_M(IJK,B,M)*&
                           VAR_S(IJKM,M)) 
		       ENDIF
                     ENDIF 
		   ENDIF
		   	 			
		 END DO
		     
   	         Do M=0,MMAX
		   
                   IF (MOMENTUM_X_EQ(M)) THEN
		     SUM_A=0.
		     SUM_B=0.0
		     do L =0,MMAX
		       SUM_A_LPRIME=0.0
		       SUM_B_LPRIME=0.0
		       DO LP=0,MMAX
			 IF ( LP .NE. M) THEN
			   SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
			   IF (LP == 0) THEN
			     SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
			   ELSE
			     SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
			   END IF
			 END IF
		       END DO
                       DEN = A(L) + SUM_A_LPRIME + F(L,M)   
		       IF ( DEN .NE. ZERO) THEN
			 SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
			 SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
			            )/DEN
		       ENDIF
		     END DO
		     A_M(IJK,0,M) = SUM_A+A(M)
		     B_M(IJK,M)  =  SUM_B+BB(M)
		     
		     DO L = 0, MMAX
		       IF(.NOT.MOMENTUM_X_EQ(L))THEN
		         IF( L /= M) THEN
			   IF(L == 0)THEN
		             A_M(IJK,0,M) = A_M(IJK,0,M) - VXF(IJK,M)
		             B_M(IJK,M)   = B_M(IJK,M)  - VXF(IJK,M) * VAR_G(IJK) 
			   ELSE
		             LM = FUNLM(L,M)	
		             IF (.NOT.IP_AT_E(IJK)) THEN 
                	       I = I_OF(IJK) 
               		       IJKE = EAST_OF(IJK) 

               		       F(M,L) = -AVG_X(F_SS(IJK,LM),F_SS(IJKE,LM),I)*VOL_U(IJK) 
		               A_M(IJK,0,M) = A_M(IJK,0,M) + F(M,L) 
		               B_M(IJK,M)   = B_M(IJK,M) + F(M,L) * VAR_S(IJK, L)   
            	             ENDIF 
                          ENDIF
			 ENDIF
		       ENDIF
		     ENDDO
		     
		   ENDIF
		 END DO
               ENDIF 
            END DO 

      RETURN  
      END SUBROUTINE PARTIAL_ELIM_U 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_V(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for Y vector quantities
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
      SUBROUTINE PARTIAL_ELIM_V(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
      USE compar 
      USE drag
      USE fldvar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, IJKN,J,LM
! 
!                      a0, b0 etc. 
!      DOUBLE PRECISION a(0:DIMENSION_M), BB(0:DIMENSION_M), F(0:DIMENSION_M,0:DIMENSION_M),&
!                       Saxf(0:DIMENSION_M) 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!	******GERA MODIFICATIONS************
      DOUBLE PRECISION SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
      INTEGER           L, M,LP 
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a(0:DIMENSION_M), BB(0:DIMENSION_M), F(0:DIMENSION_M,0:DIMENSION_M),&
                       Saxf(0:DIMENSION_M) 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
!      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

      
 
!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,  &
!$omp&  a, bb,F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, DEN) &
!$omp&  schedule(static)
		
            DO IJK = ijkstart3, ijkend3
              IF (FLOW_AT_N(IJK)) THEN 
                IMJK = IM_OF(IJK) 
                IJMK = JM_OF(IJK) 
                IPJK = IP_OF(IJK) 
                IJPK = JP_OF(IJK) 
		F = ZERO
		DO M = 0, MMAX
		
                  IF (MOMENTUM_Y_EQ(M)) THEN
		    A(M)=A_M(IJK,0,M)
		    BB(M)=B_M(IJK,M)

		    if (m .ne. 0) then
		      IF (MOMENTUM_Y_EQ(0)) F(M,0)=-VXF(IJK,M)
		      F(0,M)=F(M,0)
		    end if
		    DO L =1, MMAX
			 IF ((L .NE. M) .AND. (M .NE. 0).AND. MOMENTUM_Y_EQ(L)) THEN
			    LM = FUNLM(L,M)
			    IF (.NOT.IP_AT_N(IJK)) THEN 
			      J = J_OF(IJK) 
               		      IJKN = NORTH_OF(IJK) 
			      F(M,L)=-AVG_Y(F_SS(IJK,LM),F_SS(IJKN,LM),J)*VOL_V(IJK)
 			    END IF
			 END IF
			 F(L,M)=F(M,L)
		    END DO

		    IF (M == 0 ) THEN
		      SAXF(M) = -(A_M(IJK,E,M)*VAR_G(IPJK)+A_M(IJK,W,M)*VAR_G(IMJK&
                           )+A_M(IJK,N,M)*VAR_G(IJPK)+A_M(IJK,S,M)*VAR_G(IJMK)) 
		    ELSE
		      SAXF(M) = -(A_M(IJK,E,M)*VAR_S(IPJK,M)+A_M(IJK,W,M)*VAR_S(&
                          IMJK,M)+A_M(IJK,N,M)*VAR_S(IJPK,M)+A_M(IJK,S,M)*VAR_S(&
                          IJMK,M))
		    ENDIF 			
		  
                    IF (DO_K) THEN 
                      IJKM = KM_OF(IJK) 
                      IJKP = KP_OF(IJK) 
		      IF (M == 0) THEN
                        SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_G(IJKP)+A_M(IJK,B,M)*&
                           VAR_G(IJKM)) 
		      ELSE
                        SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_S(IJKP,M)+A_M(IJK,B,M)*&
                           VAR_S(IJKM,M))
		      ENDIF 
                    ENDIF 
		  ENDIF
		END DO
		
		Do M=0,MMAX
                  IF (MOMENTUM_Y_EQ(M)) THEN
		    SUM_A=0.
		    SUM_B=0.0
		    do L =0,MMAX
			 SUM_A_LPRIME=0.0
			 SUM_B_LPRIME=0.0
			 DO LP=0,MMAX
			   IF ( LP .NE. M) THEN
			     SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
			     IF (LP == 0) THEN
			        SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
			     ELSE
			        SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
			     END IF
			   END IF
			 END DO
			 DEN = A(L) + SUM_A_LPRIME + F(L,M)   
		         IF ( DEN .NE. ZERO) THEN
			    SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
			    SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
			            )/DEN
		         ENDIF
		    END DO
		    A_M(IJK,0,M)=SUM_A+A(M)
		    B_M(IJK,M) = SUM_B+BB(M)
		     
		     DO L = 0, MMAX
		       IF(.NOT.MOMENTUM_Y_EQ(L))THEN
		         IF( L /= M) THEN
			   IF(L == 0)THEN
		             A_M(IJK,0,M) = A_M(IJK,0,M) - VXF(IJK,M)
		             B_M(IJK,M)   = B_M(IJK,M)  - VXF(IJK,M) * VAR_G(IJK) 
			   ELSE
			     LM = FUNLM(L,M)
			     IF (.NOT.IP_AT_N(IJK)) THEN 
			       J = J_OF(IJK) 
               		       IJKN = NORTH_OF(IJK) 
			       F(M,L)=-AVG_Y(F_SS(IJK,LM),F_SS(IJKN,LM),J)*VOL_V(IJK)
		               A_M(IJK,0,M) = A_M(IJK,0,M) + F(M,L) 
		               B_M(IJK,M)   = B_M(IJK,M) + F(M,L) * VAR_S(IJK, L)   
            	             ENDIF 
                          ENDIF
			 ENDIF
		       ENDIF
		     ENDDO
		     
		  ENDIF
		END DO
              ENDIF 
            END DO 
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_V 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARTIAL_ELIM_W(Var_g, Var_s, VxF, A_m, B_m, IER)       C
!  Purpose: Do partial elimination for Z vector quantities
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
      SUBROUTINE PARTIAL_ELIM_W(VAR_G, VAR_S, VXF, A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE matrix 
      USE physprop
      USE indices
      USE run
      USE compar 
      USE drag
      USE fldvar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, IJKP,LM,K, IJKT 
! 
!                      gas phase variable 
      DOUBLE PRECISION Var_g(DIMENSION_3) 
! 
!                      solids phase variable 
      DOUBLE PRECISION Var_s(DIMENSION_3, DIMENSION_M) 
! 
!                      Volume x gas-solids transfer coefficient 
      DOUBLE PRECISION VxF(DIMENSION_3, DIMENSION_M) 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! ! !	******GERA MODIFICATIONS************
      DOUBLE PRECISION SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME, DEN
      INTEGER           L, M,LP
! 
!                      a0, b0 etc. 
      DOUBLE PRECISION a(0:DIMENSION_M), BB(0:DIMENSION_M), F(0:DIMENSION_M,0:DIMENSION_M),&
                       Saxf(0:DIMENSION_M) 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
!      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

 
!$omp  parallel do private( IMJK, IJMK, IJKM, IPJK, IJPK, IJKP, &
!$omp&  a, bb, F, Saxf,SUM_A, SUM_B, SUM_A_LPRIME,SUM_B_LPRIME,L, M,LP, DEN) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3 
               IF (FLOW_AT_T(IJK)) THEN 
                 IMJK = IM_OF(IJK) 
                 IJMK = JM_OF(IJK) 
                 IPJK = IP_OF(IJK) 
                 IJPK = JP_OF(IJK) 
		 
		 F = ZERO
		 DO M=0, MMAX
                   IF (MOMENTUM_Z_EQ(M)) THEN
		     A(M)=A_M(IJK,0,M)
		     BB(M)=B_M(IJK,M)
		   
		     if (m .ne. 0) then
			  IF (MOMENTUM_Z_EQ(0))F(M,0)=-VXF(IJK,M)
			  F(0,M)=F(M,0)
		     end if
		     DO L =1, MMAX
			IF ((L .NE. M) .AND. (M .NE. 0) .AND. MOMENTUM_Z_EQ(L) ) THEN
			  LM = FUNLM(L,M)
			  IF (.NOT.IP_AT_T(IJK)) THEN 
               		    K = K_OF(IJK) 
                            IJKT = TOP_OF(IJK) 
                            F(M,L) = -AVG_Z(F_SS(IJK,LM),F_SS(IJKT,LM),K)*VOL_W(IJK) 
            		  ENDIF 
			END IF
			F(L,M)=F(M,L)
		     END DO
		   
		     IF (M == 0 ) THEN
		       SAXF(M) = -(A_M(IJK,E,M)*VAR_G(IPJK)+A_M(IJK,W,M)*VAR_G(IMJK&
                          )+A_M(IJK,N,M)*VAR_G(IJPK)+A_M(IJK,S,M)*VAR_G(IJMK)) 
		     ELSE
		       SAXF(M) = -(A_M(IJK,E,M)*VAR_S(IPJK,M)+A_M(IJK,W,M)*VAR_S(&
                          IMJK,M)+A_M(IJK,N,M)*VAR_S(IJPK,M)+A_M(IJK,S,M)*VAR_S(&
                          IJMK,M))
		     ENDIF
		   
                     IF (DO_K) THEN 
                       IJKM = KM_OF(IJK) 
                       IJKP = KP_OF(IJK) 
		       IF (M == 0) THEN
                         SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_G(IJKP)+A_M(IJK,B,M)*&
                           VAR_G(IJKM)) 
		       ELSE
                         SAXF(M) = SAXF(M) - (A_M(IJK,T,M)*VAR_S(IJKP,M)+A_M(IJK,B,M)*&
                           VAR_S(IJKM,M)) 
		       ENDIF
                     ENDIF 
		   ENDIF			
		 END DO
		 
		 Do M=0,MMAX
                   IF (MOMENTUM_Z_EQ(M)) THEN
		     SUM_A=0.
		     SUM_B=0.0
		     do L =0,MMAX
			  SUM_A_LPRIME=0.0
			  SUM_B_LPRIME=0.0
			  DO LP=0,MMAX
			    IF ( LP .NE. M) THEN
			      SUM_A_LPRIME=SUM_A_LPRIME+F(L,LP)
			      IF (LP == 0) THEN
			        SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_G(IJK)
			      ELSE
			        SUM_B_LPRIME=SUM_B_LPRIME+F(L,LP)*VAR_S(IJK,LP)
			      END IF
			    END IF
			  END DO
			  
			  DEN = A(L) + SUM_A_LPRIME + F(L,M)   
		          IF ( DEN .NE. ZERO) THEN
			    SUM_A = SUM_A + ((F(L,M)*(A(L)+SUM_A_LPRIME))/DEN)
			    SUM_B = SUM_B + F(L,M)*(SAXF(L) + BB(L)+SUM_B_LPRIME &
			            )/DEN
		          ENDIF
		     END DO
		     A_M(IJK,0,M)=SUM_A+A(M)
		     B_M(IJK,M) = SUM_B+BB(M)
		     
		     DO L = 0, MMAX
		       IF(.NOT.MOMENTUM_Z_EQ(L))THEN
		         IF( L /= M) THEN
			   IF(L == 0)THEN
		             A_M(IJK,0,M) = A_M(IJK,0,M) - VXF(IJK,M)
		             B_M(IJK,M)   = B_M(IJK,M)  - VXF(IJK,M) * VAR_G(IJK) 
			   ELSE
			     LM = FUNLM(L,M)
			     IF (.NOT.IP_AT_T(IJK)) THEN 
               		       K = K_OF(IJK) 
                               IJKT = TOP_OF(IJK) 
                               F(M,L) = -AVG_Z(F_SS(IJK,LM),F_SS(IJKT,LM),K)*VOL_W(IJK) 
		               A_M(IJK,0,M) = A_M(IJK,0,M) + F(M,L) 
		               B_M(IJK,M)   = B_M(IJK,M) + F(M,L) * VAR_S(IJK, L)   
            	             ENDIF 
                          ENDIF
			 ENDIF
		       ENDIF
		     ENDDO
		     
		   ENDIF
		 END DO
               ENDIF 
            END DO 
      RETURN  
      END SUBROUTINE PARTIAL_ELIM_W 
      
  
