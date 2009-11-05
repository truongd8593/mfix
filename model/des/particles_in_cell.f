!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  PARTICLES_IN_CELL                                     C
!>
!!  Purpose: DES - Finding the fluid computational cell in which      
!!           a particle lies, to calculte void fraction and also       
!!           the volume averaged solids velocity of the cell            
!<
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C 
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Removed the separate volume definitions and added pic     C
!            array formulation and bed height caluclation.             C
!  Note:     PIC array is currently used only if cell linked list      C
!            search is used.                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTICLES_IN_CELL

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      USE sendrecv
      
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, I, J, K, M, MM 
      INTEGER IJK, IPJK, IJPK, IJKP
      DOUBLE PRECISION :: OVOL 
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,MMAX), OSOLVOL
      INTEGER:: npic, pos
      INTEGER, DIMENSION(DIMENSION_I,DIMENSION_J,DIMENSION_K):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Variables to calculate bed height of each solids phase
      DOUBLE PRECISION :: tmp_num(MMAX), tmp_den(MMAX), hcell 
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS1 = .TRUE.
! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
! Accounted for particles
      INTEGER PC
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'


      PINC(:) = 0
      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_LOC_DEBUG = .FALSE.

! Note : the quantities xe, zt cannot be readily replaced with the
! similar variables x_e, z_t in main mfix code as they are not the
! same.  also the variable y_n does not exist in main mfix code.
! each ijk loop starts at 2 and goes to max+2 (i.e., imin1=2,
! imax2=imax+2) 
      IF(FIRST_PASS1) THEN
         XE(1) = ZERO
         YN(1) = ZERO
         DO I = IMIN1, IMAX2
            XE(I) = XE(I-1) + DX(I)
         ENDDO
         DO J  = JMIN1, JMAX2
            YN(J) = YN(J-1) + DY(J)
         ENDDO
         IF(DIMN.EQ.3) THEN
            ZT(1) = ZERO
            DO K = KMIN1, KMAX2
               ZT(K) = ZT(K-1) + DZ(K)
            ENDDO
         ENDIF
      ENDIF
     
      PC = 1
      DO L = 1, MAX_PIS

         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(L,1)) CYCLE

         IF(FIRST_PASS1 ) THEN 
! Brute force technique to determine the particle locations in the Eulerian grid 

            DO M = 1, MMAX
               IF(ABS(2.0d0*DES_RADIUS(L)-D_P0(M)).LT.SMALL_NUMBER.AND. &
               ABS( RO_Sol(L)-RO_S(M)).LT.SMALL_NUMBER) THEN
               PIJK(L,5) = M 
               ENDIF
            ENDDO
            
            IF(PIJK(L,5).EQ.0) THEN
               IF (.NOT.DES_LOC_DEBUG) THEN
                  DES_LOC_DEBUG = .TRUE.
                  WRITE(*,1000) 
               ENDIF
               WRITE(*,'(5X,A,A,I)') &
                  'Problem determining the solids ',&
                  'association for particle: ',I
               WRITE(*,'(7X,A,(ES15.9))') &
                  'Particle position = ', DES_POS_NEW(L,:)
               WRITE(*,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle diameter = ', 2.0*DES_RADIUS(L),&
                  'and D_P0(1:MMAX)= ', D_P0(1:MMAX)
               WRITE(*,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle density = ', Ro_Sol(L), &
                  'and RO_s(1:MMAX) = ', RO_S(1:MMAX)
            ENDIF

            DO I = IMIN1, IMAX2
               IF((DES_POS_NEW(L,1).GE.XE(I-1)).AND.(DES_POS_NEW(L,1).LT.XE(I))) THEN
                  PIJK(L,1) = I
                  EXIT                  
               ENDIF
            ENDDO

            DO J = JMIN1, JMAX2
               IF((DES_POS_NEW(L,2).GE.YN(J-1)).AND.(DES_POS_NEW(L,2).LT.YN(J))) THEN
                  PIJK(L,2) = J
                  EXIT
               ENDIF
            ENDDO

            IF(DIMN.EQ.2) THEN
               PIJK(L,3)  = 1
            ELSE
               DO K = KMIN1, KMAX2
                  IF((DES_POS_NEW(L,3).GT.ZT(K-1)).AND.(DES_POS_NEW(L,3).LE.ZT(K))) THEN 
                     PIJK(L,3) = K
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

        
         ELSE     ! if not first_pass
! incremental approach to determine the new location of the particles  

            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)

            XPOS = DES_POS_NEW(L,1) 
            YPOS = DES_POS_NEW(L,2)
            IF (DIMN .EQ. 2) THEN
               ZPOS = ONE
            ELSE
               ZPOS = DES_POS_NEW(L,3)
            ENDIF

            IF( XPOS < XE(1) .OR. XPOS > XE(IMAX1)) THEN 
               IF(.NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! if particle is not new or exiting it shouldn't be in a ghost cell  
                  WRITE(*,1007) L,'I',I,'X',XPOS,'X',DES_VEL_NEW(L,1)
                  STOP
               ELSEIF(XPOS < XE(1)) THEN
                  PIJK(L,1) = 1
               ELSEIF(XPOS > XE(IMAX1)) THEN
                  PIJK(L,1) = IMAX1+1
               ENDIF
            ELSEIF(XPOS >= XE(1) .AND. XPOS < XE(IMIN1)) THEN
               PIJK(L,1) = IMIN1
            ELSEIF (XPOS >= XE(IMAX1-1) .AND. XPOS < XE(IMAX1)) THEN
               PIJK(L,1) = IMAX1
            ELSEIF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
               PIJK(L,1) = I
            ELSEIF(XPOS >= XE(I) .AND. XPOS < XE(I+1)) THEN 
               PIJK(L,1) = I+1
            ELSEIF(XPOS >= XE(I-2) .AND. XPOS < XE(I-1)) THEN 
               PIJK(L,1) = I-1
            ELSE
               WRITE(*,1008) L,'I','X',XPOS,'X',DES_VEL_NEW(L,1)
               STOP
            ENDIF

            IF(YPOS < YN(1) .OR. YPOS > YN(JMAX1)) THEN
               IF(.NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
                  WRITE(*,1007) L,'J',J,'Y',YPOS,'Y',DES_VEL_NEW(L,2)
                  STOP
               ELSEIF(YPOS < YN(1)) THEN
                  PIJK(L,2) = 1
               ELSEIF(YPOS > YN(JMAX1)) THEN
                  PIJK(L,2) = JMAX1+1
               ENDIF
            ELSEIF(YPOS >= YN(1) .AND. YPOS < YN(JMIN1)) THEN
               PIJK(L,2) = JMIN1
            ELSEIF (YPOS >= YN(JMAX1-1) .AND. YPOS < YN(JMAX1)) THEN
               PIJK(L,2) = JMAX1
            ELSEIF(YPOS >= YN(J-1) .AND. YPOS < YN(J)) THEN 
               PIJK(L,2) = J
            ELSEIF(YPOS >= YN(J) .AND. YPOS < YN(J+1)) THEN 
               PIJK(L,2) = J+1
            ELSEIF(YPOS >= YN(J-2) .AND. YPOS < YN(J-1)) THEN 
               PIJK(L,2) = J-1
            ELSE
               WRITE(*,1008) L,'J','Y',YPOS,'Y',DES_VEL_NEW(L,2)
               STOP
            ENDIF

            IF(DIMN.EQ.2) THEN
               PIJK(L,3) = 1
            ELSE
               IF(ZPOS < ZT(1) .OR. ZPOS > ZT(KMAX1)) THEN
                  IF(.NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
                     WRITE(*,1007) L,'K',K,'Z',ZPOS,'Z',DES_VEL_NEW(L,3)
                     STOP
                  ELSEIF(ZPOS < ZT(1)) THEN
                     PIJK(L,3) = 1
                  ELSEIF(ZPOS > ZT(KMAX1)) THEN
                     PIJK(L,3) = KMAX1+1
                  ENDIF
               ELSEIF(ZPOS >= ZT(1) .AND. ZPOS < ZT(KMIN1)) THEN
                  PIJK(L,3) = KMIN1
               ELSEIF(ZPOS >= ZT(KMAX1-1) .AND. ZPOS < ZT(KMAX1)) THEN
                  PIJK(L,3) = KMAX1
               ELSEIF(ZPOS >= ZT(K-1) .AND. ZPOS < ZT(K)) THEN
                  PIJK(L,3) = K
               ELSEIF(ZPOS >= ZT(K) .AND. ZPOS < ZT(K+1)) THEN
                  PIJK(L,3) = K+1
               ELSEIF(ZPOS >= ZT(K-2) .AND. ZPOS < ZT(K-1)) THEN
                  PIJK(L,3) = K-1
               ELSE
                  WRITE(*,1008) L,'K','Z',ZPOS,'Z',DES_VEL_NEW(L,3)
                  STOP
               ENDIF
            ENDIF          

         ENDIF       ! end of (if first_pass)

         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         IF(IJK.LT.0) THEN
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000) 
            ENDIF                    
            WRITE(*,'(5X,A,I,/,7X,A,(ES15.9))') &
               'WARNING IJK < 0 for particle: ', L,&
               'Particle position = ', DES_POS_NEW(L,:)            
         ENDIF
         PIJK(L,4) = IJK
         PINC(IJK) = PINC(IJK) + 1
         M = PIJK(L,5)
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) +  PVOL(L)
         DES_U_s(IJK,M) = DES_U_s(IJK,M) + PVOL(L)*DES_VEL_NEW(L,1)
         DES_V_s(IJK,M) = DES_V_s(IJK,M) + PVOL(L)*DES_VEL_NEW(L,2)
         IF(DIMN.EQ.3) DES_W_s(IJK,M) = DES_W_s(IJK,M) + PVOL(L)*DES_VEL_NEW(L,3)

         PC = PC + 1         
      ENDDO      ! end loop over L = 1,particles


! calculate the cell average solids velocity, the bulk density (if not
! des_interp_on and not first_pass1), the void fraction, and average
! height of each solids phase
! ------------------------------------------------------------
      tmp_num(:) = ZERO 
      tmp_den(:) = ZERO 
      DO IJK = IJKSTART3, IJKEND3
         J = J_OF(IJK)
         EP_G(IJK) = ONE   
         DO M = 1, MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)   
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DIMN.EQ.3) THEN
                  DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
               ENDIF
            ENDIF
            IF(VOL(IJK).GT.0) THEN 
               OVOL = ONE/(VOL(IJK))
               IF(FIRST_PASS1 .OR. &
                 ((.NOT.FIRST_PASS1).AND.(.NOT.DES_INTERP_ON))) THEN
                  ROP_S(IJK,M)  = RO_S(M)*SOLVOLINC(IJK,M)*OVOL
               ENDIF
            ENDIF
            IF(ROP_S(IJK,M) > ZERO) THEN
               EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
               IF(EP_G(IJK).LT.ZERO .AND. DES_CONTINUUM_COUPLED) THEN 
! this does not matter if pure granular flow simulation (i.e. no fluid)
                  IF (.NOT.DES_LOC_DEBUG) THEN
                     DES_LOC_DEBUG = .TRUE.
                     WRITE(*,1000)
                  ENDIF
          WRITE(*,'(5X,A,I,/,7X,A,I,2X,I,2X,A,ES15.9,/,7X,A,I)') &
                    'WARNING EP_G LT zero at IJK: ', IJK,&
                    'I,J = ', I_OF(IJK), J, ' EP_S = ', EP_S(IJK,M), & 
                    'No. of particles in cell = ', PINC(IJK)
               ENDIF 
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
! bed height calculations for each solids phase
               hcell = 0.5d0*(YN(J)+YN(J-1))
               tmp_num(M) = tmp_num(M) + EP_S(IJK,M)*hcell*VOL(IJK)
               tmp_den(M) = tmp_den(M) + EP_S(IJK,M)*VOL(IJK)
            ENDIF
         ENDDO   ! end loop over M=1,MMAX
      ENDDO     ! end loop over IJK=ijkstart3,ijkend3

! calculate avg height for each phase
      bed_height(:) = tmp_num(:)/tmp_den(:)


! in this section the variable pic(i,j,k)%p(:) is assigned; for each
! cell calculate number of particles in cell and store their id's 
! ------------------------------------------------------------

      DO IJK = ijkstart3, ijkend3
! check all cells; update entering/exiting particle regions
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         npic = PINC(IJK)

! check to see if current number of particles in the cell has changed
! and change size of p accordingly
         IF (ASSOCIATED(pic(I,J,K)%p)) THEN
            IF (npic.NE.SIZE(pic(I,J,K)%p)) THEN
               DEALLOCATE(pic(I,J,K)%p)
               IF (npic.GT.0) ALLOCATE(pic(I,J,K)%p(npic))
            ENDIF
         ELSE
            IF (npic.GT.0) ALLOCATE(pic(I,J,K)%p(npic))
         ENDIF
      ENDDO

      particle_count(:,:,:) = 1

      PC = 1
      DO L = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(L,1)) CYCLE

         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         pos = particle_count(I,J,K)
         pic(I,J,K)%p(pos) = L
         particle_count(I,J,K) = particle_count(I,J,K) + 1
         PC = PC + 1
      ENDDO


      FIRST_PASS1 = .FALSE.
      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(3X,'---------- FROM PARTICLES_IN_CELL ---------->')
 1001 FORMAT(3X,'<---------- END PARTICLES_IN_CELL ----------') 

 1007 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/&         
         ' Message: Particle ',I4,' still found in',&
         ' ghost cell with index ',A,': ',I4,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/&
         1X,70('*')/)

 1008 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/&
         ' Message: Check particle ',I4,' in cell ',A,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/&
         1X,70('*')/)         


      RETURN
      END SUBROUTINE PARTICLES_IN_CELL
