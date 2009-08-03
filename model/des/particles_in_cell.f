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
      INTEGER, DIMENSION(IMAX2,JMAX2,KMAX2):: particle_count

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
         IF(.NOT.PEA(L)) CYCLE

         IF(FIRST_PASS1) THEN 
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
               WRITE(*,*) '     Problem determining the solids',&
                  ' association for particle no: ', L
               WRITE(*,*) '     Particle diameter = ', 2.d0*DES_RADIUS(L), &
                  ' and dp =', D_P0(1:MMAX)
               WRITE(*,*) '     Particle density = ', &
                  Ro_Sol(L), 'and RO_S =', RO_S(1:MMAX)
               WRITE(*,*) '     Particle position =' , DES_POS_NEW(L,:)
            ENDIF

            DO I = IMIN1, IMAX3
               IF((DES_POS_NEW(L,1).GE.XE(I-1)).AND.(DES_POS_NEW(L,1).LT.XE(I))) THEN
                  PIJK(L,1) = I
                  EXIT                  
               ENDIF
            ENDDO

            DO J = JMIN1, JMAX3
               IF((DES_POS_NEW(L,2).GE.YN(J-1)).AND.(DES_POS_NEW(L,2).LT.YN(J))) THEN
                  PIJK(L,2) = J
                  EXIT
               ENDIF
            ENDDO

            IF(DIMN.EQ.2) THEN
               PIJK(L,3)  = 1
            ELSE
               DO K = KMIN1, KMAX3
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

            IF((DES_POS_NEW(L,1).GE.XE(I-1)).AND.(DES_POS_NEW(L,1).LT.XE(I)).OR.I.EQ.1) THEN
               
            ELSEIF((DES_POS_NEW(L,1).GE.XE(I)).AND.(DES_POS_NEW(L,1).LT.XE(I+1))) THEN 
               PIJK(L,1) = I+1
            ELSEIF((DES_POS_NEW(L,1).GE.XE(I-2)).AND.(DES_POS_NEW(L,1).LT.XE(I-1))) THEN 
               PIJK(L,1) = I-1
            ELSE
               IF (.NOT.DES_LOC_DEBUG) THEN
                  DES_LOC_DEBUG = .TRUE.
                  WRITE(*,1000)
               ENDIF                    
               WRITE(*,*) '     CHECK CELL I: ', I
               WRITE(*,*) '        Particle: ', L, ' Radius: ', &
                  DES_RADIUS(L)
               WRITE(*,*) '        X-POS: ', DES_POS_NEW(L,1), &
                  ' X-VEL: ', DES_VEL_NEW(L,1)
               STOP
            ENDIF

            IF((DES_POS_NEW(L,2).GE.YN(J-1)).AND.(DES_POS_NEW(L,2).LT.YN(J)).OR.J.EQ.1) THEN

            ELSEIF((DES_POS_NEW(L,2).GE.YN(J)).AND.(DES_POS_NEW(L,2).LT.YN(J+1))) THEN 
               PIJK(L,2) = J+1
            ELSEIF((DES_POS_NEW(L,2).GE.YN(J-2)).AND.(DES_POS_NEW(L,2).LT.YN(J-1))) THEN 
               PIJK(L,2) = J-1
            ELSE
               IF (.NOT.DES_LOC_DEBUG) THEN
                  DES_LOC_DEBUG = .TRUE.
                  WRITE(*,1000)
               ENDIF                    
               WRITE(*,*) '     CHECK CELL J: ', J
               WRITE(*,*) '        Particle: ', L, ' Radius: ', &
                  DES_RADIUS(L)
               WRITE(*,*) '        Y-POS: ', DES_POS_NEW(L,2), &
                  ' Y-VEL: ', DES_VEL_NEW(L,2)
               STOP
            ENDIF

            IF(DIMN.EQ.2) THEN
               PIJK(L,3) = 1
            ELSE
               IF((DES_POS_NEW(L,3).GE.ZT(K-1)).AND.(DES_POS_NEW(L,3).LT.ZT(K)).OR.K.EQ.1) THEN

               ELSEIF((DES_POS_NEW(L,3).GE.ZT(K)).AND.(DES_POS_NEW(L,3).LT.ZT(K+1))) THEN
                  PIJK(L,3) = K+1
               ELSEIF((DES_POS_NEW(L,3).GE.ZT(K-2)).AND.(DES_POS_NEW(L,3).LT.ZT(K-1))) THEN
                  PIJK(L,3) = K-1
               ELSE 
                  IF (.NOT.DES_LOC_DEBUG) THEN
                     DES_LOC_DEBUG = .TRUE.
                     WRITE(*,1000)
                  ENDIF                    
                  WRITE(*,*) '     CHECK CELL K: ', K
                  WRITE(*,*) '        Particle: ', L, ' Radius: ', &
                     DES_RADIUS(L)
                  WRITE(*,*) '        Z-POS: ', DES_POS_NEW(L,2), &
                     ' Z-VEL: ', DES_VEL_NEW(L,2)
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
            WRITE(*,*) '     WARNING IJK < 0  Particle: ', L,&
               ' Position: ', DES_POS_NEW(L,:)
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
! this should not matter if pure granular flow simulation (i.e. no fluid)
                  IF (.NOT.DES_LOC_DEBUG) THEN
                     DES_LOC_DEBUG = .TRUE.
                     WRITE(*,1000)
                  ENDIF
                  WRITE(*,*) '     WARNING EP_G LT zero at I,J: ',&
                     I_OF(IJK), J, '   where EP_S: ', EP_S(IJK,M)
                  !write(*,*) 'IJK = ', IJK, ' VOL(IJK) =  ', vol(IJK)
                  !write(*,*) 'Number of particles in cell = ', PINC(IJK)
                  !write(*,*) 'Eps calculated here = ', ROP_S(M)/RO_S(M)
                  !stop
               ENDIF 
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
! bed height calculations for each solids phase
               hcell = 0.5d0*(YN(J)+YN(J-1))
               tmp_num(M) = tmp_num(M) + EP_S(IJK,M)*hcell*VOL(IJK)
               tmp_den(M) = tmp_den(M) + EP_S(IJK,M)*VOL(IJK)
            ENDIF
         ENDDO   ! end loop over M=1,MMAX
      ENDDO     ! end loop over IJK=ijkstart3,ijkend3

! calculate bed height for each phase      
      bed_height(:) = tmp_num(:)/tmp_den(:)


      IF(DES_NEIGHBOR_SEARCH.EQ.4) THEN

         DO K = 1,KMAX2        !MAX(KMAX1-1,1)
            DO J = 1,JMAX2
               DO I = 1,IMAX2
                  IJK = FUNIJK(I,J,K)
                  NPIC = PINC(IJK)
               
                  IF (ASSOCIATED(pic(I,J,K)%p)) THEN
                     IF (npic.NE.SIZE(pic(I,J,K)%p)) THEN
                        DEALLOCATE(pic(I,J,K)%p)
                        IF (npic.GT.0) ALLOCATE(pic(I,J,K)%p(npic))
                     ENDIF
                  ELSE
                     IF (npic.GT.0) ALLOCATE(pic(I,J,K)%p(npic))
                  ENDIF
                  
               ENDDO
            ENDDO
         ENDDO
      
         particle_count(:,:,:) = 1
         
         PC = 1
         DO L = 1, MAX_PIS
            IF(PC .GT. PIS) EXIT
            IF(.NOT.PEA(L)) CYCLE

            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
            pos = particle_count(I,J,K)
            pic(I,J,K)%p(pos) = L
            particle_count(I,J,K) = particle_count(I,J,K) + 1

            PC = PC + 1
         ENDDO
      ENDIF

      FIRST_PASS1 = .FALSE.
      IF (DES_LOC_DEBUG) WRITE(*,1001)


 1000 FORMAT('---------- FROM PARTICLES_IN_CELL ---------->')
 1001 FORMAT('<---------- END PARTICLES_IN_CELL ----------') 

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL
