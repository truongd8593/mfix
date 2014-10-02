!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_CELL_MOVEMENT_DEM                                 !
!                                                                      !
!  Purpose: Check to see if DEM particles have moved into ghost cells. !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_DEM

      USE discretelement, only: MAX_PIP
      USE discretelement, only: PEA, PIJK
      USE discretelement, only: DES_POS_NEW
      USE discretelement, only: DES_VEL_NEW

      use mpi_utility
      USE error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER :: L
! ijk indices
      INTEGER :: I, J, K, IJK
! IER for error reporting
      INTEGER :: IER

      IER = 0

!$omp parallel default(shared)            &
!$omp private(L,I,J,K,IJK)
!$omp do reduction(+:IER) schedule (guided,50)
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1) .OR. any(PEA(L,2:4))) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1
      ENDDO
!$omp end parallel


      CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN


! Point of no return: Report errors and abort
!----------------------------------------------------------------------
      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_DEM")
      CALL OPEN_PE_LOG(IER)


      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Particles detected in a ghost cell:',/' ')
 
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1) .OR. any(PEA(L,2:4))) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(L)),'I',trim(iVal(I)),      &
               'X',DES_POS_NEW(1,L),'X',DES_VEL_NEW(1,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(L)),'J',trim(iVal(J)),      &
               'Y',DES_POS_NEW(2,L),'Y',DES_VEL_NEW(2,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF (DO_K .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(L)),'K',trim(iVal(K)),      &
               'Z',DES_POS_NEW(3,L),'Z',DES_VEL_NEW(3,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

 1101 FORMAT('Particle ',I10,' moved into cell with ',A,' index ',I6,/ &
         3x,A'-Position: ',g11.4,6x,A,'-Velocity:'g11.4,/' ')

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1102 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')


      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CHECK_CELL_MOVEMENT_DEM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_CELL_MOVEMENT_PIC                                 !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_PIC

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE cutcell
      USE mfix_pic
      USE des_rxns
      USE run
      USE error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase no.
      INTEGER M
! ijk indices
      INTEGER I, J, K, IJK, IPROC
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
! particle count in ijk fluid cell
      INTEGER, DIMENSION(DIMENSION_3):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

! count for number of particles that were found in the ghost cell and
! hence removed
      INTEGER :: PIP_DEL_COUNT
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)
! number of i,j,k nodes (parallelization)
      INTEGER :: NNODES
! IER for error reporting
      INTEGER IER
      INTEGER epg_min_loc(0:numpes-1, 4), epg_min_loc2(1)
      double precision :: epg_min(0:numpes-1), epg_min2
! Difference between a particles diameter (density) and the diameter
! (density) of a phase specified in the data file.
      DOUBLE PRECISION dDp, dRho

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

      CALL INIT_ERR_MSG("Particles_in_cell")

! following quantities are reset every call to particles_in_cell
      PIP_DEL_COUNT = 0
      PINC(:) = 0

! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index. If first_pass, also assigning PIJK(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>
!Handan Liu commented here on June 6 2013 and revised as below:
! 1)Directly adding the directives in this loop will result in the datarace;
! 2)Using Reduction-clause for PINC to avoid race condition for OpenMP
!   at 'PINC(IJK)=PINC(IJK)+1';
! 3)Using Reduction-clause for PIP_DEL_COUNT for OpenMP
!   at 'PIP_DEL_COUNT = PIP_DEL_COUNT + 1' for MPPIC.
! ----------------------------------------------------------------<<<<
!!$omp parallel do default(shared)                               &
!!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk) schedule (guided,50)
!!$      omp_start=omp_get_wtime()
!$omp parallel default(shared)            &
!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk,dDp,dRho)
!$omp do reduction(+:PINC) reduction(+:PIP_DEL_COUNT) schedule (guided,50)

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1)) CYCLE
! skipping ghost particles
         IF(PEA(L,2)) CYCLE 
         IF(PEA(L,3)) CYCLE
         IF(PEA(L,4)) CYCLE

! Assigning local aliases for particle position
         XPOS = DES_POS_NEW(1,L)
         YPOS = DES_POS_NEW(2,L)
         ZPOS = merge(0.0d0, DES_POS_NEW(3,L), NO_K)


! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

! this ijk is still an old value as it has not been updated
         IJK=PIJK(L,4)

! in MPPIC a particle can lie on the surface of the wall as only the
! centers are tracked.
         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            IF(I.EQ.IEND1+1 .AND. &
              (XPOS >= XE(IEND1-1) .AND. XPOS <= XE(IEND1)) )THEN
! This could happen if the cell adjacent to the ghost cell is a cut-cell
! and due to some tolerance issues, the particle is not detected outside
! the system. This will be a rare occurence, but it may occur and there
! is no point in stalling the simulation here.
! Implementing an easy fix for now: delete this particle.
! To add more stuff later
! 1. re distribute particle's weight among other particles in the domain
!    so that mass is conserved
! 2. rather than deactivating the particle, reflect the particle
!    inside the domain using the ghost cell bc's instead of cut-face bc

               IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'I',I,'X',&
                  XPOS,DES_POS_OLD(1,L),'X',DES_VEL_NEW(1,L)
               PIJK(L,1) = IEND1
            ELSE
               IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'I',I,'X',&
                  XPOS,DES_POS_OLD(1,L),'X',DES_VEL_NEW(1,L),&
                  DES_VEL_OLD(1,L), CUT_CELL_AT(IJK), &
                  FLUID_AT(IJK)
               PIP_DEL_COUNT = PIP_DEL_COUNT + 1
               PEA(L,1) = .false.
               CYCLE
            ENDIF
         ENDIF
         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            IF(J.EQ.JEND1+1.AND.&
              (YPOS >= YN(JEND1-1) .AND. YPOS <= YN(JEND1)) ) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'J',J,'Y',&
                  YPOS,DES_POS_OLD(2,L),'Y',DES_VEL_NEW(2,L)
               PIJK(L,2) = JEND1
            ELSE
               IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'J',J,'Y',&
                  YPOS,DES_POS_OLD(2,L),'Y',DES_VEL_NEW(2,L),&
                  DES_VEL_OLD(2,L), CUT_CELL_AT(IJK),&
                  FLUID_AT(IJK)
               PIP_DEL_COUNT = PIP_DEL_COUNT + 1
               PEA(L,1) = .false.
               CYCLE
            ENDIF
         ENDIF
         IF ( DO_K .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
            IF(K.EQ.KEND1+1.AND.&
              (ZPOS >= ZT(KEND1-1) .AND. ZPOS <= ZT(KEND1)) ) THEN
               IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'K',K,'Z',&
                  ZPOS, DES_POS_OLD(3,L),'Z',DES_VEL_NEW(3,L)
               PIJK(L,3) = KEND1
            ELSE
               IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'K',K,'Z',&
                  ZPOS,DES_POS_OLD(3,L),'Z',DES_VEL_NEW(3,L),&
                  DES_VEL_OLD(3,L), CUT_CELL_AT(IJK), &
                  FLUID_AT(IJK)
               PIP_DEL_COUNT = PIP_DEL_COUNT + 1
               PEA(L,1) = .false.
               CYCLE
            ENDIF
         ENDIF

      ENDDO   ! end loop over L = 1,particles
!!$omp end parallel do
!$omp end parallel

! ----------------------------------------------------------------<<<


      FIRST_PASS = .FALSE.
      CALL FINL_ERR_MSG

 1010     FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),/1X,A, &
         '-velocity new and old: ',2(ES17.9,4x),/&
         ' CUT_CELL and FLUID AT IJK_OLD ?', 2(L2,2x),/&
         ' Marking this particle as inactive',/&
          1X,70('*')/)

 1011 FORMAT(/1X,70('*')//,&
         ' From: PARTICLES_IN_CELL: Particle recovered',&
         ' from ghost cell -',/,&
         ' Message: Particle ',I8,' had moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,'-velocity: ',&
          ES17.9,/1X,70('*')/)


      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT_PIC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PARTICLES_IN_CELL                                       !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REPORT_PIC_STATS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE cutcell
      USE mfix_pic
      USE des_rxns
      USE run
      USE error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase no.
      INTEGER M
! ijk indices
      INTEGER I, J, K, IJK, IPROC
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
! particle count in ijk fluid cell
      INTEGER, DIMENSION(DIMENSION_3):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

! count for number of particles that were found in the ghost cell and
! hence removed
      INTEGER :: PIP_DEL_COUNT
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)
! number of i,j,k nodes (parallelization)
      INTEGER :: NNODES
! IER for error reporting
      INTEGER IER
      INTEGER epg_min_loc(0:numpes-1, 4), epg_min_loc2(1)
      double precision :: epg_min(0:numpes-1), epg_min2
! Difference between a particles diameter (density) and the diameter
! (density) of a phase specified in the data file.
      DOUBLE PRECISION dDp, dRho

!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE '../function.inc'
!-----------------------------------------------

      CALL INIT_ERR_MSG("PARTICLES_IN_CELL --> REPORT_PIC_STATS")

      PIP = PIP - PIP_DEL_COUNT
      LPIP_DEL_COUNT_ALL(:) = 0
      LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
      CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
      IF(SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
         WRITE(ERR_MSG,'(/,2x,A,2x,i10)') &
            'TOTAL NUMBER OF PARTICLES OUSIDE DOMAIN IN PIC = ', &
            SUM(LPIP_DEL_COUNT_ALL(:))
            CALL flush_err_msg(footer = .false.)

         DO IPROC = 0, NUMPES-1
            WRITE(ERR_MSG, '(/,A,i4,2x,A,2x,i5)') &
              'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,&
              ' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)

            CALL flush_err_msg(header = .false., footer = .false.)
         ENDDO

      IF(PIC_REPORT_MIN_EPG) then
         EPG_MIN(:)       = 0
         EPG_MIN_LOC(:,:) = 0

         EPG_MIN(mype)       = LARGE_NUMBER
         EPG_MIN_LOC(mype,:) = -1

         DO K = KSTART1, KEND1
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1
                  IJK = funijk(I,J,K)

                  IF(EP_G(IJK).lt.EPG_MIN(mype)) then
                     EPG_MIN_LOC(mype,1) = I
                     EPG_MIN_LOC(mype,2) = 2
                     EPG_MIN_LOC(mype,3) = 3
                     EPG_MIN_LOC(mype,4) = IJK
                     EPG_MIN(mype)       = EP_G(IJK)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         call global_all_sum(epg_min(0:numpes-1))
         call global_all_sum(epg_min_loc(0:numpes-1,1:4))

         epg_min2     = MINVAL(epg_min(0:numpes-1))
         epg_min_loc2 = MINLOC(epg_min(0:numpes-1)) - 1
         !-1, since minloc goes from 1:size of the array.
         !If not corrected by -1, then the proc id will be off by 1

         iproc = epg_min_loc2(1)

         I     = epg_min_loc(iproc, 1)
         J     = epg_min_loc(iproc, 2)
         K     = epg_min_loc(iproc, 3)
         IJK   = epg_min_loc(iproc, 4)
         WRITE(ERR_MSG,1014) EPG_MIN2, Iproc, I, J, K, IJK, &
            XE(I) - 0.5*DX(I), YN(J)-0.5*DY(J), ZT(K) - 0.5*DZ(K)

 1014       FORMAT( /, &
            &      5x,'EPGMIN                    = ', 2x,g17.8,/ &
            &      5x,'EPGMIN PROC RANK          = ', 2x, I10, / &
            &      5x,'EPGMIN (I, J, K, IJK)     = ', 3(2x,i5),2x,i10,/ &
            &      5x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8))

            call flush_err_msg(header = .false., footer = .false.)

         ENDIF
      end IF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE REPORT_PIC_STATS


