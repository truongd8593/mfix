!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_MPPIC                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MPPIC


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
!      use geometry, only: IMAX 
!      use geometry, only: JMAX
!      use geometry, only: KMAX
! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K


      USE param1
      USE geometry
      USE funits
      USE discretelement
      USE constant
      USE physprop
      USE fldvar
      USE toleranc 
      USE mfix_pic
      USE cutcell


      USE mpi_utility


! Global Parameters:
!---------------------------------------------------------------------//
!      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER :: M
      INTEGER :: I, J, K, IJK

! Number of real and comp. particles in a cell.
      DOUBLE PRECISION REAL_PARTS(DIM_M), COMP_PARTS(DIM_M)
! Volume of the cell 
      DOUBLE PRECISION :: VOLIJK, VOLIJK_UNCUT

!-----------------------------------------------
! Functions
!-----------------------------------------------
!-----------------------------------------------
! Include statement functions
!----------------------------------------------- 
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!----------------------------------------------- 



! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MPPIC")


        

      IF(MPPIC_COEFF_EN1 == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'MPPIC_COEFF_EN1'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(MPPIC_COEFF_EN1 > ONE .OR.                                &
         MPPIC_COEFF_EN1 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN1',                       &
            trim(iVal(MPPIC_COEFF_EN1))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_EN2 == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'MPPIC_COEFF_EN2'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(MPPIC_COEFF_EN2 > ONE .OR.                                &
         MPPIC_COEFF_EN2 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN2',                       &
            trim(iVal(MPPIC_COEFF_EN2))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_EN_WALL > ONE .OR.                                &
         MPPIC_COEFF_EN_WALL < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN_WALL',                   &
            trim(iVal(MPPIC_COEFF_EN_WALL))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_ET_WALL > ONE .OR.                                &
         MPPIC_COEFF_ET_WALL < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_ET_WALL',                   &
            trim(iVal(MPPIC_COEFF_ET_WALL))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')


      IF(MPPIC_CONSTANTNPC.AND.MPPIC_CONSTANTWT) THEN
         WRITE(ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: In MPPIC model, both MPPIC_CONSTANTNPC and ',&
         'MPPIC_CONSTANTWT',/'set to TRUE. Set at least one to false.',&
         ' See MFIX readme',/'Please correct the data file.')


      IF(.NOT.MPPIC_CONSTANTNPC.AND.(.NOT.MPPIC_CONSTANTWT)) THEN
         WRITE(ERR_MSG, 1101)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: In MPPIC model, both MPPIC_CONSTANTNPC and ',&
         'MPPIC_CONSTANTWT',/'are false. Set at least one to true. ',  &
         'See MFIX readme',/'Please correct the data file.')


      IF(MPPIC_CONSTANTNPC) then 
         DO M = 1, DES_MMAX
            IF(NPC_PIC(M).Eq.UNDEFINED_I) then 
               WRITE(ERR_MSG, 1000) trim(iVar('NPC_PIC',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF


      IF(MPPIC_CONSTANTWT) then 
         DO M = 1, DES_MMAX
            IF(STATWT_PIC(M).Eq.UNDEFINED) then 
               WRITE(ERR_MSG, 1000) trim(iVar('STATWT_PIC',M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDIF


! cnp_array(ijk, 0) will contain the cumulative number of real 
! particles later in the handling of inflow BC for MPPIC. See 
! the mppic_mi_bc in mppic_wallbc_mod.f 
      ALLOCATE(CNP_ARRAY(DIMENSION_3, 0:DIMENSION_M))
      CNP_ARRAY(:, :) = 0




      IF(GENER_PART_CONFIG) THEN
         ALLOCATE(RNP_PIC(DES_MMAX))
         ALLOCATE(CNP_PIC(DES_MMAX))
         RNP_PIC = ZERO
         CNP_PIC = ZERO 
         IF(DIMN.EQ.2) THEN 
! require that DZ(1)/ZLENGTH be specified for 2-dimensional case.  
! unclear why this cannot be one - other than the user may be unaware 
! that a depth has been set (a value of one implies default setting) 
            IF (DZ(1) == ONE) THEN
               WRITE(*,'(5X,A,A,/5X,A,A)') &
                  'For DIMN = 2, specify a value for DZ(1) or ',&
                  'ZLENGTH in mfix.dat which is not',&
                  'equal to one. If you want it to be one then ',&
                  'set it close to one but not exactly one'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

            IF (DZ(1) .NE. ZLENGTH) THEN 
               WRITE(ERR_MSG,'(5X,A,/5x,A,/5X,2(A20,2X,G17.8))')        &
                  'For DIMN = 2, DZ(1) and ZLENGTH are used ',         &
                  'interchangeably', ' Specify same values for ',      &
                  'DZ(1) and ZLENGTH','DZ(1) = ', DZ(1), 'ZLENGTH = ', &
                  ZLENGTH
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         DO K = KSTART1, KEND1 
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1 
                  IJK  = FUNIJK(I,J,K)
                  IF(.NOT.FLUID_AT(IJK)) CYCLE 
                  IF(EP_G(IJK).GE.1.d0-DIL_EP_s) CYCLE 
                  VOLIJK = VOL(IJK)
                  VOLIJK_UNCUT = DX(I)*DY(J)*DZ(K) 
                  DO M = 1, DES_MMAX
                     REAL_PARTS(M) = 6.d0*EP_S(IJK,M)*VOLIJK/&
                        (PI*(D_p0(M)**3.d0))
                     IF(MPPIC_CONSTANTNPC) THEN 
                        COMP_PARTS(M) = NPC_PIC(M)
                        IF(CUT_CELL_AT(IJK)) COMP_PARTS(M) = &
                         INT(VOLIJK*real(COMP_PARTS(M))/VOLIJK_UNCUT)
                     ELSEIF(MPPIC_CONSTANTWT) THEN
                        COMP_PARTS(M) = MAX(1,INT(REAL_PARTS(M)/REAL(STATWT_PIC(M))))
                     ENDIF
                        
                     RNP_PIC(M) = RNP_PIC(M) + REAL_PARTS(M)
                     CNP_PIC(M) = CNP_PIC(M) + COMP_PARTS(M)
                     CNP_ARRAY(IJK,M) = COMP_PARTS(M)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         PART_MPHASE(1:DES_MMAX) = CNP_PIC(1:DES_MMAX)
         PARTICLES = SUM(PART_MPHASE(1:DES_MMAX))
         CALL global_all_sum(PARTICLES)
      ENDIF !  end if(gener_part_config)


      CALL FINL_ERR_MSG

      RETURN  

      END SUBROUTINE CHECK_SOLIDS_MPPIC
