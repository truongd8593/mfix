!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_DES_MPPIC                                         !
!  Purpose: Check user input data                                      !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!  Comments: Moved from CHECK_DES_MPPIC.                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DES_MPPIC

!-----------------------------------------------
! Modules 
!-----------------------------------------------      
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

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------      
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


! For MPPIC, the periodicity is based on the continuum variables such 
! as cyclic_x, cyclic_x_pd, etc.
         
      IF(CYCLIC_X.OR.CYCLIC_X_PD.OR.CYCLIC_Y.OR.CYCLIC_Y_PD.OR.&
         CYCLIC_Z.OR.CYCLIC_Z_PD ) THEN 
         
         IF(myPE.eq.pe_IO) WRITE(*,2005) 
         IF(DMP_LOG) WRITE(UNIT_LOG,2005) 
         
         DES_PERIODIC_WALLS = .TRUE. 
         DES_PERIODIC_WALLS_X = CYCLIC_X.OR.CYCLIC_X_PD
         DES_PERIODIC_WALLS_Y = CYCLIC_Y.OR.CYCLIC_Y_PD 
         DES_PERIODIC_WALLS_Z = CYCLIC_Z.OR.CYCLIC_Z_PD 
            
         IF(myPE.eq.pe_IO) WRITE(*,2006)  DES_PERIODIC_WALLS, &
         DES_PERIODIC_WALLS_X, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z
         IF(DMP_LOG) WRITE(UNIT_LOG,2006)  DES_PERIODIC_WALLS, &
         DES_PERIODIC_WALLS_X, DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z     

      ENDIF
! check particle-wall normal restitution coefficient
      IF(MPPIC_COEFF_EN_WALL == UNDEFINED) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2002) 'MPPIC_COEFF_EN_WALL'
         CALL MFIX_EXIT(myPE)
      ENDIF
         
         !IF(MPPIC_COEFF_ET_WALL == UNDEFINED) THEN
         !   if(dmp_log) WRITE (UNIT_LOG, 2002) 'MPPIC_COEFF_ET_WALL'
         !   CALL MFIX_EXIT(myPE)
         !ENDIF

      IF(MPPIC_COEFF_EN1 == UNDEFINED) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2004) 'MPPIC_COEFF_EN1'
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(MPPIC_COEFF_EN2 == UNDEFINED) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2004) 'MPPIC_COEFF_EN2'
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(MPPIC_COEFF_EN_WALL > ONE .OR. MPPIC_COEFF_EN_WALL < ZERO) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN_WALL'
         CALL MFIX_EXIT(myPE)
      ENDIF
         
      IF(MPPIC_COEFF_ET_WALL > ONE .OR. MPPIC_COEFF_ET_WALL < ZERO) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_ET_WALL'
         CALL MFIX_EXIT(myPE)
      ENDIF
        
      IF(MPPIC_COEFF_EN1 > ONE .OR. MPPIC_COEFF_EN1 < ZERO) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN1'
         CALL MFIX_EXIT(myPE)
      ENDIF

      IF(MPPIC_COEFF_EN2 > ONE .OR. MPPIC_COEFF_EN2 < ZERO) THEN
         if(dmp_log) WRITE (UNIT_LOG, 2003) 'MPPIC_COEFF_EN2'
         CALL MFIX_EXIT(myPE)
      ENDIF


      IF(MPPIC_SOLID_STRESS_SNIDER) WRITE(UNIT_LOG,*) &
         'USING THE SNIDER MODEL FOR SOLID STRESS AND THE ',&
         'INTEGRATION APPROACH'

      IF(MPPIC_SOLID_STRESS_SNIDER.AND.PRINT_DES_SCREEN) &
         WRITE(*,*) 'USING THE SNIDER MODEL FOR SOLID STRESS AND ',&
         'THE INTEGRATION APPROACH'
         

      IF(MPPIC_CONSTANTNPC.AND.MPPIC_CONSTANTWT) then 
         WRITE(UNIT_LOG, 2011)
         if(PRINT_DES_SCREEN) WRITE(*, 2011)
         call  mfix_exit(mype)
      ENDIF

      IF(.not.MPPIC_CONSTANTNPC.AND.(.not.MPPIC_CONSTANTWT)) then 
         WRITE(UNIT_LOG, 2012)
         if(PRINT_DES_SCREEN) WRITE(*, 2012)
         call  mfix_exit(mype)
      ENDIF

      IF(MPPIC_CONSTANTNPC) then 
         DO M = 1, DES_MMAX
            IF(NPC_PIC(M).Eq.UNDEFINED_I) then 
               WRITE(UNIT_LOG, 2013) M
               if(PRINT_DES_SCREEN) WRITE(*, 2013) M
               call  mfix_exit(mype)
            ENDIF
         ENDDO
      ENDIF
         
      IF(MPPIC_CONSTANTWT) then 
         DO M = 1, DES_MMAX
            IF(STATWT_PIC(M).Eq.UNDEFINED) then 
               WRITE(UNIT_LOG, 2014) M
               if(PRINT_DES_SCREEN) WRITE(*, 2014) M
               call  mfix_exit(mype)
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
               IF(mype.eq.pe_IO) WRITE(*,'(5X,A,A,/5X,A,A)') &
                  'For DIMN = 2, specify a value for DZ(1) or ',&
                  'ZLENGTH in mfix.dat which is not',&
                  'equal to one. If you want it to be one then ',&
                  'set it close to one but not exactly one'

               IF(DMP_LOG) WRITE(UNIT_LOG,'(5X,A,A,/5X,A,A)') &
                  'For DIMN = 2, specify a value for DZ(1) or ',&
                  'ZLENGTH in mfix.dat which is not',&
                  'equal to one. If you want it to be one then ',&
                  'set it close to one but not exactly one'
               CALL mfix_exit(mype)
            ENDIF

            IF (DZ(1) .NE. ZLENGTH) THEN 
!this condition will probably never occur. Redundancy doesn't hurt, however!

               IF(myPE.eq.pe_IO) THEN
                  WRITE(*,'(5X,A,/5x,A)') &
                     'For DIMN = 2, DZ(1) and ZLENGTH are used ',&
                     'interchangeably', ' Specify same values for ',&
                     'DZ(1) and ZLENGTH'  
                  WRITE(*,'(5X,2(A20,2X,G17.8))') &
                     'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
               ENDIF
               IF(DMP_LOG) THEN
                  WRITE(UNIT_LOG,'(5X,A,/5x,A)') &
                     'For DIMN = 2, DZ(1) and ZLENGTH are used ',&
                     'interchangeably', ' Specify same values for ',&
                     'DZ(1) and ZLENGTH'  
                  WRITE(UNIT_LOG,'(5X,2(A20,2X,G17.8))') &
                     'DZ(1) = ', DZ(1), 'ZLENGTH = ', ZLENGTH
               ENDIF
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


      RETURN

 2001 FORMAT(/1X,70('*')//' From: CHECK_DES_MPPIC',/' Message: ',&
         'Looks like a 3-D case (IMAX, JMAX & KMAX all >1) ',&
         'but DIMN equal to ',I1,/1X,70('*')/)

 2002    FORMAT(/1X,70('*')//' From: CHECK_DES_MPPIC',/' Message: ',&
         'Particle-wall restitution coefficient:  ',  A, &
         /10X,'not specified in mfix.dat for the MPPIC model',&
         /1X,70('*')/)

 2003    FORMAT(/1X,70('*')//' From: CHECK_DES_MPPIC',/' Message: ',&
         'Unphysical ( > 1 or < 0) values of:  ', A, &
         /1X,70('*')/)

 2004    FORMAT(/1X,70('*')//' From: CHECK_DES_MPPIC',/' Message: ',&
         'Friction restitution coefficient:  ',  A, &
         /10X,'not specified in mfix.dat for the MPPIC model',&
         /1X,70('*')/)
 2005    FORMAT(/1X,70('*')//' From: CHECK_DES_MPPIC',/' Message: ', &
         'MPPIC Case: One of the cyclic_flags found to be true', /10X, & 
         'from  continuum flags, so forcing des_periodic_walls to true')

 2006    FORMAT(/10X,'Updated Periodicity related flags:', /10X,&
         'DES_PERIODIC_WALLS   :', L4, /10X, &
         'DES_PERIODIC_WALLS_X :', L4, /10X, &
         'DES_PERIODIC_WALLS_Y :', L4, /10X, &
         'DES_PERIODIC_WALLS_Z :', L4, /10X, &
         /1X,70('*')/)

 2010 FORMAT(/1X,70('*')/' From: CHECK_DES_MPPIC',/'Error 2010:',       &
         ' DES simulations cannot have defined internal obstacles.',/  &
         ' Please correct the data file.', 1X,70('*')/)
         
 2011    FORMAT(/1X,70('*')/' From: CHECK_DES_MPPIC',/' Error 2011:',       &
         ' In MPPIC model, both MPPIC_CONSTANTNPC and MPPIC_CONSTANTWT', /,  & 
         ' set to TRUE. Set at least one to false. See MFIX readme', /,  &
         ' Please correct the data file. Exiting.', /,1X,70('*')/)

 2012    FORMAT(/1X,70('*')/' From: CHECK_DES_MPPIC',/' Error 2012:',       &
         ' In MPPIC model, both MPPIC_CONSTANTNPC and MPPIC_CONSTANTWT', /,  & 
         ' set to false. Set at least one to true. See MFIX readme', /,  &
         ' Please correct the data file. Exiting.', /,1X,70('*')/)

 2013    FORMAT(/1X,70('*')/' From: CHECK_DES_MPPIC',/' Error 2013:',       &
         ' In MPPIC model, MPPIC_CONSTANTNPC specified', &
         ' but NPC_PIC undefined', /,&
         ' for phase', I4, /, &
         ' See MFIX readme', /,  &
         ' Please correct the data file. Exiting.',/, 1X,70('*')/)

 2014    FORMAT(/1X,70('*')/' From: CHECK_DES_MPPIC',/' Error 2014:',       &
         ' In MPPIC model, MPPIC_CONSTANTWT specified', &
         ' but STATWT_PIC undefined', /,&
         ' for phase', I4, /, &
         ' See MFIX readme.', /,  &
         ' Please correct the data file. Exiting.',/, 1X,70('*')/)

 2015    FORMAT( 1X,70('*')/ & 
         'From: CHECK_DES_MPPIC    ', /5x, &
         'gener_part_config specified as true', /,5x, &
         'Total IC region count with non zero solids = ', I5,2X, /5x, & 
         'Total discrete solid phases = ', I5,2X, /5x, & 
         'Particle configuration will be automatically generated: ')

 2016    FORMAT(/1X,70('-')/, 5x, &
         'IC number and Volume        = ', I4, 2x, (ES15.7,2X) )
         
 2017    FORMAT(1X,70('.'),/,5x, &
         'PHASE Index, M              =  ', I5,2X, /5x, & 
         'D_P0(M)                     =  ', (ES15.7,2X), /5x, &
         'Solids vol fraction for M   =  ', (G15.8,2X), /5x, & 
         '# of particles for phase M  = ',  (I10,2X))

 2018    FORMAT( 1X,70('*')/)

         END SUBROUTINE CHECK_DES_MPPIC
