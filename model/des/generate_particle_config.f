!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 23-Oct-08  C
!  Comments: Added a new routine for clarity of functions              C
!                                                                      C
!  Revision: Modified subroutine for parallel processing               C
!  Authour : Pradeep G                                Date:28-Feb-11   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop
      use desmpi 
      use cdist 
      use mpi_utility
      use mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L, M
      INTEGER :: lproc_parcount
      double precision lmax_dia,lfac,xp,yp,zp 
!-----------------------------------------------

        
      IF(MPPIC) THEN 
         CALL GENERATE_PARTICLE_CONFIG_MPPIC
         RETURN 
      ENDIF

      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
         '---------- START GENERATE_PARTICLE_CONFIG ---------->'

! initializing particle count      
      lproc_parcount = 0
! setting a local value of maximum diameter      
      lmax_dia = maxval(DES_D_p0(1:DES_MMAX))
      lfac =1.05 
! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      yp = lmax_dia*0.5*lfac
      xp = lmax_dia*0.5*lfac
      zp = lmax_dia*0.5*lfac

      IF (DIMN .EQ. 2) THEN

! looping of specified number of solids phases
         DO M = 1, DES_MMAX
! looping over calculated number of particles in mth solids phase 
            DO L = 1, PART_MPHASE(M)
! if current seed location within simulation domain boundaries based on
! computational FLUID grid (xe,yn,zt)
               IF (XP.GE.XE(istart1-1) .AND. XP.LT.XE(iend1) .AND.&
                   YP.GE.YN(jstart1-1) .AND. YP.LT.YN(jend1) ) THEN
                  lproc_parcount = lproc_parcount + 1
                  pea(lproc_parcount,1) = .true.
                  des_radius(lproc_parcount) = DES_D_p0(M)*half 
                  ro_sol(lproc_parcount) = DES_RO_S(M)
                  des_pos_new(lproc_parcount,1) = xp 
                  des_pos_new(lproc_parcount,2) = yp 
               ENDIF 

               xp = xp + lmax_dia*lfac 
! if current x seed position exceeds specified bounds then increment y
! seed position up and reset x seed position
               IF (xp+lmax_dia*0.5*lfac .GT. DES_EPS_XSTART) THEN
                  xp = lmax_dia*0.5*lfac
                  yp = yp + lmax_dia*lfac 
               ENDIF
            ENDDO  ! end loop over number of particles in phase M
         ENDDO   ! end loop (m=1,des_mmax)

      ELSE   ! three dimensions

! looping of specified number of solids phases              
         DO M=1,DES_MMAX
! looping over calculated number of particles in mth solids phase          
            DO L=1,PART_MPHASE(M)
! if current seed location within simulation domain boundaries based on
! computational FLUID grid (xe,yn,zt)
               IF (XP.GE.XE(istart1-1) .AND. XP.LT.XE(iend1) .AND.&
                   YP.GE.YN(jstart1-1) .AND. YP.LT.YN(jend1) .AND.&
                   ZP.GE.ZT(kstart1-1) .AND. ZP.LT.ZT(kend1)) THEN
                  lproc_parcount = lproc_parcount + 1
                  pea(lproc_parcount,1) = .true.
                  des_radius(lproc_parcount) = DES_D_p0(M)*half 
                  ro_sol(lproc_parcount) = DES_RO_S(M)
                  des_pos_new(lproc_parcount,1) = xp 
                  des_pos_new(lproc_parcount,2) = yp 
                  des_pos_new(lproc_parcount,3) = zp 
               ENDIF

               xp = xp + lmax_dia*lfac 
! if current x seed position exceeds specified bounds then increment z
! seed position up and reset x seed position
               IF (xp+lmax_dia*0.5*lfac .GT. des_eps_xstart) THEN
                  xp = lmax_dia*0.5*lfac
                  zp = zp + lmax_dia*lfac 
! if current z seed position exceeds specified bounds then increment y
! seed position up and reset z seed position
                  IF (zp+lmax_dia*0.5*lfac .GT. des_eps_zstart) THEN
                     zp = lmax_dia*0.5*lfac
                     yp = yp + lmax_dia*lfac 
                  ENDIF
               ENDIF 
            ENDDO  ! end loop over number of particles in phase M
         ENDDO   ! end loop (m=1,des_mmax)
      ENDIF   ! end if/else (dimn==2)

! setting pip to particle count
      pip = lproc_parcount 

      if(pip > 0) then
         if(maxval(des_pos_new(1:pip,2)).gt.&
         ylength-2.d0*maxval(des_radius(1:pip))) then 
            write(unit_log,1002) maxval(des_pos_new(1:pip,2)), &
               ylength-2.d0*maxval(des_radius(1:pip))
            write(*,1003)
            call des_mpi_stop
         endif
      endif
      
      IF(DMP_LOG.and.debug_des) write(UNIT_LOG,'(3x,a)') &
         '<---------- END GENERATE_PARTICLE_CONFIG ----------'
      

 1002 FORMAT(/1X,70('*')//' From: GENERATE_PARTICLE_CONFIG',/,&
         ' Message: Positive overlap with walls in y-dir. Max. ',&
         'y-position of',/10X, 'particle (=', G12.5, &
         ') > YLENGTH-DMAX = ', G12.5,/10X, 'This may occur if',&
         'starting with close packing. Increase',/10X, 'the domain ',&
         'length in the y-dir or generate the particle',/10X,&
         'configuration in a bigger box and shrink it to fit',/10X,&
         'in the desired box size',/1X,70('*')/)

 1003 FORMAT(5X,'An error has occured see the *.LOG FILE for details')

      RETURN
      END SUBROUTINE GENERATE_PARTICLE_CONFIG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         C
!  Purpose: generates particle position distribution for MPPIC         C
!  Author: Rahul Garg                                 Date: 3-May-2011 C      
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE param 
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop
      USE fldvar 
      USE indices 
      USE randomno
      USE mfix_pic 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER :: L, I, J, K, M, IDIM , IJK
      INTEGER :: PART_COUNT, CNP_CELL_COUNT, IPCOUNT
      DOUBLE PRECISION :: DOML(DIMN), CORD_START(DIMN) 
      DOUBLE PRECISION :: REAL_PARTS(DIM_M), STAT_WT
      INTEGER :: LPROC_PARCOUNT
      DOUBLE PRECISION :: VOLIJK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RANDPOS
!-----------------------------------------------  
! Include statement functions
!-----------------------------------------------  
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------  
      
      IF(dmp_log.and.debug_des) WRITE(unit_log,'(3X,A)') &
      '---------- START GENERATE_PARTICLE_CONFIG MPPIC ---------->'
      PART_COUNT = 0
      
      DO K = KSTART1, KEND1 
         DO J = JSTART1, JEND1
            DO I = ISTART1, IEND1 
               
               IJK  = FUNIJK(I,J,K)
               CORD_START(1) = XE(I-1)
               CORD_START(2) = YN(J-1)
               IF(DIMN.EQ.3) CORD_START(3) = ZT(K-1)
               DOML(1) = DX(I)
               DOML(2) = DY(J)
               IF(DIMN.EQ.3) DOML(3) = DZ(K)
               DO M = 1, DES_MMAX 
                  CNP_CELL_COUNT = CNP_ARRAY(IJK, M)
                  IF(CNP_CELL_COUNT.EQ.0) CYCLE 
                  ALLOCATE(RANDPOS(CNP_CELL_COUNT, DIMN))
                  RANDPOS = ZERO 
                  DO IDIM = 1, DIMN
                     CALL UNI_RNO(RANDPOS(1:CNP_CELL_COUNT, IDIM))
                  ENDDO

                  VOLIJK = VOL(IJK)
                 
! this is based ep_s that would be set based on IC_EPS from
! two-fluid IC specification. 
! this should be replaced with a true dem quantity                 
                  REAL_PARTS(M) = 6.d0*EP_S(IJK,M)*VOLIJK/&
                     (PI*(DES_D_p0(M)**3.d0))
                  
                  IF(MPPIC_CONSTANTNPC) THEN 
! calculate the statistical weight for CP's belonging to this solids
! phase
                     STAT_WT = REAL_PARTS(M)/REAL(CNP_CELL_COUNT)
                  ELSEIF(MPPIC_CONSTANTWT) THEN
! although the weight was specified in the input file, the integer 
! number of CP's requires recalculating statistical weight. This 
! slightly different statistical weight will ensure that the initial 
! volume fraction is as inputted. If the input statwt_pic is used,
! then the initial volume fraction might be slightly less than the 
! input initial volume fraction

                     STAT_WT = REAL_PARTS(M)/REAL(CNP_CELL_COUNT)
                  ENDIF
                  
                  DO IPCOUNT = 1, CNP_CELL_COUNT
                        DES_POS_OLD(PART_COUNT + IPCOUNT, :) = CORD_START(:) + RANDPOS(IPCOUNT, :)*DOML(:)
                     
                        DES_POS_NEW(PART_COUNT + IPCOUNT, :) =  DES_POS_OLD(PART_COUNT + IPCOUNT, :)
                        
                        DES_VEL_NEW(PART_COUNT + IPCOUNT, :) =  ZERO
                        DES_VEL_OLD(PART_COUNT + IPCOUNT, :) =  ZERO
                        DES_RADIUS(PART_COUNT + IPCOUNT) = DES_D_p0(M)*HALF
                        RO_Sol(PART_COUNT + IPCOUNT) = DES_RO_S(M)
                     
                        DES_STAT_WT(PART_COUNT + IPCOUNT) = STAT_WT

                        MARK_PART(PART_COUNT + IPCOUNT) = 1
                        
                        PIJK(PART_COUNT + IPCOUNT,1) = I
                        PIJK(PART_COUNT + IPCOUNT,2) = J
                        PIJK(PART_COUNT + IPCOUNT,3) = K
                        PIJK(PART_COUNT + IPCOUNT,4) = IJK
                        PIJK(PART_COUNT + IPCOUNT,5) = M
                        
                        IF(DES_POS_NEW(PART_COUNT + IPCOUNT,2).LE.YLENGTH/2.d0) MARK_PART(PART_COUNT + IPCOUNT) = 0
                        !MARK_PART(PART_COUNT + IPCOUNT) = myPE 

                        PEA(PART_COUNT + IPCOUNT,1) = .true.
                  ENDDO
                  PART_COUNT = PART_COUNT + CNP_CELL_COUNT 
                  DEALLOCATE(RANDPOS)
               ENDDO
! set the cnp_array to zero. It will be used for handling inflow later                
               CNP_ARRAY(IJK,:) = 0 


            end DO
         end DO
      end DO

! setting pip to part_count 
      PIP = PART_COUNT 

      IF(DMP_LOG) WRITE(UNIT_LOG,'(2X, A,/,10X, A, i5, /,2X, A, i10,/)') &
      'In GENERATE_PARTICLE_CONFIG MPPIC ',  & 
      'FROM pe =', mype, & 
      'NUMBER OF PARCELS INITIATED ON THIS PROC = ', PIP
      
      WRITE(*,'(2X, A,/,10X, A, i5, /,2X, A, i10,/)') &
      'In GENERATE_PARTICLE_CONFIG MPPIC ',  & 
      'FROM pe =', mype, & 
      'NUMBER OF PARCELS INITIATED ON THIS PROC = ', PIP
      
      IF(PART_COUNT.NE.SUM(CNP_PIC(1:DES_MMAX))) THEN 
         IF(DMP_LOG) THEN 
            WRITE(UNIT_LOG,*) 'ERROR IN GENERATE_PARTICLE_CONFIG_MPPIC'
            WRITE(UNIT_LOG,*) &
          & 'NUMBER OF PARTICLES INITIALIZED (', PART_COUNT, '), NOT &
          & EQUAL TO EARLIER CALCULATED TOTAL NUMBER OF PARTICLES (', & 
          & SUM(CNP_PIC(1:DES_MMAX)), ' )'

            WRITE(UNIT_LOG,*) 'TERMINAL ERROR: STOPPING'
         ENDIF
         CALL mfix_exit(myPE) 
         
      END IF
      
      IF(DMP_LOG.and.debug_des) WRITE(UNIT_LOG,'(3X,A)') &
      '---------- END GENERATE_PARTICLE_CONFIG MPPIC ---------->'
      !CALL mfix_exit(mypE)
      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! Subroutine : writeic
! Purpose    : write the initial position and velocity of the particles 
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE writeic
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar 
      use mpi_utility 
      USE discretelement
      use desmpi
      use cdist 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER:: NP,lunit,li,lparcnt,lcurpar
      integer :: lproc,llocalcnt,lglocnt,lgathercnts(0:numpes-1)
      character(30) :: lfilename
      double precision,dimension(:,:),allocatable :: ltemparray
!-----------------------------------------------

! based on distributed or single IO open the input file
      lunit = 25 
      if (bdist_io) then 
         write(lfilename,'("ic_",I4.4,".dat")')mype
         open(unit=lunit, file=lfilename, form="formatted")
      else 
         if(mype.eq.pe_io) then  
            lfilename= "ic.dat"
            open(unit=lunit, file=lfilename, form="formatted")
         endif 
      endif  

! Write the information      
!----------------------------------------------------------------->>>      
      if (bdist_io) then 
         lparcnt = 1 
         do lcurpar = 1,pip
            if(lparcnt.gt.pip) exit
            if(.not.pea(lcurpar,1)) cycle 
            lparcnt = lparcnt+1
            if(pea(lcurpar,4)) cycle 
            write(lunit,'(10(2x,es12.5))') (des_pos_new(lcurpar,li),li=1,dimn),&
               (des_vel_new(lcurpar,li),li=1,dimn), des_radius(lcurpar),ro_sol(lcurpar)
         enddo

!-----------------------------------------------------------------<<<
      else   ! else branch if not bdist_IO
!----------------------------------------------------------------->>>

! set parameters required for gathering info at PEIO and write as single file   
         lglocnt = 10
         llocalcnt = pip - ighost_cnt 
         call global_sum(llocalcnt,lglocnt) 
         allocate (ltemparray(lglocnt,2*dimn+2)) 
         allocate (dprocbuf(llocalcnt),drootbuf(lglocnt))
         igath_sendcnt = llocalcnt 
         lgathercnts = 0
         lgathercnts(mype) = llocalcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 
         do li = 1,dimn 
            call des_gather(des_pos_new(:,li))
            if(mype.eq.pe_io) ltemparray(1:lglocnt,li) = drootbuf(1:lglocnt)
         end do 
         do li = 1,dimn 
            call des_gather(des_vel_new(:,li))
            if(mype.eq.pe_io)ltemparray(1:lglocnt,dimn+li) = drootbuf(1:lglocnt)
         end do 
         call des_gather(des_radius)
         if(mype.eq.pe_io)ltemparray(1:lglocnt,2*dimn+1) = drootbuf(1:lglocnt)
         call des_gather(ro_sol)
         if(mype.eq.pe_io)ltemparray(1:lglocnt,2*dimn+2) = drootbuf(1:lglocnt)
         if (mype.eq.pe_io) then  
            do lcurpar = 1,lglocnt
               write(lunit,'(10(2x,es12.5))') (ltemparray(lcurpar,li),li=1,2*dimn+2)
            enddo
         end if 
         deallocate(ltemparray)
         deallocate(dprocbuf,drootbuf)

      endif   ! end if/else bdist_io
!-----------------------------------------------------------------<<<

      if(bdist_io .or. mype.eq.pe_io) close(lunit)

      RETURN
      END SUBROUTINE writeic



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE init_particles_jn
      USE randomno
      USE discretelement
      USE constant 
      USE compar 
      USE geometry
      IMPLICIT NONE 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I, J, K, L
      LOGICAL FILE_EXIST
      REAL*8 :: umf0(dimn), rsf(DIMN, DIMN)
! Quantities for LE BC      
! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      
! shear rate
      DOUBLE PRECISION SHEAR_RATE
! distance between shear boundaries      
      DOUBLE PRECISION SHEAR_LENGTH
! temporary variable - velocity of particle based on position between
! shearing boundaries and shear rate
      DOUBLE PRECISION TMP_VEL
! indices for position and velocity
      INTEGER POSI, VELI
! number of total nodes
      INTEGER NNODES      
!-----------------------------------------------      

      WRITE(*,'(3X,A)') '---------- START INIT_PARTICLES_JN ---------->'
     
      IF (.NOT.DES_LE_BC) THEN
         WRITE(*,'(5X,A)') 'Initializing normal velocity distribution:'
         WRITE(*,'(5X,A,ES17.8,A,ES17.8)') 'mean = ', pvel_mean,&
            ' and standard deviation = ', PVEL_StDev
      ELSE
         SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)
         IF(SHEAR_DIR.EQ.'DUDY') THEN
            SHEAR_LENGTH = YLENGTH
            POSI = 2
            VELI = 1
         ELSEIF(SHEAR_DIR.EQ.'DWDY') THEN
            SHEAR_LENGTH = YLENGTH
            POSI = 2
            VELI = 3
         ELSEIF(SHEAR_DIR.EQ.'DVDX') THEN
            SHEAR_LENGTH = XLENGTH
            POSI = 1
            VELI = 2
         ELSEIF(SHEAR_DIR.EQ.'DWDX') THEN
            SHEAR_LENGTH = XLENGTH
            POSI = 1
            VELI = 3
         ELSEIF(SHEAR_DIR.EQ.'DUDZ') THEN
            SHEAR_LENGTH = ZLENGTH 
            POSI = 3
            VELI = 1
         ELSEIF(SHEAR_DIR.EQ.'DVDZ') THEN
            SHEAR_LENGTH = ZLENGTH 
            POSI = 3
            VELI = 2
         ENDIF  
         SHEAR_RATE =  (2.d0*DES_LE_REL_VEL)/SHEAR_LENGTH
         pvel_mean = 0.0d0
         PVEL_StDev = DABS(SHEAR_RATE)
         WRITE(*,'(5X,A,A)') 'Setting up velocity profile consistent',&
            'with shear'
         WRITE(*,'(5X,A,ES17.8,A,ES17.8)') 'mean = ', pvel_mean,&
            ' and standard deviation = ', PVEL_StDev
      ENDIF
!-----------------------------------------------      


      DO J=1,DIMN
         umf0(j)=pvel_mean
         DO I=1,DIMN
            IF(I.EQ.J) THEN
               rsf(I,J)=PVEL_StDev
            ELSE
               rsf(I,J)=0.0
            ENDIF
         ENDDO
         WRITE(*,'(5X,A,I5,2X,A)') 'FOR DIRECTION= ', J,&
            ' where (1=X,2=Y,3=Z):   '
         CALL nor_rno(DES_VEL_OLD(1:PARTICLES,J),umf0(J),rsf(J,J))
      ENDDO


! Adjust initial condition: change position and velocity according to
! shear direction and rate
      IF (DES_LE_BC) THEN
         DO L = 1, PARTICLES
            DES_VEL_OLD(L,VELI) = SHEAR_RATE*DES_POS_OLD(L,POSI)-&
               DES_LE_REL_VEL
         ENDDO
      ENDIF
      
      DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)

! The writing of files below needs to be updated for MPI case 
! updating/writing initial particle configuration files      
      NNODES = NODESI*NODESJ*NODESK
      IF(NNODES.GT.1)   RETURN

      IF (GENER_PART_CONFIG) THEN
         INQUIRE(FILE='particle_gener_conf.dat',exist=FILE_EXIST)
         IF (FILE_EXIST) THEN
            OPEN(UNIT=24,FILE='particle_gener_conf.dat',&
                 STATUS='REPLACE')
            DO L = 1, PARTICLES
               WRITE(24,'(10(X,ES12.5))')&
                  (DES_POS_OLD(L,K),K=1,DIMN), DES_RADIUS(L),&
                  RO_Sol(L), (DES_VEL_OLD(L,K),K=1,DIMN) 
            ENDDO
            CLOSE(24)
         ENDIF
      ELSE
         OPEN(UNIT=24,FILE='particle_input2.dat',&
              STATUS='REPLACE')
         DO L = 1, PARTICLES
            WRITE(24,'(10(X,ES15.5))')&
               (DES_POS_OLD(L,K),K=1,DIMN), DES_RADIUS(L),&
               RO_Sol(L), (DES_VEL_OLD(L,K),K=1,DIMN) 
         ENDDO
         CLOSE(24)
      ENDIF

      WRITE(*,'(3X,A)') '<---------- END INIT_PARTICLES_JN ----------'

      RETURN
      END SUBROUTINE init_particles_jn
