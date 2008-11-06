!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GENERATE_PARTICLE_CONFIG                               C
!  Purpose: Generate particle configuration based on input             C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 23-Oct-08  C
!  Comments: Added a new routine for clarity of functions              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop

      IMPLICIT NONE
      
      INTEGER LN, K, M, NP
      INTEGER CHECK_MPI
      INTEGER L, I, II, PART_COUNT
      DOUBLE PRECISION DIST, R_LM, DOML(DIMN)

            PART_COUNT = 0
            DO M = 1, MMAX
               DO NP = 1, PART_MPHASE(M) 
                  PART_COUNT = PART_COUNT + 1
                  DES_RADIUS(PART_COUNT) = D_P0(M)*HALF
                  RO_Sol(PART_COUNT) = RO_S(M)
               ENDDO
            ENDDO
            
            DOML(1) = DES_EPS_XSTART
            DOML(2) = DES_EPS_YSTART
            IF(DIMN.EQ.3) THEN
               DOML(3) = DES_EPS_ZSTART
            ENDIF
            
            CALL gener_lattice_mod(PARTICLES,doml(1:DIMN)&
            &,DES_POS_OLD(1:PARTICLES,1:DIMN),&
            & DES_RADIUS(1:PARTICLES))   

            OPEN(unit = 1000, file="particle_gener_conf_out.dat", form="formatted")
            
            DO LN = 1, PARTICLES
               OMEGA_OLD(LN,:) = ZERO
               DES_POS_NEW(LN,:) = DES_POS_OLD(LN,:)
               DES_VEL_NEW(LN,:) = DES_VEL_OLD(LN,:)
               OMEGA_NEW(LN,:) = OMEGA_OLD(LN,:)
               WRITE(1000, '(10(2x,g12.5))') (DES_POS_OLD(LN,K),K=1,DIMN),DES_RADIUS(LN),RO_Sol(LN) 
            ENDDO
            CLOSE(1000, STATUS = "KEEP")

            IF(MAXVAL(DES_POS_NEW(1:PARTICLES,2)).GT.YLENGTH-2.d0*MAXVAL(DES_RADIUS(1:PARTICLES))) THEN 
               WRITE(*,*)'---------------------------------------------------'
               WRITE(*,*)'WARNING MESSAGE'
               WRITE(*,*)'MAXIMUM POSITION OF PARTICLE IN Y DIRECTION (=', MAXVAL(DES_POS_NEW(1:PARTICLES,2)), ') > YLENGTH - DMAX =  ', YLENGTH-2.d0*MAXVAL(DES_RADIUS(1:PARTICLES))
               WRITE(*,*) 'REASON: THIS WOULD PROBABLY BE HAPPENING IF YOU ARE TRYING TO START WITH CLOSE PACKING'
               WRITE(*,*)'PROBLEM: POSITIVE OVERLAPS WITH WALLS OR PARTICLES FOR PERIODIC BCs'
               WRITE(*,*) 'RESOLUTION: INCREASE THE DOMAIN LENGTH IN Y- DIRECTION IF POSSIBLE OR GENERATE THE PARTICLE CONFIGURATION IN A BIGGER BOX AND SHRINK IT TO FIT IN THE DESIRED BOX'
               WRITE(*,*) 'STOPPING THE SIMULATION'
               STOP 
            ENDIF

      END SUBROUTINE GENERATE_PARTICLE_CONFIG
      
      SUBROUTINE GENER_LATTICE_MOD(NBODY,DOMLIN,XC,DBDY)
        
      USE discretelement
      implicit none
      double precision, INTENT(IN) ::  DOMLIN(DIMN)
      double precision, INTENT(INOUT) ::  DBDY(NBODY)
      INTEGER, INTENT(IN) :: NBODY
      DOUBLE PRECISION, INTENT(OUT) :: xc(nbody,DIMN)
      double precision :: DOML(DIMN), FAC, DMAX, XP, YP, ZP
      integer :: nprob1 , i, j, k, nx, ny, nz,&
      & np1, n
      
      DOML(:) = DOMLIN(:)
      dbdy(1:NBODY) = 2.d0*dbdy(1:NBODY)
      
      fac = 1.05
      dmax =  MAXVAL(dbdy(1:nbody))
      N = 1
      yp = dmax*fac/2.d0
      
      DO While (n.lt.nbody) 
         IF(DIMN.EQ.3) THEN 
            nz = floor((real(DOML(3))/(dmax*fac)))
         ELSE 
            nz = 1
         ENDIF
         nx = floor((real(doml(1))/(dmax*fac)))
         do k = 1, nz
            zp =  dmax*half*fac +  (k-1)*dmax*fac
            do i = 1, nx
               xp = dmax*half*fac + (i-1)*dmax*fac
               XC(n,1) = xp
               XC(n,2) = yp
               IF(DIMN.EQ.3) XC(n,3) = zp
               IF(N.Eq.NBODY) GOTO 200
               n = n+1
               
            end do
            
         end do
         yp = yp + dmax*fac
         
 200     CONTINUE       
      END DO
      PRINT*,'NUMBER OF PARTICLE IN GENER LATTICE MOD = ', N
      DBDY(1:NBODY) = 0.5d0*DBDY(1:NBODY)
      END SUBROUTINE GENER_LATTICE_MOD


      SUBROUTINE SET_INITIAL_VELOCITY
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      
      IMPLICIT NONE
      
      INTEGER IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IPJPK, IPJKP, IJPKP&
           &, IPJPKP
      INTEGER L, LL, I, J, K, KK, M, MM, IJK, IER

      DOUBLE PRECISION TEMP1, TEMP2, AVG_FACTOR, AVG_FACTOR2
      DOUBLE PRECISION UGC, VGC, WGC
      DOUBLE PRECISION OEPS, OVOL

      INTEGER, DIMENSION(3):: PCELL
      INTEGER:: ONEW,PART_IJK
      INTEGER:: IB, IE, JB, JE, KB, KE, II,JJ, IJK_EAST
      DOUBLE PRECISION:: EPS_P, VCELL, TEMP_U_G(JMAX2)
      INTEGER:: NP, IJK_PER, IJKP_PER
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
      PRINT*,'SETTING THE INITIAL PARTICLE VELOCITY EQUAL TO LOCAL FLUID VELOCITY'
      
      K = 1
      J = JMAX1
      
      DO I = 1, IMAX1
         IJK = funijk(I,J,K)
         IJK_EAST = funijk(I,JMAX2,K)
         PRINT*,'UG = ', U_G(IJK), U_G(IJK_EAST)
      ENDDO
      IF(INTX_PER) THEN 
         I = 1
         K = 1
         DO J = 1, JMAX2
            
            IJK = funijk(I,J,K)
            IJK_EAST = funijk(IMAX1, J, K)
            TEMP_U_G(J) = U_G(IJK) 
            U_G(IJK) = U_G(IJK_EAST)
         end DO
      ENDIF
      CALL set_interpolation_scheme(2)
      
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)
      AVG_FACTOR2 = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)
      
      DO NP = 1, PARTICLES
          ! CALCUALTE THE DRAG FORCE ON EACH PARTICLE USING THE PARTICLE VELOCITY
          I = PIJK(NP, 1)
          J = PIJK(NP, 2)
          K = PIJK(NP, 3)
          
          IJK = PIJK(NP, 4)
          M = PIJK(NP,5)
          
          PART_IJK = PINC(IJK)
          
          PCELL(1) = I-1
          PCELL(2) = J-1
          PCELL(3) = (3-DIMN)*1+(DIMN-2)*(K-1)
          
          CALL SET_INTERPOLATION_STENCIL(PCELL, IB, IE, JB, JE, KB,&
               &KE, INTERP_SCHEME, DIMN, ORDERNEW = ONEW) 
          DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
             DO J = 1, ONEW
                DO I = 1, ONEW
                   II = IB+I-1
                   JJ = JB+J-1
                   KK = KB+K-1
                   GSTENCIL(I,J,K,1) = XE(II)
                   GSTENCIL(I,J,K,2) = YN(JJ)
                   GSTENCIL(I,J,K,3) = ZT(KK)*(DIMN-2) + DZ(1)*(3-DIMN)
                   IF(II.LT.1.and.intx_per) II = IMAX1+II
                   IF(II.GT.IMAX1.and.intx_per) II = II-IMAX1
                   IF(JJ.LT.1.and.inty_per) JJ = JMAX1+JJ
                   IF(JJ.GT.JMAX1.and.inty_per) JJ = JJ-JMAX1
                   IF(KK.LT.1.and.intz_per) KK = KMAX1+KK
                   IF(KK.GT.KMAX1.and.intz_per) KK = KK-KMAX1

                   IJK = FUNIJK(II,JJ,KK)
                   
                   IPJK = IP_OF (IJK)    
                   IJPK = JP_OF (IJK)
                   IJKP = KP_OF (IJK)
                   IJPKP = KP_OF(IJPK)
                   IPJKP = KP_OF(IPJK)
                   IPJPK = JP_OF(IPJK)
                   IPJPKP = KP_OF(IPJPK)
                   
                   
                   vstencil(i,j,k,1) = AVG_FACTOR*( u_g(ijk) + u_g(ijpk)&
                   & + (u_g(ijkp) + u_g(ijpkp))*(dimn-2) )
                   
                   
                   vstencil(i,j,k,2) = AVG_FACTOR*( v_g(ijk) + v_g(ipjk)&
                   & + (v_g(ijkp) + v_g(ipjkp))*(dimn-2) )
                   
                   
                   if(dimn.eq.3) then 
                      vstencil(i,j,k,3) = AVG_FACTOR*(w_g(ijk) +&
                      & w_g(ijpk) + w_g(ipjk) + w_g(ipjpk) ) 
                   else 
                      vstencil(i,j,k,3) = 0.d0 !doesn't matter what
                                ! ever value is put here
                   end if
                   
                   
                END DO
             END DO
          END DO
          
          !PRINT*, 'VSTENCIL = ', VSTENCIL(1:ONEW, 1:ONEW, 1:ONEW, 1:DIMN)
          if(dimn.eq.2) then 
             CALL interpolator(gstencil(1:onew,1:onew,1,1:dimn)&
                  &,vstencil(1:onew,1:onew,1,1:2),DES_POS_NEW(np,&
                  & 1:2),vel_fp(np,1:2),onew&  
                  &,Interp_scheme,weightp)
             
          else 
             CALL interpolator(gstencil(1:onew,1:onew,1:onew,1:dimn)&
                  &,vstencil(1:onew,1:onew,1:onew,1:dimn)&
                  &,DES_POS_NEW(np, 1:3),vel_fp(np,1:3),onew&
                  &,interp_scheme,weightp)   
          end if
          !Print*,'Np = ', PCELL(:), vel_fp(np,1:dimn)
          
          DES_VEL_NEW(NP,:) = VEL_FP(NP,:)
          DES_VEL_OLD(NP,:) = VEL_FP(NP,:)
          !Print*,'Np = ', PCELL(:), DES_VEL_NEW(np,1:dimn)
          
       END DO
       
       IF(INTX_PER) THEN 
          I = 1
          K = 1
          
          DO J = 1, JMAX2
             IJK = funijk(I,J,K)
             U_G(IJK) = TEMP_U_G(J)
          end DO
       ENDIF
     END SuBROUTINE SET_INITIAL_VELOCITY
       
     SUBROUTINE writeic
       USE discretelement
       Implicit None 
       INTEGER:: ng, ip   
       OPEN(1002,file='ic.dat',form='formatted')
       DO ip = 1,particles 
          WRITE(1002,32)des_pos_new(ip,1:2), DES_VEL_NEW(IP,1:2), DES_VEL_OLD(IP,1:2)
       ENDDO
       CLOSE(1002)
32     FORMAT(10(2xe17.8))

       close(1002)
     END SUBROUTINE writeic
     


      SUBROUTINE init_particles_jn
      USE randomno
      USE discretelement
      USE constant 
      IMPLICIT NONE 
      INTEGER :: i,j, k 
      REAL*8 :: umf0(dimn), rsf(DIMN, DIMN)
      
      WRITE(*,*) 'INITIALIZING NORMAL VELOCITY DISTRIBUTION'
      WRITE(*,*) 'MEAN  = ', pvel_mean, ' AND Standard Deviation = ', PVEL_StDev
      do j=1,DIMN
         umf0(j)=pvel_mean
         do i=1,dimn
            if(i.eq.j)then
               rsf(i,j)=PVEL_StDev
            else
               rsf(i,j)=0.0
            endif
         enddo
      CALL nor_rno(DES_VEL_OLD(1:PARTICLES, j),umf0(j),rsf(j,j))
      enddo

      DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)
      END SUBROUTINE init_particles_jn

