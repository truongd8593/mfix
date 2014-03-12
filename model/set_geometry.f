!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_GEOMETRY                                           !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_GEOMETRY 

      USE param 
      USE param1 
      USE run
      USE geometry
      USE compar

      use bc, only: Flux_g

      use mpi_utility

      use error_manager

      IMPLICIT NONE

! Generic loop indices 
      INTEGER :: I, J, K 
! X-direction dimension of U-momentum cell 
      DOUBLE PRECISION :: DX_E
! Y-direction dimension of V-momentum cell 
      DOUBLE PRECISION :: DY_N
! Z-direction dimension of W-momentum cell 
      DOUBLE PRECISION :: DZ_T
! Error Flag
      INTEGER :: IER

      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------


      CALL INIT_ERR_MSG("SET_GEOMETRY")


! Set DIMENSION_X variables to MAX3 values.
      DIMENSION_I   = IMAX3
      DIMENSION_J   = JMAX3
      DIMENSION_K   = KMAX3

! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      ALLOCATE( X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX    (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX_E  (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX_E (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( oDY   (0:DIMENSION_J), STAT=IER )
      ALLOCATE( oDY_N (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( Z     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( Z_T   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ_T (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FX     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FX_E     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_E_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FY_N     (0:DIMENSION_J), STAT=IER )
      ALLOCATE( FY_N_bar (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FZ_T     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( FZ_T_bar (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500


! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')


!  Determine the cyclic direction with a specified mass flux
      CYCLIC_X_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_X_PD)
      CYCLIC_Y_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_Y_PD)
      CYCLIC_Z_MF = (FLUX_G /= UNDEFINED .AND. CYCLIC_Z_PD)

! Force the cyclic flag if cyclic with pressure drop.
      IF (CYCLIC_X_PD) CYCLIC_X = .TRUE. 
      IF (CYCLIC_Y_PD) CYCLIC_Y = .TRUE. 
      IF (CYCLIC_Z_PD) CYCLIC_Z = .TRUE. 
      CYCLIC = CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z 

      IF(CYLINDRICAL .AND. COMPARE(ZLENGTH,8.D0*ATAN(ONE)) .AND. DO_K) &
         CYCLIC_Z = .TRUE. 

      IF(CYCLIC_X) THEN 
         DX(1) = DX(IMAX1) 
         DX(IMAX2) = DX(IMIN1)
         IF(NODESI.NE.1) THEN
            DX(IMIN3) = DX(IMAX1-1) 
            DX(IMAX3) = DX(IMIN1+1) 
         ENDIF
      ENDIF

      IF(CYCLIC_Y) THEN
         DY(1) = DY(JMAX1) 
         DY(JMAX2) = DY(JMIN1)
         IF(NODESJ.NE.1) THEN
            DY(JMIN3) = DY(JMAX1-1) 
            DY(JMAX3) = DY(JMIN1+1) 
         ENDIF
      ENDIF 

      IF (CYCLIC_Z) THEN 
         DZ(1) = DZ(KMAX1) 
         DZ(KMAX2) = DZ(KMIN1) 
         IF(NODESK.NE.1) THEN
            DZ(KMIN3) = DZ(KMAX1-1)
            DZ(KMAX3) = DZ(KMIN1+1)
         ENDIF
      ENDIF 


! Initialize the X-axis variables for CYLINDRICAL coordinates.
      IF(CYLINDRICAL) THEN 
         X(IMIN3:IMAX3)    = UNDEFINED 
         X_E(IMIN3:IMAX3)  = UNDEFINED 
         OX(IMIN3:IMAX3)   = UNDEFINED 
         OX_E(IMIN3:IMAX3) = UNDEFINED 
         ODX(IMIN3:IMAX3)  = UNDEFINED 

         IF(XMIN == ZERO) THEN 
            ODX(1) = ONE/DX(1) 
            OX(1) = UNDEFINED 
            OX_E(1) = UNDEFINED 
            IF (DO_I) THEN 
               X(1) = -HALF*DX(1) 
               X_E(1) = 0.0 
            ELSE 
               X(1) = HALF*DX(1) 
               X_E(1) = DX(1) 
            ENDIF 
         ELSE 
            IF (DO_I) THEN 
               X_E(1) = XMIN 
               X(1) = XMIN - HALF*DX(1) 
            ELSE 
               X_E(1) = XMIN + DX(1) 
               X(1) = XMIN + HALF*DX(1) 
            ENDIF 
            OX(1) = ONE/X(1) 
            OX_E(1) = ONE/X_E(1) 
            ODX(1) = ONE/DX(1) 
         ENDIF 

         IF (DO_I) THEN 
            DO I = IMIN1, IMAX2 
               X(I) = X(I-1) + (DX(I-1)+DX(I))/2. 
               X_E(I) = X_E(I-1) + DX(I) 
               OX(I) = ONE/X(I) 
               OX_E(I) = ONE/X_E(I) 
               ODX(I) = ONE/DX(I) 
            END DO 
         ENDIF

! Initialize the X-axis variables for CARTESIAN coordinates.
      ELSE

         X(IMIN3:IMAX3) = ONE 
         X_E(IMIN3:IMAX3) = ONE 
         OX(IMIN3:IMAX3) = ONE 
         OX_E(IMIN3:IMAX3) = ONE
         ODX(IMIN3:IMAX3) = ONE/DX(IMIN3:IMAX3) 
      ENDIF 


! Initialize the Y-Axis variables.
      ODY(JMIN3:JMAX3) = ONE/DY(JMIN3:JMAX3)


! Initialize the Z-Axis variables.
      DO K = 1, KMAX3
         IF (K == 1) THEN 
            Z(K) = ZERO - HALF*DZ(K) 
            Z_T(K) = ZERO 
            IF(NODESK.NE.1) THEN
               Z(K-1) =Z_T(K) - HALF*DZ(K-1)
               Z_T(K-1) = Z_T(K) - DZ(K-1)
            ENDIF
	      ELSE
            Z(K) = Z_T(K-1) + HALF*DZ(K) 
            Z_T(K) = Z_T(K-1) + DZ(K)
         ENDIF
         ODZ(K) = ONE/DZ(K)

         IF(NODESK.NE.1 .AND. K==1) ODZ(K-1) = ONE/DZ(K-1)
      END DO


      DX_E = HALF*(DX(1)+DX(IMIN1)) 
      DY_N = HALF*(DY(1)+DY(JMIN1)) 
      DZ_T = HALF*(DZ(1)+DZ(KMIN1)) 

      ODX_E(1) = ONE/DX_E
      ODY_N(1) = ONE/DY_N
      ODZ_T(1) = ONE/DZ_T

      FX(1) = HALF 
      FX_BAR(1) = HALF 
      FX_E(1) = HALF 
      FX_E_BAR(1) = HALF 
      FY_N(1) = HALF 
      FY_N_BAR(1) = HALF 
      FZ_T(1) = HALF 
      FZ_T_BAR(1) = HALF 

      IF(NODESI.NE.1) ODX_E(IMIN3) = ONE/DX_E
      IF(NODESJ.NE.1) ODY_N(JMIN3) = ONE/DY_N
      IF(NODESK.NE.1) ODZ_T(KMIN3) = ONE/DZ_T

      IF(NODESI.NE.1) THEN
         FX(IMIN3) = HALF
         FX_BAR(IMIN3) = HALF
         FX_E(IMIN3) = HALF
         FX_E_BAR(IMIN3) = HALF
      ENDIF

      IF(NODESJ.NE.1) THEN
         FY_N(JMIN3) = HALF
         FY_N_BAR(JMIN3) = HALF
      ENDIF

      IF(NODESK.NE.1) THEN
         FZ_T(KMIN3) = HALF
         FZ_T_BAR(KMIN3) = HALF
      ENDIF


! Look at 2 through IMAX1 U-momentum cells
      IF (DO_I) THEN 
         DO I = IMIN1, IMAX1 
            DX_E = HALF*(DX(I+1)+DX(I)) 
            ODX_E(I) = ONE/DX_E 
            FX(I) = HALF 
            FX_BAR(I) = ONE - FX(I) 
            FX_E(I) = DX(I+1)/(DX(I+1)+DX(I)) 
            FX_E_BAR(I) = ONE - FX_E(I) 
         END DO 
      ENDIF 

! Look at 2 through JMAX1 V-momentum cells
      IF (DO_J) THEN 
         DO J = JMIN1, JMAX1 
            DY_N = HALF*(DY(J+1)+DY(J)) 
            ODY_N(J) = ONE/DY_N 
            FY_N(J) = DY(J+1)/(DY(J+1)+DY(J)) 
            FY_N_BAR(J) = ONE - FY_N(J) 
         END DO 
      ENDIF 

! Look at 2 through KMAX1 W-momentum cells
      IF (DO_K) THEN 
         DO K = KMIN1, KMAX1 
            DZ_T = HALF*(DZ(K+1)+DZ(K)) 
            ODZ_T(K) = ONE/DZ_T 
            FZ_T(K) = DZ(K+1)/(DZ(K+1)+DZ(K)) 
            FZ_T_BAR(K) = ONE - FZ_T(K) 
         END DO 
      ENDIF 

! Look at last U-, V-, and W-momentum cells
      DX_E = DX(IMAX2) 
      DY_N = DY(JMAX2) 
      DZ_T = DZ(KMAX2) 
      ODX_E(IMAX2) = ONE/DX_E 
      ODY_N(JMAX2) = ONE/DY_N 
      ODZ_T(KMAX2) = ONE/DZ_T 
      FX(IMAX2) = HALF 
      FX_BAR(IMAX2) = HALF 
      FX_E(IMAX2) = HALF 
      FX_E_BAR(IMAX2) = HALF 
      
      FY_N(JMAX2) = HALF 
      FY_N_BAR(JMAX2) = HALF      

      FZ_T(KMAX2) = HALF 
      FZ_T_BAR(KMAX2) = HALF 
      FZ_T(KMAX3) = HALF 
      FZ_T_BAR(KMAX3) = HALF 

      ODX_E(IMAX3) = ONE/DX_E
      ODY_N(JMAX3) = ONE/DY_N
      ODZ_T(KMAX3) = ONE/DZ_T
      FX(IMAX3) = HALF
      FX_BAR(IMAX3) = HALF
      FX_E(IMAX3) = HALF
      FX_E_BAR(IMAX3) = HALF

      FY_N(JMAX3) = HALF
      FY_N_BAR(JMAX3) = HALF

      
      IF(CYCLIC_X) THEN
         FX_E(1) = FX_E(IMAX1)
         FX_E_BAR(1) = FX_E_BAR(IMAX1) 
          IF(NODESI.NE.1) THEN
            FX_E(IMIN3) = FX_E(IMAX1-1) 
            FX_E_BAR(IMIN3) = FX_E_BAR(IMAX1-1) 
         ENDIF
      ENDIF 
      IF (CYCLIC_Y) THEN
         FY_N(1) = FY_N(JMAX1) 
         FY_N_BAR(1) = FY_N_BAR(JMAX1) 
         IF(NODESJ.NE.1) THEN
            FY_N(JMIN3) = FY_N(JMAX1-1) 
            FY_N_BAR(JMIN3) = FY_N_BAR(JMAX1-1) 
         ENDIF
      ENDIF 
      IF (CYCLIC_Z) THEN 
         FZ_T(1) = FZ_T(KMAX1) 
         FZ_T_BAR(1) = FZ_T_BAR(KMAX1) 
         IF(NODESK.NE.1) THEN
            FZ_T(KMIN3) = FZ_T(KMAX1-1) 
            FZ_T_BAR(KMIN3) = FZ_T_BAR(KMAX1-1) 
         ENDIF
      ENDIF 

      CALL FINL_ERR_MSG

      RETURN  
      END SUBROUTINE SET_GEOMETRY
