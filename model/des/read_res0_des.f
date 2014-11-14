!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_READ_RESTART                                        !
!  Purpose : Reads either single restart file or multiple restart      !
!  fles (based on bdist_io) flag.                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_RES0_DES

      use param1
      use compar
      use discretelement
      use run
      use des_bc
      use des_rxns
      use des_thermo
      use desmpi
      use machine
      use cdist
      use mpi_utility

      use read_res1_des
      use error_manager

      implicit none

      INTEGER :: LC1, LC2
      INTEGER :: lDIMN, lNEXT_REC

      DOUBLE PRECISION :: VERSION

      lDIMN = merge(2,3,NO_K)

      CALL INIT_READ_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL READ_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL READ_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL READ_RES_DES(lNEXT_REC, DTSOLID)

! Position data is read and used to setup pARRAY reads.
      CALL READ_PAR_POS(lNEXT_REC)

      CALL READ_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      DO LC1 = 2, 4
         CALL READ_RES_pARRAY(lNEXT_REC, PEA(:,LC1))
      ENDDO

      DO LC1 = 1, lDIMN
         CALL READ_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(LC1,:))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL READ_RES_pARRAY(lNEXT_REC, OMEGA_NEW(LC1,:))
      ENDDO

      CALL READ_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL READ_RES_pARRAY(lNEXT_REC, RO_SOL)

      IF(MPPIC) CALL READ_RES_pARRAY(lNEXT_REC, DES_STAT_WT)
      IF(ENERGY_EQ) CALL READ_RES_pARRAY(lNEXT_REC, DES_T_s_NEW)

      IF(ANY_SPECIES_EQ) THEN
         DO LC1=1, DIMENSION_N_S
            CALL READ_RES_pARRAY(lNEXT_REC, DES_X_s(:,LC1))
         ENDDO
      ENDIF

! Collision/neighbor data is read and used to setup cARRAY reads.
      CALL READ_PAR_COL(lNEXT_REC)

      CALL READ_RES_cARRAY(lNEXT_REC, PV_COLL(:))
      DO LC1=1, lDIMN
         CALL READ_RES_cARRAY(lNEXT_REC, PFN_COLL(LC1,:))
         CALL READ_RES_cARRAY(lNEXT_REC, PFT_COLL(LC1,:))
      ENDDO

      CALL READ_RES_DES(lNEXT_REC, DEM_BCMI)
      DO LC1=1, DEM_BCMI
         CALL READ_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)

         LC2 = DEM_MI(LC1)%OCCUPANTS

         allocate(DEM_MI(LC1)%W(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         allocate(DEM_MI(LC1)%H(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         allocate(DEM_MI(LC1)%P(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         allocate(DEM_MI(LC1)%Q(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      CALL FINL_READ_RES_DES


      WRITE(ERR_MSG,"('DES restart file read at Time = ',g11.5)") TIME
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN
      END SUBROUTINE READ_RES0_DES
