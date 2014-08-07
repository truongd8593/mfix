!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!
!  module name: des_write_restart
!  purpose: writing des data for restart
!
!  Author : Pradeep G
!  Purpose : Reads either single restart file or multiple restart files
!            (based on bdist_io) flag
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      subroutine des_write_restart

      use param1
      use compar
      use discretelement
      use run
      use des_bc
      use des_rxns
      use des_thermo
      use mfix_pic, only: MPPIC
      use mfix_pic, only: DES_STAT_WT

      use res_des

      use mpi_utility

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      INTEGER :: LC1, LC2
      INTEGER :: lNEXT_REC
      INTEGER :: lOUT_DIMN

!-----------------------------------------------
      lOUT_DIMN = merge(2,3,NO_K)

      CALL INIT_WRITE_RES_DES(trim(RUN_NAME), lNEXT_REC)

      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, DTSOLID)


      DO LC1 = 1, lOUT_DIMN
         CALL WRITE_RES_DES(lNEXT_REC, DES_POS_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, iGLOBAL_ID)

      DO LC1 = 2, 4
         CALL WRITE_RES_DES(lNEXT_REC, PEA(:,LC1))
      ENDDO

      DO LC1 = 1, lOUT_DIMN
         CALL WRITE_RES_DES(lNEXT_REC, DES_VEL_NEW(LC1,:))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL WRITE_RES_DES(lNEXT_REC, OMEGA_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, DES_RADIUS)
      CALL WRITE_RES_DES(lNEXT_REC, RO_SOL)

      IF(MPPIC) &
         CALL WRITE_RES_DES(lNEXT_REC, DES_STAT_WT)

      IF(ENERGY_EQ) &
         CALL WRITE_RES_DES(lNEXT_REC, DES_T_s_NEW)

      IF(ANY_SPECIES_EQ) THEN
         DO LC1=1, DIMENSION_N_S 
            CALL WRITE_RES_DES(lNEXT_REC, DES_X_s(:,LC1))
         ENDDO
      ENDIF

      CALL WRITE_RES_DES(lNEXT_REC, NEIGHBOURS(:,1))
      CALL WRITE_RES_DES(lNEXT_REC, PN(1,:))

      DO LC1=2, MAXNEIGHBORS
         CALL WRITE_RES_DES(lNEXT_REC, NEIGHBOURS(:,LC1), pLOC2GLB=.TRUE.)
         CALL WRITE_RES_DES(lNEXT_REC, PN(LC1,:), pLOC2GLB=.TRUE.)
         CALL WRITE_RES_DES(lNEXT_REC, PV(LC1,:))

         DO LC2=1, lOUT_DIMN
            CALL WRITE_RES_DES(lNEXT_REC, PFT(:,LC1,LC2))
         ENDDO
      ENDDO

      CALL FINL_WRITE_RES_DES

      RETURN
      END SUBROUTINE DES_WRITE_RESTART
