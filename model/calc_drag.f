!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: CALC_DRAG                                               C
!  Purpose: Calculate the gas solids and solids-solids drag terms if   C
!           directed to do so by the corresponding flags               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Modifications for variable grid size capability,           C
!           logic for volume-weighted averaging                        C
!  Author: W. Rogers                                  Date: 20-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: MFIX 2.0 mods                                              C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To call solids drag Drag_SS only when not using DES        C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number 4                                  Date: 2-July-07  C
!  Author: Rahul Garg                                                  C
!  Purpose: If DES_INTERP_ON is true, then additional call is made     C
!           to des/drag_fgs.f in order to obtain cell corner drag      C
!           related AM and BM coeff's                                  C
!                                                                      C
!  Revision Number: 5                                                  C
!  Purpose: To incorporate calls and checks for QMOM for               C
!           solution of the particle phase kinetic equation            C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
     
      SUBROUTINE CALC_DRAG(DRAGD, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE drag
      USE compar
      USE discretelement
      USE qmom_kinetic_equation      
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Error index
      INTEGER :: IER
! Flag for exchange functions
      LOGICAL :: DRAGD(0:DIMENSION_M, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Solids phase indices
      INTEGER :: M, L, DM
! Flag only used when the hybrid model is invoked and notifies the
! routine that the solid phase index M refers to the indice of a
! discrete 'phase' not a continuous phase so that the appropriate
! variables are referenced.  This is currently required until moves 
! are made to fully separate use of F_GS from DEM (i.e. strictly 
! assign F_GS for use with gas-continuum solids drag and assign a 
! new variable, like F_GDS, when refering to gas-discrete solids drag)
      LOGICAL :: DISCRETE_FLAG 
!-----------------------------------------------
    
! Alberto Passalacqua:  QMOMB       
      IF (QMOMK) RETURN


! calculate drag between continuum phases (gas-solid & solids-solids)      
      IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
         DISCRETE_FLAG = .FALSE.   ! only matters if des_continuum_hybrid   
         DO M = 1, SMAX 
            IF (DRAGD(0,M) .AND. RO_G0/=ZERO) THEN 
               CALL DRAG_GS (M, DISCRETE_FLAG, IER)
            ENDIF
            IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN
               DO L = 1, SMAX
                  IF (DRAGD(L,M)) CALL CALC_IA_NONEP_DRAG_SS (L,M,IER)
               ENDDO
            ELSEIF (TRIM(KT_TYPE) /= 'GHD') THEN  ! do nothing for GHD theory
               DO L = 1, M - 1 
                  IF (DRAGD(L,M)) CALL DRAG_SS (L, M, IER)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

! calculate drag between continuum phases and discrete particles
! (gas-particle & solids-particle)
      IF (DES_CONTINUUM_COUPLED) THEN
         CALL DES_DRAG_GS

! calculate drag between continuum solids and discrete solids
! for now non-interpolated solid-solid version is only option
         IF(DES_CONTINUUM_HYBRID) THEN   
            CALL DES_DRAG_SS
         ENDIF
      ENDIF         
 

      RETURN  
      END SUBROUTINE CALC_DRAG 
