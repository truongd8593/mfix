!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DISCRETIZE                                             C
!  Purpose: A collection of functions for higher-order discretization  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
!
      DOUBLE PRECISION FUNCTION SUPERBEE (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: TWO = 2. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TH 
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region 
         TH = PHI_C/(ONE - PHI_C) 
         SUPERBEE = HALF*MAX(ZERO,MIN(ONE,TWO*TH),MIN(TWO,TH)) 
      ELSE IF (PHI_C == ONE) THEN 
         SUPERBEE = ONE 
      ELSE                                       !first order upwinding 
         SUPERBEE = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION SUPERBEE 
!
!
      DOUBLE PRECISION FUNCTION SMART (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TH 
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region 
         TH = PHI_C/(ONE - PHI_C) 
         SMART = HALF*MAX(ZERO,MIN(4.*TH,0.75 + 0.25*TH,2.)) 
      ELSE IF (PHI_C == ONE) THEN 
         SMART = ONE 
      ELSE                                       !first order upwinding 
         SMART = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION SMART 
!
!
!
      DOUBLE PRECISION FUNCTION ULTRA_QUICK (PHI_C, CF) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C, CF 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: FIVEOSIX = 5./6. 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TH, OCF 
!-----------------------------------------------
!
      OCF = MAX(ONE,ONE/MAX(1.D-2,CF)) 
      IF (PHI_C > ONE) THEN 
         ULTRA_QUICK = HALF 
      ELSE IF (PHI_C > FIVEOSIX) THEN 
         ULTRA_QUICK = ONE 
      ELSE IF (PHI_C > 3./(8.*OCF - 6.)) THEN 
         TH = PHI_C/(ONE - PHI_C) 
         ULTRA_QUICK = HALF - 0.125*(ONE - TH) 
      ELSE IF (PHI_C > ZERO) THEN 
         TH = PHI_C/(ONE - PHI_C) 
         ULTRA_QUICK = (OCF - ONE)*TH 
      ELSE 
         TH = PHI_C/(ONE - PHI_C) 
         ULTRA_QUICK = HALF*TH 
      ENDIF 
!
      RETURN  
      END FUNCTION ULTRA_QUICK 
!
!
!
      DOUBLE PRECISION FUNCTION QUICKEST (PHI_C, CF, ODXC, ODXUC, ODXCD) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C, CF, ODXC, ODXUC, ODXCD 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: FCF, TH 
!-----------------------------------------------
!
      IF (PHI_C>ZERO .AND. PHI_C<ONE) THEN       !monotonic region 
         FCF = -(ONE - CF*CF)/3. 
         TH = PHI_C/(ONE - PHI_C + SMALL_NUMBER) 
         QUICKEST = HALF*(ONE - CF) + FCF*(ODXC/ODXCD - ODXC*ODXUC*TH/ODXCD**2) 
         IF(PHI_C<CF)QUICKEST=MIN(QUICKEST,(ONE/CF-ONE)*PHI_C/(ONE-PHI_C)) 
         QUICKEST = MAX(ZERO,MIN(ONE,QUICKEST)) 
      ELSE                                       !first order upwinding 
         QUICKEST = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION QUICKEST 
!
!
      DOUBLE PRECISION FUNCTION MUSCL (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TH 
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region 
         TH = PHI_C/(ONE - PHI_C) 
         MUSCL = HALF*MAX(ZERO,MIN(2.*TH,(ONE + TH)/2.,2.)) 
      ELSE IF (PHI_C == ONE) THEN 
         MUSCL = ONE 
      ELSE                                       !first order upwinding 
         MUSCL = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION MUSCL 
!
!
!
      DOUBLE PRECISION FUNCTION VANLEER (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region 
         VANLEER = PHI_C 
      ELSE                                       !first order upwinding 
         VANLEER = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION VANLEER 
!
!
!
      DOUBLE PRECISION FUNCTION MINMOD (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<=ONE) THEN     !monotonic region 
         IF (PHI_C >= HALF) THEN                 !central differencing 
            MINMOD = HALF 
         ELSE                                    !second order upwinding 
            MINMOD = HALF*PHI_C/(ONE - PHI_C) 
         ENDIF 
      ELSE                                       !first order upwinding 
         MINMOD = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION MINMOD 
!
!
      DOUBLE PRECISION FUNCTION UMIST (PHI_C) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: TH 
!-----------------------------------------------
!
      IF (PHI_C>=ZERO .AND. PHI_C<ONE) THEN      !monotonic region 
         TH = PHI_C/(ONE - PHI_C) 
         UMIST = HALF*MAX(ZERO,MIN(2.*TH,0.75 + 0.25*TH,2.)) 
      ELSE IF (PHI_C == ONE) THEN 
         UMIST = ONE 
      ELSE                                       !first order upwinding 
         UMIST = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION UMIST 
!
!----------------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION XSI (V, DWF) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION V, DWF 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      IF (V >= ZERO) THEN 
         XSI = DWF 
      ELSE 
         XSI = ONE - DWF 
      ENDIF 
!
      RETURN  
      END FUNCTION XSI 
!
      DOUBLE PRECISION FUNCTION PHI_C_OF (PHI_U, PHI_C, PHI_D) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_U, PHI_C, PHI_D 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: DEN 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
!
      IF (COMPARE(PHI_D,PHI_U)) THEN 
         PHI_C_OF = ZERO 
      ELSE 
         DEN = PHI_D - PHI_U 
         PHI_C_OF = (PHI_C - PHI_U)/DEN 
      ENDIF 
!
      RETURN  
      END FUNCTION PHI_C_OF 
