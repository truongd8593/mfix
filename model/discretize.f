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
      DOUBLE PRECISION, PARAMETER :: TWO = 2.D0 
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
         SMART = HALF*MAX(ZERO,MIN(4.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0)) 
      ELSE IF (PHI_C == ONE) THEN 
         SMART = ONE 
      ELSE                                       !first order upwinding 
         SMART = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION SMART 

      DOUBLE PRECISION FUNCTION Chi_SMART (PHI_C, Chi) 
!	calculate DWF from Chi_SMART scheme
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C, Chi 

      if(PHI_C < ONE)then
        Chi_SMART = Chi * (3.D0/8.D0 - PHI_C/4.d0) / (ONE-PHI_C)
      elseif(Chi == zero)then ! insures that all species equations uses the same chi
        Chi_SMART = zero
      else
      Chi_SMART = ONE
      endif

      RETURN  
      END FUNCTION Chi_SMART 
      
      DOUBLE PRECISION FUNCTION Chi4SMART (PHI_C, PHIU, PHIC, PHID) 
!	calculate CHI for SMART scheme
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C , PHIU, PHIC, PHID

      IF(PHI_C > ZERO .AND. PHI_C <= (1.D0/6.D0))THEN
        Chi4SMART = 16.D0 * PHI_C/(3.D0 - 2.D0 * PHI_C)
      ELSEIF(PHI_C > (1.D0/6.D0) .AND. PHI_C <= (5.D0/6.D0))THEN
        Chi4SMART = ONE
      ELSEIF(PHI_C > (5.D0/6.D0) .AND. PHI_C <= ONE)THEN
        Chi4SMART = 8.D0*(ONE-PHI_C)/(3.D0 - 2.D0 * PHI_C)
      ELSEIF( (PHIU < 1d-15) .AND. (PHIC < 1d-15) .AND. (PHID < 1d-15) )THEN 
      ! if a species Xg is less than machine precision, do not use its chi.
      ! this will avoid calculating chi = 0 for a vanishing species.
        Chi4SMART = LARGE_NUMBER
      ELSE
        Chi4SMART = ZERO
      ENDIF

      RETURN  
      END FUNCTION Chi4SMART 
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
      DOUBLE PRECISION, PARAMETER :: FIVEOSIX = 5.D0/6.D0 
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
      ELSE IF (PHI_C > 3.D0/(8.D0*OCF - 6.D0)) THEN 
         TH = PHI_C/(ONE - PHI_C) 
         ULTRA_QUICK = HALF - 0.125D0*(ONE - TH) 
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
         FCF = -(ONE - CF*CF)/3.D0 
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
         MUSCL = HALF*MAX(ZERO,MIN(2.0D0*TH,(ONE + TH)/2.0D0,2.0D0)) 
      ELSE IF (PHI_C == ONE) THEN 
         MUSCL = ONE 
      ELSE                                       !first order upwinding 
         MUSCL = ZERO 
      ENDIF 
!
      RETURN  
      END FUNCTION MUSCL 

      DOUBLE PRECISION FUNCTION Chi_MUSCL (PHI_C, Chi) 
!	calculate DWF from Chi_MUSCL scheme
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C, Chi 

      if(PHI_C < ONE)then
        Chi_MUSCL = Chi /(4.D0 * (ONE - PHI_C))
      elseif(Chi == zero)then ! insures that all species equations uses the same chi
        Chi_MUSCL = zero
      else
        Chi_MUSCL = ONE
      endif

      RETURN  
      END FUNCTION Chi_MUSCL 
      
      DOUBLE PRECISION FUNCTION Chi4MUSCL (PHI_C, PHIU, PHIC, PHID) 
!	calculate CHI for MUSCL scheme
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_C, PHIU, PHIC, PHID

      IF(PHI_C > ZERO .AND. PHI_C <= 0.25d0)THEN
        Chi4MUSCL = 4.D0 * PHI_C
      ELSEIF(PHI_C > 0.25d0 .AND. PHI_C <= 0.75d0)THEN
        Chi4MUSCL = ONE
      ELSEIF(PHI_C > 0.75d0 .AND. PHI_C <= ONE)THEN
        Chi4MUSCL = 4.D0*(ONE - PHI_C)
      ELSEIF( (PHIU < 1d-15) .AND. (PHIC < 1d-15) .AND. (PHID < 1d-15) )THEN 
      ! if a species Xg is less than machine precision, do not use its chi.
      ! this will avoid calculating chi = 0 for a vanishing species.
        Chi4MUSCL = LARGE_NUMBER
      ELSE
        Chi4MUSCL = ZERO
      ENDIF

      RETURN  
      END FUNCTION Chi4MUSCL 
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




! Central scheme for MMS cases.
      DOUBLE PRECISION FUNCTION CENTRAL_SCHEME (PHI_C)
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
!
      CENTRAL_SCHEME = HALF

      RETURN  
      END FUNCTION CENTRAL_SCHEME

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
         UMIST = HALF*MAX(ZERO,MIN(2.0D0*TH,0.75D0 + 0.25D0*TH,2.0D0)) 
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
!
!
      DOUBLE PRECISION FUNCTION FPFOI_OF (PHI_D, PHI_C, &
                       PHI_U, PHI_UU) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_D, PHI_C, PHI_U, PHI_UU
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: UNIV_LIMITER_OF
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
      FPFOI_OF = PHI_C + (5.0D0/16.0D0)*(PHI_D - PHI_C) + &
                         (1.0D0/4.0D0)*(PHI_C - PHI_U) - &
                         (1.0D0/16.0D0)*(PHI_U - PHI_UU)
!
!	LIMIT THE HIGH ORDER INTERPOLATION
!
      FPFOI_OF = UNIV_LIMITER_OF(FPFOI_OF , PHI_D, PHI_C, PHI_U)
      RETURN  
      END FUNCTION FPFOI_OF 
!
!
      DOUBLE PRECISION FUNCTION UNIV_LIMITER_OF (PHI_TEMP, PHI_D, PHI_C, &
                                                   PHI_U)
! 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1
      USE run
!      USE fldvar
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION PHI_TEMP, PHI_D, PHI_C, PHI_U
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
	DOUBLE PRECISION PHI_REF, DEL, CURV
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
      	DEL = PHI_D - PHI_U
      	CURV = PHI_D - 2.0*PHI_C + PHI_U
      IF (ABS (CURV) >= ABS (DEL)) THEN     ! NON-MONOTONIC
      	UNIV_LIMITER_OF = PHI_C
      ELSE
        PHI_REF = PHI_U + (PHI_C-PHI_U)/C_FAC
       IF (DEL > ZERO) THEN
         IF (PHI_TEMP < PHI_C)UNIV_LIMITER_OF = PHI_C
         IF (PHI_TEMP > MIN(PHI_REF,PHI_D)) &
                      UNIV_LIMITER_OF = MIN(PHI_REF,PHI_D)
         IF (PHI_TEMP >= PHI_C .AND. PHI_TEMP <= MIN(PHI_REF,PHI_D)) &
                      UNIV_LIMITER_OF = PHI_TEMP
       ELSE
         IF (PHI_TEMP > PHI_C)UNIV_LIMITER_OF = PHI_C
         IF (PHI_TEMP < MAX(PHI_REF,PHI_D)) &
                      UNIV_LIMITER_OF = MAX(PHI_REF,PHI_D)
         IF (PHI_TEMP >= MAX(PHI_REF,PHI_D) .AND. PHI_TEMP <= PHI_C) &
                      UNIV_LIMITER_OF = PHI_TEMP
       ENDIF
      ENDIF
      RETURN  
      END FUNCTION UNIV_LIMITER_OF 
