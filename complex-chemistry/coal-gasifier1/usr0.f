!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!  Author: S. Venkatesan, M. Syamlal                  Date: 21-JUN-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE USR0
!
      Use param
      Use param1
!
!  Use files defining common blocks here
!
      Use physprop
      Use funits
      Use usr
      
      IMPLICIT NONE
      
      INCLUDE 'usrnlst.inc'
!
!  Function subroutines
!
      LOGICAL COMPARE
!
!  Define local variables here
!
!
!     Type of Coal
      CHARACTER*20     SCOALNAM(5)
!
!                      Temperatures, Pressures, Proximate Analysis,
!                      Rate constants, activation energies,
!                      reaction rates
      DOUBLE PRECISION SAK2(5), SAE2(5), SAK5(5), SAE5(5), SAKM(5),&
                       SAEM(5), SAKD(5), SAED(5), SAKC(5), SAEC(5),&
                       SWG3(5), SUM
!
!  KINETIC CONSTANTS FOR VARIOUS COALS -- WEN ET AL. (1982)
!  These constants are acccessed through parameter COAL:
!    COAL     Coal type
!      1      Pittsburgh No.8 (Bituminous)
!      2      Arkwright Pittsburgh (Bituminous)
!      3      Illinois No. 6 (Bituminous)
!      4      Rosebud (Subbituminous)
!      5      North Dakota Lignite
!
!                   1       2      3       4        5
!
      DATA SAK2 /   930.,   600.,  2250.,    70.,    1.1/
      DATA SAE2 / 45000., 45000., 42000., 30000., 30000./
      DATA SAK5 /   930.,   600.,  2250.,    70.,    1.1/
      DATA SAE5 / 45000., 45000., 42000., 30000., 30000./
      DATA SAKM /  1.1E5,  1.1E5,  1.1E5,  7.5E4,  5.1E4/
      DATA SAEM /  21200., 21200., 21200., 18700., 16200./
      DATA SAKD /  1.1E5,  1.1E5,  1.1E5,  7.5E4,  5.1E4/
      DATA SAED / 21200., 21200., 21200., 18700., 16200./
      DATA SAKC /  2.5E7,  2.5E7,  2.5E7,  9.0E7, 2.1E8 /
      DATA SAEC / 29000., 29000., 29000., 27750., 26500./
      DATA SWG3 / 0.0068, 0.0068, 0.0155,  0.014,    0.5/
      DATA SCOALNAM/'Pittsburgh No. 8',&
                    'Arkwright',&
                    'Illinois No. 6',&
                    'Rosebud',&
                    'N. Dakota Lignite'/
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!devol
      Allocate( VMSTAR(DIMENSION_3) )
      Allocate( N_sh (DIMENSION_3, DIMENSION_M) )
!
! CHECK COAL type, proximate, and ultimate analyses
!
      IF (USE_DEVOL) THEN
!
!      the next statement just reminds you to use the correct version
!          of RRATES rather than the original DUMMY version if DEVOL
!          is set to .TRUE.
!
       CALL ERROR_ROUTINE ('USR0',&
         'Be sure to use the correct ver. of RRATES - not DUMMY ver',&
          0, 2)
!
        IF(COAL .EQ. UNDEFINED_I .OR. COAL .LE. 0&
           .OR. COAL .GT. 4) THEN
         CALL ERROR_ROUTINE ('USR0',&
           'Coal Type not specified or illegal', 1, 1)
        ELSE IF&
        (PAFC .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'PAFC not specified', 1, 1)
        ELSE IF&
        (PAVM .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'PAVM not specified', 1, 1)
        ELSE IF&
        (PAM .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'PAM not specified', 1, 1)
        ELSE IF&
        (UAC .EQ. UNDEFINED )THEN
          CALL ERROR_ROUTINE ('USR0', 'UAC not specified', 1, 1)
        ELSE IF&
        (UAH .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'UAH not specified ', 1, 1)
        ELSE IF&
        (UAO .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'UAO not specified ', 1, 1)
        ELSE IF&
        (UAN .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'UAN not specified', 1, 1)
        ELSE IF&
        (UAS .EQ. UNDEFINED ) THEN
          CALL ERROR_ROUTINE ('USR0', 'UAS not specified', 1, 1)
        END IF
        IF(PAA .NE. UNDEFINED)THEN
          SUM = PAFC + PAVM + PAM + PAA
          IF( .NOT.COMPARE(ONE,SUM) )THEN
            WRITE(UNIT_LOG,'(A,F10.5/A)')&
              ' *** PAFC + PAVM + PAM + PAA = ',SUM,&
              '       It should be equal to 1.0'
            CALL MFIX_EXIT
          ENDIF
        ELSE
          PAA = 1.0 - (PAFC + PAVM + PAM)
        ENDIF
        SUM = UAC + UAH + UAO + UAN + UAS + PAM + PAA
        IF(.NOT.COMPARE(ONE,SUM))THEN
          WRITE(UNIT_LOG,'(A,F10.5/A)')&
            ' *** UAC + UAH + UAO + UAN + UAS + PAM + PAA = ',SUM,&
            '       It should be equal to 1.0'
          CALL MFIX_EXIT
        ENDIF
      END IF
!
!  Set constants needed for MGAS chemistry
!
!
      COALNAM = SCOALNAM(COAL)
      AK2 = SAK2(COAL)
      AE2 = SAE2(COAL)
      AK5 = SAK5(COAL)
      AE5 = SAE5(COAL)
      AKM = SAKM(COAL)
      AEM = SAEM(COAL)
      AKD = SAKD(COAL)
      AED = SAED(COAL)
      AKC = SAKC(COAL)
      AEC = SAEC(COAL)
      WG3 = SWG3(COAL)
!
!  DENSITY OF DRY, ASH-FREE, COAL
!
!      DAFC = RO_s(1) * (PAFC + PAVM)
      DAFC = (PAFC + PAVM)            !for constant density
!
!  Void fraction of ash layer
      EP_A = 0.25 + 0.75 * ( 1.0 - PAA )
      f_EP_A = EP_A**2.5
!
!  Set constants for devolatilization calculations
!
      CALL SET_DEV_CONS
!
      RETURN
      END SUBROUTINE USR0
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: set_dev_cons.for                                       C
!  Purpose: set devolatilization constants                             C
!                                                                      C
!  Author: S. Venkatesan                              Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Implement modified MGAS chemistry                          C
!  Author: M. Syamlal                                 Date: 22-JUN-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: MGAS code, Wen, et al. (1982)       C
!                                                                      C
!  Variables referenced: COAL, IJK, T_g, T_s1, UAC, UAH, UAO, UAS,     C
!                        UAN, PAFC, PAVM, PAM, TOLFC, RO_s
!
!  Variables modified: FVC, FVH, FVO, FVN, FVS, ALPHAD, ALPHAC, BETAD, C
!                      BETAC, HHVC, HHVT, HEATD, HEATC, DAFC, AK2,     C
!                      AE2, AK5, AE5, AKM, AEM, AKD, AED, AKC, AEC,    C
!                      WG3, FTC, FTH, FTO, FTN, FTS, DOCO, DOCO2,      C
!                      DOH2O, DHH2, DHCH4, DHC2H6, DHC2H4, DHC3H8,     C
!                      DHC6H6, COCO, COCO2, COH2O, CHH2, CHCH4,        C
!                      CHC2H6, CHC2H4, CHC3H8, CHC6H6                  C
!                                                                      C
!  Local variables: PAFCN, SAK2, SAE2, SAK5, SAE5, SAKM, SAEM, SAKD,   C
!                   SAED, SAKC, SAEC, SWG3                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_DEV_CONS
      
      Use param
      Use param1
      Use fldvar
      Use geometry
      Use run
      Use indices
      Use physprop
      Use constant
      Use funits
      Use usr
      Use compar
      
      IMPLICIT NONE
      
      INCLUDE 'usrnlst.inc'
!
!  FUnction subroutines
      LOGICAL COMPARE
!               needed for function.inc to work
      INTEGER  IJK
!
!
!  Local variables
!
      DOUBLE PRECISION PAFCN, H1, H2, H3, H4, H5, H6, H7, H8,&
                       YYY, SUM
!
      INCLUDE 'function.inc'
      
    

!
      IF(PAVM .EQ. ZERO) RETURN
      CALL START_LOG
!
!  Error check the data
!
!
      SUM = FTC + FTH + FTO + FTN + FTS
      IF(.NOT.COMPARE(ONE,SUM))THEN
        WRITE(UNIT_LOG,'(A,F10.5/A)')&
          ' *** FTC + FTH + FTO + FTN + FTS = ',SUM,&
          '       It should be equal to 1.0'
        CALL MFIX_EXIT
      ENDIF
!
      SUM = DOCO + DOCO2 + DOH2O
      IF(.NOT.COMPARE(ONE,SUM))THEN
        WRITE(UNIT_LOG,'(A,F10.5/A)')&
          ' *** DOCO + DOCO2 + DOH2O = ',SUM,&
          '       It should be equal to 1.0'
        CALL MFIX_EXIT
      ENDIF
!
      SUM = DHH2 + DHCH4 + DHC2H6 + DHC2H4 + DHC3H8 + DHC6H6
      IF(.NOT.COMPARE(ONE,SUM))THEN
        WRITE(UNIT_LOG,'(A,F10.5/A)')&
          ' *** DHH2 + DHCH4 + DHC2H6 + DHC2H4 + DHC3H8 + DHC6H6 = ',&
          SUM,'       It should be equal to 1.0'
        CALL MFIX_EXIT
      ENDIF
      
!
      SUM = COCO + COCO2 + COH2O
      IF(.NOT.COMPARE(ONE,SUM))THEN
        WRITE(UNIT_LOG,'(A,F10.5/A)')&
          ' *** COCO + COCO2 + COH2O = ',SUM,&
          '       It should be equal to 1.0'
        CALL MFIX_EXIT
      ENDIF
!
      SUM = CHH2 + CHCH4 + CHC2H6 + CHC2H4 + CHC3H8 + CHC6H6
      IF(.NOT.COMPARE(ONE,SUM))THEN
        WRITE(UNIT_LOG,'(A,F10.5/A)')&
          ' *** CHH2 + CHCH4 + CHC2H6 + CHC2H4 + CHC3H8 + CHC6H6 = ',&
          SUM,'       It should be equal to 1.0'
        CALL MFIX_EXIT
      ENDIF
!
!  Set constants for tar combustion rxn
!    Tar + f3_1 O2 --> f3_3 CO2 + f3_6 H2O    (HEATF3)
!
      F3_1 = MW_g(8) * ( FTC/12. + FTH/4. - FTO/32.)
      F3_3 = MW_g(8) * FTC/12.
      F3_6 = MW_g(8) * FTH/2.
      HEATF3 = MW_g(8) * ( (FTC/12.) * (-94052.) + (FTH/2.) * (-57798.))
!
!    Determine the composition of VM
!
        PAFCN = UAC + UAH + UAO + UAS + UAN - PAVM
        IF(.NOT.COMPARE(PAFCN, PAFC) )THEN
          IF(PAFCN .GT. PAFC )THEN
            YYY = (PAFCN - PAFC)*100/PAFC
            WRITE(UNIT_LOG,*)&
              ' *** Fixed Carbon (FC) increased by ',YYY,' %'
          ELSE
            YYY = (PAFC - PAFCN)*100/PAFC
            WRITE(UNIT_LOG,*)&
              ' *** Fixed Carbon (FC) decreased by ',YYY,' %'
          ENDIF
          PAFC = PAFCN
        ENDIF
        FVC = (UAC - PAFC)/PAVM
        FVH = UAH/PAVM
        FVO = UAO/PAVM
        FVN = UAN/PAVM
        FVS = UAS/PAVM
        IF(FVN .NE. ZERO .OR. FVS .NE. ZERO) THEN
          WRITE(UNIT_LOG,1000)
          CALL MFIX_EXIT
        ENDIF
!
!      Tar fraction in the devolatilization reaction
!
        ALPHAD = ((42*DHCH4+168*DHC6H6+63*DHC3H8+56*DHC2H6+84*DHC2H4)*FVS&
         +((84*DHCH4+336*DHC6H6+126*DHC3H8+112*DHC2H6+168*DHC2H4)*DOH2O&
         -84*DOCO2-168*DOCO)*FVO+(144*DHCH4+576*DHC6H6+216*DHC3H8+192*DHC2H6&
         +288*DHC2H4)*FVN+(-672*DHCH4-2688*DHC6H6-1008*DHC3H8-896*DHC2H6&
         -1344*DHC2H4)*FVH+224*FVC)/((42*DHCH4+168*DHC6H6+63*DHC3H8+&
         56*DHC2H6+84*DHC2H4)*FTS+((84*DHCH4+336*DHC6H6+126*DHC3H8+112*DHC2H6&
         +168*DHC2H4)*DOH2O-84*DOCO2-168*DOCO)*FTO+(144*DHCH4+576*DHC6H6&
         +216*DHC3H8+192*DHC2H6+288*DHC2H4)*FTN+(-672*DHCH4-2688*DHC6H6&
         -1008*DHC3H8-896*DHC2H6-1344*DHC2H4)*FTH+224*FTC)
        IF(ALPHAD .GT. ONE .OR. ALPHAD .LT. ZERO)THEN
          WRITE(UNIT_LOG,*)&
            ' *** Error: Tar fraction (ALPHAD) = ', ALPHAD
          CALL ERROR_ROUTINE('SET_DEV_CONS',&
            'Invalid Tar Frac.  (ALPHAD)', 1, 2)
        ENDIF
!
!     Char fraction in the cracking reaction
!
        ALPHAC =  FTC - (FTH - FTS * 2/32 - FTN * 3/14 - FTO * COH2O *&
                2/16) *(CHCH4 * 12/4 + CHC2H6 * 24/6 + CHC2H4 * 24/4 +&
                CHC3H8 * 36/8 + CHC6H6 * 72/6) - FTO * ( COCO * 12/16 +&
                COCO2 * 12/32)
        IF(ALPHAC .GT. ONE .OR. ALPHAC .LT. ZERO)THEN
          WRITE(UNIT_LOG,*)&
            ' *** Error: Char fraction (ALPHAC) = ', ALPHAC
          CALL ERROR_ROUTINE('SET_DEV_CONS',&
            'Invalid Char Frac.  (ALPHAC)', 1, 2)
        ENDIF
!
!   Devolatilization reaction
!
!      H consumed in H2S formation,
!
        H1 = (FVS - ALPHAD*FTS) * 2/32
        IF(H1 .LT. ZERO)THEN
          WRITE(UNIT_LOG,*)&
            ' *** Not enough H for H2S during devolatilization'
          CALL ERROR_ROUTINE('SET_DEV_CONS', 'Invalid H1', 1, 2)
        ENDIF
!
!      H consumed in NH3 formation,
!
        H2 = (FVN - ALPHAD * FTN) * 3/14
        IF(H2 .LT. ZERO) THEN
          WRITE(UNIT_LOG,*)&
            ' *** Not enough H for NH3 during devolatilization'
          CALL ERROR_ROUTINE('SET_DEV_CONS', 'Invalid H2', 1, 2)
        ENDIF
!
!      H consumed in H2O formation,
!
        H3 = (FVO - ALPHAD * FTO) * DOH2O * 2/16
        IF(H3 .LT. ZERO) THEN
          WRITE(UNIT_LOG,*)&
            ' *** Not enough H for H2O during devolatilization'
          CALL ERROR_ROUTINE('SET_DEV_CONS', 'Invalid H3', 1, 2)
        ENDIF
!
!      H remaining,
!
        H4 = FVH - ALPHAD * FTH - H1 - H2 - H3
        IF(H4 .LT. ZERO) THEN
          WRITE(UNIT_LOG,*)&
            ' *** Not enough H for H2S, NH3, and H2O formation',&
            'during devolatilization'
          CALL ERROR_ROUTINE('SET_DEV_CONS', 'Invalid H4', 1, 2)
        ENDIF
!
!   Cracking reaction
!
!      H consumed in H2S formation,
!
        H5 = FTS * 2/32
!
!      H consumed in NH3 formation,
!
        H6 = FTN * 3/14
!
!      H consumed in H2O formation,
!
        H7 = FTO * COH2O * 2/16
!
!      H remaining,
!
        H8 = FTH - H5 - H6 - H7
        IF(H8 .LT. ZERO) THEN
          WRITE(UNIT_LOG,*)&
            ' *** Not enough H for H2S, NH3, and H2O formation',&
            'during tar cracking'
          CALL ERROR_ROUTINE('SET_DEV_CONS', 'Invalid H8', 1, 2)
        ENDIF
!
!     Coefficients for devolatilization reaction
!      there is no BETAD(1) since that would correspond to O2
!
!      CO
        BETAD(2) = (FVO - ALPHAD * FTO) * DOCO * 28/16
!      CO2
        BETAD(3) = (FVO - ALPHAD * FTO) * DOCO2 * 44/32
!      CH4
        BETAD(4) = H4 * DHCH4 * 16/4
!      C2H4
        BETAD(11) = H4 * DHC2H4 * 28/4
!      C2H6
        BETAD(12) = H4 * DHC2H6 * 30/6
!      C3H8
        BETAD(13) = H4 * DHC3H8 * 44/8
!      C6H6
        BETAD(14) = H4 * DHC6H6 * 78/6
!      H2
        BETAD(5) = H4 * DHH2
!      H2O
        BETAD(6) = H3 * 18/2
!      H2S
        BETAD(7) = H1 * 34/2
!      NH3
        BETAD(9) = H2 * 17/3
!
!  Coefficients for cracking reaction
!      again, there would be no BETAC(1) since that would correspond
!      to O2
!
!      CO
        BETAC(2) = FTO * COCO * 28/16
!      CO2
        BETAC(3) = FTO * COCO2 * 44/32
!      CH4
        BETAC(4) = H8 * CHCH4 * 16/4
!      C2H4
        BETAC(11) = H8 * CHC2H4 * 28/4
!      C2H6
        BETAC(12) = H8 * CHC2H6 * 30/6
!      C3H8
        BETAC(13) = H8 * CHC3H8 * 44/8
!      C6H6
        BETAC(14) = H8 * CHC6H6 * 78/6
!      H2
        BETAC(5) = H8 * CHH2
!      H2O
        BETAC(6) = H7 * 18/2
!      H2S
        BETAC(7) = H5 * 34/2
!      NH3
        BETAC(9) = H6 * 17/3
!
!  CALCULATE THE HEATING VALUES OF COAL & TAR, IF NOT GIVEN, FROM
!     THE DULONG FORMULA
!        IF(HHVC .EQ. UNDEFINED)THEN
!          HHVC = 8080. * UAC + 34444.4 * (UAH - UAO/8.) +
!     &           2277.8 * UAS
!        ENDIF
!        IF(HHVT .EQ. UNDEFINED)THEN
!          HHVT = 8080. * FTC + 34444.4 * (FTH - FTO/8.) + 2277.8 * FTS
!        ENDIF
         HHVC = ZERO
         HHVT = ZERO
!
!  HEAT OF DEVOLATILIZATION REACTION (cal/g-VM)
!
        IF( HHVC .GT. ZERO)THEN
          HEATD = (-HHVC -&
                  PAVM * ( (-HHVT) * ALPHAD +&
                              ( -2415.6) * BETAD(2) +&
                              (-13300.0) * BETAD(4) +&
                              (-12043.9) * BETAD(11) +&
                              (-12427.3) * BETAD(12) +&
                              (-12059.1) * BETAD(13) +&
                              (-10012.6) * BETAD(14) +&
                              (-34158.5) * BETAD(5) +&
                              (  -584.4) * BETAD(6) +&
                              ( -3956.5) * BETAD(7) +&
                              ( -5394.2) * BETAD(9)   ) -&
                  PAFC * (-7837.7)                    )/PAVM
!
!  HEAT OF CRACKING REACTION (cal/g-Tar)
!
          HEATC = -HHVT - (&
                  (-7837.7) * ALPHAC +&
                  ( -2415.6) * BETAC(2) +&
                  (-13300.0) * BETAC(4) +&
                  (-12043.9) * BETAC(11) +&
                  (-12427.3) * BETAC(12) +&
                  (-12059.1) * BETAC(13) +&
                  (-10012.6) * BETAC(14) +&
                  (-34158.5) * BETAC(5) +&
                  (  -584.4) * BETAC(6) +&
                  ( -3956.5) * BETAC(7) +&
                  ( -5394.2) * BETAC(9)   )
        ELSE
          HEATD = 0.0
          HEATC = 0.0
        ENDIF
!
      CALL END_LOG
      RETURN
1000  FORMAT(/1X,70('*')//' From: SET_DEV_CONS',&
      /' Message: This version of Chemistry cannot account for '&
      /' sulfur or nitrogen in coal'&
       ,/1X, 70('*')/)
      END SUBROUTINE SET_DEV_CONS
