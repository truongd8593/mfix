!      Program Test; CALL Read_Therm_tester; END Program Test
      
      SUBROUTINE READ_Therm_tester 
      Implicit none
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom, MW
      DOUBLE PRECISION Cp1, Cp2, h1, h2, T, Hf298oR
      DOUBLE PRECISION, EXTERNAL :: calc_CpoR, calc_H0oR
      CHARACTER PATH*132, SPECIES*18
      integer i, funit, IER
      CHARACTER(len=142) FILENAME
      CHARACTER(len=10) :: THERM = 'BURCAT.THR'
      LOGICAL LocalCopy
      
      SPECIES = 'CH4'
      PATH = '.'
      funit = 5
      
      INQUIRE(FILE=TRIM(THERM),EXIST=LocalCopy)
      IF(LocalCopy)Then
        OPEN(UNIT=funit,FILE=TRIM(THERM))
      ELSE
        FILENAME = TRIM(PATH) // '/'  // TRIM(THERM)
        OPEN(UNIT=funit,FILE=TRIM(FILENAME), ERR=500)
      ENDIF
 
 !      Call Read_Therm(PATH, 'N2', Thigh, Tlow, Tcom, Ahigh, Alow, Hf298oR)
      Call Read_Therm(funit, SPECIES, Thigh, Tlow, Tcom, MW, Ahigh, Alow, Hf298oR, IER)
      IF(IER /= 0) GOTO 200
      
      print *, SPECIES
      print *, Thigh, Tlow, Tcom, MW, Hf298oR*1.987207
      
!      print *, Hf298oR
!      T = 300
!      DO i = 1, 12
!        Cp1 = calc_CpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
!        T = T + 100
!        print *, T, Cp1
!      ENDDO
      
      Cp1 = calc_CpoR(8D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
      h1 = calc_H0oR(4D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
      h2 = calc_H0oR(12D2, Thigh, Tlow, Tcom, Ahigh, Alow)*1.987207
      print *, Cp1, h1, h2
      CLOSE(UNIT=funit)
      STOP
200   PRINT *, 'READ_Therm_tester: Species ', TRIM(SPECIES), ' not found in Database!'
      STOP
500   PRINT *, 'READ_Therm_tester: Cannot Open file ', TRIM(THERM), '!'
      PRINT *, 'Check path or copy mfix/model/thermochemical/', TRIM(THERM), ' into run directory'
      STOP       
      END Subroutine READ_Therm_tester

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv  C
!                                                                        C
!     Module name: READ_Therm()                                          C
!     Purpose: Read Thermo coefficients from Burcat and Ruscic (2005)    C
!     Author: M. Syamlal                                 Date: 30-SEP-05 C
!                                                                        C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  C
!     
      SUBROUTINE READ_Therm(funit, Sp, Thigh, Tlow, Tcom, MW, Ahigh, Alow, Hf298oR, IER) 
!     
!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      
      IMPLICIT NONE
      
      CHARACTER(*) SP
      
!     holds one line in the input file
      CHARACTER(len=80) :: LINE_STRING
      
      CHARACTER(len=18) :: SPECIES, ss
      INTEGER i, funit, IER
      DOUBLE PRECISION Ahigh(7), Alow(7), Hf298oR
      DOUBLE PRECISION Thigh, Tlow, Tcom, MW


!     Assoiciate a simple species name with that in BURCAT.THR file
      INTEGER, PARAMETER :: Max = 3
      CHARACTER, DIMENSION(2,Max) :: SPECIES_ALIAS*18
      SPECIES_ALIAS = RESHAPE ( (/ &
!        common name             BURCAT.THR name
!        123456789012345678     123456789012345678
        'O2                ',  'O2 REF ELEMENT    ', &
        'N2                ',  'N2  REF ELEMENT   ', &
        'CH4               ',  'CH4   ANHARMONIC  ' &
       /), (/2,Max/))

      
      IER = 0
      SPECIES = SP
      DO i = 1, MAX
        IF(TRIM(SP) == TRIM(SPECIES_ALIAS(1,I)))THEN
          SPECIES = SPECIES_ALIAS(2,I)
        ENDIF
      ENDDO
       
     
      LINE_STRING = '                '
      DO WHILE(LINE_STRING(1:11) /= 'THERMO DATA')
        READ(UNIT=funit, FMT='(A)',ERR=100,END=100)LINE_STRING
      END DO
      
      ss = '                 '
      call trimTab(SPECIES)
      DO WHILE(TRIM(ss) /= TRIM(SPECIES))
        READ(UNIT=funit, FMT='(A)',ERR=100,END=100)LINE_STRING
        ss = LINE_STRING(1:18)
        call trimTab(ss)
      END DO
!      print *, LINE_STRING
       
      call get_values(LINE_STRING, Tlow, Thigh, MW)
      Tcom = 1000.D0
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Ahigh(1:5)
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Ahigh(6:7), Alow(1:3)
      READ(UNIT=funit, FMT='(5E15.0)',ERR=300,END=300)Alow(4:7), Hf298oR
      
      RETURN
      ! species not found or THERMO DATA not found!
100   IER = 1
      RETURN
            
300   PRINT *, 'READ_Therm: Error reading coefficients for Species ', TRIM(LINE_STRING(1:18))
      STOP       
     
      END SUBROUTINE READ_Therm 
      
      
      Double Precision Function calc_CpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow) 
      !Cp/R
      
      IMPLICIT NONE
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom
      DOUBLE PRECISION T
      
      If(T > Thigh .or. T < Tlow)then
        print *, 'Calc_Cp: Temperature ', T, ' not in the range: ', Tlow, Thigh
        STOP
      elseif (T < Tcom) then
        calc_CpoR = (((Alow(5)*T +Alow(4))*T + Alow(3))*T + Alow(2))*T + Alow(1)
      else
        calc_CpoR = (((Ahigh(5)*T +Ahigh(4))*T + Ahigh(3))*T + Ahigh(2))*T + Ahigh(1)
      endif
      
      RETURN
               
      END Function calc_CpoR 
      
      Double Precision Function calc_ICpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow) 
      ! integral_0_to_T(Cp/R*dT)
      
      IMPLICIT NONE
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom
      DOUBLE PRECISION T, ICpoRT
      
      If(T > Thigh .or. T < Tlow)then
        print *, 'Calc_ICp: Temperature ', T, ' not in the range: ', Tlow, Thigh
        STOP
      elseif (T < Tcom) then
        ICpoRT = (((Alow(5)*T/5.D0 +Alow(4))*T/4.D0 + Alow(3))*T/3.D0 + Alow(2))*T/2.D0 + Alow(1)
      else
        ICpoRT = (((Ahigh(5)*T/5.D0 +Ahigh(4))*T/4.D0 + Ahigh(3))*T/3.D0 + Ahigh(2))*T/2.D0 + Ahigh(1)
      endif
      
      calc_ICpoR = T * ICpoRT
      
      RETURN
      END Function calc_ICpoR 
      
      
      
      Double Precision Function calc_H0oR(T, Thigh, Tlow, Tcom, Ahigh, Alow) 
      ! H0/R
      
      IMPLICIT NONE
      DOUBLE PRECISION Ahigh(7), Alow(7)
      DOUBLE PRECISION Thigh, Tlow, Tcom
      DOUBLE PRECISION T, H0oT, ICp
      DOUBLE PRECISION calc_ICpoR
      
      ICp = calc_ICpoR(T, Thigh, Tlow, Tcom, Ahigh, Alow)
      If (T < Tcom) then
        calc_H0oR = ICp + Alow(6)
      else
        calc_H0oR = ICp + Ahigh(6)
      endif
      
      RETURN
      END Function calc_H0oR 
     
      subroutine replaceTab(C)
      IMPLICIT NONE
      CHARACTER(*) C
      INTEGER I
      DO I = 1, len(C)
        IF(C(I:I) == '	')C(I:I)=' '
      ENDDO
      RETURN
      END subroutine replaceTab
      
      subroutine trimTab(C)
      IMPLICIT NONE
      CHARACTER(*) C
      INTEGER I
      LOGICAL tabFound
      tabFound = .FALSE.
      DO I = 1, len(C)
        IF(C(I:I) == '	')tabFound = .TRUE.
        if(tabFound) C(I:I)=' '
      ENDDO
      RETURN
      END subroutine trimTab
      
      subroutine stripTab(C)
      IMPLICIT NONE
      CHARACTER(*) C
      CHARACTER C1*20
      INTEGER I, I1
      I1 = 1
      DO I = 1, len(C)
        IF(C(I:I) /= '	')C1(I1:I1)=C(I:I)
        I1 = I1+1
      ENDDO
      DO I = 1, len(C)
        IF(I < I1) THEN
          C(I:I)=C1(I:I)
        ELSE
          C(I:I)=' '
        ENDIF
      ENDDO
      RETURN
      END subroutine stripTab
       
