!
!VD$G NOVECTOR
!VD$G NOCONCUR
!deck xerabt
      subroutine xerabt(messg, nmessg) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg 
      character messg*(*) 
!-----------------------------------------------
!***begin prologue  xerabt
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  abort program execution and print error message.
!***description
!
!     abstract
!        ***note*** machine dependent routine
!        xerabt aborts the execution of the program.
!        the error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     description of parameters
!        messg and nmessg are as in xerror, except that nmessg may
!        be zero, in which case no message is being supplied.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  1 august 1982
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  (none)
!***end prologue  xerabt
!***first executable statement  xerabt
      call exit (1) 
      end subroutine xerabt 
!deck xerctl
      subroutine xerctl(messg1, nmessg, nerr, level, kontrl) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg, nerr, level, kontrl 
      character messg1*20 
!-----------------------------------------------
!***begin prologue  xerctl
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  allow user control over handling of errors.
!***description
!
!     abstract
!        allows user control over handling of individual errors.
!        just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to xerctl.
!        if the user has provided his own version of xerctl, he
!        can then override the value of kontrol used in processing
!        this message by redefining its value.
!        kontrl may be set to any value from -2 to 2.
!        the meanings for kontrl are the same as in xsetf, except
!        that the value of kontrl changes only for this message.
!        if kontrl is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     description of parameters
!
!      --input--
!        messg1 - the first word (only) of the error message.
!        nmessg - same as in the call to xerror or xerrwv.
!        nerr   - same as in the call to xerror or xerrwv.
!        level  - same as in the call to xerror or xerrwv.
!        kontrl - the current value of the control flag as set
!                 by a call to xsetf.
!
!      --output--
!        kontrl - the new value of kontrl.  if kontrl is not
!                 defined, it will remain at its original value.
!                 this changed value of control affects only
!                 the current occurrence of the current message.
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  (none)
!***end prologue  xerctl
!***first executable statement  xerctl
      return  
      end subroutine xerctl 
!deck xerprt
      subroutine xerprt(messg, nmessg) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg 
      character messg*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(5) :: lun 
      integer :: nunit, lenmes, kunit, iunit, ichar, last 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: i1mach 
!-----------------------------------------------
!***begin prologue  xerprt
!***date written   790801   (yymmdd)
!***revision date  851213   (yymmdd)
!***category no.  r3
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  print error messages.
!***description
!
!     abstract
!        print the hollerith message in messg, of length nmessg,
!        on each file indicated by xgetua.
!     latest revision ---  1 august 1985
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  i1mach,xgetua
!***end prologue  xerprt
!     obtain unit numbers and write line to each unit
!***first executable statement  xerprt
      call xgetua (lun, nunit) 
      lenmes = len(messg) 
      do kunit = 1, nunit 
         iunit = lun(kunit) 
         if (iunit == 0) iunit = i1mach(4) 
         do ichar = 1, lenmes, 72 
            last = min0(ichar + 71,lenmes) 
            write (iunit, '(1x,a)') messg(ichar:last) 
         end do 
      end do 
      return  
      end subroutine xerprt 
!deck xerror
      subroutine xerror(messg, nmessg, nerr, level) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg, nerr, level 
      character messg*(*) 
!-----------------------------------------------
!***begin prologue  xerror
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  process an error (diagnostic) message.
!***description
!
!     abstract
!        xerror processes a diagnostic message, in a manner
!        determined by the value of level and the current value
!        of the library error control flag, kontrl.
!        (see subroutine xsetf for details.)
!
!     description of parameters
!      --input--
!        messg - the hollerith message to be processed, containing
!                no more than 72 characters.
!        nmessg- the actual number of characters in messg.
!        nerr  - the error number associated with this message.
!                nerr must not be zero.
!        level - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (i.e., it is
!                   non-fatal if xsetf has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!
!     examples
!        call xerror('smooth -- num was zero.',23,1,2)
!        call xerror('integ  -- less than full accuracy achieved.',
!    1                43,2,1)
!        call xerror('rooter -- actual zero of f found before interval f
!    1ully collapsed.',65,3,0)
!        call xerror('exp    -- underflows being set to zero.',39,1,-1)
!
!     written by ron jones, with slatec common math library subcommittee
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  xerrwv
!***end prologue  xerror
!***first executable statement  xerror
      call xerrwv (messg, nmessg, nerr, level, 0, 0, 0, 0, 0., 0.) 
      return  
      end subroutine xerror 
!deck xerrwv
      subroutine xerrwv(messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg, nerr, level, ni, i1, i2, nr 
      real r1, r2 
      character messg*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(5) :: lun 
      integer :: lkntrl, maxmes, kdummy, junk, kount, lmessg, lerr, llevel, &
         mkntrl, nunit, isizei, isizef, kunit, iunit, i, ifatal 
      character :: lfirst*20, form*37 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save, i1mach 
!-----------------------------------------------
!***begin prologue  xerrwv
!***date written   800319   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  process an error message allowing 2 integer and 2 real
!            values to be included in the message.
!***description
!
!     abstract
!        xerrwv processes a diagnostic message, in a manner
!        determined by the value of level and the current value
!        of the library error control flag, kontrl.
!        (see subroutine xsetf for details.)
!        in addition, up to two integer values and two real
!        values may be printed along with the message.
!
!     description of parameters
!      --input--
!        messg - the hollerith message to be processed.
!        nmessg- the actual number of characters in messg.
!        nerr  - the error number associated with this message.
!                nerr must not be zero.
!        level - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (i.e., it is
!                   non-fatal if xsetf has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!        ni    - number of integer values to be printed. (0 to 2)
!        i1    - first integer value.
!        i2    - second integer value.
!        nr    - number of real values to be printed. (0 to 2)
!        r1    - first real value.
!        r2    - second real value.
!
!     examples
!        call xerrwv('smooth -- num (=i1) was zero.',29,1,2,
!    1   1,num,0,0,0.,0.)
!        call xerrwv('quadxy -- requested error (r1) less than minimum (
!    1r2).,54,77,1,0,0,0,2,errreq,errmin)
!
!     latest revision ---  1 august 1985
!     written by ron jones, with slatec common math library subcommittee
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  fdump,i1mach,j4save,xerabt,xerctl,xerprt,xersav,
!                    xgetua
!***end prologue  xerrwv
!     get flags
!***first executable statement  xerrwv
      lkntrl = j4save(2,0,.FALSE.) 
      maxmes = j4save(4,0,.FALSE.) 
!     check for valid input
      if (.not.(nmessg>0 .and. nerr/=0 .and. level>=(-1) .and. level<=2)) then 
         if (lkntrl > 0) call xerprt ('fatal error in...', 17) 
         call xerprt ('xerror -- invalid input', 23) 
!        if (lkntrl.gt.0) call fdump
         if (lkntrl > 0) call xerprt ('job abort due to fatal error.', 29) 
         if (lkntrl > 0) call xersav (' ', 0, 0, 0, kdummy) 
         call xerabt ('xerror -- invalid input', 23) 
         return  
      endif 
      junk = j4save(1,nerr,.TRUE.) 
      call xersav (messg, nmessg, nerr, level, kount) 
!     let user override
      lfirst = messg 
      lmessg = nmessg 
      lerr = nerr 
      llevel = level 
      call xerctl (lfirst, lmessg, lerr, llevel, lkntrl) 
!     reset to original values
      lmessg = nmessg 
      lerr = nerr 
      llevel = level 
      lkntrl = max0(-2,min0(2,lkntrl)) 
      mkntrl = iabs(lkntrl) 
!     decide whether to print message
      if (llevel>=2 .or. lkntrl/=0) then 
         if (.not.(llevel==(-1) .and. kount>min0(1,maxmes) .or. llevel==0&
             .and. kount>maxmes .or. llevel==1 .and. kount>maxmes .and. mkntrl&
            ==1 .or. llevel==2 .and. kount>max0(1,maxmes))) then 
            if (lkntrl > 0) then 
               call xerprt (' ', 1) 
!           introduction
               if (llevel == (-1)) call xerprt (&
                  'warning message...this message will only be printed once.', &
                  57) 
               if (llevel == 0) call xerprt ('warning in...', 13) 
               if (llevel == 1) call xerprt ('recoverable error in...', 23) 
               if (llevel == 2) call xerprt ('fatal error in...', 17) 
            endif 
            call xerprt (messg, lmessg) 
            call xgetua (lun, nunit) 
            isizei = log10(float(i1mach(9))) + 1.0 
            isizef = log10(float(i1mach(10))**i1mach(11)) + 1.0 
            do kunit = 1, nunit 
               iunit = lun(kunit) 
               if (iunit == 0) iunit = i1mach(4) 
               do i = 1, min(ni,2) 
                  write (form, 21) i, isizei 
   21             format('(11x,21hin above message, i',i1,'=,i',i2,')   ') 
                  if (i == 1) write (iunit, form) i1 
                  if (i == 2) write (iunit, form) i2 
               end do 
               do i = 1, min(nr,2) 
                  write (form, 23) i, isizef + 10, isizef 
   23             format('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')') 
                  if (i == 1) write (iunit, form) r1 
                  if (i == 2) write (iunit, form) r2 
               end do 
               if (lkntrl > 0) then 
!              error number
                  write (iunit, 30) lerr 
   30             format(' error number =',i10) 
               endif 
            end do 
         endif 
      endif 
      ifatal = 0 
      if (llevel==2 .or. llevel==1 .and. mkntrl==2) ifatal = 1 
!     quit here if message is not fatal
      if (ifatal <= 0) return  
      if (lkntrl>0 .and. kount<=max0(1,maxmes)) then 
!        print reason for abort
         if (llevel == 1) call xerprt ('job abort due to unrecovered error.', &
            35) 
         if (llevel == 2) call xerprt ('job abort due to fatal error.', 29) 
!        print error summary
         call xersav (' ', -1, 0, 0, kdummy) 
      endif 
      if (llevel==2 .and. kount>max0(1,maxmes)) lmessg = 0 
      call xerabt (messg, lmessg) 
      return  
      end subroutine xerrwv 
!deck xersav
      subroutine xersav(messg, nmessg, nerr, level, icount) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer nmessg, nerr, level, icount 
      character messg*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(5) :: lun 
      integer, dimension(10) :: nertab, levtab, kount 
      integer :: kountx, nunit, kunit, iunit, i, ii 
      character, dimension(10) :: mestab*20 
      character :: mes*20 

      save mestab, nertab, levtab, kount, kountx 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: i1mach 
!-----------------------------------------------
!***begin prologue  xersav
!***date written   800319   (yymmdd)
!***revision date  851213   (yymmdd)
!***category no.  r3
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  record that an error has occurred.
!***description
!
!     abstract
!        record that this error occurred.
!
!     description of parameters
!     --input--
!       messg, nmessg, nerr, level are as in xerror,
!       except that when nmessg=0 the tables will be
!       dumped and cleared, and when nmessg is less than zero the
!       tables will be dumped and not cleared.
!     --output--
!       icount will be the number of times this message has
!       been seen, or zero if the table has overflowed and
!       does not contain this message specifically.
!       when nmessg=0, icount will not be altered.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  1 august 1985
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  i1mach,xgetua
!***end prologue  xersav
!     next two data statements are necessary to provide a blank
!     error table initially
      data kount(1), kount(2), kount(3), kount(4), kount(5), kount(6), kount(7)&
         , kount(8), kount(9), kount(10)/0, 0, 0, 0, 0, 0, 0, 0, 0, 0/ 
      data kountx/0/ 
!***first executable statement  xersav
      if (nmessg <= 0) then 
!     dump the table
         if (kount(1) == 0) return  
!        print to each unit
         call xgetua (lun, nunit) 
         do kunit = 1, nunit 
            iunit = lun(kunit) 
            if (iunit == 0) iunit = i1mach(4) 
!           print table header
            write (iunit, 10) 
   10       format('0          error message summary'/&
               ' message start             nerr     level     count') 
!           print body of table
            do i = 1, 10 
               if (kount(i) == 0) exit  
               write (iunit, 15) mestab(i), nertab(i), levtab(i), kount(i) 
   15          format(1x,a20,3i10) 
            end do 
            if (kountx /= 0) write (iunit, 40) kountx 
   40       format('0other errors not individually tabulated=',i10) 
            write (iunit, 50) 
   50       format(1x) 
         end do 
         if (nmessg < 0) return  
!        clear the error tables
         kount = 0 
         i = 11 
         kountx = 0 
         return  
      endif 
      mes = messg 
      do i = 1, 10 
         ii = i 
         if (kount(i) == 0) go to 110 
         if (mes /= mestab(i)) cycle  
         if (nerr /= nertab(i)) cycle  
         if (level /= levtab(i)) cycle  
         go to 100 
      end do 
      kountx = kountx + 1 
      icount = 1 
      return  
!     message found in table
  100 continue 
      kount(ii) = kount(ii) + 1 
      icount = kount(ii) 
      return  
!     empty slot found for new message
  110 continue 
      mestab(ii) = mes 
      nertab(ii) = nerr 
      levtab(ii) = level 
      kount(ii) = 1 
      icount = 1 
      return  
      end subroutine xersav 
      subroutine xgetf(kontrl) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer kontrl 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xgetf
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  return the current value of the error control flag.
!***description
!
!   abstract
!        xgetf returns the current value of the error control flag
!        in kontrl.  see subroutine xsetf for flag value meanings.
!        (kontrl is an output parameter only.)
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xgetf
!***first executable statement  xgetf
      kontrl = j4save(2,0,.FALSE.) 
      return  
      end subroutine xgetf 
!deck xgetua
      subroutine xgetua(iunita, n) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n 
      integer, dimension(5) :: iunita 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, index 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xgetua
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  return unit number(s) to which error messages are being
!            sent.
!***description
!
!     abstract
!        xgetua may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        these unit numbers may have been set by a call to xsetun,
!        or a call to xsetua, or may be a default value.
!
!     description of parameters
!      --output--
!        iunit - an array of one to five unit numbers, depending
!                on the value of n.  a value of zero refers to the
!                default unit, as defined by the i1mach machine
!                constant routine.  only iunit(1),...,iunit(n) are
!                defined by xgetua.  the values of iunit(n+1),...,
!                iunit(5) are not defined (for n .lt. 5) or altered
!                in any way by xgetua.
!        n     - the number of units to which copies of the
!                error messages are being sent.  n will be in the
!                range from 1 to 5.
!
!     latest revision ---  19 mar 1980
!     written by ron jones, with slatec common math library subcommittee
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xgetua
!***first executable statement  xgetua
      n = j4save(5,0,.FALSE.) 
      do i = 1, n 
         index = i + 4 
         if (i == 1) index = 3 
         iunita(i) = j4save(index,0,.FALSE.) 
      end do 
      return  
      end subroutine xgetua 
!deck j4save
      integer function j4save (iwhich, ivalue, iset) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer iwhich, ivalue 
      logical iset 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(9) :: iparam 

      save iparam 
!-----------------------------------------------
!***begin prologue  j4save
!***refer to  xerror
!***routines called  (none)
!***description
!
!     abstract
!        j4save saves and recalls several global variables needed
!        by the library error handling routines.
!
!     description of parameters
!      --input--
!        iwhich - index of item desired.
!                = 1 refers to current error number.
!                = 2 refers to current error control flag.
!                 = 3 refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                 = 4 refers to the maximum number of times any
!                     message is to be printed (as set by xermax).
!                 = 5 refers to the total number of units to which
!                     each error message is to be written.
!                 = 6 refers to the 2nd unit for error messages
!                 = 7 refers to the 3rd unit for error messages
!                 = 8 refers to the 4th unit for error messages
!                 = 9 refers to the 5th unit for error messages
!        ivalue - the value to be set for the iwhich-th parameter,
!                 if iset is .true. .
!        iset   - if iset=.true., the iwhich-th parameter will be
!                 given the value, ivalue.  if iset=.false., the
!                 iwhich-th parameter will be unchanged, and ivalue
!                 is a dummy parameter.
!      --output--
!        the (old) value of the iwhich-th parameter will be returned
!        in the function value, j4save.
!
!     written by ron jones, with slatec common math library subcommittee
!    adapted from bell laboratories port library error handler
!     latest revision ---  1 august 1985
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***end prologue  j4save
      data iparam(1), iparam(2), iparam(3), iparam(4)/0, 2, 0, 10/ 
      data iparam(5)/1/ 
      data iparam(6), iparam(7), iparam(8), iparam(9)/0, 0, 0, 0/ 
!***first executable statement  j4save
      j4save = iparam(iwhich) 
      if (iset) iparam(iwhich) = ivalue 
      return  
      end function j4save 
!deck xerclr
      subroutine xerclr 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: junk 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xerclr
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  reset current error number to zero.
!***description
!
!     abstract
!        this routine simply resets the current error number to zero.
!        this may be necessary to do in order to determine that
!        a certain error has occurred again since the last time
!        numxer was referenced.
!
!     written by ron jones, with slatec common math library subcommittee
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xerclr
!***first executable statement  xerclr
      junk = j4save(1,0,.TRUE.) 
      return  
      end subroutine xerclr 
      subroutine xerdmp 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kount 
!-----------------------------------------------
!***begin prologue  xerdmp
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  print the error tables and then clear them.
!***description
!
!     abstract
!        xerdmp prints the error tables, then clears them.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  xersav
!***end prologue  xerdmp
!***first executable statement  xerdmp
      call xersav (' ', 0, 0, 0, kount) 
      return  
      end subroutine xerdmp 
      subroutine xermax(max) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer max 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: junk 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xermax
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  set maximum number of times any error message is to be
!            printed.
!***description
!
!     abstract
!        xermax sets the maximum number of times any message
!        is to be printed.  that is, non-fatal messages are
!        not to be printed after they have occured max times.
!        such non-fatal messages may be printed less than
!        max times even if they occur max times, if error
!        suppression mode (kontrl=0) is ever in effect.
!
!     description of parameter
!      --input--
!        max - the maximum number of times any one message
!              is to be printed.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xermax
!***first executable statement  xermax
      junk = j4save(4,max,.TRUE.) 
      return  
      end subroutine xermax 
      subroutine xgetun(iunit) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer iunit 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xgetun
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3c
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  return the (first) output file to which error messages
!            are being sent.
!***description
!
!     abstract
!        xgetun gets the (first) output file to which error messages
!        are being sent.  to find out if more than one file is being
!        used, one must use the xgetua routine.
!
!     description of parameter
!      --output--
!        iunit - the logical unit number of the  (first) unit to
!                which error messages are being sent.
!                a value of zero means that the default file, as
!                defined by the i1mach routine, is being used.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision --- 23 may 1979
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xgetun
!***first executable statement  xgetun
      iunit = j4save(3,0,.FALSE.) 
      return  
      end subroutine xgetun 
      subroutine xsetf(kontrl) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer kontrl 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: junk 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xsetf
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3a
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  set the error control flag.
!***description
!
!     abstract
!        xsetf sets the error control flag value to kontrl.
!        (kontrl is an input parameter only.)
!        the following table shows how each message is treated,
!        depending on the values of kontrl and level.  (see xerror
!        for description of level.)
!
!        if kontrl is zero or negative, no information other than the
!        message itself (including numeric values, if any) will be
!        printed.  if kontrl is positive, introductory messages,
!        trace-backs, etc., will be printed in addition to the message.
!
!              iabs(kontrl)
!        level        0              1              2
!        value
!          2        fatal          fatal          fatal
!
!          1     not printed      printed         fatal
!
!          0     not printed      printed        printed
!
!         -1     not printed      printed        printed
!                                  only           only
!                                  once           once
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save,xerrwv
!***end prologue  xsetf
!***first executable statement  xsetf
      if (kontrl<(-2) .or. kontrl>2) then 
         call xerrwv ('xsetf  -- invalid value of kontrl (i1).', 33, 1, 2, 1, &
            kontrl, 0, 0, 0., 0.) 
         return  
      endif 
      junk = j4save(2,kontrl,.TRUE.) 
      return  
      end subroutine xsetf 
      subroutine xsetua(iunita, n) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer n 
      integer, dimension(5) :: iunita 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, index, junk 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xsetua
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3b
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  set logical unit numbers (up to 5) to which error
!            messages are to be sent.
!***description
!
!     abstract
!        xsetua may be called to declare a list of up to five
!        logical units, each of which is to receive a copy of
!        each error message processed by this package.
!        the purpose of xsetua is to allow simultaneous printing
!        of each error message on, say, a main output file,
!        an interactive terminal, and other files such as graphics
!        communication files.
!
!     description of parameters
!      --input--
!        iunit - an array of up to five unit numbers.
!                normally these numbers should all be different
!                (but duplicates are not prohibited.)
!        n     - the number of unit numbers provided in iunit
!                must have 1 .le. n .le. 5.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save,xerrwv
!***end prologue  xsetua
!***first executable statement  xsetua
      if (n<1 .or. n>5) then 
         call xerrwv ('xsetua -- invalid value of n (i1).', 34, 1, 2, 1, n, 0, &
            0, 0., 0.) 
         return  
      endif 
      do i = 1, n 
         index = i + 4 
         if (i == 1) index = 3 
         junk = j4save(index,iunita(i),.TRUE.) 
      end do 
      junk = j4save(5,n,.TRUE.) 
      return  
      end subroutine xsetua 
      subroutine xsetun(iunit) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer iunit 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: junk 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      integer , EXTERNAL :: j4save 
!-----------------------------------------------
!***begin prologue  xsetun
!***date written   790801   (yymmdd)
!***revision date  851111   (yymmdd)
!***category no.  r3b
!***keywords  error,xerror package
!***author  jones, r. e., (snla)
!***purpose  set output file to which error messages are to be sent.
!***description
!
!     abstract
!        xsetun sets the output file to which error messages are to
!        be sent.  only one file will be used.  see xsetua for
!        how to declare more than one file.
!
!     description of parameter
!      --input--
!        iunit - an input parameter giving the logical unit number
!                to which error messages are to be sent.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
!                 handling package', sand82-0800, sandia laboratories,
!                 1982.
!***routines called  j4save
!***end prologue  xsetun
!***first executable statement  xsetun
      junk = j4save(3,iunit,.TRUE.) 
      junk = j4save(5,1,.TRUE.) 
      return  
      end subroutine xsetun 
      REAL FUNCTION RAND (R) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  15:55:30  12/05/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL R 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA1, IA0, IA1MA0, IC, IX1, IX0, IY0, IY1 

      SAVE IA1, IA0, IA1MA0, IC, IX1, IX0 
!-----------------------------------------------
!***BEGIN PROLOGUE  RAND
!***DATE WRITTEN   770401   (YYMMDD)
!***REVISION DATE  861211   (YYMMDD)
!***CATEGORY NO.  L6A21
!***KEYWORDS  LIBRARY=SLATEC(FNLIB),TYPE=SINGLE PRECISION(RAND-S),
!             RANDOM NUMBER,SPECIAL FUNCTIONS,UNIFORM
!***AUTHOR  FULLERTON, W., (LANL)
!***PURPOSE  Generates a uniformly distributed random number.
!***DESCRIPTION
!
!      This pseudo-random number generator is portable among a wide
! variety of computers.  RAND(R) undoubtedly is not as good as many
! readily available installation dependent versions, and so this
! routine is not recommended for widespread usage.  Its redeeming
! feature is that the exact same random numbers (to within final round-
! off error) can be generated from machine to machine.  Thus, programs
! that make use of random numbers can be easily transported to and
! checked in a new environment.
!      The random numbers are generated by the linear congruential
! method described, e.g., by Knuth in Seminumerical Methods (p.9),
! Addison-Wesley, 1969.  Given the I-th number of a pseudo-random
! sequence, the I+1 -st number is generated from
!             X(I+1) = (A*X(I) + C) MOD M,
! where here M = 2**22 = 4194304, C = 1731 and several suitable values
! of the multiplier A are discussed below.  Both the multiplier A and
! random number X are represented in double precision as two 11-bit
! words.  The constants are chosen so that the period is the maximum
! possible, 4194304.
!      In order that the same numbers be generated from machine to
! machine, it is necessary that 23-bit integers be reducible modulo
! 2**11 exactly, that 23-bit integers be added exactly, and that 11-bit
! integers be multiplied exactly.  Furthermore, if the restart option
! is used (where R is between 0 and 1), then the product R*2**22 =
! R*4194304 must be correct to the nearest integer.
!      The first four random numbers should be .0004127026,
! .6750836372, .1614754200, and .9086198807.  The tenth random number
! is .5527787209, and the hundredth is .3600893021 .  The thousandth
! number should be .2176990509 .
!      In order to generate several effectively independent sequences
! with the same generator, it is necessary to know the random number
! for several widely spaced calls.  The I-th random number times 2**22,
! where I=K*P/8 and P is the period of the sequence (P = 2**22), is
! still of the form L*P/8.  In particular we find the I-th random
! number multiplied by 2**22 is given by
! I   =  0  1*P/8  2*P/8  3*P/8  4*P/8  5*P/8  6*P/8  7*P/8  8*P/8
! RAND=  0  5*P/8  2*P/8  7*P/8  4*P/8  1*P/8  6*P/8  3*P/8  0
! Thus the 4*P/8 = 2097152 random number is 2097152/2**22.
!      Several multipliers have been subjected to the spectral test
! (see Knuth, p. 82).  Four suitable multipliers roughly in order of
! goodness according to the spectral test are
!    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5
!    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5
!    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5
!    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1
!
!      In the table below LOG10(NU(I)) gives roughly the number of
! random decimal digits in the random numbers considered I at a time.
! C is the primary measure of goodness.  In both cases bigger is better.
!
!                   LOG10 NU(I)              C(I)
!       A       I=2  I=3  I=4  I=5    I=2  I=3  I=4  I=5
!
!    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
!    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
!    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
!    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
!   Best
!    Possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
!
!             Input Argument --
! R      If R=0., the next random number of the sequence is generated.
!        If R .LT. 0., the last generated number will be returned for
!          possible use in a restart procedure.
!        If R .GT. 0., the sequence of random numbers will start with
!          the seed R mod 1.  This seed is also returned as the value of
!          RAND provided the arithmetic is done exactly.
!
!             Output Value --
! RAND   a pseudo-random number between 0. and 1.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  RAND
      DATA IA1, IA0, IA1MA0/1536, 1029, 507/ 
      DATA IC/1731/ 
      DATA IX1, IX0/0, 0/ 
!***FIRST EXECUTABLE STATEMENT  RAND
      IF (R >= 0.) THEN 
         IF (R > 0.) GO TO 20 
!
!           A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
!                   + IA0*IX0) + IA0*IX0
!
         IY0 = IA0*IX0 
         IY1 = IA1*IX1 + IA1MA0*(IX0 - IX1) + IY0 
         IY0 = IY0 + IC 
         IX0 = MOD(IY0,2048) 
         IY1 = IY1 + (IY0 - IX0)/2048 
         IX1 = MOD(IY1,2048) 
!
      ENDIF 
   10 CONTINUE 
      RAND = IX1*2048 + IX0 
      RAND = RAND/4194304. 
      RETURN  
!
   20 CONTINUE 
      IX1 = AMOD(R,1.)*4194304. + 0.5 
      IX0 = MOD(IX1,2048) 
      IX1 = (IX1 - IX0)/2048 
      GO TO 10 
!
      END FUNCTION RAND 
