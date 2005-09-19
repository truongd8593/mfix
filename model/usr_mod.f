

      MODULE usr


        Use param
        Use param1


!
!       Include user-defined variables in this module.  To access the variables
!       from a subroutine add the statement "Use usr".  If allocatable arrays
!       are defined in this module allocate them in usr0.  To turn on the
!       user defined subroutines (usr0, usr1, and usr2) set call_usr to true in
!       mfix.dat.
!
!                        a dummy variable to keep the compiler happy.                     
        DOUBLE PRECISION usr_dummy

      END MODULE usr                                                                             
