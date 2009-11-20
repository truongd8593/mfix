! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Random Number Generation Utilities                     C
!>  Purpose: Removed from interpolation mod and added built-in random   
!>           number routines instead of Pope's                          
!                                                                      C
!                                                                      C
!  Author: Sreekanth Pannala and Rahul Garg           Date: 23-Oct-08  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE randomno
      USE constant
      Implicit None
      Private
      Public :: uni_rno, nor_rno
Contains

      SUBROUTINE UNI_RNO(Y)

            implicit none
            double precision, intent(out), dimension(:) :: y
            double precision rmean, variance, sigma
            integer i, nsize

            nsize = size(y(:))

            call init_random_seed
            call random_number(y)

            rmean = sum(y(:))/nsize

            write(*,*) 'Generating Uniform Random Variables for size', nsize
            write(*,*) 'mean', rmean

            variance = 0.0
            do i = 1, nsize
               write(20,*) i, y(i)
               variance = variance + (y(i)-rmean)**2
            end do

            close(20)

            variance = variance/nsize
            sigma = sqrt(variance)

            write(*,*) 'sigma', sigma

            RETURN
      END SUBROUTINE UNI_RNO



      SUBROUTINE NOR_RNO(Y, mean, sigma)

            implicit none
            double precision, intent(out), dimension(:) :: y
            double precision mean, sigma

            double precision lmean, lvariance, lsigma
            double precision x(2), w
            integer i, nsize, n
            ! no. of times this routine has been called (should = dimn)
            integer, save :: COUNTER = 0
            ! so all components are written
            integer fileunit 

            COUNTER = COUNTER + 1
            fileunit = 20 + COUNTER

            nsize = size(y(:))

            call init_random_seed

            do i = 1, ceiling(real(nsize/2))
               do n = 1,100000
                  call random_number(x)
                  x = 2.0 * x - 1.0
                  w = x(1)**2 + x(2)**2
                  if(w.lt.1.0) exit
               end do

               w = sqrt( (-2.0 * log( w ) ) / w )
               y(2*i-1) = x(1) * w * sigma + mean
               if(2*i.lt.nsize) y(2*i) = x(2) * w * sigma + mean
            end do

            lmean = sum(y(:))/nsize

            write(*,'(7X,A)') 'Generating Normal Random Variables'
            write(*,'(7X,A,2X,ES)') 'specified mean =', mean
            write(*,'(7X,A,2X,ES)') 'computed mean =', lmean

            write(fileunit,'(A)') 'FROM NOR_RNO'
! specific to the call from init_particles_jn            
            write(fileunit,'(A,I5,A)') 'FOR DIRECTION = ', &
               COUNTER, ' where (1=X,2=Y,3=Z)'
            write(fileunit,'(5X,A,5X,A)') 'particle no.', 'velocity component'

            lvariance = 0.0
            do i = 1, nsize
               write(fileunit,'(I,5X,F)') i, y(i)
               lvariance = lvariance + (y(i)-lmean)**2
            end do

            close(fileunit)

            lvariance = lvariance/nsize
            lsigma = sqrt(lvariance)

            write(*,'(7X,A,2X,ES)') 'specified sigma =', sigma
            write(*,'(7X,A,2X,ES)') 'computed sigma =', lsigma

            RETURN
      END SUBROUTINE NOR_RNO

      SUBROUTINE init_random_seed
            INTEGER              :: isize,idate(8)
            INTEGER,ALLOCATABLE  :: iseed(:)

            CALL DATE_AND_TIME(VALUES=idate)
            CALL RANDOM_SEED(SIZE=isize)
            ALLOCATE( iseed(isize) )
            CALL RANDOM_SEED(GET=iseed)
            iseed = iseed * (idate(8)-500) ! idate(8) contains millisecond
            CALL RANDOM_SEED(PUT=iseed)

            DEALLOCATE( iseed )

      END SUBROUTINE init_random_seed

End MODULE randomno

