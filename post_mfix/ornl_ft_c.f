!
! *********************************************************************
! ******************************* do_power_spectrum ************************
! *********************************************************************
!
      subroutine do_power_spectrum(time_series,nt)
      use usr_input

      implicit none
      

      integer*4 :: nu, nb, nrm, ibeg, iend
      integer :: nt,k 

      real*8  :: sr, ol, hng
      real*8    :: time_series(*)

      real*8, allocatable    :: f(:), p(:)

      allocate (f(nt))
      allocate (p(nt))

      ibeg = 1
      iend = nt

      write(*,*) 'Data sampling rate [samples/time unit]-suggested=2'
      read(*,*) sr
      write(*,*) 'binary power of FFT block size (block=2**nu) '
      write(*,*) '(typically, this depends on the desired frequency '
      write(*,*) 'resolution, but 10-13 are useful values; frequency '
      write(*,*) 'resolution (PSD bin width) is  df = sr / 2**(nu-1))'
      read(*,*) nu
      write(*,*) 'number of averaging blocks (1 is minimum, 15 is a useful max.; the constraint is nb <= fix(length(TS)/(2**nu)).)'
      read(*,*) nb
      write(*,*) 'percentage overlap in blocks [0,10]; typically, use 0'
      read(*,*) ol
      write(*,*) 'Hanning-window power (0 -> no windowing, use 2-5 )'
      read(*,*) hng
      write(*,*) 'normalization of power (1=sum of power, 2=data variance)'
      read(*,*) nrm

      if(nb*(2**nu).gt.nt) then
         write(*,*) 'TS is not long enough for nb*n (n=2**nu)'
         write(*,*) 'nb =', nb
         write(*,*) '2**nu =', 2**nu
         write(*,*) 'length(TS) =', nt
         return
      endif


      call psd(time_series, sr, nu, nb, ol, hng, nrm, f, p)

 !
      if (usr_fname(1:1) .eq. '*') then
         write (*,*) ' '
         write (*,*) ' doing power spectral density:'
         write (*,*)  ' '
      else
         write (40,*) ' '
         write (40,*) ' doing power spectral density:'
         write (40,*)  ' '
      end if

      call usr_write_input(1)

      do k = 1, 2** (nu-1)
         if (usr_fname(1:1) .eq. '*') then
            write (*,*) f(k), p(k)
         else
            write (40,*) f(k), p(k)
         end if
      end do ! k
!
      deallocate (f)
      deallocate (p)
 !
      return
      end


