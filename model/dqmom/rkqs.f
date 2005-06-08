       SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext)
       IMPLICIT NONE
       INTEGER n,NMAX
       double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
       PARAMETER (NMAX=50)
       INTEGER i
       double precision errmax,h,htemp,xnew
       double precision yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,&
       PSHRNK,ERRCON
       PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
       h=htry
       do 
         call rkck(y,dydx,n,x,h,ytemp,yerr)
         errmax=0.
         do i=1,n
          errmax=max(errmax,abs(yerr(i)/yscal(i)))
         end do
         errmax=errmax/eps
         if(errmax<=1.0) exit
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1*abs(h)),h)
         xnew=x+h
         if(xnew==x)pause 'stepsize underflow in rkqs'
       end do 
     
        if(errmax>ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.0*h
        endif
        hdid=h
        x=x+h
          y(:)=ytemp(:)
        return
      END
