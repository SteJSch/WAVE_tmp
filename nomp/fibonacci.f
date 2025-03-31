*CMZ :  4.00/04 17/05/2019  14.22.20  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.48.55  by  Michael Scheer
*CMZ :  2.63/04 22/07/2009  07.39.04  by  Michael Scheer
*CMZ :  2.63/03 15/05/2009  12.38.54  by  Michael Scheer
*CMZ :  2.63/02 05/03/2008  17.56.28  by  Michael Scheer
*-- Author :    Michael Scheer   30/01/2008
      subroutine fibonacci(npoles,xpol,eta,icutfibo)
*KEEP,gplhint.
*KEND.

c calculates Fibonacci-series for quasi-periodic undulator according
c to program qpu_field of Sasaki and Bahrdt

c xpol contains on return poles to be shimmed

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      integer npoles,nmax,m,icuts,ii,iisum,k,j,i,ifibo1,ifibo2,ifibo,
     &  icutfibo,ifind,kcut

      double precision eta
      real xpol(npoles)

      integer, dimension(:), allocatable :: iflag,iiflag

      nmax=10*npoles

      allocate(iflag(nmax))
      allocate(iiflag(nmax))

      xpol=0.0

      open(unit=99,file='fibonacci_available_cuts.dat',status='unknown')
      write(99,*)'* ',icode
      write(99,*)'* ',code

      ifibo1=0
      do m=2, nmax
        ifibo2=INT(m/(1.0d0+eta))+1
        ifibo=ifibo2-ifibo1
        if(ifibo.lt.1) then
          iflag(m-1)=0
        else
          iflag(m-1)=(-1)**(m-1)
        endif
        ifibo1=ifibo2
      enddo

c------------------------ start loop over quais periodc stuctures ------
      icuts=nmax-npoles+1
      ifind=0
      kcut=0

      do i=1,icuts

        ii=0
        iisum=0

        do j=i,i+npoles-1
          if(abs(iflag(j)).ge.0.5)then
            ii=ii+1
            iisum=iisum+iflag(j)
            iiflag(ii)=(j-i+1)*iflag(j)
          endif
        enddo

        if(iisum.eq.0.and.mod(ii,2).eq.0.and.ii.gt.0) then

          if (i.ge.icutfibo.and.kcut.eq.0) then
            kcut=1
            icutfibo=i
            ifind=1
            do k=1,ii
              xpol(abs(iiflag(k)))=iiflag(k)
            enddo
          endif

          write(99,*)i

        endif

      enddo    ! i1=1,icuts

9999  if (ifind.eq.0) then

        icutfibo=-1

*KEEP,QPU.
*KEND.
        if (iuout.ne.0) then
          call uout
        endif

        write(6,*)
        write(6,*)'*** Warning in FIBONACCI:'
        write(6,*)'*************************************'
        write(6,*)' did not find quasiperiodic structure'
        write(6,*)'*************************************'
        write(6,*)

        write(lungfo,*)
        write(lungfo,*)'*** Warning in FIBONACCI:'
        write(lungfo,*)'*************************************'
        write(lungfo,*)' did not find quasiperiodic structure'
        write(lungfo,*)'*************************************'
        write(lungfo,*)

        xpol=0.

        if (ibatch.eq.0) then
          call sleep(3)
        endif

      endif

      close(99)

      open(unit=99,file='fibonacci_used_cut.dat',status='unknown')
      write(99,*)'* ',icode
      write(99,*)'* ',code
      do k=1,npoles
        write(99,*)k,xpol(k)
      enddo
      close(99)

      deallocate(iflag)
      deallocate(iiflag)

      return
      end
