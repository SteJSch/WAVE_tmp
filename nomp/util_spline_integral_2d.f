*CMZ :          16/08/2024  10.03.52  by  Michael Scheer
*CMZ :  4.01/03 17/05/2023  11.24.58  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  11.49.33  by  Michael Scheer
*CMZ : 00.00/16 21/11/2014  14.53.59  by  Michael Scheer
*-- Author :    Michael Scheer   21/11/2014
      subroutine util_spline_integral_2d(nx,ny,x,y,f,result,istat,kalloc)
*KEEP,gplhint.
*KEND.

      implicit none

      double precision x(nx),y(ny),f(nx,ny),result
      integer :: istat,nx,ny,ix,iy,kstat,kalloc,nxyo=0,kallo=0

      double precision, allocatable :: fb(:),fb2(:),coef(:),
     &  w1(:),w2(:),w3(:),w4(:)

      save

      if (kalloc.eq.0) then
        if (kallo.eq.0) then
          kalloc=1
        else
          if (max(nx,ny).gt.nxyo) then
            deallocate(fb)
            deallocate(fb2)
            deallocate(coef)
            deallocate(w1)
            deallocate(w2)
            deallocate(w3)
            deallocate(w4)
          endif
          kalloc=1
        endif
      endif

      if (kalloc.gt.0) then
        allocate(fb(max(nx,ny)))
        allocate(fb2(max(nx,ny)))
        allocate(coef(max(nx,ny)))
        allocate(w1(max(nx,ny)))
        allocate(w2(max(nx,ny)))
        allocate(w3(max(nx,ny)))
        allocate(w4(max(nx,ny)))
        kallo=1
        kalloc=0
      else if (kalloc.lt.0) then
        deallocate(fb)
        deallocate(fb2)
        deallocate(coef)
        deallocate(w1)
        deallocate(w2)
        deallocate(w3)
        deallocate(w4)
        return
      endif

      kstat=0

      if (ny.gt.nx) then
        do ix=1,nx
          fb(1:ny)=f(ix,1:ny)
          call util_spline_integral_stat(y,fb,ny,fb2(ix)
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
        call util_spline_integral_stat(x,fb2,nx,result
     &    ,coef,w1,w2,w3,w4,istat)
        kstat=kstat+istat
      else !nx.gt.ny?
        do iy=1,ny
          fb(1:nx)=f(1:nx,iy)
          call util_spline_integral_stat(x,fb,nx,fb2(iy)
     &      ,coef,w1,w2,w3,w4,istat)
          kstat=kstat+istat
        enddo
        call util_spline_integral_stat(y,fb2,ny,result
     &    ,coef,w1,w2,w3,w4,istat)
        kstat=kstat+istat
      endif !nx.gt.ny

      nxyo=max(nx,ny)

      return
      end
