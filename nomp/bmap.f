*CMZ :  4.00/16 09/08/2022  09.07.08  by  Michael Scheer
*CMZ :  4.00/07 07/06/2020  15.15.28  by  Michael Scheer
*CMZ :  3.05/05 13/07/2018  11.51.31  by  Michael Scheer
*CMZ :  3.03/04 29/11/2017  10.21.49  by  Michael Scheer
*CMZ :  3.02/00 10/09/2014  14.10.09  by  Michael Scheer
*CMZ :  3.01/04 26/05/2014  16.15.21  by  Michael Scheer
*CMZ :  3.01/00 06/05/2013  09.13.42  by  Michael Scheer
*CMZ :  2.68/05 01/10/2012  14.11.38  by  Michael Scheer
*CMZ :  2.68/04 04/09/2012  09.22.40  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.01.42  by  Michael Scheer
*CMZ :  2.68/02 02/07/2012  13.51.47  by  Michael Scheer
*-- Author :    Michael Scheer   18/06/2012
C**********************************************************************
      subroutine bmap(xin,yin,zin,bxout,byout,bzout)
C**********************************************************************
*KEEP,gplhint.
*KEND.

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      double precision, dimension (:,:), allocatable :: bmappe

      double precision xin,yin,zin,x,y,z,bx,by,bz,step,stepy,xold,yold,
     &  bxout,byout,bzout,
c     &  bmxmin,bmxmax,bmymin,bmymax,bmzmin,bmzmax,
c     &  bmbxmin,bmbxmax,bmbymin,bmbymax,bmbzmin,bmbzmax,
     &  offsetx,offsety,offsetz,offsetbx,offsetby,offsetbz,
     &  scalex,scaley,scalez,scalebx,scaleby,scalebz,
     &  bmapdy,bmapdz,x1,x2,a3(3),x3(3),b3(3),
     &  b11,b12,b13,b21,b22,b23,b31,b32,b33,bb1,bb2,bb3,dy,dz

      double precision b111(3),b211(3), b121(3),b221(3)
      double precision b112(3),b212(3), b122(3),b222(3)
      double precision b112111(3),b212211(3), b122121(3),b222221(3)
      double precision blow(3),bhig(3),b(3),dxx,dyy,dzz,xx,yy,zz,y1,z1

      real*8 :: eps=1.0d-6

      integer ical,iwarnx,iwarny,iwarnz,nx,ny,nz,last,ntot,ianf,
     &  kd,ix1,ix2,iy1,iy2,iz1,iz2,nyz,i,k,
     &  kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3

      character(2048) cline

      data ical/0/
      data iwarnx/0/
      data iwarny/0/
      data iwarnz/0/
      data scalex,scaley,scalez/1.0d0,1.0d0,1.0d0/
      data scalebx,scaleby,scalebz/1.0d0,1.0d0,1.0d0/
      data xold/1.0d30/
      data yold/1.0d30/

      save

      iwarnbmap=0

      if (ical.eq.0) then

        write(lungfo,*)
        write(lungfo,*)'      Subroutine BMAP:'
        write(lungfo,*)
        write(lungfo,*)'      Reading file'
        write(lungfo,*)'      ',fileb0
        write(lungfo,*)

        bmxmin= 1.0d99
        bmxmax=-1.0d99
        bmymin= 1.0d99
        bmymax=-1.0d99
        bmzmin= 1.0d99
        bmzmax=-1.0d99

        bmbxmin= 1.0d99
        bmbxmax=-1.0d99
        bmbymin= 1.0d99
        bmbymax=-1.0d99
        bmbzmin= 1.0d99
        bmbzmax=-1.0d99

        offsetx=0.0d0
        offsety=0.0d0
        offsetz=0.0d0

        open(unit=lunb0,file=fileb0,status='old')
        ntot=0
        nx=-1
        ny=-1
        nz=0
 1      read(lunb0,'(a)',end=9) cline
        last=len_trim(cline)

        if (last.le.1.or.
     &      cline(1:1).eq.'%'.or.
     &      cline(1:1).eq.'!'.or.
     &      cline(1:1).eq.'#'.or.
     &      cline(1:1).eq.'*'.or.
     &      cline(1:1).eq.'@'.or.
     &      cline(1:2).eq.' %'.or.
     &      cline(1:2).eq.' !'.or.
     &      cline(1:2).eq.' #'.or.
     &      cline(1:2).eq.' *'.or.
     &      cline(1:2).eq.' @'
     &      ) then

          write(lungfo,*) cline(1:last)

          ianf=index(cline,"scaling")
          if (ianf.gt.0) then
            ianf=index(cline,"=")
            read(cline(ianf+1:last),*)scalex,scaley,scalez,
     &        scalebx,scaleby,scalebz
          endif

          ianf=index(cline,"offset")
          if (ianf.gt.0) then
            ianf=index(cline,"=")
            read(cline(ianf+1:last),*)offsetx,offsety,offsetz,
     &        offsetbx,offsetby,offsetbz
          endif

        else

          ntot=ntot+1
          read(cline(1:last),*)x,y,z,bx,by,bz

          x=x*scalex+offsetx
          y=y*scaley+offsety
          z=z*scalez+offsetz
          bx=bx*scalebx+offsetbx
          by=by*scaleby+offsetby
          bz=bz*scalebz+offsetbz

          if (x.lt.bmxmin) bmxmin=x
          if (x.gt.bmxmax) bmxmax=x
          if (y.lt.bmymin) bmymin=y
          if (y.gt.bmymax) bmymax=y
          if (z.lt.bmzmin) bmzmin=z
          if (z.gt.bmzmax) bmzmax=z
          if (bx.lt.bmbxmin) bmbxmin=bx
          if (bx.gt.bmbxmax) bmbxmax=bx
          if (by.lt.bmbymin) bmbymin=by
          if (by.gt.bmbymax) bmbymax=by
          if (bz.lt.bmbzmin) bmbzmin=bz
          if (bz.gt.bmbzmax) bmbzmax=bz

          if (ntot.eq.1) then
            xold=x
            yold=y
          endif

          if (ntot.eq.2.and.abs(x-xold).gt.eps) then
            write(6,*)'*** Error in BMAP: Bad file format! ***'
            write(6,*)'*** x must run latest! ***'
            write(lungfo,*)'*** Error in BMAP: Bad file format! ***'
            write(lungfo,*)'*** x must run latest! ***'
            stop '*** WAVE aborted ***'
          endif

          if (nx.eq.-1) then
            if (ny.eq.-1) then
              if (abs(y-yold).le.eps) then
                nz=nz+1
              else
                ny=nz+1
              endif
            else if (abs(x-xold).le.eps) then
              ny=ny+1
            else
              ny=ny/nz
              nx=0
            endif !ny
          endif !nx

        endif !line type

        goto 1
 9      rewind(lunb0)

        nx=ntot/(ny*nz)

        if (nx.lt.2.and.irfilb0.eq.6.or.nx.lt.3.and.irfilb0.eq.-6) then
          write(lungfo,*)'*** Error in BMAP: Too few data for field map on'
          write(lungfo,*)fileb0
          write(lungfo,*)'*** Program WAVE aborted ***'
          write(6,*)'*** Error in BMAP: Too few data for field map on'
          write(6,*)fileb0
          write(6,*)'*** Program WAVE aborted ***'
          stop
        else
          allocate(bmappe(6,ntot))
        endif

        ntot=0
 11     read(lunb0,'(a)',end=99) cline

        last=len_trim(cline)

        if (last.le.1.or.
     &      cline(1:1).eq.'%'.or.
     &      cline(1:1).eq.'!'.or.
     &      cline(1:1).eq.'#'.or.
     &      cline(1:1).eq.'*'.or.
     &      cline(1:1).eq.'@'.or.
     &      cline(1:2).eq.' %'.or.
     &      cline(1:2).eq.' !'.or.
     &      cline(1:2).eq.' #'.or.
     &      cline(1:2).eq.' *'.or.
     &      cline(1:2).eq.' @'
     &      ) then
        else
          ntot=ntot+1
          read(cline(1:last),*) x,y,z,bx,by,bz
          x=x*scalex+offsetx
          y=y*scaley+offsety
          z=z*scalez+offsetz
          bx=bx*scalebx+offsetbx
          by=by*scaleby+offsetby
          bz=bz*scalebz+offsetbz
          bmappe(1,ntot)=x
          bmappe(2,ntot)=y
          bmappe(3,ntot)=z
          bmappe(4,ntot)=bx
          bmappe(5,ntot)=by
          bmappe(6,ntot)=bz
        endif
        goto 11
 99     close(lunb0)

        step=1.0d0/myinum
        stepy=step/10.0d0
        bmapdy=1.0d0
        if (ny.gt.1) bmapdy=(bmymax-bmymin)/(ny-1)
        bmapdz=1.0d0
        if (nz.gt.1) bmapdz=(bmzmax-bmzmin)/(nz-1)

        ix1=1
        ix2=ntot
        nyz=ny*nz

        write(lungfo,*)
        write(lungfo,*)' nx, ny, nz, and number of data lines:',
     &    nx,ny,nz,ntot
        write(lungfo,*)' xmin, xmax                   :',sngl(bmxmin),sngl(bmxmax)
        write(lungfo,*)' ymin, ymax , step size of map:',sngl(bmymin),sngl(bmymax),sngl(bmapdy)
        write(lungfo,*)' zmin, zmax , step size of map:',sngl(bmzmin),sngl(bmzmax),sngl(bmapdz)
        write(lungfo,*)' Bxmin, Bxmax of map:',sngl(bmbxmin),sngl(bmbxmax)
        write(lungfo,*)' Bymin, Bymax of map:',sngl(bmbymin),sngl(bmbymax)
        write(lungfo,*)' Bzmin, Bzmax of map:',sngl(bmbzmin),sngl(bmbzmax)
        write(lungfo,*)

        ical=1
      endif

      if (bxout.eq.-9999.0d0) then
        xin=bmxmin
        return
      else if (bxout.eq.9999.0d0) then
        xin=bmxmax
        return
      endif

      x=xin
      y=yin
      z=zin

      IF (IWARNX.EQ.0.AND.(X.LT.BMXMIN-STEP.OR.X.GT.BMXMAX+STEP)) THEN
        WRITE(6,*)'*** WARNING: IN BMAP: X OUT OF RANGE'
        WRITE(6,*)'X:',X
        WRITE(6,*)'Y:',Y
        WRITE(6,*)'Z:',Z
        WRITE(LUNGFO,*)'*** WARNING IN BMAP: X OUT OF RANGE'
        WRITE(LUNGFO,*)'X:',X
        WRITE(LUNGFO,*)'Y:',Y
        WRITE(LUNGFO,*)'Z:',Z
        IWARNX=1
      ENDIF

      IF (IWARNX.NE.0.AND.(X.LT.BMXMIN-STEP.OR.X.GT.BMXMAX+STEP)) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        iwarnbmap=1
        RETURN
      ENDIF

      IF (ny.gt.1.and.(Y.LT.BMYMIN-STEPy.OR.Y.GT.BMYMAX+STEPy)) THEN
        if (iwarny.eq.0) then
          WRITE(6,*)'*** ERROR IN BMAP: Y OUT OF RANGE'
          WRITE(6,*)'X:',X
          WRITE(6,*)'Y:',Y
          WRITE(6,*)'Z:',Z
          WRITE(LUNGFO,*)'*** ERROR IN BMAP: Y OUT OF RANGE'
          WRITE(LUNGFO,*)'X:',X
          WRITE(LUNGFO,*)'Y:',Y
          WRITE(LUNGFO,*)'Z:',Z
          iwarny=1
        endif
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        iwarnbmap=1
        RETURN
c        STOP
      ENDIF

      IF (nz.gt.1.and.(Z.LT.BMZMIN-STEP.OR.Z.GT.BMZMAX+STEP)) THEN
        if (iwarnz.eq.0) then
          WRITE(6,*)'*** ERROR IN BMAP: Z OUT OF RANGE'
          WRITE(6,*)'X:',X
          WRITE(6,*)'Y:',Y
          WRITE(6,*)'Z:',Z
          WRITE(LUNGFO,*)'*** ERROR IN BMAP: Z OUT OF RANGE'
          WRITE(LUNGFO,*)'X:',X
          WRITE(LUNGFO,*)'Y:',Y
          WRITE(LUNGFO,*)'Z:',Z
          iwarnz=1
        endif
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        iwarnbmap=1
        RETURN
c        STOP
      ENDIF

      if (x.lt.bmxmin) then
        ix1=1
        ix2=2
      else if (x.gt.bmxmax) then
        ix1=nx-1
        ix2=nx
      else

        if (x.ge.bmappe(1,(ix1-1)*nyz+1)) then
c hunt up
          kd=1
111       ix2=min(ix1+kd,nx)
          if (x.gt.bmappe(1,nyz*(ix2-1)+1)) then
            kd=2*kd
            ix1=ix2
            goto 111
          endif
        else    !(x.ge.bmappe(1,ix1))
c hunt down
          kd=1
          ix2=ix1
22        ix1=max(ix2-kd,1)
          if (x.lt.bmappe(1,(ix1-1)*nyz+1)) then
            kd=2*kd
            ix2=ix1
            goto 22
          endif
        endif

1111    if (ix2-ix1.gt.1) then
          k=(ix2+ix1)/2
          if(bmappe(1,(k-1)*nyz+1).gt.x)then
            ix2=k
          else
            ix1=k
          endif
          goto 1111
        endif
      endif

      x1=bmappe(1,(ix1-1)*nyz+1)
      x2=bmappe(1,(ix2-1)*nyz+1)
      dxx=(x-x1)/(x2-x1)

      if (ny.gt.1) then
        if (y.lt.bmymin) then
          iy1=1
        else
          iy1=int((y-bmymin)/bmapdy)+1
          iy1=max(1,iy1)
        endif
        iy2=iy1+1
        if (iy2.gt.ny) then
          iy2=ny
          iy1=iy2-1
        endif
      else
        iy1=1
        iy2=1
      endif

      if (nz.gt.1) then
        if (z.lt.bmzmin) then
          iz1=1
        else
          iz1=int((z-bmzmin)/bmapdz)+1
          iz1=max(1,iz1)
        endif
        iz2=iz1+1
        if (iz2.gt.nz) then
          iz2=nz
          iz1=iz2-1
        endif
      else
        iz1=1
        iz2=1
      endif

      if (irfilb0.eq.6) then

        dyy=(y-(bmymin+(iy1-1)*bmapdy))/bmapdy
        dzz=(z-(bmzmin+(iz1-1)*bmapdz))/bmapdz

        do i=1,3

          b111(i)=bmappe(3+i,iz1+(iy1-1)*nz+(ix1-1)*nyz)
          b211(i)=bmappe(3+i,iz2+(iy1-1)*nz+(ix1-1)*nyz)
          b121(i)=bmappe(3+i,iz1+(iy2-1)*nz+(ix1-1)*nyz)
          b221(i)=bmappe(3+i,iz2+(iy2-1)*nz+(ix1-1)*nyz)
          b112(i)=bmappe(3+i,iz1+(iy1-1)*nz+(ix2-1)*nyz)
          b212(i)=bmappe(3+i,iz2+(iy1-1)*nz+(ix2-1)*nyz)
          b122(i)=bmappe(3+i,iz1+(iy2-1)*nz+(ix2-1)*nyz)
          b222(i)=bmappe(3+i,iz2+(iy2-1)*nz+(ix2-1)*nyz)

          b112111(i)=b111(i)+(b112(i)-b111(i))*dxx
          b122121(i)=b121(i)+(b122(i)-b121(i))*dxx
          b212211(i)=b211(i)+(b212(i)-b211(i))*dxx
          b222221(i)=b221(i)+(b222(i)-b221(i))*dxx

          blow(i)=b112111(i)+(b212211(i)-b112111(i))*dzz
          bhig(i)=b122121(i)+(b222221(i)-b122121(i))*dzz

          b(i)=blow(i)+(bhig(i)-blow(i))*dyy

        enddo

      else !ifilb0

        if (ix1.gt.1) then
          kx1=ix1-1
          kx2=ix1
          kx3=ix1+1
        else
          kx1=ix1
          kx2=ix1+1
          kx3=ix1+2
        endif

        if (ny.gt.1) then
          if (iy1.gt.1) then
            ky1=iy1-1
            ky2=iy1
            ky3=iy1+1
          else
            ky1=iy1
            ky2=iy1+1
            ky3=iy1+2
          endif
        else
          ky1=1
          ky2=1
          ky3=1
        endif

        if (nz.gt.1) then
          if (iz1.gt.1) then
            kz1=iz1-1
            kz2=iz1
            kz3=iz1+1
          else
            kz1=iz1
            kz2=iz1+1
            kz3=iz1+2
          endif
        else
          kz1=1
          kz2=1
          kz3=1
        endif

        x3(1)=0.0d0

        x1=bmappe(1,kz1+(ky1-1)*nz+(kx1-1)*nyz)
        xx=x-x1

        y1=bmappe(2,kz1+(ky1-1)*nz+(kx1-1)*nyz)
        yy=y-y1
        dy=bmappe(2,kz1+(ky2-1)*nz+(kx1-1)*nyz)-y1

        z1=bmappe(3,kz1+(ky1-1)*nz+(kx1-1)*nyz)
        zz=z-z1
        dz=bmappe(3,kz2+(ky1-1)*nz+(kx1-1)*nyz)-z1

        do i=1,3

          x3(2)=bmappe(1,kz1+(ky1-1)*nz+(kx2-1)*nyz)-x1
          x3(3)=bmappe(1,kz1+(ky1-1)*nz+(kx3-1)*nyz)-x1

          b3(1)=bmappe(3+i,kz1+(ky1-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz1+(ky1-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz1+(ky1-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b11=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz1+(ky2-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz1+(ky2-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz1+(ky2-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b21=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz1+(ky3-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz1+(ky3-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz1+(ky3-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b31=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz2+(ky1-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz2+(ky1-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz2+(ky1-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b12=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz2+(ky2-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz2+(ky2-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz2+(ky2-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b22=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz2+(ky3-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz2+(ky3-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz2+(ky3-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b32=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz3+(ky1-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz3+(ky1-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz3+(ky1-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b13=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz3+(ky2-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz3+(ky2-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz3+(ky2-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b23=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=bmappe(3+i,kz3+(ky3-1)*nz+(kx1-1)*nyz)
          b3(2)=bmappe(3+i,kz3+(ky3-1)*nz+(kx2-1)*nyz)
          b3(3)=bmappe(3+i,kz3+(ky3-1)*nz+(kx3-1)*nyz)
          call parabel_short(x3,b3,a3)
          b33=a3(1)+(a3(2)+a3(3)*xx)*xx

          x3(2)=dy
          x3(3)=x3(2)+dy

          if (ny.gt.1) then

            b3(1)=b11
            b3(2)=b21
            b3(3)=b31
            call parabel_short(x3,b3,a3)
            bb1=a3(1)+(a3(2)+a3(3)*yy)*yy

            b3(1)=b12
            b3(2)=b22
            b3(3)=b32
            call parabel_short(x3,b3,a3)
            bb2=a3(1)+(a3(2)+a3(3)*yy)*yy

            b3(1)=b13
            b3(2)=b23
            b3(3)=b33
            call parabel_short(x3,b3,a3)
            bb3=a3(1)+(a3(2)+a3(3)*yy)*yy

          else
            bb3=b11
          endif

          if (nz.gt.1) then
            x3(2)=dz
            x3(3)=x3(2)+dz
            b3(1)=bb1
            b3(2)=bb2
            b3(3)=bb3
            call parabel_short(x3,b3,a3)
            b(i)=a3(1)+(a3(2)+a3(3)*zz)*zz
          else
            b(i)=bb1
          endif

        enddo !i=1,3

      endif !irfilb0

      bxout=b(1)
      byout=b(2)
      bzout=b(3)

      return
      end
