*CMZ :  4.00/15 24/03/2022  07.43.56  by  Michael Scheer
*CMZ :          24/03/2022  07.14.23  by  Michael Scheer
      subroutine bhalbasy2arg(xin,yin,zin,xcen,b0halbasy,ahwpol,
     &  zlhalbasy,xlhalbasy,bxout,byout,bzout,axout,ayout,azout)

*KEEP,gplhint.
*KEND.

c subroutine calculates magnetic field and vector potential for simple
c wavelength shifter model with end poles. the field of the single poles
c corresponds to halbachs formula.
c input and output correspond to lab.-system, where x is coordinate on
c longitudinal axis

      implicit none

      double precision xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  xkx,yky,zkz,dsnxkx,dcsxkx,dshyky,dchyky,dsnzkz,dcszkz
     &  ,bxh,byh,bzh,axh,ayh,azh,ahwmod,totlen,totlen2,ahwpol,
     &  x,x2,xcen,b0halbasy,
     &  zlhalbasy,ylhalbasy,xlhalbasy,
     &  zkhalbasy,ykhalbasy,xkhalbasy,
     &  zlhalbasy2,xlhalbasy2,
     &  zkhalbasy2,ykhalbasy2,xkhalbasy2

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

c--- k-values

      xkhalbasy=0.0d0
      ykhalbasy=0.0d0
      zkhalbasy=0.0d0

      if (zlhalbasy.ne.0.0d0) zkhalbasy=2.0d0*pi1/zlhalbasy
      if (ylhalbasy.ne.0.0d0) ykhalbasy=2.0d0*pi1/ylhalbasy
      if (xlhalbasy.ne.0.0d0) xkhalbasy=2.0d0*pi1/xlhalbasy

c--- adjust k-values

      ykhalbasy=dsqrt(zkhalbasy**2+xkhalbasy**2)
      ylhalbasy=2.0d0*pi1/ykhalbasy

c--- bending radius and device length

      totlen=zlhalbasy*((ahwpol-1.0d0)/2.0d0+1.0d0)
      totlen2=totlen/2.0d0

      totlen=zlhalbasy*((ahwpol-1.0d0)/2.0d0+1.0d0)
      totlen2=totlen/2.0d0

      x=xin-xcen
      if (dabs(xin).gt.totlen2) then
        bxout=0.0
        byout=0.0
        bzout=0.0
        axout=0.0
        ayout=0.0
        azout=0.0
        return
      endif

      if (dabs(x).le.totlen2-zlhalbasy/2.0d0) then

        xkx=xkhalbasy*(-zin)
        yky=ykhalbasy*yin
        zkz=zkhalbasy*x

        dsnxkx=dsin(xkx)
        dcsxkx=dcos(xkx)
        dshyky=dsinh(yky)
        dchyky=dsqrt(1.0d0+dshyky*dshyky)
        dsnzkz=dsin(zkz)
        dcszkz=dcos(zkz)

        bxh=-xkhalbasy/ykhalbasy*b0halbasy*dsnxkx*dshyky*dcszkz
        byh=                     b0halbasy*dcsxkx*dchyky*dcszkz
        bzh=-zkhalbasy/ykhalbasy*b0halbasy*dcsxkx*dshyky*dsnzkz

        axh=b0halbasy/zkhalbasy*                    dcsxkx*dchyky*dsnzkz
        ayh=b0halbasy/zkhalbasy*xkhalbasy/ykhalbasy*dsnxkx*dshyky*dsnzkz
        azh=0.0

        bzout=-bxh
        byout=byh
        bxout=bzh

        azout=-axh
        ayout=ayh
        axout=azh

      else

        ahwmod=-isign(1,-(mod(nint(ahwpol),4)-2))/2.0d0

        xkhalbasy2=xkhalbasy

        zkhalbasy2=zkhalbasy
        zlhalbasy2=2.0d0*pi1/zkhalbasy2
        ykhalbasy2=dsqrt(zkhalbasy2**2+xkhalbasy2**2)
        ylhalbasy=2.0d0*pi1/ykhalbasy2

        x2=x+totlen2+zlhalbasy/2.0d0

        xkx=xkhalbasy2*(-zin)
        yky=ykhalbasy2*yin
        zkz=zkhalbasy2*(x2)

        dsnxkx=dsin(xkx)
        dcsxkx=dcos(xkx)
        dshyky=dsinh(yky)
        dchyky=dsqrt(1.0d0+dshyky*dshyky)
        dsnzkz=dsin(zkz)
        dcszkz=dcos(zkz)

        bxh=-xkhalbasy2/ykhalbasy2*b0halbasy*dsnxkx*dshyky*dcszkz
        byh=                      b0halbasy*dcsxkx*dchyky*dcszkz
        bzh=-zkhalbasy2/ykhalbasy2*b0halbasy*dcsxkx*dshyky*dsnzkz

        axh=b0halbasy/zkhalbasy2*                    dcsxkx*dchyky*dsnzkz
        ayh=b0halbasy/zkhalbasy2*xkhalbasy2/ykhalbasy2*dsnxkx*dshyky*dsnzkz
        azh=0.0

        zkhalbasy2=zkhalbasy*2.0d0
        zlhalbasy2=2.0d0*pi1/zkhalbasy2
        ykhalbasy2=dsqrt(zkhalbasy2**2+xkhalbasy2**2)
        ylhalbasy=2.0d0*pi1/ykhalbasy2

        xkx=xkhalbasy2*(-zin)
        yky=ykhalbasy2*yin
        zkz=zkhalbasy2*(x2)

        dsnxkx=dsin(xkx)
        dcsxkx=dcos(xkx)
        dshyky=dsinh(yky)
        dchyky=dsqrt(1.0d0+dshyky*dshyky)
        dsnzkz=dsin(zkz)
        dcszkz=dcos(zkz)

        bxh=bxh-xkhalbasy2/ykhalbasy2*b0halbasy*dsnxkx*dshyky*dcszkz
        byh=byh+                      b0halbasy*dcsxkx*dchyky*dcszkz
        bzh=bzh-zkhalbasy2/ykhalbasy2*b0halbasy*dcsxkx*dshyky*dsnzkz

        axh=axh+b0halbasy/zkhalbasy2*                    dcsxkx*dchyky*dsnzkz
        ayh=ayh+b0halbasy/zkhalbasy2*xkhalbasy2/ykhalbasy2*dsnxkx*dshyky*dsnzkz
        azh=0.0

        bzout=bxh*ahwmod
        byout=-byh*ahwmod
        bxout=-bzh*ahwmod

        azout=axh*ahwmod
        ayout=-ayh*ahwmod
        axout=-azh*ahwmod

      endif

      return
      end
