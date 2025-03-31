*CMZ :  4.01/04 15/11/2023  12.38.14  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.10.03  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/03 12/05/2010  13.34.28  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2009
      subroutine aphase(isour)
*KEEP,gplhint.
*KEND.

c calculates phase of field amplitudes afreq and afreqrphi

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,spectf90.
      include 'spectf90.cmn'
*KEND.

      double precision cenxexi,dphase,phiy,phiz,phiy1,phiz1
     &  ,ddist,h2,dist0,dist02,dphi

      integer iphi,ir,isour,ifreq,iobrp,idphi2pi

      cenxexi=(min(sourceeo(1,1,isour),xiend)
     &  +max(sourceao(1,1,isour),xianf))/2.d0

      if (ipin.ne.0) then
        dist0=pincen(1)-cenxexi
      else
        dist0=obs1x-cenxexi
      endif

      dist02=dist0**2

      do ifreq=1,nfreq

        if (mpinr.ne.0) then

          iobrp=0

          do iphi=1,nobsvphi

            phiy1=0.0d0
            phiz1=0.0d0

            do ir=1,nobsvr

              iobrp=iobrp+1
              iobfr=iobrp+nobsvrphi*(ifreq-1)

              h2=(obsvr(ir)/dist0)**2
              if (h2.lt.0.01) then
                ddist=dist0*(h2/2.0d0-h2**2/8.0d0)
              else
                ddist=dist0*(sqrt(1.0d0+h2)-1.0d0)
              endif

              dphase=ddist/freq(ifreq)*wtoe1*1.0d9*twopi1
              dphase=0.0d0

c              phiy=atan2(reaimarphi(2,2,iobfr),reaimarphi(2,1,iobfr))
c              phiz=atan2(reaimarphi(3,2,iobfr),reaimarphi(3,1,iobfr))
              if (reaimarphi(2,1,iobfr).ne.0.0d0) then
                phiy=atan(reaimarphi(2,2,iobfr)/reaimarphi(2,1,iobfr))
              else
                phiy=sign(halfpi1,reaimarphi(2,2,iobfr))
              endif

              if (reaimarphi(3,1,iobfr).ne.0.0d0) then
                phiz=atan(reaimarphi(3,2,iobfr)/reaimarphi(3,1,iobfr))
              else
                phiz=sign(halfpi1,reaimarphi(3,2,iobfr))
              endif

              reaimarphi(4,1,iobfr)=
     &          sqrt(reaimarphi(2,1,iobfr)**2+reaimarphi(2,2,iobfr)**2)
              reaimarphi(4,2,iobfr)=phiy

              reaimarphi(5,1,iobfr)=
     &          sqrt(reaimarphi(3,1,iobfr)**2+reaimarphi(3,2,iobfr)**2)
              reaimarphi(5,2,iobfr)=phiz

              phiz1=phiz

            enddo !nobsvr

          enddo !iphi

        endif !(mpinr.ne.0) then

      enddo !nfreq

      return
      end
